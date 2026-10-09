# Frontier locality scaling benchmark

Measured with R 4.5.1 on Apple Silicon. The main table reports native-engine
latencies (`execution_time_us`). Each cell is the
distribution over 30 paired frontier/full perturbations at `dt = 0.05` and
horizon 20. Random effect strengths prevent repeated identical requests from
reusing the result cache. The directed active ring is the perturbed SCC and
causal cone; a disconnected background ring independently varies the largest
SCC in the network. All cases had maximum frontier/full delta error 0.

| Total nodes | Largest SCC | Affected SCC / cone | Frontier p50/p95/p99 (µs) | Full p50/p95/p99 (µs) |
|---:|---:|---:|---:|---:|
| 1,000 | 39 | 1 | 0.167 / 1.056 / 4.917 | 14.105 / 16.669 / 16.917 |
| 10,000 | 39 | 1 | 0.209 / 2.369 / 4.459 | 217.792 / 225.122 / 226.292 |
| 1,000 | 39 | 39 | 8.792 / 10.425 / 10.708 | 260.667 / 268.875 / 269.583 |
| 10,000 | 39 | 39 | 19.188 / 22.647 / 27.500 | 3,822.938 / 3,864.260 / 3,900.916 |
| 1,000 | 717 | 1 | 0.167 / 0.263 / 0.334 | 10.167 / 10.652 / 11.042 |
| 10,000 | 717 | 1 | 0.208 / 1.864 / 2.042 | 205.584 / 213.886 / 217.250 |
| 1,000 | 717 | 39 | 8.792 / 9.063 / 9.417 | 193.625 / 199.186 / 204.250 |
| 10,000 | 717 | 39 | 22.688 / 26.475 / 27.041 | 3,902.583 / 3,968.265 / 3,969.292 |
| 100,000 | 717 | 1 | 1.230 / 8.788 / 18.917 | 2,103.792 / 2,167.594 / 2,168.834 |
| 100,000 | 717 | 39 | 24.813 / 30.897 / 35.042 | 45,953.792 / 47,663.386 / 49,737.916 |
| 100,000 | 717 | 717 | 299.417 / 321.902 / 331.500 | 41,536.625 / 45,088.083 / 45,265.166 |
| 1,000,000 | 7,171 | 1 | 1.938 / 7.087 / 18.916 | 22,113.771 / 24,641.145 / 31,674.541 |
| 1,000,000 | 7,171 | 39 | 23.417 / 27.667 / 37.583 | 409,762.750 / 442,951.983 / 459,601.500 |

Raw measurements: [`benchmarks/frontier_scaling_1k_10k.csv`](benchmarks/frontier_scaling_1k_10k.csv)
[`benchmarks/frontier_scaling_100k.csv`](benchmarks/frontier_scaling_100k.csv),
and [`benchmarks/frontier_scaling_1m.csv`](benchmarks/frontier_scaling_1m.csv).
All 30-repetition rows had maximum frontier/full delta error 0. At 1M nodes,
the 7,171-node SCC was also tested as an unaffected background component.
Making it the affected cone completed once with 3,774 µs frontier, 434,616 µs
full, and zero delta error; that single run is not a percentile estimate. The
30-repetition version of this largest-cone case was stopped after 14 minutes.
Its one-run result is preserved separately at
[`benchmarks/frontier_scaling_1m_worstcase_single.csv`](benchmarks/frontier_scaling_1m_worstcase_single.csv).
The benchmark runner writes each completed scenario incrementally. Run
`benchmarks/run_frontier_scaling.R` with larger total-node and SCC/cone sizes
before making claims at those scales. Large SCCs that are themselves affected
remain frontier's worst case, as documented in the README.

## R wrapper wall time after native causal-distance traversal

The previous `.new_result()` calculated causal distances in R by repeatedly
scanning the complete edge table once per BFS depth. It now consumes distances
computed in C++ from the compiled dependency CSR. The C++ handle reuses its
visited-generation arrays and BFS queue, and stops once all returned changed
nodes have been found. This removes the `BFS depth × total edges` wrapper path.

On the 1M-node network (largest SCC 7,171; perturbed cone 39), 30 paired calls
reported frontier native-engine p50/p95/p99 of 25.3/33.0/40.1 µs, while the
measured end-to-end R `perturb()` wall time was 1,000/1,150/2,000 µs. Full
strategy measured 415,876/818,879/862,241 µs natively and
416,500/819,650/863,000 µs end-to-end. The R wall clock has millisecond-scale
quantization at this duration, and includes data-frame creation and attribute
handling. Frontier/full maximum delta error was zero. Thus native compute is
not the entire R call cost, but the wrapper no longer repeatedly scans one
million edges per BFS level.

The 1-node cone rows are self-loop/timer-floor controls, not representative
frontier performance claims; use cone sizes 39 and above to assess useful
workloads. The benchmark currently uses synthetic rings. The repository has
the small A375 model but no large curated signaling-network topology suitable
for a representative cone-distribution benchmark, so that biological-topology
measurement remains outstanding.

Raw end-to-end measurements:
[`benchmarks/frontier_scaling_1m_walltime.csv`](benchmarks/frontier_scaling_1m_walltime.csv).
