# Frontier locality smoke benchmark

Measured with R 4.5.1 on Apple Silicon. These are native-engine-reported
latencies (`execution_time_us`), not R wrapper wall times. Each cell is the
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

Raw measurements: [`benchmarks/frontier_scaling_1k_10k.csv`](benchmarks/frontier_scaling_1k_10k.csv)
and [`benchmarks/frontier_scaling_100k.csv`](benchmarks/frontier_scaling_100k.csv).
All reported rows had maximum frontier/full delta error 0. The 1M grid has not
been run yet. Run
`benchmarks/run_frontier_scaling.R` with larger total-node and SCC/cone sizes
before making claims at those scales. Large SCCs that are themselves affected
remain frontier's worst case, as documented in the README.
