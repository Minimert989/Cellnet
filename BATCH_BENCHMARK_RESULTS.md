# cellnetR batch benchmark

Measured on the local Apple Silicon host with R 4.5.1. The workload has a
100-node star network and one changed node per perturbation. Each row is one
wall-clock run of `perturb_batch()` and includes R argument/result conversion,
the native C++ batch call, and materialization of returned `CellNetResult`
frames. It is an end-to-end wrapper measurement, not a pure solver timing.

| Perturbations | Elapsed (s) | Per perturbation (µs) | Returned changed rows |
|---:|---:|---:|---:|
| 100 | 0.014 | 140.00 | 100 |
| 1,000 | 0.125 | 125.00 | 1,000 |
| 10,000 | 1.413 | 141.30 | 10,000 |
| 100,000 | 13.573 | 135.73 | 100,000 |

Raw values: [`benchmarks/batch_results.csv`](benchmarks/batch_results.csv).
This is one run per batch size, so it does not establish p50/p95 or a
native-vs-Python-vs-R comparison. The three optional data-source packages
OmniPathR, CellNOptR, and CARNIVAL plus DoRothEA/decoupleR were not installed on
this host; their conversion functions have fixture-level tests, not live-dataset
integration tests.
