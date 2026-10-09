# cellnetR

`cellnetR` is an R package for signed causal-network perturbation analysis,
backed by a vendored C++ execution engine. The repository is a standalone R
package checkout; install it directly from the repository root.

## Install from this checkout

After publication, install directly from GitHub with:

```r
install.packages("remotes")
remotes::install_github("Minimert989/Cellnet")
```

To install from a local source checkout:

```sh
R -q -e 'install.packages("Rcpp")'
R CMD INSTALL .
```

Optional adapters require their corresponding R/Bioconductor packages. The
package itself has only Rcpp as a required R dependency.

## Minimal example

```r
library(cellnetR)
edges <- data.frame(source = c("EGFR", "MAP2K1"),
                    target = c("MAP2K1", "MAPK1"),
                    sign = c(1, 1), mechanism = c("activation", "phosphorylation"),
                    source_db = "example", confidence = 1)
net <- as_cellnet(edges)
net <- set_initial_state(net, c(EGFR = 0, MAP2K1 = 0, MAPK1 = 0))
compiled <- cellnet_compile(net, dt = 0.05, horizon = 100)
res <- perturb(compiled, "EGFR", effect = 1)
summary(res)
batch <- perturb_batch(compiled, c("EGFR", "MAP2K1"), effect = -1)
edge_inhibition <- perturb_reaction(compiled, reaction_id = 1, effect = -1)
```

`effect < 0` means activity inhibition and `effect > 0` means activity
activation. Sign is the causal direction; `weight` controls reaction rate and
confidence is retained as model confidence/provenance. `perturb_batch()` sends
all requested targets through one native API call.
`perturb_reaction()` applies inhibition or activation to one edge reaction without
clamping its target node.

## Conversion boundary

- Edge tables require `source`, `target`, `sign`; optional `weight`, `mechanism`,
  `source_db`, `references`, and `confidence` are retained.
- `igraph`, `Matrix::dgCMatrix` (source rows, target columns), and `graphNEL`.
- OmniPath interactions/PTM, CellNOptR/CARNIVAL/graphite tables, and DoRothEA
  regulons through explicit adapters.
- Activity vectors or selected `SummarizedExperiment`/
  `SingleCellExperiment` columns are context, never topology.

Context values are interpreted as node multipliers in `[0,1]`; values outside
that range are rejected rather than silently rescaled. Applying a context
invalidates and rebuilds the engine baseline before the next perturbation.

PTM membership alone does not imply that phosphorylation activates a substrate.
The PTM adapter requires a signed effect in the source data or an explicit sign
mapping. Identifier conversion (HGNC/UniProt/Ensembl/Entrez) is never implicit.
Categorical confidence tiers are retained as labels and are not converted to
arbitrary numeric kinetic multipliers.

## Scope and limitations

The first binding lowers signed node-level edges to CellNet relaxation
reactions. It preserves edge provenance in R but does not yet encode arbitrary
Boolean AND gates, scheduled reaction-level protocol events, or full omics
objects inside C++. CellNOptR fitting remains outside CellNet.

The `frontier` strategy schedules work at SCC granularity. It is most useful
when a perturbation affects a small causal region of the SCC condensation DAG.
If the affected region is a large strongly connected component (SCC), the
engine must still evaluate much of that component; frontier may then approach
or be slower than `full`. This is an expected limitation of SCC-granularity
locality, not a guarantee of network-size-independent latency. Benchmarks
should report SCC sizes and affected-cone size alongside timings.

The source snapshot is vendored so the package builds independently. Core
maintainers should synchronize it from the matching engine source and run the
package and C++ test suites before publishing a release.

Optional adapters are tested only when their packages are installed. See
`vignettes/bioconductor-workflows.Rmd` and `tests/testthat/`.

The CellNet R package and its vendored CellNet core snapshot are distributed
under the MIT License. The full terms are in `COPYING`; `LICENSE` also carries
the package copyright metadata.

The native batch scaling harness is `benchmarks/run_batch_benchmark.R`; it
measures 100, 1k, 10k, and 100k perturbations on a 100-node network with a
one-node affected cone per target. It is a wrapper/batch benchmark, not a
locality scaling benchmark.

`benchmarks/run_frontier_scaling.R` compares frontier and full p50/p95/p99
latency while varying total node count and a directed-ring SCC that is also
the causal cone. Its output records total nodes, largest SCC, cone size,
reaction/SCC evaluation counts, and maximum frontier/full delta error. For a
small smoke run:

```sh
Rscript benchmarks/run_frontier_scaling.R /tmp/frontier-smoke.csv 100,500 1,4,16 30
```

The default grid includes 1k, 10k, 100k, and 1M total nodes and cone/SCC sizes
from 1 through 10k. An optional fifth argument sets the size of a disconnected
background SCC (for example
`Rscript benchmarks/run_frontier_scaling.R /tmp/out.csv 10000 1,100 30 39,7171`),
so largest-SCC size can be varied independently of the perturbed SCC/cone. The
background SCC is not in the perturbation's causal cone. Large grids can
require substantial time and memory.

A 1k/10k, selected 100k, and selected 1M, 30-repetition results are in
`FRONTIER_SCALING_RESULTS.md` with raw CSVs under `benchmarks/`. The 1M grid
is partial: a 7,171-node affected SCC did not finish 30 paired full/frontier
repetitions within 14 minutes and is called out in the report.
