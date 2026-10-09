# cellnetR

This is a separate R package directory for the CellNet C++ causal execution
backend. It is intentionally kept under `integrations/cellnetR/`; the C++ core
and public headers are not modified by this integration. A source snapshot of
the core is vendored under `src/` and can be refreshed from the monorepo with
`sh integrations/cellnetR/tools/sync-core.sh` from the repository root.

## Install from this checkout

After publication, install directly from GitHub with:

```r
install.packages("remotes")
remotes::install_github("Minimert989/Cellnet")
```

To install from a local source checkout:

```r
install.packages("Rcpp")
Rcpp::compileAttributes("integrations/cellnetR")
R CMD INSTALL integrations/cellnetR
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
objects inside C++. CellNOptR fitting remains outside CellNet. Solver semantics and reference
implementation are unchanged. The package includes the vendored snapshot so
it can build from a source tarball without reaching outside this directory.
Refresh that snapshot before release and review the copied C++ diff.

Optional adapters are tested only when their packages are installed. See
`vignettes/bioconductor-workflows.Rmd` and `tests/testthat/`.

The CellNet R package and its vendored CellNet core snapshot are distributed
under the MIT License. The full terms are in `COPYING`; `LICENSE` also carries
the package copyright metadata.

The native batch scaling harness is `benchmarks/run_batch_benchmark.R`; it
measures 100, 1k, 10k, and 100k perturbations on a 100-node network with a
one-node affected cone per target.
