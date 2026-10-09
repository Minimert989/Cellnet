\name{cellnetR}
\alias{cellnetR}
\alias{as_cellnet}
\alias{as_cellnet.data.frame}
\alias{as_cellnet.igraph}
\alias{as_cellnet.dgCMatrix}
\alias{as_cellnet.graphNEL}
\alias{as.data.frame.CellNetNetwork}
\alias{cellnet_compile}
\alias{perturb}
\alias{perturb_reaction}
\alias{perturb_batch}
\alias{run_protocol}
\alias{causal_cone}
\alias{affected_nodes}
\alias{as_igraph}
\alias{as_sparse_matrix}
\alias{set_initial_state}
\alias{set_context}
\alias{cellnet_context}
\alias{set_activity_context}
\alias{cellnet_from_omnipath}
\alias{cellnet_from_omnipath_ptm}
\alias{cellnet_from_cellnoptr}
\alias{cellnet_from_carnival}
\alias{cellnet_from_graphite}
\alias{cellnet_from_dorothea}
\alias{print.CellNetResult}
\alias{summary.CellNetResult}
\alias{as.data.frame.CellNetResult}
\alias{plot.CellNetResult}
\title{R interface to the CellNet causal execution engine}
\description{
Compile signed causal networks once, retain the compiled C++ engine behind an
external pointer, then execute single or batch perturbations from R. Input
conversion helpers support edge tables and optional graph packages. Provenance
columns remain available in the R-side network object.
}
\details{
The minimum edge schema is \code{source}, \code{target}, and signed direction
\code{sign} (\code{1} activation, \code{-1} inhibition). Optional \code{weight},
\code{mechanism}, \code{source_db}, \code{references}, and \code{confidence}
are retained. \code{perturb()} uses negative effects for activity inhibition
and positive effects for activation. \code{perturb_batch()} passes its target
set through a single native call. PTM enzyme-substrate relationships require an
explicit signed effect when direction is not available in the source data.

\code{cellnet_compile()} currently lowers signed node-level edges to native
CellNet relaxation reactions; it does not model arbitrary logical gates or
perform CellNOptR fitting. Identifier mapping is not implicit. See the package
README and vignette for integration boundaries.
}
\usage{
as_cellnet(x, ...)
cellnet_compile(network, dt = 0.05, horizon = 100L, epsilon = 1e-06,
  baseline_max_steps = 2000L)
perturb(net, target, effect = -1, strategy = c("frontier", "full"), context = NULL)
perturb_reaction(net, reaction_id, effect = -1, strategy = c("frontier", "full"))
perturb_batch(net, targets, effect = -1, strategy = c("frontier", "full"))
run_protocol(net, schedule, end_time, sample_times, strategy = c("frontier", "full"))
causal_cone(network, seeds, include_self = TRUE)
affected_nodes(network, target)
as_igraph(x)
as_sparse_matrix(x)
set_initial_state(network, values)
set_context(net, context)
cellnet_context(x, cell = NULL, assay = 1L, transform = c("identity", "positive"), condition = NULL)
set_activity_context(net, activities, condition = NULL)
cellnet_from_omnipath(interactions)
cellnet_from_omnipath_ptm(interactions, sign = NULL)
cellnet_from_cellnoptr(pkn)
cellnet_from_carnival(result)
cellnet_from_graphite(pathway, sign = NULL)
cellnet_from_dorothea(regulon)
}
\arguments{
\item{x, network, net, pathway, interactions, pkn, result, regulon}{A supported
network/input object or a network returned by \code{as_cellnet()}.}
\item{...}{Additional format-specific options.}
\item{dt}{Positive solver step in model time units.}
\item{horizon}{Positive number of solver steps per perturbation.}
\item{epsilon}{Positive convergence tolerance. Baseline convergence requires
the maximum per-step state delta to be at most \code{epsilon * min(dt, 1)}.}
\item{baseline_max_steps}{Maximum number of full-network integration steps
allowed for baseline equilibrium construction.}
\item{target, targets, seeds}{One or more node identifiers.}
\item{reaction_id}{1-based row index into the compiled edge table.}
\item{effect}{Nonzero signed perturbation strength in \code{[-1,1]}.}
\item{strategy}{\code{"frontier"} or \code{"full"} execution strategy.}
\item{context, activities}{Named, finite, nonnegative node multipliers. A
decoupleR-style table may provide \code{source}, \code{score}, and optionally
\code{condition}.}
\item{schedule}{Data frame with \code{time}, \code{target}, and \code{effect}.}
\item{end_time}{Protocol end time.}
\item{sample_times}{Times at which full state vectors are returned.}
\item{include_self}{Whether seed nodes are included in the causal cone.}
\item{values}{Named initial-state values in \code{[0,1]}.}
\item{cell}{Selected sample/cell for matrix or assay input.}
\item{condition}{Condition selector for activity result tables.}
\item{assay}{Assay name or index in a SummarizedExperiment.}
\item{transform}{Activity transform; positive clips negative values to zero.}
\item{sign}{Explicit activation/inhibition mapping for unsigned pathway/PTM input.}
}
\value{
Conversion functions return a \code{CellNetNetwork}; \code{cellnet_compile()}
returns a \code{CellNetCompiled} object with a persistent native pointer.
Perturbation functions return \code{CellNetResult} data frames. Protocol runs
return sample times and a time-by-node state matrix.
}
\seealso{\code{README.md}, \code{vignettes/bioconductor-workflows.Rmd}}
\keyword{models}
