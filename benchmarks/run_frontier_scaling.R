#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[[1L]] else "frontier_scaling.csv"
parse_csv_ints <- function(text, default) {
    if (is.null(text) || !nzchar(text)) return(default)
    values <- suppressWarnings(as.integer(strsplit(text, ",", fixed = TRUE)[[1L]]))
    if (!length(values) || anyNA(values) || any(values < 1L))
        stop("expected comma-separated positive integers")
    values
}
total_sizes <- parse_csv_ints(if (length(args) >= 2L) args[[2L]] else NULL,
                              c(1000L, 10000L, 100000L, 1000000L))
cone_sizes <- parse_csv_ints(if (length(args) >= 3L) args[[3L]] else NULL,
                             c(1L, 4L, 8L, 16L, 32L, 64L, 128L, 256L,
                               512L, 1000L, 2000L, 5000L, 10000L))
repetitions <- parse_csv_ints(if (length(args) >= 4L) args[[4L]] else NULL, 30L)[[1L]]
background_scc_sizes <- parse_csv_ints(if (length(args) >= 5L) args[[5L]] else NULL,
                                       cone_sizes)
if (repetitions < 30L)
    warning("fewer than 30 repetitions requested; tail quantiles will be noisy")

suppressPackageStartupMessages(library(cellnetR))
set.seed(20261010)

measure <- function(total_nodes, cone_nodes, background_scc_nodes) {
    nodes <- sprintf("N%07d", seq_len(total_nodes))
    cone <- nodes[seq_len(cone_nodes)]
    # The active directed ring is both the causal cone and the SCC reached by
    # the perturbation. An optional disconnected ring varies the network's
    # largest SCC independently, while remaining nodes are singleton SCCs.
    background_begin <- cone_nodes + 1L
    background_end <- if (background_scc_nodes > cone_nodes)
        background_begin + background_scc_nodes - 1L else cone_nodes
    background_ring <- if (background_end > cone_nodes)
        nodes[seq.int(background_begin, background_end)] else character()
    singleton_nodes <- if (background_end < total_nodes)
        nodes[seq.int(background_end + 1L, total_nodes)] else character()
    ring_edges <- function(ring) {
        if (!length(ring)) return(data.frame(source = character(), target = character()))
        data.frame(source = ring, target = c(ring[-1L], ring[[1L]]),
                   stringsAsFactors = FALSE)
    }
    edges <- rbind(ring_edges(cone), ring_edges(background_ring),
                   data.frame(source = singleton_nodes, target = singleton_nodes,
                              stringsAsFactors = FALSE))
    edges$sign <- 1
    edges$weight <- .1
    edges$confidence <- 1
    initial <- stats::setNames(rep(.5, total_nodes), nodes)
    network <- set_initial_state(as_cellnet(edges), initial)
    compiled <- cellnet_compile(network, dt = .05, horizon = 20L,
                                baseline_max_steps = 100L)

    frontier_us <- full_us <- numeric(repetitions)
    frontier_reactions <- full_reactions <- numeric(repetitions)
    frontier_sccs <- full_sccs <- numeric(repetitions)
    maximum_error <- 0
    for (i in seq_len(repetitions)) {
        target <- cone[[(i - 1L) %% cone_nodes + 1L]]
        effect <- runif(1L, .2, 1)
        frontier <- perturb(compiled, target, effect = effect, strategy = "frontier")
        full <- perturb(compiled, target, effect = effect, strategy = "full")
        frontier_us[[i]] <- attr(frontier, "execution_time_us")
        full_us[[i]] <- attr(full, "execution_time_us")
        frontier_reactions[[i]] <- attr(frontier, "reaction_evaluations")
        full_reactions[[i]] <- attr(full, "reaction_evaluations")
        frontier_sccs[[i]] <- attr(frontier, "scc_evaluations")
        full_sccs[[i]] <- attr(full, "scc_evaluations")
        fd <- stats::setNames(frontier$delta, frontier$node)
        rd <- stats::setNames(full$delta, full$node)
        ids <- union(names(fd), names(rd))
        fv <- fd[ids]
        rv <- rd[ids]
        fv[is.na(fv)] <- 0
        rv[is.na(rv)] <- 0
        maximum_error <- max(maximum_error, max(abs(fv - rv)))
    }
    quantile_row <- function(x, prefix) {
        q <- stats::quantile(x, c(.5, .95, .99), names = FALSE, type = 8)
        stats::setNames(as.list(q), paste0(prefix, c("_p50_us", "_p95_us", "_p99_us")))
    }
    c(list(total_nodes = total_nodes,
           largest_scc_nodes = max(cone_nodes, background_scc_nodes),
           affected_scc_nodes = cone_nodes,
           causal_cone_nodes = cone_nodes,
           repetitions = repetitions),
      quantile_row(frontier_us, "frontier"),
      quantile_row(full_us, "full"),
      list(frontier_reaction_evaluations_p50 = stats::median(frontier_reactions),
           full_reaction_evaluations_p50 = stats::median(full_reactions),
           frontier_scc_evaluations_p50 = stats::median(frontier_sccs),
           full_scc_evaluations_p50 = stats::median(full_sccs),
           maximum_frontier_full_delta_error = maximum_error))
}

grid <- expand.grid(total_nodes = total_sizes, cone_nodes = cone_sizes,
                    background_scc_nodes = background_scc_sizes,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$required_nodes <- ifelse(grid$background_scc_nodes > grid$cone_nodes,
                              grid$cone_nodes + grid$background_scc_nodes,
                              grid$cone_nodes)
grid <- grid[grid$required_nodes <= grid$total_nodes, , drop = FALSE]
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
    total <- grid$total_nodes[[i]]
    cone <- grid$cone_nodes[[i]]
    background_scc <- grid$background_scc_nodes[[i]]
    cat(sprintf("total=%d largest-SCC=%d cone=%d (%d/%d)\n", total,
                max(cone, background_scc), cone, i, nrow(grid)))
    rows[[i]] <- as.data.frame(measure(total, cone, background_scc), check.names = FALSE)
    if (rows[[i]]$maximum_frontier_full_delta_error >= 1e-5)
        stop("frontier/full accuracy threshold exceeded")
    # Persist each completed scenario so a long 1M-node grid can be stopped
    # without losing prior measurements.
    utils::write.csv(do.call(rbind, rows[seq_len(i)]), out_path, row.names = FALSE)
}

cat("Wrote", out_path, "\n")
