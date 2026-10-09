#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args)) args[[1]] else "integrations/cellnetR/benchmarks/batch_results.csv"
suppressPackageStartupMessages(library(cellnetR))

n_nodes <- 100L
node_ids <- sprintf("N%03d", seq_len(n_nodes))
edges <- data.frame(source = node_ids[[1]], target = node_ids[-1], sign = 1,
                    weight = 1, confidence = 1, stringsAsFactors = FALSE)
network <- set_initial_state(as_cellnet(edges), stats::setNames(rep(1, n_nodes), node_ids))
compiled <- cellnet_compile(network, dt = .05, horizon = 10L)
targets <- node_ids[-1]
counts <- c(100L, 1000L, 10000L, 100000L)
rows <- vector("list", length(counts))

for (i in seq_along(counts)) {
    n <- counts[[i]]
    batch_targets <- rep(targets, length.out = n)
    invisible(gc())
    timing <- system.time(result <- perturb_batch(compiled, batch_targets, effect = -1))
    rows[[i]] <- data.frame(
        api = "R_to_native_perturb_batch",
        total_nodes = n_nodes,
        causal_cone_nodes_per_target = 1L,
        perturbations = n,
        elapsed_seconds = unname(timing[["elapsed"]]),
        microseconds_per_perturbation = unname(timing[["elapsed"]]) * 1e6 / n,
        changed_result_rows = sum(vapply(result, nrow, integer(1))),
        stringsAsFactors = FALSE
    )
    rm(result)
    invisible(gc())
    cat(sprintf("%d perturbations: %.3f s (%.2f us/perturbation)\n",
                n, rows[[i]]$elapsed_seconds,
                rows[[i]]$microseconds_per_perturbation))
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(do.call(rbind, rows), out_path, row.names = FALSE)
cat("Wrote", out_path, "\n")
