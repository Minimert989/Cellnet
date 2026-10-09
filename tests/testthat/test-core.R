test_that("edge-list conversion preserves signs and provenance", {
    edges <- data.frame(source = c("A", "B"), target = c("B", "C"),
                        sign = c(1, -1), source_db = "fixture",
                        mechanism = c("activation", "inhibition"), confidence = c(.9, .8))
    net <- as_cellnet(edges)
    expect_equal(net$edges$sign, c(1, -1))
    expect_equal(net$edges$source_db, c("fixture", "fixture"))
    expect_equal(net$nodes, c("A", "B", "C"))
})

test_that("compiled engine runs perturbations and native batches", {
    edges <- data.frame(source = c("A", "B"), target = c("B", "C"),
                        sign = c(1, -1), weight = 1)
    network <- set_initial_state(as_cellnet(edges), c(A = 0, B = 0, C = 1))
    compiled <- cellnet_compile(network, dt = .05, horizon = 100L)
    result <- perturb(compiled, "A", effect = 1)
    expect_s3_class(result, "CellNetResult")
    expect_true(any(result$node == "B" & result$delta > 0))
    expect_true(any(result$node == "C" & result$delta < 0))
    expect_equal(result$causal_distance[match(c("B", "C"), result$node)], c(1L, 2L))
    edge_result <- perturb_reaction(compiled, reaction_id = 1, effect = -1)
    expect_s3_class(edge_result, "CellNetResult")
    batch <- perturb_batch(compiled, c("A", "B"), effect = -1)
    expect_length(batch, 2)
    expect_s3_class(batch[[1]], "CellNetResult")
    expect_error(perturb(compiled, "missing"), "unknown target node")
    expect_error(perturb(compiled, c("A", "B")), "target")
    expect_error(perturb(compiled, NA_character_), "target")
})

test_that("native causal-distance traversal handles cyclic dependency graphs", {
    edges <- data.frame(source = c("A", "B", "C"),
                        target = c("B", "C", "A"), sign = 1)
    compiled <- cellnet_compile(set_initial_state(as_cellnet(edges),
                                                  c(A = 0, B = 0, C = 0)))
    result <- perturb(compiled, "A", effect = 1)
    expect_equal(result$causal_distance[match(c("A", "B", "C"), result$node)],
                 c(0L, 1L, 2L))
})

test_that("reaction inhibition is not a node clamp", {
    edges <- data.frame(source = "A", target = "B", sign = 1)
    compiled <- cellnet_compile(set_initial_state(as_cellnet(edges), c(A = 1, B = 1)))
    result <- perturb_reaction(compiled, 1, effect = -1)
    expect_true(any(result$node == "B" & result$delta < 0))
    expect_false(any(result$node == "A" & abs(result$delta) > 1e-7))
})

test_that("causal cone and sparse round trip retain direction", {
    edges <- data.frame(source = c("A", "B"), target = c("B", "C"), sign = c(1, -1))
    network <- as_cellnet(edges)
    expect_equal(causal_cone(network, "A"), c("A", "B", "C"))
    if (requireNamespace("Matrix", quietly = TRUE)) {
        restored <- as_cellnet(as_sparse_matrix(network))
        expect_setequal(paste(restored$edges$source, restored$edges$target, restored$edges$sign),
                        paste(edges$source, edges$target, edges$sign))
    }
})

test_that("PTM adapter refuses unsigned enzyme-substrate edges", {
    ptm <- data.frame(enzyme = "MAP2K1", substrate = "MAPK1", modification = "phosphorylation")
    expect_error(cellnet_from_omnipath_ptm(ptm), "explicit sign mapping")
    net <- cellnet_from_omnipath_ptm(ptm, sign = 1)
    expect_equal(net$edges$mechanism, "phosphorylation")
})

test_that("OmniPath direction and protocol sampling are explicit", {
    op <- data.frame(source = c("A", "B"), target = c("B", "C"),
                     is_stimulation = c(TRUE, FALSE), is_inhibition = c(FALSE, TRUE),
                     sources = c("SIGNOR", "SIGNOR"), references = c("1", "2"))
    net <- cellnet_from_omnipath(op)
    expect_equal(net$edges$sign, c(1, -1))
    expect_equal(net$edges$references, c("1", "2"))
    compiled <- cellnet_compile(set_initial_state(net, c(A = 0, B = 0, C = 1)))
    trace <- run_protocol(compiled, data.frame(time = 0, target = "A", effect = 1),
                          end_time = 1, sample_times = c(0, .5, 1))
    expect_equal(dim(trace$states), c(3, 3))
    expect_equal(trace$times, c(0, .5, 1), tolerance = 1e-7)
})

test_that("igraph and graphNEL adapters work when available", {
    if (requireNamespace("igraph", quietly = TRUE)) {
        g <- igraph::graph_from_data_frame(data.frame(from = c("A", "B"),
                                                       to = c("B", "C"), sign = c(1, -1)))
        expect_equal(as_cellnet(g)$edges$sign, c(1, -1))
    }
    if (requireNamespace("graph", quietly = TRUE)) {
        g <- graph::graphNEL(nodes = c("A", "B"), edgeL = list(A = "B", B = character()),
                             edgemode = "directed")
        expect_error(as_cellnet(g), "causal sign")
        expect_equal(nrow(as_cellnet(g, sign = 1)$edges), 1)
    }
})

test_that("curated network adapters retain signed direction", {
    pkn <- data.frame(Node1 = "A", Node2 = "B", interaction = -1)
    expect_equal(cellnet_from_cellnoptr(pkn)$edges$sign, -1)
    carnival <- data.frame(Node1 = "A", Node2 = "B", Sign = -1)
    expect_equal(cellnet_from_carnival(carnival)$edges$sign, -1)
    pathway <- data.frame(source = "A", target = "B", type = "binding")
    expect_error(cellnet_from_graphite(pathway), "unsigned")
    expect_equal(cellnet_from_graphite(pathway, sign = 1)$edges$sign, 1)
    regulon <- data.frame(tf = "A", target = "B", mor = -1, confidence = "B")
    converted <- cellnet_from_dorothea(regulon)
    expect_equal(converted$edges$sign, -1)
    expect_equal(converted$edges$confidence, 1)
    expect_equal(converted$edges$confidence_label, "B")
})

test_that("decoupleR-style activities map into selected node context", {
    activities <- data.frame(source = c("A", "A", "B"),
                             condition = c("control", "treated", "treated"),
                             score = c(0.4, 0.8, 0.5))
    ctx <- cellnet_context(activities, condition = "treated", transform = "positive")
    expect_equal(unname(ctx), c(.8, .5))
    expect_equal(names(ctx), c("A", "B"))
    expect_error(cellnet_context(c(A = 2), transform = "positive"), "within")
})

test_that("context resets baseline before perturbation and enforces multiplier range", {
    edges <- data.frame(source = "A", target = "B", sign = 1)
    compiled <- cellnet_compile(set_initial_state(as_cellnet(edges), c(A = 0, B = 0)))
    ctx <- cellnet_context(c(A = 0.5), transform = "identity")
    expect_silent(set_context(compiled, ctx))
    expect_s3_class(perturb(compiled, "A", effect = 1), "CellNetResult")
    expect_error(set_context(compiled, c(A = 2)), "within")
    expect_error(cellnet_context(c(A = 2), transform = "positive"), "within")
})

test_that("baseline convergence controls are exposed and failures are diagnostic", {
    edges <- data.frame(
        source = c("n0", "n0", "n1"),
        target = c("n1", "n1", "n0"),
        sign = 1, weight = c(.5, .5, 1), confidence = c(1, 1, .8)
    )
    network <- set_initial_state(as_cellnet(edges), c(n0 = .789, n1 = .065))
    expect_error(
        cellnet_compile(network, baseline_max_steps = 1L),
        "baseline did not converge after 1/1 steps.*last max delta.*required"
    )
    compiled <- cellnet_compile(network, baseline_max_steps = 10000L)
    expect_s3_class(compiled, "CellNetCompiled")
    expect_error(cellnet_compile(network, epsilon = 0), "epsilon must be finite and positive")
    expect_error(cellnet_compile(network, baseline_max_steps = 1.5), "positive integer")
})

test_that("frontier and full agree on the reported cyclic multi-input fixture", {
    edges <- data.frame(
        source = c("n1", "n0", "n2", "n2"),
        target = c("n0", "n1", "n1", "n0"),
        sign = c(-1, 1, 1, -1), weight = 1
    )
    network <- set_initial_state(as_cellnet(edges), c(n0 = .25, n1 = .75, n2 = .75))
    compiled <- cellnet_compile(network, dt = .05, horizon = 100L)
    frontier <- perturb(compiled, "n2", effect = 1, strategy = "frontier")
    full <- perturb(compiled, "n2", effect = 1, strategy = "full")
    frontier_delta <- stats::setNames(frontier$delta, frontier$node)
    full_delta <- stats::setNames(full$delta, full$node)
    all_nodes <- union(names(frontier_delta), names(full_delta))
    frontier_delta <- frontier_delta[all_nodes]
    full_delta <- full_delta[all_nodes]
    frontier_delta[is.na(frontier_delta)] <- 0
    full_delta[is.na(full_delta)] <- 0
    expect_lt(max(abs(frontier_delta - full_delta)), 1e-5)
    expect_lt(abs(full$perturbed[match("n0", full$node)] - .2), 1e-5)
    expect_lt(abs(full$perturbed[match("n1", full$node)] - .6), 1e-5)
})

test_that("frontier and full agree across seeded random cyclic multi-input networks", {
    set.seed(20261009)
    for (case in seq_len(12L)) {
        node_count <- sample(2:10, 1L)
        nodes <- paste0("v", case, "_", seq_len(node_count))
        ring_sources <- nodes
        ring_targets <- c(nodes[-1L], nodes[[1L]])
        extra_sources <- vapply(seq_along(nodes), function(i) {
            sample(setdiff(nodes, nodes[[i]]), 1L)
        }, character(1))
        edges <- data.frame(
            source = c(ring_sources, extra_sources),
            target = c(ring_targets, nodes),
            sign = sample(c(-1, 1), 2L * node_count, replace = TRUE),
            weight = runif(2L * node_count, .2, .4),
            confidence = runif(2L * node_count, .5, 1)
        )
        edges <- edges[sample(seq_len(nrow(edges))), , drop = FALSE]
        # Uniform half-activity is an equilibrium for signed Modulate terms,
        # avoiding baseline-transient failures in this scheduler regression.
        initial <- stats::setNames(rep(.5, node_count), nodes)
        compiled <- cellnet_compile(
            set_initial_state(as_cellnet(edges), initial), dt = .05, horizon = 100L
        )
        target <- sample(nodes, 1L)
        frontier <- perturb(compiled, target, effect = 1, strategy = "frontier")
        full <- perturb(compiled, target, effect = 1, strategy = "full")
        frontier_delta <- stats::setNames(frontier$delta, frontier$node)
        full_delta <- stats::setNames(full$delta, full$node)
        all_nodes <- union(names(frontier_delta), names(full_delta))
        frontier_delta <- frontier_delta[all_nodes]
        full_delta <- full_delta[all_nodes]
        frontier_delta[is.na(frontier_delta)] <- 0
        full_delta[is.na(full_delta)] <- 0
        expect_lt(max(abs(frontier_delta - full_delta)), 1e-5)
    }
})
