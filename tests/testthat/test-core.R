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
