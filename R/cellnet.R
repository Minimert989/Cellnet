#' Create an R-side CellNet causal network
#' @param x edge table, igraph, Matrix::dgCMatrix, or graphNEL.
#' @param ... format-specific options.
#' @return A CellNetNetwork containing signed edges and preserved provenance.
#' @export
as_cellnet <- function(x, ...) UseMethod("as_cellnet")

#' @export
as_cellnet.data.frame <- function(x, ...) {
    edges <- .normalize_edges(x)
    structure(list(edges = edges,
                   nodes = unique(c(edges$source, edges$target)),
                   initial_state = stats::setNames(rep(0, length(unique(c(edges$source, edges$target)))),
                                                   unique(c(edges$source, edges$target))),
                   metadata = list(input_class = class(x))),
              class = "CellNetNetwork")
}

#' @export
as_cellnet.igraph <- function(x, ...) {
    .require_namespace("igraph")
    vertices <- igraph::as_data_frame(x, what = "vertices")
    edges <- igraph::as_data_frame(x, what = "edges")
    names(edges)[seq_len(2)] <- c("source", "target")
    edges$source <- as.character(edges$source)
    edges$target <- as.character(edges$target)
    if (!"sign" %in% names(edges) && "interaction" %in% names(edges)) edges$sign <- edges$interaction
    if (!"sign" %in% names(edges)) stop("igraph has no signed edge attribute; add sign before conversion")
    if (!"weight" %in% names(edges)) edges$weight <- 1
    if ("name" %in% names(vertices)) attr(edges, "vertex_names") <- as.character(vertices$name)
    as_cellnet.data.frame(edges)
}

#' @export
as_cellnet.dgCMatrix <- function(x, ...) {
    .require_namespace("Matrix")
    if (is.null(rownames(x)) || is.null(colnames(x)))
        stop("dgCMatrix input must have row/column names (source rows, target columns)")
    triplets <- Matrix::summary(x)
    edges <- data.frame(source = rownames(x)[triplets$i],
                        target = colnames(x)[triplets$j],
                        sign = sign(triplets$x), weight = abs(triplets$x),
                        stringsAsFactors = FALSE)
    as_cellnet.data.frame(edges)
}

#' @export
as_cellnet.graphNEL <- function(x, sign = NULL, ...) {
    .require_namespace("graph")
    node_names <- graph::nodes(x)
    adjacency <- graph::edges(x)
    rows <- lapply(names(adjacency), function(source) {
        targets <- adjacency[[source]]
        if (!length(targets)) return(NULL)
        data.frame(source = source, target = targets, sign = 1,
                   stringsAsFactors = FALSE)
    })
    edges <- do.call(rbind, rows)
    if (is.null(edges)) edges <- data.frame(source = character(), target = character(), sign = numeric())
    if (is.null(sign)) stop("graphNEL topology has no causal sign; provide sign = 1 or -1 explicitly")
    edges$sign <- rep_len(sign, nrow(edges))
    out <- as_cellnet.data.frame(edges)
    out$nodes <- unique(c(node_names, out$nodes))
    out$initial_state <- stats::setNames(rep(0, length(out$nodes)), out$nodes)
    out$metadata$input_class <- "graphNEL"
    out
}

.normalize_edges <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    nms <- names(x)
    if (!all(c("source", "target") %in% nms)) stop("edge table requires source and target columns")
    x$source <- as.character(x$source)
    x$target <- as.character(x$target)
    if ("sign" %in% nms) {
        raw_sign <- x$sign
    } else if (all(c("is_stimulation", "is_inhibition") %in% nms)) {
        raw_sign <- ifelse(as.logical(x$is_inhibition), -1, 1)
        raw_sign[is.na(x$is_stimulation) & is.na(x$is_inhibition)] <- NA
    } else {
        stop("edge table requires sign, or Omnipath is_stimulation/is_inhibition columns")
    }
    if (is.numeric(raw_sign)) {
        x$sign <- sign(as.numeric(raw_sign))
    } else {
        label <- tolower(as.character(raw_sign))
        x$sign <- ifelse(label %in% c("+", "+1", "1", "activation", "stimulation", "activating"), 1,
                         ifelse(label %in% c("-", "-1", "inhibition", "inhibitory", "inhibiting"), -1, NA_real_))
    }
    if (anyNA(x$sign) || any(x$sign == 0)) stop("sign must encode activation (+1) or inhibition (-1)")
    if (!"weight" %in% names(x)) x$weight <- 1
    x$weight <- abs(as.numeric(x$weight))
    if (any(!is.finite(x$weight))) stop("weight must be finite")
    if (!"confidence" %in% names(x)) x$confidence <- 1
    if (is.character(x$confidence)) {
        original_confidence <- x$confidence
        numeric_confidence <- suppressWarnings(as.numeric(original_confidence))
        if (all(is.finite(numeric_confidence))) {
            x$confidence <- numeric_confidence
        } else {
            x$confidence_label <- original_confidence
            x$confidence <- rep(1, nrow(x))
        }
    }
    x$confidence <- as.numeric(x$confidence)
    if (anyNA(x$confidence) || any(!is.finite(x$confidence)) || any(x$confidence < 0 | x$confidence > 1))
        stop("numeric confidence must be finite and in [0,1]")
    for (field in c("mechanism", "source_db", "references", "original_source", "original_target"))
        if (!field %in% names(x)) x[[field]] <- NA_character_
    if (any(!nzchar(x$source)) || any(!nzchar(x$target))) stop("node identifiers cannot be empty")
    x
}

#' @export
as.data.frame.CellNetNetwork <- function(x, row.names = NULL, optional = FALSE, ...) x$edges

#' Compile a signed network once into a persistent native CellNet engine
#' @param network CellNetNetwork or a supported network object.
#' @param dt timestep in model time units.
#' @param horizon number of solver steps per perturbation.
#' @param epsilon positive convergence tolerance; baseline convergence uses
#'   epsilon * min(dt, 1) as its maximum per-step state delta.
#' @param baseline_max_steps maximum full-network steps used to construct the
#'   baseline equilibrium.
#' @export
cellnet_compile <- function(network, dt = 0.05, horizon = 100L,
                            epsilon = 1e-6, baseline_max_steps = 2000L) {
    if (!inherits(network, "CellNetNetwork")) network <- as_cellnet(network)
    if (!is.numeric(dt) || length(dt) != 1 || !is.finite(dt) || dt <= 0) stop("dt must be positive")
    if (length(horizon) != 1 || is.na(horizon) || horizon < 1) stop("horizon must be positive")
    if (!is.numeric(epsilon) || length(epsilon) != 1 || !is.finite(epsilon) || epsilon <= 0)
        stop("epsilon must be finite and positive")
    if (!is.numeric(baseline_max_steps) || length(baseline_max_steps) != 1 ||
        !is.finite(baseline_max_steps) || baseline_max_steps < 1 ||
        baseline_max_steps > .Machine$integer.max || baseline_max_steps != floor(baseline_max_steps))
        stop("baseline_max_steps must be a positive integer no larger than .Machine$integer.max")
    edges <- network$edges
    nodes <- network$nodes
    initial <- unname(network$initial_state[nodes])
    ptr <- cn_compile(edges, nodes, as.numeric(initial), as.numeric(dt), as.integer(horizon),
                      as.numeric(epsilon), as.integer(baseline_max_steps))
    structure(list(ptr = ptr, nodes = nodes, edges = edges,
                   provenance = edges[, intersect(c("source", "target", "mechanism", "source_db",
                                                    "references", "confidence", "original_source",
                                                    "original_target"), names(edges)), drop = FALSE],
                   initial_state = network$initial_state, dt = dt, horizon = as.integer(horizon),
                   epsilon = epsilon, baseline_max_steps = as.integer(baseline_max_steps),
                   metadata = network$metadata), class = "CellNetCompiled")
}

#' Execute a single causal perturbation
#' @param net compiled CellNet network.
#' @param target node identifier.
#' @param effect signed effect in [-1,1]. Negative inhibits; positive activates.
#' @param strategy frontier or full reference execution.
#' @export
perturb <- function(net, target, effect = -1, strategy = c("frontier", "full"), context = NULL) {
    stopifnot(inherits(net, "CellNetCompiled"))
    strategy <- match.arg(strategy)
    if (!is.null(context)) set_context(net, context)
    if (length(target) != 1 || is.na(target)) stop("target must be one node identifier")
    target <- as.character(target)
    if (length(effect) != 1 || !is.finite(effect) || effect == 0 || abs(effect) > 1)
        stop("effect must be nonzero and within [-1,1]")
    .new_result(cn_perturb(net$ptr, target, effect, strategy), net)
}

#' Inhibit or activate one specific edge reaction
#' @param reaction_id 1-based row index in the compiled edge table.
#' @export
perturb_reaction <- function(net, reaction_id, effect = -1,
                             strategy = c("frontier", "full")) {
    stopifnot(inherits(net, "CellNetCompiled"))
    strategy <- match.arg(strategy)
    if (length(reaction_id) != 1 || is.na(reaction_id) || reaction_id < 1L ||
        reaction_id > nrow(net$edges)) stop("reaction_id must index a compiled edge")
    if (length(effect) != 1 || !is.finite(effect) || effect == 0 || abs(effect) > 1)
        stop("effect must be nonzero and within [-1,1]")
    .new_result(cn_perturb_reaction(net$ptr, as.integer(reaction_id), effect, strategy), net)
}

#' Execute many perturbations through one R-to-native batch call
#' @export
perturb_batch <- function(net, targets, effect = -1, strategy = c("frontier", "full")) {
    stopifnot(inherits(net, "CellNetCompiled"))
    strategy <- match.arg(strategy)
    targets <- as.character(targets)
    if (any(!targets %in% net$nodes)) stop("all targets must be known nodes")
    raw <- cn_perturb_batch(net$ptr, targets, effect, strategy)
    Map(function(data) .new_result(data, net), raw)
}

.new_result <- function(data, net) {
    # Native materialization computes distances with the compiled dependency CSR.
    # Keep the R wrapper O(number of returned changes), not O(depth * all edges).
    data$affected <- !is.na(data$causal_distance)
    class(data) <- c("CellNetResult", class(data))
    attr(data, "network") <- net
    data
}

#' @export
print.CellNetResult <- function(x, ...) {
    cat(sprintf("CellNetResult: %d changed nodes; %d reaction evaluations; %s us\n",
                nrow(x), attr(x, "reaction_evaluations"),
                format(attr(x, "execution_time_us"), digits = 4)))
    invisible(x)
}

#' @export
summary.CellNetResult <- function(object, ...) {
    data.frame(changed_nodes = nrow(object),
               affected_nodes = sum(object$affected),
               reaction_evaluations = attr(object, "reaction_evaluations"),
               scc_evaluations = attr(object, "scc_evaluations"),
               execution_time_us = attr(object, "execution_time_us"))
}

#' @export
as.data.frame.CellNetResult <- function(x, row.names = NULL, optional = FALSE, ...) {
    x[, c("node", "baseline", "perturbed", "delta", "affected", "causal_distance"), drop = FALSE]
}

#' @export
plot.CellNetResult <- function(x, top = 20L, ...) {
    if (!nrow(x)) return(invisible(graphics::plot.new()))
    x <- x[order(abs(x$delta), decreasing = TRUE), , drop = FALSE]
    x <- utils::head(x, top)
    graphics::barplot(x$delta, names.arg = x$node, las = 2,
                      col = ifelse(x$delta < 0, "#4477AA", "#CC6677"), ...)
}

#' Return downstream causal cone from one or more seed nodes
#' @export
causal_cone <- function(network, seeds, include_self = TRUE) {
    edges <- if (inherits(network, "CellNetCompiled") || inherits(network, "CellNetNetwork")) network$edges else as_cellnet(network)$edges
    seeds <- as.character(seeds)
    known_nodes <- if (inherits(network, "CellNetCompiled") || inherits(network, "CellNetNetwork"))
        network$nodes else unique(c(edges$source, edges$target))
    unknown <- setdiff(seeds, known_nodes)
    if (length(unknown)) stop("unknown seed: ", unknown[[1]])
    seen <- seeds
    frontier <- seeds
    while (length(frontier)) {
        next_nodes <- unique(edges$target[edges$source %in% frontier])
        frontier <- setdiff(next_nodes, seen)
        seen <- c(seen, frontier)
    }
    if (include_self) seen else setdiff(seen, seeds)
}

#' @export
affected_nodes <- function(network, target) causal_cone(network, target)

#' Convert CellNet network to igraph
#' @export
as_igraph <- function(x) {
    .require_namespace("igraph")
    edges <- if (inherits(x, "CellNetCompiled") || inherits(x, "CellNetNetwork")) x$edges else as_cellnet(x)$edges
    igraph::graph_from_data_frame(edges, directed = TRUE)
}

#' Convert CellNet network to Matrix::dgCMatrix signed adjacency
#' @export
as_sparse_matrix <- function(x) {
    .require_namespace("Matrix")
    edges <- if (inherits(x, "CellNetCompiled") || inherits(x, "CellNetNetwork")) x$edges else as_cellnet(x)$edges
    nodes <- if (inherits(x, "CellNetCompiled") || inherits(x, "CellNetNetwork")) x$nodes else unique(c(edges$source, edges$target))
    Matrix::sparseMatrix(i = match(edges$source, nodes), j = match(edges$target, nodes),
                         x = edges$sign * edges$weight, dims = c(length(nodes), length(nodes)),
                         dimnames = list(nodes, nodes))
}

#' @export
set_initial_state <- function(network, values) {
    if (!inherits(network, "CellNetNetwork")) network <- as_cellnet(network)
    if (is.null(names(values)) || any(!names(values) %in% network$nodes))
        stop("initial state must be a named vector of known node IDs")
    if (any(!is.finite(values)) || any(values < 0 | values > 1)) stop("initial values must be in [0,1]")
    network$initial_state[names(values)] <- values
    network
}

#' Apply nonnegative node-context multipliers to the compiled engine
#' @export
set_context <- function(net, context) {
    stopifnot(inherits(net, "CellNetCompiled"))
    if (is.data.frame(context)) context <- stats::setNames(context$multiplier, context$node)
    if (is.null(names(context)) || any(!names(context) %in% net$nodes)) stop("context must be a named node multiplier vector")
    if (any(!is.finite(context)) || any(context < 0 | context > 1))
        stop("context multipliers must be finite and within [0,1]")
    cn_set_context(net$ptr, names(context), as.numeric(context))
    invisible(net)
}

#' Convert a named activity vector or omics sample to a context vector
#' @export
cellnet_context <- function(x, cell = NULL, assay = 1L, transform = c("identity", "positive"),
                            condition = NULL) {
    transform <- match.arg(transform)
    if (is.data.frame(x) && all(c("source", "score") %in% names(x))) {
        if ("condition" %in% names(x)) {
            conditions <- unique(as.character(x$condition))
            if (length(conditions) > 1 && is.null(condition))
                stop("select one decoupleR condition explicitly")
            if (!is.null(condition)) x <- x[as.character(x$condition) == condition, , drop = FALSE]
        }
        if (anyDuplicated(x$source)) stop("activity table has repeated source IDs after condition selection")
        x <- stats::setNames(x$score, x$source)
    } else if (inherits(x, "SummarizedExperiment")) {
        if (is.null(cell)) stop("cell/sample must be selected for SummarizedExperiment input")
        mat <- SummarizedExperiment::assay(x, assay)
        x <- mat[, cell]
    } else if (is.matrix(x) || inherits(x, "Matrix")) {
        if (is.null(cell)) stop("column/cell must be selected for matrix input")
        x <- x[, cell]
    }
    if (is.null(names(x))) stop("context input must have node identifiers as names")
    if (transform == "positive") {
        node_names <- names(x)
        x <- pmax(0, x)
        names(x) <- node_names
    }
    if (any(!is.finite(x)) || any(x < 0 | x > 1))
        stop("engine context accepts finite multipliers within [0,1]")
    x
}

#' Supply activity inference output as node context (not topology)
#' @export
set_activity_context <- function(net, activities, condition = NULL) {
    set_context(net, cellnet_context(activities, transform = "positive", condition = condition))
}

#' Execute scheduled signed interventions and sample full network state
#' @export
run_protocol <- function(net, schedule, end_time, sample_times, strategy = c("frontier", "full")) {
    stopifnot(inherits(net, "CellNetCompiled"))
    strategy <- match.arg(strategy)
    schedule <- as.data.frame(schedule, stringsAsFactors = FALSE)
    if (!all(c("time", "target", "effect") %in% names(schedule))) stop("schedule needs time, target, effect")
    out <- cn_run_protocol(net$ptr, schedule, end_time, as.numeric(sample_times), strategy)
    structure(out, class = "CellNetProtocolResult")
}

.require_namespace <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) stop("package '", package, "' is required for this adapter")
}
