.adapter_table <- function(x, source_aliases, target_aliases, sign_aliases = character()) {
    if (inherits(x, "CellNetNetwork")) return(x$edges)
    if (inherits(x, "igraph") || inherits(x, "graphNEL") || inherits(x, "dgCMatrix"))
        return(as_cellnet(x)$edges)
    if (isS4(x)) {
        for (slot_name in c("interMat", "interactions", "edges", "data")) {
            if (slot_name %in% methods::slotNames(x)) {
                candidate <- methods::slot(x, slot_name)
                if (is.data.frame(candidate)) { x <- candidate; break }
            }
        }
    }
    if (is.list(x) && !is.data.frame(x)) {
        for (key in c("interactions", "edges", "network", "results"))
            if (!is.null(x[[key]])) { x <- x[[key]]; break }
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    find_name <- function(candidates) {
        hit <- candidates[candidates %in% names(x)]
        if (length(hit)) hit[[1]] else NULL
    }
    source_col <- find_name(source_aliases)
    target_col <- find_name(target_aliases)
    sign_col <- find_name(sign_aliases)
    if (is.null(source_col) || is.null(target_col))
        stop("could not identify source/target columns in adapter input")
    if (!"sign" %in% names(x) && all(c("is_stimulation", "is_inhibition") %in% names(x))) {
        x$sign <- ifelse(as.logical(x$is_inhibition), -1,
                         ifelse(as.logical(x$is_stimulation), 1, NA_real_))
    }
    names(x)[match(source_col, names(x))] <- "source"
    names(x)[match(target_col, names(x))] <- "target"
    if (!"sign" %in% names(x) && !is.null(sign_col)) names(x)[match(sign_col, names(x))] <- "sign"
    if (!"sign" %in% names(x)) x$sign <- 1
    if (!"confidence" %in% names(x)) {
        confidence_col <- find_name(c("confidence", "dorothea_confidence"))
        if (!is.null(confidence_col)) x$confidence <- x[[confidence_col]]
    }
    if (!"source_db" %in% names(x)) {
        db_col <- find_name(c("sources", "databases", "source_database"))
        if (!is.null(db_col)) x$source_db <- as.character(x[[db_col]])
    }
    if (!"references" %in% names(x)) {
        ref_col <- find_name(c("references", "pmids", "pubmed", "reference"))
        if (!is.null(ref_col)) x$references <- as.character(x[[ref_col]])
    }
    if (!"weight" %in% names(x)) x$weight <- 1
    as_cellnet.data.frame(x)$edges
}

#' Import OmniPath signed interactions while retaining source/reference columns
#' @export
cellnet_from_omnipath <- function(interactions) {
    if (missing(interactions)) {
        .require_namespace("OmnipathR")
        interactions <- OmnipathR::import_all_interactions()
    }
    edges <- .adapter_table(interactions,
                            c("source", "source_genesymbol", "source_uniprot"),
                            c("target", "target_genesymbol", "target_uniprot"),
                            c("sign", "is_stimulation", "is_inhibition"))
    if (all(c("is_stimulation", "is_inhibition") %in% names(interactions))) {
        edges$sign <- ifelse(as.logical(interactions$is_inhibition), -1, 1)
        edges$sign[is.na(interactions$is_stimulation) & is.na(interactions$is_inhibition)] <- NA
        if (anyNA(edges$sign)) stop("OmniPath contains unsigned interactions; filter or annotate them before compilation")
    }
    as_cellnet.data.frame(edges)
}

#' Import OmniPath PTM enzyme-substrate relations
#' @param interactions PTM table, usually OmniPathR::import_post_translational()
#' @param sign Optional explicit signed effect, scalar or row-wise. Required when
#'   the table does not provide stimulation/inhibition direction.
#' @export
cellnet_from_omnipath_ptm <- function(interactions, sign = NULL) {
    x <- as.data.frame(interactions, stringsAsFactors = FALSE)
    source_col <- intersect(c("enzyme", "source", "enzyme_genesymbol", "enzyme_uniprot"), names(x))
    target_col <- intersect(c("substrate", "target", "substrate_genesymbol", "substrate_uniprot"), names(x))
    if (!length(source_col) || !length(target_col)) stop("PTM table needs enzyme and substrate identifiers")
    out <- x
    out$source <- as.character(x[[source_col[[1]]]])
    out$target <- as.character(x[[target_col[[1]]]])
    if (is.null(sign)) {
        if (all(c("is_stimulation", "is_inhibition") %in% names(x))) {
            out$sign <- ifelse(x$is_inhibition, -1, ifelse(x$is_stimulation, 1, NA_real_))
            if (anyNA(out$sign)) stop("unsigned PTM rows require an explicit sign mapping")
        } else stop("PTM does not imply activating direction; provide an explicit sign mapping")
    } else out$sign <- rep_len(sign, nrow(out))
    if (!"mechanism" %in% names(out) && "modification" %in% names(out)) out$mechanism <- out$modification
    as_cellnet.data.frame(out)
}

#' Import CellNOptR prior-knowledge network (does not perform model fitting)
#' @export
cellnet_from_cellnoptr <- function(pkn) {
    edges <- .adapter_table(pkn, c("source", "from", "Node1", "node1"),
                            c("target", "to", "Node2", "node2"),
                            c("sign", "interaction", "signal"))
    if ("interaction" %in% names(edges)) edges$sign <- ifelse(as.numeric(edges$interaction) < 0, -1, 1)
    as_cellnet.data.frame(edges)
}

#' Import signed CARNIVAL network inference results
#' @export
cellnet_from_carnival <- function(result) {
    edges <- .adapter_table(result, c("source", "Node1", "node1", "from"),
                            c("target", "Node2", "node2", "to"),
                            c("sign", "Sign", "signaling", "interaction"))
    as_cellnet.data.frame(edges)
}

#' Convert a graphite pathway object to causal edges
#' @export
cellnet_from_graphite <- function(pathway, sign = NULL) {
    if (inherits(pathway, "graphNEL")) return(as_cellnet.graphNEL(pathway, sign = sign))
    if (inherits(pathway, "igraph")) {
        if (is.null(sign)) return(as_cellnet(pathway))
        edges <- igraph::as_data_frame(pathway, what = "edges")
        names(edges)[seq_len(2)] <- c("source", "target")
        edges$sign <- rep_len(sign, nrow(edges))
        return(as_cellnet.data.frame(edges))
    }
    if (!is.data.frame(pathway) && !is.list(pathway) && requireNamespace("graphite", quietly = TRUE)) {
        pathway_graph <- graphite::pathwayGraph(pathway)
        if (is.null(sign)) stop("graphite pathway topology is unsigned here; provide an explicit sign mapping")
        return(as_cellnet.graphNEL(pathway_graph, sign = sign))
    }
    pathway_table <- as.data.frame(pathway)
    type_col <- intersect(c("sign", "type", "interaction"), names(pathway_table))
    encoded_direction <- FALSE
    if (length(type_col)) {
        values <- pathway_table[[type_col[[1]]]]
        if (is.numeric(values)) {
            encoded_direction <- all(is.finite(values) & values != 0)
        } else {
            encoded_direction <- all(tolower(as.character(values)) %in%
                                         c("+", "+1", "1", "activation", "stimulation", "activating",
                                           "-", "-1", "inhibition", "inhibitory", "inhibiting"))
        }
    }
    if (is.null(sign) && !encoded_direction)
        stop("graphite pathway is unsigned; provide an explicit sign mapping")
    edges <- .adapter_table(pathway, c("source", "from", "src", "src_name"),
                            c("target", "to", "dest", "dest_name"),
                            if (is.null(sign)) c("sign", "type", "interaction") else "sign")
    if (!is.null(sign)) edges$sign <- rep_len(sign, nrow(edges))
    if (anyNA(edges$sign)) stop("graphite edge direction/sign is ambiguous; provide an explicit mapping")
    as_cellnet.data.frame(edges)
}

#' Import DoRothEA TF-target regulons using signed mode of regulation
#' @export
cellnet_from_dorothea <- function(regulon) {
    x <- as.data.frame(regulon, stringsAsFactors = FALSE)
    tf <- intersect(c("tf", "source", "TF"), names(x))
    target <- intersect(c("target", "target_gene", "Target"), names(x))
    mor <- intersect(c("mor", "mode_of_regulation", "sign"), names(x))
    if (!length(tf) || !length(target) || !length(mor)) stop("DoRothEA table needs tf, target, and mor")
    out <- x
    out$source <- as.character(x[[tf[[1]]]])
    out$target <- as.character(x[[target[[1]]]])
    out$sign <- sign(as.numeric(x[[mor[[1]]]]))
    confidence <- intersect(c("confidence", "dorothea_confidence", "confidence_level"), names(x))
    if (length(confidence)) out$confidence <- x[[confidence[[1]]]]
    as_cellnet.data.frame(out)
}
