#'
#' Plot geneset enrichment from CREATE dds
#' 
#' @import DESeq2
#' @import dplyr
#' @import fgsea
#' @import ggplot2
#' 
#' @param dds A DESeqDataSet object from DESeqDataSetFromCreate
#' @param gene_ids A character vector of gene IDs or a logical vector indicating which genes to plot
#' @param counts Either "raw" or "normalized" counts to plot
#' @param contrast A vector specifying the contrast for DESeq2 results, or a single factor to plot all levels
#' @param color_by The column in colData to color the boxplots by - defaults to gene IDs
#' @param show_ns Whether to show non-significant comparisons [default = TRUE]
#' @param step The step size for positioning significance annotations [default = 0.1]
#' @param nrow The number of rows in the facet wrap - defaults to a square
#' @param size Base font size for the plot
#' 
#' @return A ggplot2 object
#' @export
plot_expression <- function(dds, gene_ids, counts="normalized", contrast = NA, color_by=NA, show_ns = TRUE, step = 0.1, nrow = NA, size = 8, labeller = NA) {
    if (class(gene_ids) == "logical") {
        if (length(gene_ids) != nrow(dds)) {
            stop("If gene_ids is a logical vector, it must be the same length as the number of rows in dds")
        }
        gene_ids <- rownames(dds)[gene_ids]
    } else if (class(gene_ids) != "character") {
        stop("gene_ids must be a character vector or a logical vector")
    }

    if (is.na(nrow)) {
        nrow <- ceiling(sqrt(length(gene_ids)))
    }

    if (!(counts %in% c("raw", "normalized"))) {
        stop("counts must be either 'raw' or 'normalized'")
    }
    
    message("Extracting expression data...")
    if (counts == 'raw') {
        expr <- SummarizedExperiment::assays(dds)$counts[gene_ids, , drop = FALSE]
    } else {
        expr <- DESeq2::counts(dds, normalized=TRUE)[gene_ids, , drop = FALSE]
    }

    expr <- as.data.frame(t(expr))
    expr[["sample"]] <- rownames(expr)
    expr <- tidyr::pivot_longer(expr, cols = -sample, names_to = "gene_id", values_to = "count")
    expr <- dplyr::left_join(expr, as.data.frame(SummarizedExperiment::colData(dds)), by = c("sample" = "names"))
    expr <- dplyr::left_join(expr, as.data.frame(SummarizedExperiment::rowData(dds)), by = "gene_id")

    if (counts == "normalized") {
        message("Calculating significance...")
        if (length(contrast) == 1) {
            comps <- expand.grid(levels(expr[[contrast[1]]]), levels(expr[[contrast[1]]]))
            comps <- dplyr::filter(comps, as.numeric(Var1) > as.numeric(Var2))
            message(paste0("...found ", nrow(comps), " comparisons to test"))
            signif <- apply(comps, 1, function(comp) {
                res <- DESeq2::results(dds, contrast = c(contrast[1], comp[1], comp[2]))
                res <- as.data.frame(res)
                res <- dplyr::mutate(res, gene_id = rownames(res))
                res <- dplyr::filter(res, gene_id %in% gene_ids)
                res <- dplyr::mutate(res, 
                    numer = factor(comp[1], levels = levels(expr[[contrast[1]]])), 
                    denom = factor(comp[2], levels = levels(expr[[contrast[1]]]))
                )
                return(res)
            })
            signif <- do.call(rbind, signif)
            if (!show_ns) {
                signif <- dplyr::filter(signif, padj < 0.05)
            }
            signif <- dplyr::left_join(signif, expr[,c('gene_id', 'count')], by='gene_id', relationship='many-to-many')
            signif <- dplyr::summarize(signif, ymax = max(count), .by=colnames(signif)[colnames(signif) != "count"])
            signif <- dplyr::group_by(signif, gene_id)
            signif <- dplyr::mutate(signif, ymax = max(ymax))
            signif <- dplyr::arrange(signif, dplyr::desc(padj))
            signif <- dplyr::mutate(signif, ypos = ymax * (1 + step * (dplyr::row_number())))
        } else if (length(contrast) == 3) {
            message(paste0("Comparing ", contrast[2], " vs ", contrast[3], " for factor ", contrast[1]))
            expr <- dplyr::filter(expr, .data[[contrast[1]]] %in% c(contrast[2], contrast[3]))
            signif <- DESeq2::results(dds, contrast = contrast)
            signif <- as.data.frame(signif)
            signif <- dplyr::filter(signif, rownames(signif) %in% gene_ids)
            signif <- dplyr::mutate(signif, gene_id = rownames(signif))
            if (!show_ns) {
                signif <- dplyr::filter(signif, padj < 0.05)
            }
            signif <- dplyr::mutate(
                signif,
                numer = factor(contrast[2], levels=levels(expr[[contrast[1]]])),
                denom = factor(contrast[3], levels = levels(expr[[contrast[1]]]))
            )
            signif <- dplyr::left_join(signif, expr[,c('gene_id', 'count')], by='gene_id', relationship='many-to-many')
            signif <- dplyr::summarize(signif, ymax = max(count), .by=colnames(signif)[colnames(signif) != "count"])
            signif <- dplyr::group_by(signif, gene_id)
            signif <- dplyr::mutate(signif, ymax = max(ymax))
            signif <- dplyr::arrange(signif, dplyr::desc(padj))
            signif <- dplyr::mutate(signif, ypos = ymax * 1.05)
        } else {
            stop("contrast must be either of length 1 or length 3")
        }
        signif <- dplyr::ungroup(signif)
        signif <- dplyr::mutate(signif, signif = ifelse(is.na(padj), 'ns', ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "ns")))))
        signif <- dplyr::left_join(signif, as.data.frame(SummarizedExperiment::rowData(dds))[,c('gene_id', 'gene_name', color_by)], by = "gene_id")
        signif <- dplyr::mutate(signif, xstart = as.numeric(numer), xend = as.numeric(denom))
    } else {
        contrast <- c('gene_name')
    }

    if (is.na(color_by)) {
        color_by <- contrast[1]
    }
    
    message("Creating plot...")
    expr_plot <- ggplot2::ggplot(expr, ggplot2::aes(y = count, x = .data[[contrast[1]]], fill=.data[[color_by]])) +
        ggplot2::geom_boxplot() +
        ggplot2::ylab(ifelse(counts == "raw", "Raw Counts", "Normalized Counts")) +
        ggplot2::xlab("Sample") +
        theme_create(base_size = size)

    if (counts == "normalized") {
        message("...adding significance annotations")
        expr_plot <- expr_plot +
            ggplot2::facet_wrap(~gene_name, scales = "free_y", axes="all") +
            ggplot2::geom_segment(
                data = signif, ggplot2::aes(
                    x = xstart, xend = xend,
                    y = ypos, yend = ypos
                ), colour='black'
            ) +
            ggplot2::geom_text(
                data = signif, ggplot2::aes(
                    label = signif, y = ypos,
                    x = (xstart + xend) / 2,
                ),
                colour='black',
                vjust = -0.15, hjust = 0.5
            ) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(step, step)))
    }

    message("Done!")
    return(expr_plot)
}
