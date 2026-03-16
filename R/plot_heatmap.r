#'
#' Plot a clustered heatmap of expression data
#'
#' @param dds A DESeqDataSet object from DESeqDataSetFromCreate
#' @param genes A character vector of gene IDs or a logical vector indicating which genes to plot
#' @param base_size Base font size for the plot [default=8]
#' @param samples A character or logical vector of samples to include [default=NA, all samples]
#' @param gene_label Column in rowData to use as row labels [default="gene_id"]
#' @param show_rownames Whether to show row (gene) labels [default=TRUE]
#' @param show_colnames Whether to show column (sample) labels [default=FALSE]
#' @param colors A length-3 vector of colors for the gradient (low, mid, high) [default=c("red","white","blue")]
#' @param scale Scaling method: "log" (log10+1), "logz" (log10+1 then z-score), or "z" (z-score) [default="logz"]
#' @param col_annotations A data.frame of column annotations, rownames must match sample names [default=NA]
#' @param row_annotations A data.frame of row annotations, rownames must match gene IDs [default=NA]
#' @param annotation_colors A named list of color vectors for each annotation column [default=NA]
#' @param row_annotation_scale Relative width of row annotation tiles [default=0.0625]
#' @param col_annotation_scale Relative height of column annotation tiles [default=0.0625]
#' @param row_dendro_scale Relative width of the row dendrogram [default=0.125]
#' @param col_dendro_scale Relative height of the column dendrogram [default=0.125]
#' @param labeller A named character vector mapping gene IDs to display labels [default=NA]
#' @param filename Path to save the plot; if NA the plot is returned only [default=NA]
#' @param assemble If TRUE (default), return the assembled patchwork. If FALSE, return a list
#'   with `$plots`, `$design`, `$widths`, and `$heights` so the panels can be integrated into
#'   an outer patchwork at the same level (required for `guides = 'collect'` to work across
#'   the boundary).
#'
#' @return A patchwork plot object, or a list of components when `assemble = FALSE`
#'
#' @export
plot_heatmap <- function(
    dds, genes, base_size = 8,
    samples = NA, gene_label = 'gene_id',
    show_rownames = TRUE, show_colnames = FALSE,
    colors = c('red', 'white', 'blue'), scale = 'logz',
    col_annotations = NA, row_annotations = NA, annotation_colors = NA,
    row_annotation_scale = 0.0625, col_annotation_scale = 0.0625,
    row_dendro_scale = 0.125, col_dendro_scale = 0.125,
    labeller = NA, filename = NA, assemble = TRUE
) {
    ## Order leaves intelligently
    treeorder <- function(hclust, values) {
        reorder_node <- function(node, hclust, values) {
            if (node < 0) return(node)

            left_child <- hclust$merge[node, 1]
            left_child <- reorder_node(left_child, hclust, values)
            left_leaves <- get_leaves(left_child, hclust)
            left_mean <- mean(values[left_leaves])

            right_child <- hclust$merge[node, 2]
            right_child <- reorder_node(right_child, hclust, values)
            right_leaves <- get_leaves(right_child, hclust)
            right_mean <- mean(values[right_leaves])

            if (left_mean > right_mean) {
                hclust$merge[node, 1] <<- right_child
                hclust$merge[node, 2] <<- left_child
            } else {
                hclust$merge[node, 1] <<- left_child
                hclust$merge[node, 2] <<- right_child
            }

            return(node)
        }

        get_leaves <- function(node, hclust) {
            if (node < 0) return(-node)

            c(
                get_leaves(hclust$merge[node, 1], hclust),
                get_leaves(hclust$merge[node, 2], hclust)
            )
        }

        reorder_node(nrow(hclust$merge), hclust, values)
        hclust$order <- get_leaves(nrow(hclust$merge), hclust)
        return(hclust)
    }

    if (class(genes) != 'logical') {
        genes <- rownames(dds) %in% genes
    }

    ## Make primary data
    mat <- DESeq2::counts(dds, normalized = TRUE)
    mat <- as.data.frame(mat)
    mat <- dplyr::filter(mat, genes)
    if (!any(is.na(samples))) {
        if (!is.logical(samples)) {
            samples <- colnames(mat) %in% samples
        }
        mat <- mat[, samples, drop = FALSE]
    }
    if (scale == 'log') {
        mat <- mat + 1
        mat <- log10(mat)
    } else if (scale == 'logz') {
        mat <- mat + 1
        mat <- log10(mat)
        mat <- sweep(mat, 1, rowMeans(mat), '-')
        mat <- sweep(mat, 1, matrixStats::rowSds(as.matrix(mat)), '/')
    } else if (scale == 'z') {
        mat <- sweep(mat, 1, rowMeans(mat), '-')
        mat <- sweep(mat, 1, matrixStats::rowSds(as.matrix(mat)), '/')
    } else {
        stop("Invalid scale option")
    }
    mat <- as.matrix(mat)

    ## Make column dendrograms
    col_hclust <- stats::hclust(stats::dist(t(mat)))
    col_hclust <- treeorder(col_hclust, matrixStats::colMedians(mat))
    col_order <- col_hclust$order
    col_hclust <- stats::as.dendrogram(col_hclust)
    col_hclust <- ggdendro::dendro_data(col_hclust)

    ## Make row dendrograms
    row_hclust <- stats::hclust(stats::dist(mat))
    row_hclust <- treeorder(row_hclust, matrixStats::rowMedians(mat))
    row_order <- row_hclust$order
    row_hclust <- stats::as.dendrogram(row_hclust)
    row_hclust <- ggdendro::dendro_data(row_hclust)

    ## Reorder matrix for plotting
    mat <- mat[row_order, col_order]

    mat_melt <- reshape2::melt(as.matrix(mat))
    colnames(mat_melt) <- c('Row', 'Col', 'Value')
    mat_melt$Row <- factor(mat_melt$Row, levels = rownames(mat))
    mat_melt$Col <- factor(mat_melt$Col, levels = colnames(mat))

    ## Plot main heatmap
    heatmap <- ggplot2::ggplot(mat_melt, ggplot2::aes(x = Col, y = Row, fill = Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3], na.value = 'lightgray') +
        theme_crate(base_size = base_size) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
            axis.title = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.x = if (show_colnames) ggplot2::element_text(angle = 90, hjust = 1) else ggplot2::element_blank(),
            axis.text.y = if (show_rownames & any(is.na(row_annotations))) ggplot2::element_text(size = 6, hjust = 0) else ggplot2::element_blank(),
            plot.tag = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 1)
        ) +
        ggplot2::coord_cartesian(clip = 'off')
    ## Plot dendrograms
    col_dendro <- ggplot2::ggplot(ggdendro::segment(col_hclust)) +
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggplot2::scale_x_continuous(expand = c(0, 0.5)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        theme_crate(base_size = base_size) +
        ggdendro::theme_dendro() +
        ggplot2::theme(
            plot.tag = ggplot2::element_blank(),
            panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
            plot.margin = ggplot2::margin(t = 0, r = 0, b = 1, l = 0),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_cartesian(clip = 'off')
    row_dendro <- ggplot2::ggplot(ggdendro::segment(row_hclust)) +
        ggplot2::geom_segment(ggplot2::aes(x = y, y = x, xend = yend, yend = xend)) +
        ggplot2::scale_y_continuous(expand = c(0, 0.5)) +
        ggplot2::scale_x_reverse(expand = c(0, 0)) +
        theme_crate(base_size = base_size) +
        ggdendro::theme_dendro() +
        ggplot2::theme(
            plot.tag = ggplot2::element_blank(),
            panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
            plot.margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_cartesian(clip = 'off')

    ## Plot annotations
    if (!any(is.na(col_annotations))) {
        col_annotation_titles <- patchwork::wrap_plots(
            lapply(colnames(col_annotations), function(feature) {
                ggplot2::ggplot() +
                    ggplot2::geom_point(x=1, y=1) +
                    theme_crate(base_size = base_size) +
                    ggplot2::scale_y_continuous(
                        breaks = c(1), labels = feature, expand = c(0, 0)
                    ) +
                    ggplot2::theme(
                        axis.title.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.line.x = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_text(
                            face='italic', angle = 0, vjust = 0.5, hjust = 1
                        ),
                        panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.tag = ggplot2::element_blank()
                    ) +
                    ggplot2::coord_cartesian(clip = 'off')
            }),
            ncol = 1
        )
        col_annotation_tiles <- patchwork::wrap_plots(
            lapply(colnames(col_annotations), function(feature) {
                df <- col_annotations[rownames(col_annotations) %in% colnames(mat), , drop = FALSE]
                df[["sample"]] <- factor(rownames(df), levels = colnames(mat))
                ggplot2::ggplot(df, ggplot2::aes(x = sample, y = 1, fill = .data[[feature]])) +
                    ggplot2::geom_tile() +
                    ggplot2::scale_fill_manual(values = annotation_colors[[feature]]) +
                    ggplot2::scale_x_discrete(expand = c(0, 0)) +
                    ggplot2::scale_y_continuous(expand = c(0, 0)) +
                    theme_crate(base_size = base_size) +
                    ggplot2::theme(
                        axis.title.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.line.x = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.tag = ggplot2::element_blank()
                    ) +
                    ggplot2::coord_cartesian(clip = 'off')
            }),
            ncol = 1
        )
    }
    if (!any(is.na(row_annotations))) {
        rd <- as.data.frame(SummarizedExperiment::rowData(dds[rownames(mat), ]))
        gene_labels <- rd[[gene_label]]
        names(gene_labels) <- rownames(mat)
        row_annotation_tiles <- patchwork::wrap_plots(
            lapply(colnames(row_annotations), function(feature) {
                df <- row_annotations[rownames(row_annotations) %in% rownames(mat), , drop = FALSE]
                df[["gene"]] <- factor(rownames(df), levels = rownames(mat))
                ggplot2::ggplot(df, ggplot2::aes(x = 1, y = gene, fill = .data[[feature]])) +
                    ggplot2::geom_tile() +
                    ggplot2::scale_fill_manual(values = annotation_colors[[feature]]) +
                    ggplot2::scale_y_discrete(
                        expand = c(0, 0),
                        labels = gene_labels
                    ) +
                    ggplot2::scale_x_continuous(
                        expand = c(0, 0), breaks = c(1), labels = feature
                    ) +
                    theme_crate(base_size = base_size) +
                    ggplot2::theme(
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.title.x = ggplot2::element_blank(),
                        axis.line.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank(),
                        axis.text.y = if (show_rownames) ggplot2::element_text(size = 6, hjust = 0) else ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_text(face = 'italic', angle = 60, vjust = 1.05, hjust = 1.1),
                        panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.tag = ggplot2::element_blank()
                    ) +
                    ggplot2::coord_cartesian(clip = 'off')
            }),
            nrow = 1
        )
    }

    tag = ggplot2::ggplot() +
        theme_crate(base_size = base_size) +
        ggplot2::theme(
            panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
            plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
        )
    spacer <- tag + ggplot2::theme(plot.tag = ggplot2::element_blank())
    if (any(is.na(col_annotations)) & any(is.na(row_annotations))) {
        panels  <- list(tag, col_dendro, row_dendro, heatmap)
        design  <- "AB\nCD"
        widths  <- c(row_dendro_scale, 1)
        heights <- c(col_dendro_scale, 1)
    } else if (any(is.na(col_annotations))) {
        panels  <- list(tag, col_dendro, row_dendro, row_annotation_tiles, heatmap)
        design  <- "AAB\nCDE"
        widths  <- c(row_dendro_scale * ncol(mat), row_annotation_scale * ncol(mat) * ncol(row_annotations), ncol(mat))
        heights <- c(col_dendro_scale, 1)
    } else if (any(is.na(row_annotations))) {
        panels  <- list(tag, col_dendro, col_annotation_titles, col_annotation_tiles, row_dendro, heatmap)
        design  <- "AB\nCD\nEF"
        widths  <- c(row_dendro_scale, 1)
        heights <- c(col_dendro_scale * nrow(mat), col_annotation_scale * nrow(mat) * ncol(col_annotations), nrow(mat))
    } else {
        panels  <- list(tag, col_dendro, spacer, col_annotation_titles, col_annotation_tiles, row_dendro, row_annotation_tiles, heatmap)
        design  <- "AAB\nCDE\nFGH"
        widths  <- c(row_dendro_scale * ncol(mat), row_annotation_scale * ncol(mat) * ncol(row_annotations), ncol(mat))
        heights <- c(col_dendro_scale * nrow(mat), col_annotation_scale * nrow(mat) * ncol(col_annotations), nrow(mat))
    }

    if (!assemble) {
        return(list(plots = panels, design = design, widths = widths, heights = heights))
    }

    plot <- patchwork::wrap_plots(panels) +
        patchwork::plot_layout(
            design = design, widths = widths, heights = heights,
            guides = 'collect'
        )
    if (!is.na(filename)) {
        ggplot2::ggsave(filename = filename, plot = plot)
    }
    return(plot)
}
