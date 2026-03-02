plot_heatmap <- function(
    dds, genes, base_size = 8,
    samples = NA, gene_label = 'gene_id',
    show_rownames = TRUE, show_colnames = FALSE,
    colors = c('red', 'white', 'blue'), scale = 'logz',
    col_annotations = NA, row_annotations = NA, annotation_colors = NA,
    row_annotation_scale = 0.0625, col_annotation_scale = 0.0625,
    row_dendro_scale = 0.125, col_dendro_scale = 0.125,
    labeller = NA, filename = NA
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
    if (!is.na(samples)) {
        if (class(samples) != 'logical') {
            samples <- colnames(mat) %in% samples
        }
        mat <- dplyr::select(mat, dplyr::all_of(samples))
    }
    if (scale == 'log') {
        mat <- magrittr::add(mat, 1)
        mat <- log10(mat)
    } else if (scale == 'logz') {
        mat <- magrittr::add(mat, 1)
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
    col_hclust <- hclust(dist(t(mat)))
    col_hclust <- treeorder(col_hclust, matrixStats::colMedians(mat))
    col_order <- col_hclust$order
    col_hclust <- as.dendrogram(col_hclust)
    col_hclust <- ggdendro::dendro_data(col_hclust)

    ## Make row dendrograms
    row_hclust <- hclust(dist(mat))
    row_hclust <- treeorder(row_hclust, matrixStats::rowMedians(mat))
    row_order <- row_hclust$order
    row_hclust <- as.dendrogram(row_hclust)
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
        creater::theme_create(base_size = base_size) +
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
        creater::theme_create(base_size = base_size) +
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
        creater::theme_create(base_size = base_size) +
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
        col_annotation_titles <- lapply(colnames(col_annotations), function(feature) {
            ggplot2::ggplot() +
                ggplot2::geom_point(x=1, y=1) +
                creater::theme_create(base_size = base_size) +
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
        }) %>%
        patchwork::wrap_plots(ncol = 1)
        col_annotation_tiles <- lapply(colnames(col_annotations), function(feature) {
            col_annotations %>%
                dplyr::filter(rownames(.) %in% colnames(mat)) %>%
                dplyr::mutate(sample = factor(rownames(.), levels = colnames(mat))) %>%
                dplyr::select(sample, feature) %>%
                ggplot2::ggplot(ggplot2::aes(x = sample, y = 1, fill = .data[[feature]])) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_manual(values = annotation_colors[[feature]]) +
                ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::scale_y_continuous(expand = c(0, 0)) +
                theme_create(base_size = base_size) +
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
        }) %>%
        patchwork::wrap_plots(ncol = 1)
    }
    if (!any(is.na(row_annotations))) {
        row_annotation_tiles <- lapply(colnames(row_annotations), function(feature) {
            row_annotations %>%
                dplyr::filter(rownames(.) %in% rownames(mat)) %>%
                dplyr::mutate(gene = factor(rownames(.), levels = rownames(mat))) %>%
                dplyr::select(gene, feature) %>%
                ggplot2::ggplot(ggplot2::aes(x = 1, y = gene, fill = .data[[feature]])) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_manual(values = annotation_colors[[feature]]) +
                ggplot2::scale_y_discrete(
                    expand = c(0, 0),
                    labels = dds[rownames(mat),] %>%
                        SummarizedExperiment::rowData() %>%
                        as.data.frame() %>%
                        {
                            res <- .[[gene_label]]
                            names(res) <- rownames(mat)
                            res
                        } %>%
                        ggplot2::labeller()
                ) +
                ggplot2::scale_x_continuous(
                    expand = c(0, 0), breaks = c(1), labels = feature
                ) +
                theme_create(base_size = base_size) +
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
        }) %>%
        patchwork::wrap_plots(nrow = 1)
    }

    tag = ggplot2::ggplot() +
        theme_create(base_size = base_size) +
        ggplot2::theme(
            panel.spacing = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
            plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
        )
    spacer = tag + ggplot2::theme(plot.tag = ggplot2::element_blank())
    if (any(is.na(col_annotations)) & any(is.na(row_annotations))) {
        plot <- patchwork::wrap_plots(list(tag, col_dendro, row_dendro, heatmap)) +
            patchwork::plot_layout(
                design = "AB\nCD",
                widths = c(row_dendro_scale, 1),
                heights = c(col_dendro_scale, 1),
                guides = 'collect'
            )
    } else if (any(is.na(col_annotations))) {
        plot <- patchwork::wrap_plots(list(tag, col_dendro, row_dendro, row_annotation_tiles, heatmap)) +
            patchwork::plot_layout(
                design = "AAB\nCDE",
                widths = c(row_dendro_scale*ncol(mat), row_annotation_scale*ncol(mat)*ncol(row_annotations), ncol(mat)),
                heights = c(col_dendro_scale, 1),
                guides = 'collect'
            )
    } else if (any(is.na(row_annotations))) {
        plot <- patchwork::wrap_plots(list(
                tag, col_dendro,
                col_annotation_titles, col_annotation_tiles,
                row_dendro, heatmap
            )) +
            patchwork::plot_layout(
                design = "AB\nCD\nEF",
                widths = c(row_dendro_scale, 1),
                heights = c(col_dendro_scale*nrow(mat), col_annotation_scale*nrow(mat)*ncol(col_annotations), nrow(mat)),
                guides = 'collect'
            )
    } else {
        plot <- patchwork::wrap_plots(list(
                tag, col_dendro,
                spacer, col_annotation_titles, col_annotation_tiles,
                row_dendro, row_annotation_tiles, heatmap
            )) +
            patchwork::plot_layout(
                design = "AAB\nCDE\nFGH",
                widths = c(row_dendro_scale*ncol(mat), row_annotation_scale*ncol(mat)*ncol(row_annotations), ncol(mat)),
                heights = c(col_dendro_scale*nrow(mat), col_annotation_scale*nrow(mat)*ncol(col_annotations), nrow(mat)),
                guides = 'collect'
            )
    }
    if (!is.na(filename)) {
        ggplot2::ggsave(
            filename = filename,
            plot = plot
        )
    }
    return(plot)
}