#'
#' Plot geneset enrichment from CRATE dds
#' 
#' @param dds A DESeqDataSet object from DESeqDataSetFromCrate
#' @param geneset A list of gene sets for enrichment analysis
#' @param top The number of top pathways to display
#' @param hide_insig Whether to hide insignificant pathways
#' @param signif_level The significance level for filtering pathways
#' @param colors A vector of colors for the plot
#' @param labeller A function for labeling the pathways
#' 
#' @return A ggplot2 object
#' @export
plot_enrichment <- function(dds, geneset, top=10, hide_insig = TRUE, signif_level = 0.01, colors = c("#50FA7B", "#44475A"), labeller = NA, rm_dupes = FALSE, size = 10, ...) {
    gene_order <- DESeq2::results(dds, ...)
    gene_order <- as.data.frame(gene_order)
    gene_order[["gene_id"]] <- rownames(gene_order)
    gene_order <- dplyr::filter(gene_order, grepl("^ENSG", gene_id))
    gene_names[is.na(gene_names)] = gene_order$gene_id[is.na(gene_names)]
    gene_order <- dplyr::left_join(gene_order, as.data.frame(SummarizedExperiment::rowData(dds)), by = "gene_id")
    gene_order <- dplyr::arrange(gene_order, desc(stat))
    if (rm_dupes) { 
        gene_order <- gene_order[!duplicated(gene_order$gene_name), ]
        gene_names <- gene_order$gene_name
    }
    else { gene_names = make.unique(gene_order$gene_name) }
    gene_order = gene_order$stat
    names(gene_order) = gene_names
    gene_order = na.omit(gene_order)

    enrichment = fgsea::fgsea(
        pathways = geneset,
        stats = gene_order
    )

    if (hide_insig) enrichment = dplyr::filter(enrichment$padj < signif_level)
    enrichment = dplyr::mutate(signif = ifelse(padj < signif_level, TRUE, FALSE))
    enrichment = dplyr::mutate(direction = ifelse(NES > 0, "UP", "DOWN"))
    enrichment = dplyr::group_by(enrichment, direction)
    enrichment = dplyr::slice_max(enrichment, order_by = abs(NES), n = top, with_ties = FALSE)

    plot = ggplot2::ggplot(enrichment, ggplot2::aes(x = reorder(pathway, NES), y = NES, color = signif)) +
        ggplot2::geom_segment(ggplot2::aes(yend = 0)) +
        ggplot2::geom_point(ggplot2::aes(size = -log2(padj))) +
        ggplot2::labs(x = "Pathway", y = "Normalized Enrichment Score (NES)") +
        ggplot2::scale_color_manual(values = c("TRUE" = colors[1], "FALSE" = colors[2])) +
        theme_crate(base_size=size) +
        ggplot2::theme(legend.position = 'none') +
        ggplot2::coord_flip()
    
    if (!is.na(labeller)) plot = plot + ggplot2::scale_x_discrete(labels = labeller)

    return(plot)
}
