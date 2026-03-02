#' 
#' Plot Differential expression from DESeq2 results as a volcano plot
#' 
#' @param dds A DESeqDataSet object after running DESeq()
#' @param contrast A character vector of length 3 specifying the contrast to plot, the numerator and the denominator
#' @param color_by A column name in the rowData of the dds to color points by. Default is "gene_biotype"
#' @param colors A named vector of colors to use for the different levels of color_by. Default is biotype_colors (assume biotypes are summarized)
#' @param p_threshold Adjusted p-value threshold for significance. Default is 0.01
#' @param lfc_threshold Log2 fold change threshold for significance. Default is 1
#' @param size Base size for the plot text and points. Default is 8
#' 
#' @return A ggplot2 object representing the volcano plot
#' 
#' @import ggplot2
#' @import dplyr
#' @import SummarizedExperiment
#' @import DESeq2
#' 
#' @export
plot_volcano <- function(dds, contrast, color_by = "gene_biotype", colors = biotype_colors, p_threshold=0.01, lfc_threshold=1, size = 8) {
    message("Gathering DESeq2 results...")
    res <- DESeq2::results(dds, contrast=contrast)
    res <- as.data.frame(res)
    res <- dplyr::filter(res, !is.na(padj))
    res <- dplyr::mutate(res, gene_id = rownames(res), negLogPadj = -log10(padj))
    res <- dplyr::left_join(res, as.data.frame(SummarizedExperiment::rowData(dds)), by="gene_id")
    res <- dplyr::mutate(res, signif = dplyr::case_when(
        (padj < p_threshold) & (abs(log2FoldChange) >= lfc_threshold) ~ "Significant",
        TRUE ~ "Not Significant"
    ))
    res <- dplyr::mutate(res, tmp_color_signif = dplyr::case_when(
        signif == "Significant" ~ res[[color_by]],
        TRUE ~ "Not Significant"
    ))
    res <- dplyr::mutate(res, tmp_color_signif = factor(tmp_color_signif, levels=c(names(colors), "Not Significant")))
    
    missing_levels <- setdiff(levels(res$tmp_color_signif), levels(factor(res$tmp_color_signif)))
    if (length(missing_levels) > 0) {
        missing_rows <- lapply(missing_levels, function(lvl) {
            row <- c(rep(NA, ncol(res)))
            names(row) <- colnames(res)
            row[["tmp_color_signif"]] <- lvl
            return(row)
        })
        missing_rows <- do.call(rbind, missing_rows)
        res <- rbind(res, missing_rows)
    }

    message("Creating volcano plot...")
    volcano_plot <- ggplot2::ggplot(res, ggplot2::aes(x=as.numeric(log2FoldChange), y=as.numeric(negLogPadj))) +
        ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype="dashed", color="grey") +
        ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype="dashed", color="grey") +
        ggplot2::geom_point(ggplot2::aes(color=tmp_color_signif), size=size/10, na.rm = TRUE) +
        ggplot2::xlab("Log2 Fold Change") +
        ggplot2::ylab("-Log10 Adjusted P-value") +
        ggplot2::scale_color_manual(
            values = c("Not Significant" = "lightgrey", colors),
            breaks = names(colors), limits = levels(res$tmp_color_signif)
        ) +
        ggplot2::scale_y_continuous(expand = c(0.01, 0)) +
        ggplot2::guides(color=ggplot2::guide_legend(override.aes = list(shape = 15, size = size/2))) +
        ggplot2::labs(
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value",
            color = stringr::str_to_title(gsub("_", " ", color_by))
        ) +
        theme_create(base_size = size)
    return(volcano_plot)
}
