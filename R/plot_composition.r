#'
#' Make a plot of the composition of each sample
#' 
#' @param dds A DESeqDataSet object from DESeqDataSetFromCreate
#' @param position Position adjustment for geom_bar [default="fill"]
#' @param fill Variable to fill bars with [default=gene_biotype]
#' @param facet Variables to facet the plot with [default=NA]
#' @param subset Samples to include in the plot [default=NA]
#' 
#' @return A ggplot2 object
#'
#' @export
plot_composition <- function(dds, position="fill", fill="gene_biotype", facet=NA, subset=NA, nested=FALSE, add_counts=FALSE, base_size = 10) {
    counts = DESeq2::counts(dds, normalized=TRUE)
    counts = as.data.frame(counts)

    if (!is.na(subset)) {
        counts = dplyr::select(counts, dplyr::all_of(subset))
    }

    biotypes = SummarizedExperiment::rowData(dds)
    biotypes = as.data.frame(biotypes)
    samples = colnames(counts)
    counts = dplyr::mutate(counts, gene_id = rownames(counts))
    counts = dplyr::left_join(counts, biotypes, by="gene_id")
    counts = tidyr::pivot_longer(
        counts, dplyr::all_of(samples),
        names_to="names", values_to="count"
    )

    groups = c(facet, fill, "names")
    groups = groups[!is.na(groups)]
    counts = dplyr::summarise(counts, count = sum(count), .by = dplyr::all_of(groups))

    counts = dplyr::left_join(counts, as.data.frame(SummarizedExperiment::colData(dds)), by = "names")
    counts = dplyr::mutate(counts, dplyr::across(dplyr::all_of(fill), function(bt){factor(bt, levels = names(biotype_colors))}))

    plot = ggplot2::ggplot(counts, ggplot2::aes(x = names, y = count)) +
        ggplot2::geom_bar(ggplot2::aes(fill = .data[[fill]]), stat = "identity", position = position)

    if (!is.na(facet) && length(facet) > 0) {
        if (length(facet) == 2 && 'list' %in% class(facet[[1]])) {
            if (nested) {
                plot = plot + ggh4x::facet_nested(rows = facet[[1]], cols = facet[[2]], scales = "free", space= "free")
            } else {
                plot = plot + ggplot2::facet_grid(rows = facet[[1]], cols = facet[[2]], scales = "free", space= "free")
            }
        } else {
            if (nested) {
                plot = plot + ggh4x::facet_nested(~facet, scales = "free", space= "free")
            } else {
                plot = plot + ggplot2::facet_grid(~facet, scales = "free", space= "free")
            }
        }
    }

    plot = plot +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_fill_manual(values = biotype_colors) +
        ggplot2::labs(x = "Sample", y = ifelse(position == "fill", "Proportion of counts", "Counts"), fill = "Biotypes") +
        theme_crate(base_size = base_size)

    if (add_counts) {
        sample_counts = dplyr::group_by(counts, names)
        sample_counts = dplyr::summarize(sample_counts, total = sum(count))
        plot = plot + ggplot2::geom_text(
            data = sample_counts,
            ggplot2::aes(x = names, y = Inf, label = scales::comma(round(total))),
            inherit.aes = FALSE,
            vjust = 1.2, size = base_size / ggplot2::.pt
        )
    }

    return(plot)
}
