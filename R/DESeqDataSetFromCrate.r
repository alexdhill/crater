#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#'
#' @param create_se A H5SummarizedExperiment object created by the CRATE pipeline
#' @param design A formula specifying the model design
#'
#' @return A DESeqDataSet object
#' @export
DESeqDataSetFromCrate <- function(
    create_se,
    design,
    min_expr = 3,
    min_frac = 0.25,
    min_lib_size = 100000
) {
    dds <- create_se

    expr_filter <- SummarizedExperiment::assays(dds)
    expr_filter <- expr_filter$count
    expr_filter <- expr_filter > min_expr
    expr_filter <- matrixStats::rowSums2(expr_filter)
    expr_filter <- expr_filter > (ncol(dds) * min_frac)

    lib_filter <- SummarizedExperiment::assays(dds)
    lib_filter <- lib_filter$counts
    lib_filter <- matrixStats::colSums2(lib_filter)
    lib_filter <- lib_filter > min_lib_size

    dds <- dds[expr_filter, lib_filter]

    SummarizedExperiment::assays(dds) <- lapply(
        SummarizedExperiment::assays(dds),
        as.matrix
    )
    SummarizedExperiment::assays(dds) <- lapply(
        SummarizedExperiment::assays(dds),
        floor
    )
    dds <- DESeq2::DESeqDataSet(dds, design = design)
    dds
}
