#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#'
#' @param create_se A H5SummarizedExperiment object created by the CRATE pipeline
#' @param design A formula specifying the model design
#'
#' @return A DESeqDataSet object
#' @export
DESeqDataSetFromCrate <- function(create_se, design) {
    dds <- create_se
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
