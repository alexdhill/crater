#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#' 
#' @import SummarizedExperiment
#' @import DelayedArray
#' @import DESeq2
#'
#' @param create A H5SummarizedExperiment object created by the CREATE pipeline
#' @param design A formula specifying the model design
#' 
#' @return A DESeqDataSet object
#' @export
DESeqDataSetFromCreate <- function(create, design) {
    dds <- create
    SummarizedExperiment::assays(dds) <- lapply(SummarizedExperiment::assays(dds), as.matrix)
    SummarizedExperiment::assays(dds) <- lapply(SummarizedExperiment::assays(dds), floor)
    dds <- DESeq2::DESeqDataSet(dds, design = design)
    return(dds)
}