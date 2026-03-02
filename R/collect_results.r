#'
#' Collect the results from a DESeq2 run and add in gene biotype information
#'
#' @param dds A DESeqDataSet object from DESeqDataSetFromCreate
#' @param ... Additional arguments to pass to DESeq2::results
#' 
#' @return A data.frame of results with gene biotype information
#' @export
collect_results <- function(dds, ...) {
    res <- DESeq2::results(dds, ...)
    res <- as.data.frame(res)
    res[["gene_id"]] <- rownames(res)
    res <- dplyr::left_join(res, as.data.frame(SummarizedExperiment::rowData(dds)), by = "gene_id")
    return(res)
}
