#'
#' Collect the results from a DESeq2 run and add in gene biotype information
#'
#' @param dds A DESeqDataSet object from DESeqDataSetFromCreate
#' @param cols The columns to include in the results
#' @param ... Additional arguments to pass to DESeq2::results
#'
#' @return A data.frame of results with gene biotype information
#' @export
collect_results <- function(dds, cols = NA, ...) {
    res <- DESeq2::results(dds, ...)
    res <- as.data.frame(res)
    res[["gene_id"]] <- rownames(res)
    if (is.na(cols)) {
        res <- dplyr::left_join(
            res,
            as.data.frame(SummarizedExperiment::rowData(dds)),
            by = "gene_id"
        )
    } else {
        res <- dplyr::left_join(
            res,
            dplyr::select(
                as.data.frame(
                    SummarizedExperiment::rowData(
                        dds
                    )
                ),
                dplyr::all_of(cols)
            ),
            by = "gene_id"
        )
    }
    res
}
