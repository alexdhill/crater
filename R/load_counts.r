#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#'
#' @param counts A H5SummarizedExperiment directory created by the CREATE pipeline
#'
#' @return A H5-backed SummarizedExperiment object
#'
#' @export
load_counts <- function(counts, summarize_biotypes = TRUE) {
    se <- HDF5Array::loadHDF5SummarizedExperiment(counts)
    if (summarize_biotypes) {
        row_dat <- SummarizedExperiment::rowData(se)
        row_dat <- as.data.frame(row_dat)
        row_dat <- summarize_biotypes(row_dat)
        row_dat <- S4Vectors::DataFrame(row_dat)
        SummarizedExperiment::rowData(se) <- row_dat
    }
    se
}
