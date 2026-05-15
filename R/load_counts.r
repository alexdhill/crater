#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#'
#' @param counts A H5SummarizedExperiment object created by the CREATE pipeline
#'
#' @return A H5-backed SummarizedExperiment object
#'
#' @export
load_counts <- function(counts) {
    se <- HDF5Array::loadHDF5SummarizedExperiment(counts)
    se
}
