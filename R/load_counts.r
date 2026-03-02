#'
#' Make a DDS object from a CREATE H5SummarizedExperiment
#' 
#' @param create A H5SummarizedExperiment object created by the CREATE pipeline
#' 
#' @return A H5-backed SummarizedExperiment object
#'
#' @export
load_counts <- function(countsDir) {
    se <- HDF5Array::loadHDF5SummarizedExperiment(countsDir)
    return(se)
}
