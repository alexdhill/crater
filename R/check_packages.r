check_packages <- function(libs) {
    return(unlist(lapply(libs, function(lib) {lib %in% rownames(installed.packages())})))
}