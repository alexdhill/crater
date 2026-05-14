check_packages <- function(libs) {
    unlist(lapply(libs, function(lib) {
        lib %in% rownames(installed.packages())
    }))
}

