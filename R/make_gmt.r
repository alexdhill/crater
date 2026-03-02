#'
#' Load in gene lists as a 'GMT' object
#' 
#' @param genes A list of gene sets (a list of list of character vectors)
#' @param names The name of each gene set
#' 
#' @return A GMT object
#' @export
make_gmt <- function(genes, names, save_to=NA) {
    if (length(genes) != length(names)) {
        stop("Length of genes and names must be the same")
    }

    outfile = ifelse(is.na(save_to), file.path("tmp" + as.numeric(Sys.time()) + ".gmt"), save_to)
    
    dat = data.frame(Name=names, Description=rep("", length(names)), Genes=I(list(genes)))
    genesets = apply(dat, 1, function(row) {
        paste(row["Name"], row["Description"], row["Genes"], sep="\t")
    })
    writeLines(genesets, outfile)

    res = fgsea::gmtPathways(outfile)
    if (is.na(save_to)) file.remove(outfile)

    return(res)
}
