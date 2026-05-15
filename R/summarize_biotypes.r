#'
#' Collect the GENCODE/RepeatMasker biotypes into relevant subsets
#'
#' @param data A data frame containing a column with biotype information
#' @param biotypes A string specifying the column name with biotype information [default="gene_biotype"]
#'
#' @return The data frame with the summarized biotypes
#'
#' @export
summarize_biotypes <- function(
    data,
    biotypes = "gene_biotype",
    gene_ids = 'gene_id',
    gene_names = "gene_name",
    biotypes_to = NA
) {
    # Gather columns
    summarized_biotypes <- data.frame(bt = data[[biotypes]])
    if (gene_ids %in% colnames(data)) {
        summarized_biotypes$gi <- data[[gene_ids]]
    } else {
        stop("Need gene_ids to summarize biotypes")
    }
    if (gene_names %in% colnames(data)) {
        summarized_biotypes$gn <- data[[gene_names]]
    } else {
        stop("Need gene_names to summarize biotypes")
    }

    summarized_biotypes <- dplyr::mutate(
        summarized_biotypes,
        biotype_class = dplyr::case_when(
            stringr::str_starts(gi, "ENSG") ~ "gene",
            TRUE ~ "repeat"
        )
    )

    summarized_biotypes <- dplyr::mutate(
        summarized_biotypes,
        bt <- dplyr::case_when(
            startsWith(bt, "Mt_") | startsWith(gn, "MT-") ~ "Mitochondrial",
            bt == "protein_coding" ~ "Coding",
            bt %in% c("lncRNA", "miRNA", "LINE", "SINE", "LTR", "DNA") ~ bt,
            bt == "Simple_repeat" ~ "Microsatellite",
            bt == "Satellite" ~ "Human satellite",
            biotype_class == "gene" ~ "Other gene",
            TRUE ~ "Other repeat"
        )
    )

    summarized_biotypes <- dplyr::pull(summarized_biotypes, bt)
    if (is.na(biotypes_to)) {
        biotype_col <- biotypes
    } else {
        biotype_col <- biotypes_to
    }
    data[[biotype_col]] <- summarized_biotypes

    data
}
