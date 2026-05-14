#'
#' Plot geneset enrichment from CRATE dds
#'
#' @param dds A DESeqDataSet object from **after** running DESeq()
#' @param geneset A list of gene sets for enrichment analysis
#' @param top The number of top pathways to display
#' @param hide_insig Whether to hide insignificant pathways
#' @param signif_level The significance level for filtering pathways
#' @param colors A vector of colors for the plot
#' @param labeller A function for labeling the pathways
#'
#' @return A ggplot2 object
#'
#' @export
plot_enrichment <- function(
    dds,
    gene_list,
    ontology = "BP",
    contrast = NA,
    pval = 0.01,
    fold_change = 1,
    signif = 0.01,
    top = 10,
    hide_insig = TRUE,
    colors = c("#50FA7B", "#44475A"),
    labeller = NA,
    ...
) {
    libraries <- c("dplyr", "DESeq2", "topGO", "org.Hs.eg.db")
    has_libs <- check_packages(libraries)
    if (!all(has_libs)) {
        stop(paste0(
            "Missing required packages:\n",
            paste(libraries[!has_libs], collapse = "\n")
        ))
    }

    if (is.na(contrast)) {
        de_genes <- DESeq2::results(dds)
    } else {
        de_genes <- DESeq2::results(dds, contrast = contrast)
    }
    de_genes <- as.data.frame(de_genes)
    de_genes <- dplyr::filter(de_genes, !is.na(padj))
    de_genes <- dplyr::mutate(de_genes, gene_id = rownames(.))
    de_genes <- dplyr::mutate(
        de_genes,
        signif = dplyr::case_when(
            gene_id %in% gene_list ~ TRUE,
            TRUE ~ FALSE
        )
    )
    all_genes <- de_genes$signif
    names(all_genes) <- de_genes$gene_id

    go_data <- methods::new(
        "topGOdata",
        ontology = ontology,
        allGenes = all_genes,
        geneSel = ~ . == TRUE,
        annot = topGO::annFUN.org,
        mapping = "org.Hs.eg.db",
        ID = "ensembl",
        ... = ...
    )

    go_fisher_classic <- topGO::runTest(
        go_data,
        algorithm = "classic",
        statistic = "fisher"
    )
    go_fisher_elim <- topGO::runTest(
        go_data,
        algorithm = "elim",
        statistic = "fisher"
    )
    go_kruskal_classic <- topGO::runTest(
        go_data,
        algorithm = "classic",
        statistic = "KruskalWallis"
    )
    go_kruskal_elim <- topGO::runTest(
        go_data,
        algorithm = "elim",
        statistic = "KruskalWallis"
    )

    go_results <- topGO::GenTable(
        go_data,
        classicFisher = go_fisher_classic,
        elimFisher = go_fisher_elim,
        classicKruskal = go_kruskal_classic,
        elimKruskal = go_kruskal_elim,
        orderBy = "classicFisher",
        ranksOf = "classicFisher",
        topNodes = length(topGO::usedGO(go_data))
    )

    return(go_results)
}
