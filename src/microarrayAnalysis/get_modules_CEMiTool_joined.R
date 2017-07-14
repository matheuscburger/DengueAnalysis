#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

"Get modules from CEMiTool joined results

Usage: get_modules_CEMiTool_joined.R --input=<file> --tsv=<file> --graphml=<v> --gmt=<file> --min-studies=<v> --algorithms=<val>... --chosen-method=<val> 


Options:
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             input file in TSV format
  --tsv=<file>               output file in TSV format
  --graphml=<file>           output file in GraphML format
  --gmt=<file>               output file in GMT format
  --min-studies=<v>          minimum number of studies to consider [default: 1]
  --algorithms=<val>         algorithms to use fast_greedy, label_prop, leading_eigen, louvain, walktrap
  --chosen-method=<val>      chosen method to use in gmt file

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

write_gmt <- function(gene_list, filename){
    genes <- lapply(names(gene_list), function(x) c(x, "", gene_list[[x]]))
    genes_pasted <- sapply(genes, paste, collapse="\t")
    writeLines(genes_pasted, filename)
}

if (!interactive() && !exists('SOURCE')) {
    # Get and check arguments.
    suppressMessages(library("docopt"))
    arg <- docopt(doc, version="0.0.1\n", strict=T)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))

    library("igraph")
    library("readr")
    library("dplyr")
    library("parallel")

    input <- arg[["input"]]
    output_tsv <- arg[["tsv"]]
    output_graphml <- arg[["graphml"]]
    min_studies <- as.numeric(arg[["min_studies"]])
    #power_beta <- 1
    mods_to_run <- arg[["algorithms"]]
    chosen_method <- arg[["chosen_method"]]
    gmt <- arg[["gmt"]]

    in_tib <- read_tsv(input) %>%
        dplyr::filter(Sum >= min_studies) %>%
        mutate(powered=Sum**2/3**2)

    g <- graph_from_data_frame(in_tib, directed=FALSE)
    edge_attr(g, "weight") <- E(g)$"powered"

    cluster_methods <- list()
    cluster_methods[["edge_betweeness"]] <- cluster_edge_betweenness  # time
    cluster_methods[["fast_greedy"]] <- cluster_fast_greedy
    cluster_methods[["label_prop"]] <- cluster_label_prop
    cluster_methods[["leading_eigen"]] <- cluster_leading_eigen
    cluster_methods[["louvain"]] <- cluster_louvain
    cluster_methods[["optimal"]] <- cluster_optimal  # memory
    cluster_methods[["spinglass"]] <- cluster_spinglass  # time
    cluster_methods[["walktrap"]] <- cluster_walktrap

    for(mod_name in mods_to_run) {
        mods <- (cluster_methods[[mod_name]])(g)
        vertex_attr(g, mod_name) <- mods$membership
    }


    vertex_df <- do.call(cbind.data.frame, vertex_attr(g))
    write_tsv(vertex_df, output_tsv)
    write_graph(g, output_graphml, "graphml")

    modules <- split(as.character(vertex_df[["name"]]), vertex_df[[chosen_method]])
    modules <- modules[order(sapply(modules, length), decreasing=T)]
    names(modules) <- paste0("Mod", 1:length(modules))
    write_gmt(modules, gmt)

}
