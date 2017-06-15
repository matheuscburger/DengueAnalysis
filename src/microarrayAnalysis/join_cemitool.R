#!/usr/bin/env Rscript

"Join CEMiTool results (modules file)

Usage: join_cemitool.R (--input=<file>...) --output=<file> [(--names=<v>...) --modules-column=<v> --genes-column=<v> --not-correlated=<v>] 

Options:
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             input file name (at least 2)
  --output=<file>            output file name
  --names=<v>                name of each input file (same order that input files)
  --modules-column=<v>       name of the column containing the modules [default: modules]
  --genes-column=<v>         name of the columns containing the genes [default: genes]
  --not-correlated=<v>       name of the module of genes not correlated [default: Not.Correlated]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc


get_adjacency_tib <- function(s1_tib, study_name="Present", 
                              modules_cname="modules",
                              genes_cname="genes",
                              not_correlated="Not.Correlated"){
    adj_tib <- s1_tib %>%
        dplyr::filter_(paste0(modules_cname, "!=\"", not_correlated, "\"")) %>%
        with(., split(get(genes_cname), get(modules_cname))) %>%
        lapply(., sort) %>%
        lapply(., combn, 2) %>%
        lapply(., t) %>%
        do.call(rbind, .) %>%
        as_data_frame %>%
        mutate(., Present=T) %>%
        setNames(., c("Gene1", "Gene2", study_name))
    return(adj_tib)
}

main <- function(studies_files, studies_names,
                 modules_cname, genes_cname, not_correlated){
    all_tibs <- lapply(studies_files, read_tsv) %>%
        setNames(., studies_names)
    all_adjs <- lapply(studies_names, function(x) 
                       get_adjacency_tib(all_tibs[[x]], x,
                                         modules_cname,
                                         genes_cname,
                                         not_correlated)) %>%
        Reduce(function(t1, t2) 
               full_join(t1, t2, by=c("Gene1", "Gene2")), .) %>%
        mutate(Sum=rowSums(.[studies_names], na.rm=T)) %>%
        arrange(desc(Sum))
    return(all_adjs)
}

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

    suppressMessages(library("readr"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("tidyr"))

    modules_cname <- parameters[["modules_column"]]
    genes_cname <- parameters[["genes_column"]]
    not_correlated <- parameters[["genes_column"]]
    out_file <- parameters[["output"]]
    studies_files <- parameters[["input"]]
    if(!"names" %in% names(parameters)) {
        studies_names <- paste0("S", 1:length(studies_files))
    } else {
        studies_names <- parameters[["names"]]
    }
    all_adjs <- main(studies_files, studies_names,
         modules_cname, genes_cname, not_correlated)

    write_tsv(all_adjs, out_file)
}
