#!/usr/bin/Rscript

"Integrate Results from CEMiTool

 WARNING: Script uses 4 to 8 cores. 
 
 Script derives a PPI from CEMiTool results taking as inputs 2 or more GenesInModules files
 	generated by CEMiTool. 
 All genes from the same module will be connected within themselfs and integrated between
       files. For each file submitted to the script, a column containing logicals
       will be created. In each column there will be ones and zeros showing whether the 
       interaction described in the row is present in any module of the file.
 
 Output default filename is hardcoded to IntersecMods.tsv
 Any questions, ask Diogenes


Usage: integrateCemitResult.R [--threads=<val>] [--output=<file>] INPUT (INPUT...)

Input:
  INPUT                      the input file containing 

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --threads=<val>            number of threads to use
  --output=<file>            output file name [default: IntersecMods.tsv]

Authors:
  Diogenes Saulo de Lima - s.lima.diogenes at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc


suppressPackageStartupMessages({
  require(data.table)
  require("parallel")
})

# Parallel apply as implemented by Matheus Burger.
mc.apply <- function(in.df, fun, n.cores){
    if(missing(n.cores)){
        n.cores <- getOption("mc.cores", 2L)
    }
    n.cores <- min(n.cores, nrow(in.df))
    split.df <- split(in.df, cut(1:nrow(in.df), n.cores))
    res.list <- mclapply(split.df, function(x) {
                             sub.res <- apply(x, 1, fun)
                             if(is.vector(sub.res)){
                                sub.res <- as.matrix(sub.res)
                             } else {
                                sub.res <- t(sub.res)
                             }
                             return(sub.res)
                }, mc.cores=n.cores)
    names(res.list) <- NULL
    res.df <- do.call(rbind.data.frame, res.list)
    return(res.df)
}

# Function returns all combinations in one data frame
makeintact <- function(x, ncores){
  gn_in <- split(x, x$mods)
  tti <- mclapply(gn_in, function(oline){
    cb <- combn(rownames(oline), 2, simplify = F)
    cb 
  }, mc.cores = ncores)
  ppi_from_mods <- lapply(tti, function(mod) do.call(rbind, mod))
  ppi_from_mods <- as.data.frame(do.call(rbind, ppi_from_mods))
  ppi_from_mods <- mc.apply(ppi_from_mods, sort, n.cores = ncores)
  unique_intacts <- paste(ppi_from_mods[,1], ppi_from_mods[,2], 
			  sep = '//')
  ppi_from_mods <- data.frame(unique_intacts,1)
  ppi_from_mods
}

# Function opens, checks file consistence and calls for the function to
# 	perform interactions 
run_intact <- function(filepath, rname, ncores){
  #   x = valid file path
  #   rname = unique identifier of result
  #   ncores = number of cores to use
  k <- read.table(filepath)		
  if(ncol(k) != 5){stop( 'Inconsistent file dimensions' )}
  colnames(k) <- c('colors', 'WGCNA', 'color2', 'mods', 'membership')
  k <- subset(k, !is.na(mods))
  k <- makeintact(k, ncores = ncores)
  colnames(k) <- c('ppi', rname) 
  pp <- data.table(k)
  pp
}

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	## RUN
	# Start running
	out_filename <- parameters$output
	arg <- parameters$input
	if(!all(file.exists(arg))){ stop('Please submit valid files') }

	avCores <- detectCores()
	if("threads" %in% names(parameters)) {
		cores2use <- as.numeric(parameters$threads)
	} else {
		if(avCores <= 2){
		  cores2use <- 1
		}else if(avCores %in% c(4:16)){
		  cores2use <- avCores/2 
		}else if(avCores > 16){
		  cores2use <- avCores/4 
		}
	}

	results <- lapply(arg, function(filepath){
	  # x <- gsub('^results', '', filepath)
	  x <- gsub('_.*', '', filepath)
	  x <- gsub('.*/', '', x)
	  if(x == 'GenesInModules.txt'){ stop('All files must have unique names') }
	  m <- run_intact(filepath, rname = x, ncores = cores2use)
	  m <- as.data.frame(m)
	  m
	})

	message(paste0('\n\t', cores2use, ' cores are being used.'))
	message(paste0('\tCross-referencing ', length(arg), ' results. This may take a while.'))

	out <- Reduce(function(...) {merge(..., all = T, by = 'ppi')}, results) 
	out <- as.data.frame(out)
	out[is.na(out)] <- 0
	output <- out
	rownames(output) <- output$ppi
	output <- output[,!colnames(output) == 'ppi']

	# Make character vector describing which files row interaction is present
	presentIn <- mc.apply(output, function(oneline){
	  colns <- names(oneline) 
	  oneline <- as.logical(as.numeric(oneline))
	  g <- paste0(colns[oneline], collapse = '_')
	  g
	}, n.cores = cores2use)

	presentIn <- as.character(presentIn[,1])

	summed <- mc.apply(output, sum, n.cores = cores2use)
	summed <- as.numeric(summed[,1])

	G1 <- sapply(strsplit(as.character(out$ppi), '//'), '[', 1) 
	G2 <- sapply(strsplit(as.character(out$ppi), '//'), '[', 2) 
	mods_tti <- data.frame(G1 = G1, G2 = G2, output, presentIn = presentIn, sumOfPairs = summed)
	rownames(mods_tti) <- NULL
	mods_tti <- mods_tti[order(mods_tti$sumOfPairs, decreasing = T), ]

	write.table(mods_tti, out_filename, sep = '\t', quote = F, row.names = F)
	message('\tOutput was written on ', out_filename, ' \n')
}
