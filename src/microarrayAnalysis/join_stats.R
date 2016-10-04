#!/usr/bin/env Rscript

"Join Statistics

Usage: join_stats.R OUTPUT [--comp-col=<value>] [--probe-col=<value>] [--gene-col=<value>] [--lfc-col=<value>] [--pv-col=<value>] [--padj-col=<value>] [--cores=<value>] [--key-col=<value>...] [--maintain-col=<value>] (--input=<file>...)


Options:
  OUTPUT                     output file
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             The input file containing the statistics values (Log2FC, p-value, ...)
  --comp-col=<value>         Column containing the comparison [default: title] 
  --probe-col=<value>        Column containing the probe name
  --gene-col=<value>         Column containing the gene symbol
  --lfc-col=<value>          Column containing Log2(fold-change) [default: Log2FC]
  --pv-col=<value>           Column containing raw p-values [default: LIMMA.rawp]
  --padj-col=<value>         Column containing adjusted p-values [default: LIMMA.adjp]
  --cores=<value>            Number of cores/threads
  --key-col=<value>          Column to use as key on full join operation
  --maintain-col=<value>     Column to maintain after full join operation

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc


options(stringsAsFactors=FALSE)


mc.apply <- function(in.df, fun, n.cores, ...){
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
                }, mc.cores=n.cores, ...)
    names(res.list) <- NULL
    res.df <- do.call(rbind.data.frame, res.list)
    return(res.df)
}

# parameters



# Join tables with statistics
join.stats <- function(filenames, maintain.cols, key.cols, probe.col){
    # extract basename without extension
    bnames <- str_match(basename(filenames), "(\\S+)\\.\\S+")[,2]
    names(bnames) <- filenames	

    res <- data.table()
	for(k in key.cols){
		res[[k]] <- character(0)
	}

    for( fn in filenames ){
		# read file
        fa <- fread(fn)
        message("First 6 rows in ", fn, " ...")
        write.table(head(fa), file=stderr(), sep="\t", quote=FALSE, row.names=F)
		# remove _PM do probename
		if(!missing(probe.col)){
			fa[[probe.col]] <- sub("_PM", "", fa[[probe.col]])
		}
		# keep only columns that will be used
        fa <- fa[, union(maintain.cols, key.cols), with=FALSE]
		# get the positions of columns that should be maintained
        pos <- which(colnames(fa) %in% setdiff(maintain.cols, key.cols))
		# include filename on column name
        colnames(fa)[pos] <- paste(colnames(fa)[pos], bnames[fn], sep=".")
		# merge tables using key columns
        res <- merge(res, fa, by=union(comp.col, key.cols), all=T)
    }
    return(res)
}

# for a numeric vector with p-values return combined p-values
get_metap <- function(x, methods=c("logitp", "meanp", "minimump", "sumlog", "sump", "sumz", "votep"), calc.min=T, calc.max=T){
    res <- numeric(length=length(methods)+sum(calc.min, calc.max))
    names.res <- paste0(methods, ".metap")
    if(calc.min){
        names.res <- c(names.res, "minp")
    }
    if(calc.max){
        names.res <- c(names.res, "maxp")
    }
    names(res) <- names.res
    for(m.name in methods){
        m <- get(m.name)
        res[paste0(m.name, ".metap")] <- tryCatch(m(x[which(!is.na(x))])$p, error=function(e) NA)
    }
    if(calc.min){
        res["minp"] <- min(x, na.rm=T)
    }
    if(calc.max){
        res["maxp"] <- max(x, na.rm=T)
    }
    return(res)
}

combine.stats <- function(in.df, pv.col, padj.col, lfc.col, n.cores){
    ## comparison between parallel and serial
    #ps <- res %>% select(starts_with(pv.col))
    #system.time(mps <- t(apply(ps, 1, get_metap)))
    #system.time(mc.mps <- mc.apply(ps, get_metap, 10))
	res <- in.df
	if(!missing(pv.col)){
		ps <- res %>% select(starts_with(pv.col)) %>%
			mc.apply(get_metap, n.cores)
		colnames(ps) <- paste0("rawp.", colnames(ps))
		res <- cbind(res, ps)
	}
	if(!missing(padj.col)){
		adjps <- res %>% select(starts_with(padj.col)) %>%
			mc.apply(get_metap, n.cores)
		colnames(adjps) <- paste0("adjp.", colnames(adjps))
		res <- cbind(res, adjps)
	}
	if(!missing(lfc.col)){
		lfc <- res %>% select(starts_with(lfc.col))
		res$mean.lfc <- mc.apply(lfc, function(x) mean(x, na.rm=T), n.cores=n.cores)
		res$min.lfc <- mc.apply(lfc, function(x) x[which.min(abs(x))], n.cores=n.cores)
	}
    return(res)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	comp.col <- arg[["comp_col"]]
	probe.col <- arg[["probe_col"]]
	lfc.col <- arg[["lfc_col"]]
	pv.col <- arg[["pv_col"]]
	padj.col <- arg[["padj_col"]]
	filenames <- arg[["input"]]
	n.cores <- ifelse("cores" %in% names(arg), as.numeric(arg[["cores"]]), getOption("mc.cores", 2L))
	maintain.cols <- c(lfc.col, pv.col, padj.col)
	if("maintain_col" %in% names(arg)){
		maintain.cols <- c(maintain.cols, arg[["maintain_col"]])
	}
	key.cols <- c(comp.col, probe.col)
	if("probe_col" %in% names(arg)){
		key.cols <- union(key.cols, arg[["probe_col"]])
	}
	if("key_col" %in% names(arg)){
		key.cols <- union(key.cols, arg[["key_col"]])
	}
	if("gene_col" %in% names(arg)){
		key.cols <- union(key.cols, arg[["gene_col"]])
	}

	message("Using as key following columns: ", paste(key.cols, collapse=", "))
	message("Will maintain following columns: ", paste(maintain.cols, collapse=", "))

	## RUN
	suppressMessages(library("stringr"))
	suppressMessages(library("readr"))
	suppressMessages(library("tidyr"))
	suppressMessages(library("dplyr"))
	suppressMessages(library("data.table"))
	suppressMessages(library("metap"))
	suppressMessages(library("parallel"))

	joined <- join.stats(filenames=filenames, maintain.cols=maintain.cols, 
						 key.cols=key.cols, probe.col=probe.col)
	res <- combine.stats(joined, pv.col=pv.col, padj.col=padj.col, 
						 lfc.col=lfc.col, n.cores=n.cores)
	write.table(res, arg[["output"]], quote=FALSE, sep="\t", row.names=FALSE)
}
