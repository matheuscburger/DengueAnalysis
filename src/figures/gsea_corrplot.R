library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("forcats")
library("scales")
library("stringr")

fgsea_dir <- file.path("results", "fgsea", "Log2FC")

my_squish <- function(...){
    return(scales::squish(..., only.finite=FALSE))
}

# FOR BTMs
for(db in list.files(fgsea_dir)) {
    db_dir <- file.path(fgsea_dir, db)
    gsea_files <- list.files(db_dir, full.names=TRUE)
    study <- gsub(".tsv", "", basename(gsea_files))

    gsea <- lapply(gsea_files, read_tsv) %>% 
        lapply(., function(x) dplyr::select(x, pathway, ends_with("_padj"), ends_with("_NES"))) %>%
        lapply(., function(x) gather(x, "Name", "Value", -pathway) ) %>%
            lapply(., function(x) separate(x, Name, c("Comparison", "Stat"), "_"))
    names(gsea) <- study
    gsea  <- lapply(study, function(x) mutate(gsea[[x]], Study=x))
    gsea_tib <- do.call(rbind, gsea) %>%
        spread(key=Stat, value=Value) %>%
        filter(padj <= 0.1) %>%
        mutate(pathway=str_wrap(pathway, width=40))

    custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                    "#D6604D", "#B2182B", "#67001F")
    custom_pal <- colorRampPalette(custom_pal)(200)

    output_dir <- file.path("figures", "GSEA", "Log2FC", db)
    comparisons <- unique(gsea_tib[["Comparison"]])
    for(comp in comparisons) {
        nes <- gsea_tib %>%
            filter(Comparison == comp) %>%
            group_by(pathway, Comparison) %>%
            mutate(sum_nes = sum(NES),
                   n=n()) %>%
            ungroup() %>%
            arrange(sum_nes) %>%
            mutate(pathway=fct_inorder(pathway)) %>%
            filter(n >= max(n))

        res <- ggplot(nes, aes(x=Study, y=pathway, size=abs(NES), fill=NES)) + 
            geom_point(color = "white", shape=21) +
                scale_fill_gradientn(colours=custom_pal, space = "Lab", 
                                     limits=c(-2, 2),
                                     oob=my_squish) +
            scale_size(range=c(0,12)) +
            guides(size="none") +
            theme_minimal() +
            theme(panel.grid.major = element_blank()) +
            scale_x_discrete(position = "top")
        h <- min(max(length(nes[["pathway"]])*.15, 7), 49)
        message(paste("comp ", comp))
        message(paste("h ", h))
        message(paste("l ", length(nes[["pathway"]])))
        ggsave(res, file=file.path(output_dir, paste0(comp, ".pdf")), height=h)

    }
}
