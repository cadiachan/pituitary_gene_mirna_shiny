# load("data/pituitary_mirna_counts.rda")
# load("data/pituitary_utr_counts.rda")
# 
set.seed(2)
# Libraries
library(dplyr)
library(stringr)
library(edgeR)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(scales)
library(colorspace)
library(EDASeq)
library(biomaRt)
library(gridExtra)
library(grid)
library(cowplot)
library(UpSetR)

# Load in UTR and miRNA counts
utr_obj <- readRDS("data/pit_utr_2019__RUV_k2_set1_2019-07-03.rds")
colnames(utr_obj) <- substr(colnames(utr_obj),8,50)
utr_normcounts <- normCounts(utr_obj)
utr_log2 <- log2(normCounts(utr_obj) + 1)

mirna_obj <- readRDS("data/RUV_corrected_pit_srna-seq_counts_combined.rds")
mirna_normcounts <- normCounts(mirna_obj)
mirna_log2 <- log2(normCounts(mirna_obj) + 1)

# Load in UTR and miRNA DE comparisons
utr_de <- readRDS("data/pit_utr_2019_de_result_list_2019-07-03.rds")
mirna_de <- readRDS("data/mirna_pit_all_edgeR.RDS")

# Reformat DE tables
utr_de_table <- bind_rows(lapply(names(utr_de), function(x) mutate(utr_de[[x]], comparison = x)))
utr_de_table <- utr_de_table[,c(7,1,2,6,8)] %>%
  filter(., abs(logFC) > log2(1.5) & FDR < 0.05) %>%
  mutate(FDR = -log10(FDR)) %>%
  dplyr::rename(ensembl_id = genes, ID = genename, log2FC = logFC, `-log10(FDR)` = FDR)

mirna_de_table <- bind_rows(lapply(names(mirna_de), function(x) mutate(mirna_de[[x]], comparison = x)))
mirna_de_table <- mirna_de_table[,c(1,2,6,ncol(mirna_de_table))] %>%
  filter(., abs(logFC) > log2(1.5) & FDR < 0.05) %>%
  mutate(FDR = -log10(FDR)) %>%
  dplyr::rename(ID = genes, log2FC = logFC, `-log10(FDR)` = FDR)


# Load in GWAS and disease gene lists
pub_genes <- read.table("data/pituitary_puberty_genes_combined_2020-09-28.txt",
                        sep = "\t", header = T)
pub_genes <- dplyr::rename(pub_genes, gene_symbol = MGI.symbol)
splitted <- strsplit(as.character(pub_genes$source), ";")
pub_genes_split <- data.frame(ID = rep.int(pub_genes$gene_symbol, sapply(splitted, length)),
           source = unlist(splitted))



# Parse input gene list
# Precompute ensembl gene ids from gene symbols present as rownames in utr object using biomart.
# Script only loads in the dataframe

# mart <- useMart("ensembl", "mmusculus_gene_ensembl")
# ensemble2gene <- getBM(attributes=c("mgi_symbol","ensembl_gene_id","ensembl_gene_id_version"),
#                        filters = "mgi_symbol",
#                        values = rownames(utr_obj), 
#                        mart = mart)
# saveRDS(ensemble2gene, "data/20211015_ensembl_gene_id_mgi_biomart_conversion.rds")

gene_ensembl_convert <- readRDS("data/20211015_ensembl_gene_id_mgi_biomart_conversion.rds")
# genelist <- "let-7a-5p,mmu-let-7e-5p,ENSMUSG00000027120.7,test"
# genelist <- "let-7a-5p,mmu-let-7e-5p,test"
# genelist <- "ENSMUSG00000027120.7,test" # need to fix
genelist <- "test"
# test_file <- read.table("data/example_input_genes.txt", header = F)
parse_list <- function(genelist, type) {
  if(type == "text") {
    genelist <- unlist(strsplit(genelist, ","))
  }
  
  if(type == "file") {
    genelist <- genelist[,1]
  }

  
  # 1: match gene symbols with utr row names
  keep_genes <- genelist[which(genelist %in% rownames(utr_obj))]
  if(length(keep_genes) > 0) {
    genelist <- genelist[-which(genelist %in% keep_genes)]
  } 
  
  
  # 2: try to convert non-matches from ensembl ID to gene symbol
  # Parse string for "ENSMUSG" to find a mouse ensembl genes
  ens_match <- genelist[grep("ENSMUSG", genelist)]
  ensid_match <- ens_match[grep(".", ens_match, fixed = T)]
  ens_match <- ens_match[-grep(".", ens_match, fixed = T)]
  
  use_gene_ensembl_convert <- filter(gene_ensembl_convert,
                                     ensembl_gene_id %in% ens_match |
                                       ensembl_gene_id_version %in% ensid_match)$mgi_symbol
  keep_genes <- c(keep_genes, use_gene_ensembl_convert) # this is a valid gene list
  
  # 3A: match non-matches with mirna row names
  if(length(ens_match) | length(ensid_match) > 0) {
    genelist <- genelist[-c(which(genelist %in% ens_match), which(genelist %in% ensid_match))]
  }
  
  keep_mirnas <- genelist[which(genelist %in% rownames(mirna_obj))]
  
  if(length(keep_mirnas) > 0) {
    genelist <- genelist[-which(genelist %in% keep_mirnas)]
  }
  
  # 3B: match non-matches with the addition of mmu- prefix

  genelist <- paste0("mmu-", genelist)
  mmu_match <- genelist[which(genelist %in% rownames(mirna_obj))]
  keep_mirnas <- c(keep_mirnas, mmu_match) # this is a valid miRNA list
  
  if(length(mmu_match) > 0) {
    genelist <- genelist[-which(genelist %in% mmu_match)]
  }
  no_match <- gsub("mmu-", "", genelist)
  
  print(keep_genes)
  print(keep_mirnas)
  print(no_match)
  return(list("genes" = unique(keep_genes),
              "mirnas" = unique(keep_mirnas),
              "invalid" = unique(no_match)))
}

# Return expression table
expr_table <- function(genelist, count_data) {
  genes_not_found <- which(!(genelist %in% rownames(count_data)))
  if(length(genes_not_found) > 0) {
    # error <- (paste0(paste0(genelist[genes_not_found],collapse=","), " not found in dataset."))
    use_genelist <- genelist[-genes_not_found]
  } else {  use_genelist <- genelist
  # error <- NULL
  }
  if(length(use_genelist) > 0) {
    return(t(count_data[use_genelist,,drop = F]))
  }
}

# Plot gene/miRNA expression
exprplot_hhtheme <- function(genelist,
                             count_data,
                             metadata,
                             counttype,
                             type = "log2", #Default
                             pal_cols = c("F" = "tomato", "M" = "steelblue")) {
  
  if(type == "log2") {
    ylabel = "log2(normCounts)"
  } else { ylabel = "Normalized Counts"}
  
  use_genelist <- genelist[which(genelist %in% rownames(count_data))]
  if(length(genelist %in% rownames(count_data)) > 0) {
    
    genecounts <- as.data.frame(t(count_data[which(rownames(count_data) %in% use_genelist),, drop = F]))
    # colnames(genecounts) <- use_genelist
    genecounts$sample <- rownames(genecounts)
    
    genecounts <- melt(genecounts)
    
    genecounts <- cbind(genecounts, "sex" = as.character(metadata$sex), "time" = gsub("d", "", metadata$age))
    
    genecounts$time <- as.numeric(as.character(genecounts$time))
    
    genecounts$reps <- paste0(metadata$sex, metadata$rep)
    # print(genecounts)
    gene_med <- aggregate(genecounts$value, by=list(genecounts$sex, genecounts$time, genecounts$variable), median)
    # print(gene_med)
    colnames(gene_med) <- c("sex", "time", "variable", "median")
    if(counttype == "mirnas") {
      x_breaks <- c(12, 22, 27, 32)
    }
    if(counttype == "genes") {
      x_breaks <- c(12, 22, 27, 32, 37)
    }
    
    p <- ggplot() +
      geom_point(data = gene_med, aes(x = time, y = median, fill = sex, shape = sex, group = sex, color = sex), size = 4, alpha = 0.8) +
      geom_jitter(data = genecounts, aes(x = time, y = value, shape = sex, color = sex), width=0.2, height = 0) + 
      geom_line(data=gene_med, aes(x = time, y = median, color = sex, group = sex))  +
      xlab("Age (postnatal days)") +
      ylab(ylabel) +
      scale_fill_manual(values = pal_cols, name = "Sex") +
      scale_color_manual(values = pal_cols, name = "Sex") +
      scale_shape_manual(values = c(21, 24)) + 
      theme_light() +
      facet_wrap( .~variable, scales = "free", ncol = 1) +
      theme(strip.text = element_text(size=18),
            axis.text = element_text(size=18 - 2, color="black"),
            axis.title = element_text(size=18),
            # strip.background = element_rect(color = "gray50") ,
            legend.position = "top", text = element_text(size = 18))+
      # strip.background = element_rect(colour = "red", fill = alpha("blue",0.2) )) +
      scale_x_continuous(breaks = x_breaks) +
      guides(shape = F)
    # print(p)
  } else {
    # print(paste0("Genes/miRNAs not found in count data", genelist))
    p <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label=paste0("No ", counttype, " inputted found in count data.")), size = 5) +
      xlab(NULL) +
      theme(legend.position = "none")
  }
  return(p)
}

print_de_table <- function(genelist,
                           de_table,
                           counttype) {
  if(length(genelist) > 0) {
    dtable <- filter(de_table, ID %in% genelist) %>%
      arrange(desc(abs(log2FC)))
  }
  else {
    dtable <- data.frame(paste0("No ", counttype, " inputted."))
    colnames(dtable) <- ""
  }
  
  if(nrow(dtable) == 0) {
    dtable <- data.frame(paste0(genelist, " is not DE."))
    colnames(dtable) <- ""
  }
  # print(dtable)
  return(dtable)
}

add_study <- function(use_table, study) {
  pub_genes_use <- filter(pub_genes_split, source %in% study)
  utr_de_merge <- left_join(use_table, pub_genes_use, by = "ID")
  utr_de_merge <- utr_de_merge %>% group_by(ID, ensembl_id, log2FC, `-log10(FDR)`, comparison) %>%
    summarise(source  = toString(source)) %>%
    arrange(desc(abs(log2FC)))
  return(utr_de_merge)
}