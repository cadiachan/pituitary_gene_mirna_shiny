# Libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)

# Load in UTR and miRNA counts
# utr_obj <- readRDS("data/pit_utr_2019__RUV_k2_set1_2019-07-03.rds")
# colnames(utr_obj) <- substr(colnames(utr_obj),8,50)
# saveRDS(list(norm = normCounts(utr_obj),
#              meta = pData(utr_obj)), "data/pit_utr_2019__RUV_k2_set1_2019-07-03_counts.rds")
utr_obj <- readRDS("data/pit_utr_2019__RUV_k2_set1_2019-07-03_counts.rds")
utr_normcounts <- utr_obj[["norm"]]
utr_log2 <- log2(utr_normcounts + 1)
utr_meta <- utr_obj[["meta"]]

# mirna_obj <- readRDS("data/RUV_corrected_pit_srna-seq_counts_combined.rds")
# saveRDS(list(norm = normCounts(mirna_obj),
#              meta = pData(mirna_obj)), "data/RUV_corrected_pit_srna-seq_counts_combined_counts.rds")
mirna_obj <- readRDS("data/RUV_corrected_pit_srna-seq_counts_combined_counts.rds")
mirna_normcounts <- mirna_obj[["norm"]]
mirna_log2 <- log2(mirna_normcounts + 1)
mirna_meta <- mirna_obj[["meta"]]

# Load in UTR and miRNA DE comparisons
utr_de <- readRDS("data/pit_utr_2019_de_result_list_2019-07-03.rds")
mirna_de <- readRDS("data/mirna_pit_all_edgeR.RDS")

# Rename comparisons in utr_de so that they match comparisons in mirna_de
names(utr_de) <- c(names(utr_de)[1:5],
                   "d22_12_M",
                   "d27_22_M",
                   "d32_27_M",
                   "d37_32_M",
                   "d22_12_F",
                   "d27_22_F",
                   "d32_27_F",
                   "d37_32_F",
                   "d37_12_M",
                   "d37_12_F",
                   "d37_22_M",
                   "d37_22_F")

# Reformat DE tables
utr_de_table <- bind_rows(lapply(names(utr_de), function(x) mutate(utr_de[[x]], comparison = x)))
utr_de_all_table <- utr_de_table[,c(7,1,2,6,8)] %>%
  mutate(FDR = -log10(FDR)) %>%
  dplyr::rename(ensembl_id = genes, ID = genename, log2FC = logFC, `-log10(FDR)` = FDR)
utr_de_table <- utr_de_table[,c(7,1,2,6,8)] %>%
  filter(., abs(logFC) > log2(1.5) & FDR < 0.05) %>%
  mutate(FDR = -log10(FDR)) %>%
  dplyr::rename(ensembl_id = genes, ID = genename, log2FC = logFC, `-log10(FDR)` = FDR)

mirna_de_table <- bind_rows(lapply(names(mirna_de), function(x) mutate(mirna_de[[x]], comparison = x)))
mirna_de_all_table <-  mirna_de_table[,c(1,2,6,ncol(mirna_de_table))] %>%
  mutate(FDR = -log10(FDR)) %>%
  dplyr::rename(ID = genes, log2FC = logFC, `-log10(FDR)` = FDR)
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
# genelist <- "let-7a-5p,mmu-let-7e-5p,test,fshb,ENsmUSG00000027120.7,mmu-miR-224-5p,mir-224-5p,MIR-383-5P"
# genelist <- "let-7a-5p,mmu-let-7e-5p,test"
# genelist <- "ENSMUSG00000027120.7,test" # need to fix
# genelist <- "test"
# test_file <- read.table("data/example_input_genes.txt", header = F)
parse_list <- function(genelist, type) {
  if(type == "text") {
    genelist <- unlist(strsplit(genelist, ","))
    genelist <- gsub(" ", "", genelist) # In case user inputs comma separated with space
  }
  
  if(type == "file") {
    genelist <- genelist[,1]
  }
  
  
  # 1: match gene symbols with utr row names
  keep_genes <- str_to_sentence(genelist[which(tolower(genelist) %in% tolower(rownames(utr_normcounts)))])
  if(length(keep_genes) > 0) {
    genelist <- genelist[-which(tolower(genelist) %in% tolower(keep_genes))]
  } 
  
  
  # 2: try to convert non-matches from ensembl ID to gene symbol
  # Parse string for "ENSMUSG" to find a mouse ensembl genes
  ens_match <- genelist[grep("ENSMUSG", genelist, ignore.case = T)] 
  ensid_match <- ens_match[grep(".", ens_match, fixed = T)]
  ens_match <- ens_match[-grep(".", ens_match, fixed = T)]
  
  
  use_gene_ensembl_convert <- filter(gene_ensembl_convert,
                                     ensembl_gene_id %in% toupper(ens_match) |   # Ensure ENS matches are upper case
                                       ensembl_gene_id_version %in% toupper(ensid_match))$mgi_symbol   # Ensure ENS matches are upper case
  keep_genes <- c(keep_genes, use_gene_ensembl_convert) # this is a valid gene list
  
  # 3A: match non-matches with mirna row names
  if(length(ens_match) | length(ensid_match) > 0) {
    genelist <- genelist[-c(which(genelist %in% ens_match), which(genelist %in% ensid_match))]
  }
  
  genelist <- gsub("mir", "miR", tolower(genelist), ignore.case = T) # Allows users to enter mirna ids without case sensitive format
  keep_mirnas <- genelist[which(genelist %in% rownames(mirna_normcounts))]
  
  if(length(keep_mirnas) > 0) {
    genelist <- genelist[-which(genelist %in% keep_mirnas)]
  }
  
  # 3B: match non-matches with the addition of mmu- prefix
  
  genelist <- paste0("mmu-", genelist)
  mmu_match <- genelist[which(genelist %in% rownames(mirna_normcounts))]
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
    if(length(use_genelist) > 50) {
      genesub <- order(colSums(genecounts), decreasing = T)[1:50]
      genecounts <- genecounts[,genesub]
    }
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
      study_lab <- unique(left_join(dplyr::select(gene_med, variable),
                             dplyr::rename(pub_genes_split, variable = ID),
                             by = "variable")) %>%
        group_by(variable) %>%
        summarise(study_id  = toString(study_id))
      study_lab$study_id <- gsub(" ", "", study_lab$study_id)
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
      theme(strip.text = element_text(size=16),
            axis.text = element_text(size=16 - 2, color="black"),
            axis.title = element_text(size=16),
            legend.position = "top", text = element_text(size = 16))+
      scale_x_continuous(breaks = x_breaks) +
      guides(shape = F)
    
    if(counttype == "genes") {
      p <- p + 
        geom_label(data = study_lab, aes(x = -Inf, y = Inf, label = paste0("StudyID: ", study_id)),
                  hjust = -0.1, vjust = 1.2, size = 5, alpha = 0.5,
                label.size = 0.005)
    }
  } else {
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
  } else {
    dtable <- data.frame(paste0("No ", counttype, " inputted."))
    colnames(dtable) <- ""
  }
  
  if(nrow(dtable) == 0) {
    dtable <- data.frame(paste0(genelist, " does not pass set cutoff."))
    colnames(dtable) <- ""
  } 
  # else {
  #   if(counttype == "genes") {
  #     pub_genes_use <- filter(pub_genes_split, source %in% unique(pub_genes_split$source))
  #     dtable <- left_join(dtable, pub_genes_use, by = "ID")
  #     dtable <- dtable %>% group_by(ID, ensembl_id, log2FC, `-log10(FDR)`, comparison) %>%
  #       summarise(source  = toString(source)) %>%
  #       arrange(desc(abs(log2FC)))
  #   }
  # }
  return(dtable)
}

pub_genes_key <- as.data.frame(matrix(c(1:9, c("Perry2014", "Day2015_GWAS_VB",
                                               "Day2017_nearest", "Hollis2020_GWAS_VB_FH",
                                               "Ye2015_PA_GWAS", "Fang2016_CPHD",
                                               "Hauser2019_PA", "Kurtoglu2019_hypopituitarism",
                                               "IHH/Kallmann")), ncol = 2,
                                      dimnames = list(NULL,c("study_id", "source"))))
pub_genes_split <- left_join(pub_genes_split, pub_genes_key, by = "source")

add_study <- function(use_table, study=unique(pub_genes_split$study_id)) {
  if(ncol(use_table) > 1) { # Check to see that dataframe has DE info and is not just an error message
    pub_genes_use <- filter(pub_genes_split, study_id %in% study)
    utr_de_merge <- left_join(use_table, pub_genes_use, by = "ID")
    utr_de_merge <- utr_de_merge %>% group_by(ID, ensembl_id, log2FC, `-log10(FDR)`, comparison) %>%
      summarise(study_id = toString(study_id)) %>%
      arrange(desc(abs(log2FC)))
    return(utr_de_merge)
  } else { return(use_table)}
}

add_study_corr <- function(use_table, study=unique(pub_genes_split$study_id),
                           de_filt = T) {
  pub_genes_use <- filter(pub_genes_split, study_id %in% study)
  corr_merge <- dplyr::rename(use_table, ID = gene)
  corr_merge <- dplyr::left_join(corr_merge, pub_genes_use, by = "ID")
  if(de_filt) {
    corr_merge <- corr_merge %>% group_by(pair, rho, fdr,mirna,ID, database,gene_comparison, mirna_comparison) %>%
      summarise(study_id  = toString(study_id)) %>%
      arrange(desc(abs(rho)))
  }
  else {
    corr_merge <- corr_merge %>% group_by(pair, rho, fdr,mirna,ID, database) %>%
      summarise(study_id  = toString(study_id)) %>%
      arrange(desc(abs(rho)))
  }
  return(corr_merge)
}

corr <- readRDS("data/pit_mirna_mrna_pairs_w_correlation.rds")
corr <- corr[, -grep("Pd", colnames(corr))]
corr <- corr[, -grep("pval", colnames(corr))]
corr <- filter(corr, rho < 0 & fdr < 0.1)
corr <- mutate(corr, fdr = -log10(fdr)) %>%
  dplyr::rename(`-log10(FDR)` = fdr)
utr_de_corr <- mutate(utr_de_table, comparison_dir = ifelse(log2FC > 0,
                                                            "up",
                                                            "down")) %>%
  dplyr::select(ID, comparison, comparison_dir) %>%
  dplyr::rename(gene = ID, gene_comparison = comparison, gene_dir = comparison_dir)
mirna_de_corr <- mutate(mirna_de_table, comparison_dir = ifelse(log2FC > 0,
                                                                "up",
                                                                "down")) %>% 
  dplyr::select(ID, comparison, comparison_dir) %>%
  dplyr::rename(mirna = ID, mirna_comparison = comparison, mirna_dir = comparison_dir)
# mirnalist <- c("mmu-miR-224-5p", "mmu-miR-383-5p", "mmu-let-7a-5p",
#                sample(rownames(mirna_log2), 20))
# genelist <- c("Fshb", "Lhb", "Ammecr1",
#               sample(rownames(utr_log2), 20))
print_corr_table <- function(mirnalist = NULL,
                             genelist = NULL,
                             de_filt = T) {
  if(de_filt) {
    if(!is.null(mirnalist)) {
      if(length(which(mirnalist %in% corr$mirna)) > 0) {
        test_de <- length(which(unique(mirnalist) %in% mirna_de_table$ID))
        if(test_de == 0) {
          cortable <- data.frame(paste0("No DE miRNAs inputted."))
          colnames(cortable) <- ""
          return(cortable)
        }
        else {
          cortable <- filter(corr, mirna %in% mirnalist)
        }
      } else {
        cortable <- data.frame(paste0(mirnalist, " not found in correlation table."))
        colnames(cortable) <- ""
        return(cortable)
      }
    }
    else if(!is.null(genelist)) {
      if(length(which(genelist %in% corr$gene)) > 0) {
        test_de <- length(which(unique(genelist) %in% utr_de_table$ID))
        if(test_de == 0) {
          cortable <- data.frame(paste0("No DE genes inputted."))
          colnames(cortable) <- ""
          return(cortable)
        }
        else {
          cortable <- filter(corr, gene %in% genelist)
        }
      } else {
        cortable <- data.frame(paste0(genelist, " not found in correlation table."))
        colnames(cortable) <- ""
        return(cortable)
      }
      
    }
    else{ 
      cortable <- data.frame(paste0("Please input mirna list or gene list"))
      colnames(cortable) <- ""
      return(cortable)
    }
    cortable <- left_join(cortable, utr_de_corr, by = "gene") %>%
      left_join(., mirna_de_corr, by = "mirna") 
    sextable <- filter(cortable, grepl("sex", gene_comparison) & grepl("sex", mirna_comparison)) %>%
      filter(gene_dir != mirna_dir)
    agetable <- filter(cortable,!grepl("sex", gene_comparison)) %>%
      mutate(gene_comparison_nosex = substr(gene_comparison, 1,6),
             mirna_comparison_nosex = substr(mirna_comparison, 1,6)) %>%
      filter(gene_comparison_nosex == mirna_comparison_nosex & gene_dir != mirna_dir) %>%
      dplyr::select(-gene_comparison_nosex, -mirna_comparison_nosex)
    cortable <- bind_rows(sextable, agetable) %>%
      mutate(gene_comparison = paste0(gene_comparison, "_", gene_dir),
             mirna_comparison = paste0(mirna_comparison, "_", mirna_dir)) %>%
      dplyr::select(-gene_dir,-mirna_dir)
    cortable_aggr <- aggregate(cortable$gene_comparison,
                               by = list(cortable$pair),
                               function(x) toString(unique(x))) %>%
      dplyr::rename(pair = Group.1, gene_comparison = x)
    cortable_aggr_mirna <- aggregate(cortable$mirna_comparison,
                                     by = list(cortable$pair),
                                     function(x) toString(unique(x))) %>%
      dplyr::rename(pair = Group.1, mirna_comparison = x)
    cortable <- unique(inner_join(dplyr::select(cortable, -gene_comparison, -mirna_comparison),
                                  cortable_aggr, by = "pair") %>%
                         inner_join(., cortable_aggr_mirna, by = "pair"))
  }
  else {
    if(!is.null(mirnalist)) {
      if(length(which(mirnalist %in% corr$mirna)) > 0) {
        cortable <- filter(corr, mirna %in% mirnalist)
      } else {
        cortable <- data.frame(paste0(mirnalist, " not found in correlation table."))
        colnames(cortable) <- ""
        return(cortable)
      }
    }
    else if(!is.null(genelist)) {
      if(length(which(genelist %in% corr$gene)) > 0) {
        cortable <- filter(corr, gene %in% genelist)
      } else {
        cortable <- data.frame(paste0(genelist, " not found in correlation table."))
        colnames(cortable) <- ""
        return(cortable)
      }
    }
    else{ 
      cortable <- data.frame(paste0("Please input mirna list or gene list"))
      colnames(cortable) <- ""
      return(cortable)
    }
    cortable <- left_join(cortable, utr_de_corr, by = "gene") %>%
      left_join(., mirna_de_corr, by = "mirna")
    cortable <- mutate(cortable, gene_comparison = paste0(gene_comparison, "_", gene_dir),
                       mirna_comparison = paste0(mirna_comparison, "_", mirna_dir)) %>%
      dplyr::select(-gene_dir,-mirna_dir)
    cortable_aggr <- aggregate(cortable$gene_comparison,
                               by = list(cortable$pair),
                               function(x) toString(unique(x))) %>%
      dplyr::rename(pair = Group.1, gene_comparison = x)
    cortable_aggr_mirna <- aggregate(cortable$mirna_comparison,
                                     by = list(cortable$pair),
                                     function(x) toString(unique(x))) %>%
      dplyr::rename(pair = Group.1, mirna_comparison = x)
    cortable <- unique(inner_join(dplyr::select(cortable, -gene_comparison, -mirna_comparison),
                                  cortable_aggr, by = "pair") %>%
                         inner_join(., cortable_aggr_mirna, by = "pair"))
    cortable <- mutate(cortable, gene_comparison = ifelse(gene_comparison == "NA_NA", NA, gene_comparison)) %>%
      mutate(mirna_comparison = ifelse(mirna_comparison == "NA_NA", NA, mirna_comparison))
    print(cortable)
  }
  return(cortable)
}



