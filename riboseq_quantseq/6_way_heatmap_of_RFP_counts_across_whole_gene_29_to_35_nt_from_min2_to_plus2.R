library(tidyverse)
library(DESeq2)
library(here)

setwd("A:/home/users/wilkino/ALL_RIBOSEQ/Roberto_Simone/RS_mapped/results/genome_dedup/")

samples <- c("A1",
                      "A2",
                      "A3",
                      "B1",
                      "B2",
                      "B3",
                      "C1",
                      "C2",
                      "C3",
                      "D1",
                      "D2",
                      "D3",
                      "E1",
                      "E2",
                      "E3",
                      "G1",
                      "G2")

for(i in samples){
  message(i)
  
  filename <- paste0(i, ".Aligned.sortedByCoord.out.bam29to35nt.sam.bam.rfp_counts.gz")
  filename_with_dir <- filename
                   
  trimmed <- str_sub(filename, 1, (str_length(filename)-14))
  this_data <- read.csv(filename_with_dir, sep="", head=T, skip=1, row.names = "Geneid") %>%
    select(!!i := trimmed) %>%
    rownames_to_column(var = "gene_id")
  
  if(i==samples[1]){
    combined <- this_data
  }else{
    combined <- combined %>%
      full_join(this_data, by = "gene_id")
  }
  
}

col_data <- data.frame(samples = colnames(combined)[2:18]) %>%
  mutate(condition = str_sub(samples,1,1))

#ENSG00000186868.15 is mapt

matrix <- as.matrix(combined %>% column_to_rownames(var = "gene_id"))  

dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = col_data,
                              design = ~ condition)
# Do all the comparisons

sample_names <- c("A", "B", "C", "D", "E", "G")
all_comparisons <- data.frame(comparison = c(), reference = c(), pvalue = c(), log2fold = c())

for(i in sample_names){
  dds$condition <- relevel(dds$condition, ref = i)
  dds <- DESeq(dds)
  
  for(j in sample_names){
    if (j != i){
      comparison_name <- paste0("condition_", j, "_vs_", i)
      
      res <- results(dds, name=comparison_name)
      
      mapt <- data.frame(res) %>% 
        rownames_to_column(var = "gene_id") %>%
        filter(gene_id == "ENSG00000186868.15")
      
        this_comparison = data.frame(comparison = j, 
                                     reference = i, 
                                     pvalue = mapt$pvalue, 
                                     log2fold = mapt$log2FoldChange)
        
        all_comparisons <- bind_rows(all_comparisons, this_comparison)
        
    }
  }
}

square_pvalue <- all_comparisons %>% 
  dplyr::select(-log2fold) %>%
  pivot_wider(names_from = "reference", values_from = "pvalue") %>%
  arrange(comparison)

write_csv(square_pvalue, paste0(here(), "/data/square_pvalues.csv"))

square_lg2 <- all_comparisons %>% 
  dplyr::select(-pvalue) %>%
  pivot_wider(names_from = "reference", values_from = "log2fold") %>%
  arrange(comparison)

write_csv(square_lg2, paste0(here(), "/data/square_lg2.csv"))

## Save 6 Way heatmap ##

ggplot(all_comparisons, aes(x = reference, y = comparison)) + 
  geom_tile(aes(fill=log2fold)) +
  scale_fill_gradient2(low="blue", mid= "white" ,high="red", midpoint=0,
                       limits = c(-2,2), breaks =c(-2, -1, 0, 1, 2))  +
  theme_minimal() +
  ggtitle("Pair-wise comparison of number of RFPs aligning to MAPT gene") +
  xlab("Reference Sample") +
  ylab("Comparison Sample")
