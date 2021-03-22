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

min_counts_per_sample_average <- 10
min_counts <- 10*17

filtered <- combined %>% 
  filter(A1 +A2+A3+B1+B2+B3+C1+C2+C3+D1+D2+D3+E1+E2+E3+G1+G2 >= min_counts) # verified that this gives identical results to dds[keep,] method on DESeq2 vignette 

matrix <- as.matrix(filtered %>% column_to_rownames(var = "gene_id"))  

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



## RFP volcano plot

matrix <- as.matrix(filtered %>% column_to_rownames(var = "gene_id"))  

dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = col_data,
                              design = ~ condition)

# state reference
i <- "A"

# state comparison
j <- "B"

dds$condition <- relevel(dds$condition, ref = i)
dds <- DESeq(dds)
comparison_name <- paste0("condition_", j, "_vs_", i)
res <- results(dds, name=comparison_name)

info <- read_tsv(paste0(here(), "/data/longest_proteincoding_transcript_hs_details.txt"))

res_data <- data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(info, by="gene_id")%>%
  mutate(minus_log10_p = -log10(pvalue)) %>%
  mutate(is_MAPT = gene_id == "ENSG00000186868.15") %>%
  mutate(color = ifelse(is_MAPT, "red", "grey60")) %>%
  mutate(name = ifelse(abs(log2FoldChange) >2 | pvalue < 0.001 | gene_name == "MAPT", gene_name, "")) %>%
  arrange(log2FoldChange)

ggplot(res_data, aes(log2FoldChange, minus_log10_p, colour = color, name=name)) +
  geom_point(size=1, alpha = 1) +
  ggtitle(paste0("Volcano plot of condition B versus A with minimum \n", min_counts, " counts total across 17 samples")) +
  scale_colour_identity() +
  theme_minimal() +
  geom_text(aes(label=name),hjust=0.6, vjust=1.2, size =3)

write_csv(res_data, paste0(here(), "/data/RFP_B_vs_A_volcano_data.csv"))

