library(tidyverse)
library(DESeq2)
library(here)

setwd("A:/home/users/wilkino/ALL_RIBOSEQ/Roberto_Simone/quantseq/Mapped")

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
             "G2",
             "G3")

for(i in samples){
  print(i)
  
  filename <- paste0(i, ".Aligned.sortedByCoord.out.bam.cds_counts.gz")
  filename_with_dir <- paste0(filename)
  
  trimmed <- str_sub(filename, 1, (str_length(filename)-14))
  this_data <- read.csv(filename_with_dir, sep="", head=T, skip=1, row.names = "Geneid") %>%
    dplyr::select(!!i := trimmed) %>%
    rownames_to_column(var = "gene_id")
  
  if(i==samples[1]){
    combined <- this_data
  }else{
    combined <- combined %>%
      full_join(this_data, by = "gene_id")
  }
  
}

col_data <- data.frame(samples = colnames(combined)[2:19]) %>%
  mutate(condition = str_sub(samples,1,1))

### Make Volcano plot ####
min_counts <- 1000

col_data <- data.frame(samples = samples) %>%
  mutate(condition = str_sub(samples,1,1))

filtered <- combined %>% 
  filter(A1 +A2+A3+B1+B2+B3+C1+C2+C3+D1+D2+D3+E1+E2+E3+G1+G2+G3 >= min_counts)

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
  mutate(colour = ifelse(is_MAPT, "red", "grey60"),
         alpha =ifelse(is_MAPT, 1, 0.2),
         colour2 = ifelse(padj<0.1, "blue", colour)) %>%
  mutate(name = ifelse(abs(log2FoldChange) >2.5 | pvalue < 0.0001 | gene_name == "MAPT", gene_name, ""))

ggplot(res_data, aes(log2FoldChange, minus_log10_p, color = colour2, alpha =1, name = name)) +
  geom_point() +
  ggtitle(comparison_name) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_text(aes(label=name),hjust=0.6, vjust=1.2, size =3) +
  ggtitle(paste0("Quantseq B versus A for min total \ncounts of ", min_counts, " across 18 samples")) +
  theme_minimal()

write_csv(res_data, paste0(here(), "/data/quantseq_B_vs_A_volcano_data.csv"))







