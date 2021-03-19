library(tidyverse)
library(Rsubread)
library(Biostrings)
library(DESeq2)

setwd("//data.thecrick.org/lab-ulej/working/Oscar/Roberto_Simone/RS_mapped/results/genome_dedup")
longest_transcript_name <- "D:/transcriptome_dedup/bedfiles/longest_proteincoding_transcript_hs_details.txt"
info <- read_tsv(longest_transcript_name)

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

col_data <- data.frame(samples = samples) %>%
  mutate(condition = str_sub(samples,1,1))

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



### Make Volcano plot ####
min_counts <- 170

filtered <- combined %>% 
  filter(A1 +A2+A3+B1+B2+B3+C1+C2+C3+D1+D2+D3+E1+E2+E3+G1+G2 >= min_counts)


matrix <- as.matrix(filtered %>% column_to_rownames(var = "gene_id"))  

dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = col_data,
                              design = ~ condition)



# state reference
i <- "A"

# state comparison
j <- "B"


## add lnc info

#lnc <- read_csv("MIRlnc S_AS gene targets list.csv",
#                col_names = c("gene_id2", "gene_name", "lncRNA", "exon_type",
#                             "orientation_of_lnc", "orientation_of_repeat"))
  
dds$condition <- relevel(dds$condition, ref = i)
dds <- DESeq(dds)
comparison_name <- paste0("condition_", j, "_vs_", i)
res <- results(dds, name=comparison_name)

res_data <- data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(info, by="gene_id")%>%
  mutate(minus_log10_p = -log10(pvalue)) %>%
  mutate(is_MAPT = gene_id == "ENSG00000186868.15") %>%
  mutate(color = ifelse(is_MAPT, "red", "grey60")) %>%
  mutate(name = ifelse(abs(log2FoldChange) >2 | pvalue < 0.001 | gene_name == "MAPT", gene_name, ""))

#ggplot(res_data, aes(log2FoldChange, minus_log10_p)) +
#  geom_point(aes(colour = is_MAPT), size=1, alpha = 1) +
#  ggtitle(comparison_name)

ggplot(res_data, aes(log2FoldChange, minus_log10_p, colour = color, name=name)) +
  geom_point(size=1, alpha = 1) +
  ggtitle(paste0("Volcano plot of condition B versus A with minimum \n", min_counts, " counts total across 17 samples")) +
  scale_colour_identity() +
  theme_minimal() +
  geom_text(aes(label=name),hjust=0.6, vjust=1.2, size =3)

write_csv(res_data, "C:/Users/ogw/Google Drive/UCL PhD/Year 1/Roberto Simone/Figures/code/RFP_A_vs_B_volcano_data.csv")

images_dir <- "//data.thecrick.org/lab-ulej/working/Oscar/Roberto_Simone/figures/image_files/"

ggsave(paste0(images_dir, "RFP volcano plot.svg"))
# svg might be difficult because of huge number of points
ggsave(paste0(images_dir, "RFP volcano plot.png"))





















res_data %>% mutate(target = !is.na(gene_id2)) %>%
  group_by(target) %>%
  mutate(average_P = mean(pvalue)) %>%
  mutate(average_log2fc = mean(log2FoldChange)) %>%
  ungroup() %>% select(target, average_P, average_log2fc) %>% 
  distinct()


mapt_pvalue <- filter(res_data, gene_id == "ENSG00000186868.15")$pvalue
mapt_FC <- filter(res_data, gene_id == "ENSG00000186868.15")$log2FoldChange

nrow(res_data)
nrow(dplyr::filter(res_data, pvalue <= mapt_pvalue & log2FoldChange < 1))
nrow(dplyr::filter(res_data, pvalue < 0.05, log2FoldChange < mapt_FC))
nrow(dplyr::filter(res_data, pvalue < 0.01))



yo <- dplyr::filter(res_data, pvalue < 0.05, log2FoldChange < -1)


### gather to make tidy dataframe ###

tidy <- combined %>%
  gather("Sample", "Counts", -gene_id) %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(Counts)) %>%
  ungroup() %>%
  mutate(normalised_counts = 1000000*Counts/sample_total)

tidy_filtered <- filter(tidy, gene_id == "ENSG00000186868.15") %>%
  mutate(condition = str_sub(Sample,1,1)) %>%
  group_by(Sample) #%>%
  #mutate(sample_mean = mean(Counts), sample_sd = sd(counts))

ggplot(tidy_filtered, aes(condition, log(normalised_counts))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center')


filter(res_data, pvalue < 0.01) %>%
  left_join(info, by="gene_id") %>%
  select(gene_name)


#### Extract fasta for motif analysis #######

nt_seq_data <- data.frame(seq = readDNAStringSet("D:/transcriptome_dedup/bedfiles/gencodev29transcriptome_fixed_longestonly.fasta")) %>%
  rownames_to_column(var = "transcript_id") %>%
  left_join(info, by = "transcript_id") %>%
  mutate(fivep_utr_seq = str_sub(seq,1,cds_start)) %>%
  select(gene_id, fivep_utr_seq)

# filter for sig genes

sig_res <- res_data %>% filter(pvalue < 0.05 & log2FoldChange < 0)

fivep_utr_sig <- filter(nt_seq_data,
                   gene_id %in% sig_res$gene_id) %>%
  select(name = gene_id, seq = fivep_utr_seq)

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

writeFasta(fivep_utr_sig, "sig_down_5p_utr.fa")

fivep_utr_notsig <- filter(nt_seq_data,
                        !(gene_id %in% sig_res$gene_id)) %>%
  select(name = gene_id, seq = fivep_utr_seq) %>%
  sample_n(500)

writeFasta(fivep_utr_notsig, "not_sig_5p_utr.fa")

sig_ss <- readDNAStringSet("sig_down_5p_utr.fa")
notsig_ss <- readDNAStringSet("not_sig_5p_utr.fa")


library(motifRG)
category <- c(rep(1,nrow(data.frame(sig_ss))), rep(0, nrow(data.frame(notsig_ss))))
motifs <- motifRG::findMotif(append(sig_ss, notsig_ss),category = category,
                         max.motif=3,enriched=T, min.ratio=1)

motifRG::summaryMotif(motifs$motifs, motifs$category)

### search for specific sequences ####

sequences = c()
nt_list = c("A", "C", "G", "T")

for(i in nt_list){
  for(j in nt_list){
    for(k in nt_list){
      for(l in nt_list){
        for(m in nt_list){
          for(n in nt_list){
            for(o in nt_list){
          
          this_seq <- paste0(i,j,k,l,m,n,o)
          sequences <- c(sequences,this_seq)
            }
          }
        }
      }
    }
  }
}

motifs <- data.frame("fivemers" = sequences, stringsAsFactors = F) %>%
  mutate(background = 0, sig = 0)

# Now, for each motif, test how many times it appears

for(i in 1:nrow(motifs)){
  this_motif <- motifs$fivemers[i]
  background_counts <- sum(str_count(fivep_utr_notsig$seq, this_motif))
  sig_counts <- sum(str_count(fivep_utr_sig$seq, this_motif))
  
  motifs$background[i] <- background_counts
  motifs$sig[i] <- sig_counts
  
}

motifs2 <- motifs %>%
  mutate(background_normalised = background/sum(background),
         sig_normalised = sig/sum(sig)) %>%
  filter(background > 0) %>%
  mutate(sig_over_background = sig_normalised/background_normalised) %>%
  mutate(av_no = sig_normalised + background_normalised) %>%
  filter(background > mean(background) & sig > mean(sig))









