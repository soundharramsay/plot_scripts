########################### MA plot and tric color MA-plot (search blue for tricolor)
### stringties counts
stringtie_count <- read.csv("./stringtie_gene_count_matrix.csv",header = T) 
colnames(stringtie_count) <- paste(colnames(stringtie_count),"count",sep = "_")
colnames(stringtie_count)[1] <- "gene_id"
###


#### deseq normalized 
## 
setwd("./")
str.LT_D12_nt_vs_e13_deseq.norm <- read.table("./stringtieall.normalised_counts.tsv",header = T)
str.LT_D12_nt_vs_e13_deseq.res <- read.table("./stringtiecondition_control_treated.deseq2.results.tsv",header=T)
##### marking colnames with id before merging 
colnames(str.LT_D12_nt_vs_e13_deseq.norm) <- paste(colnames(str.LT_D12_nt_vs_e13_deseq.norm),"str_LT_nt_vs_e13_deseq.norm",sep="-")
colnames(str.LT_D12_nt_vs_e13_deseq.res) <- paste(colnames(str.LT_D12_nt_vs_e13_deseq.res),"str_LT_nt_vs_e13_deseq",sep="-")
colnames(str.LT_D12_nt_vs_e13_deseq.res)[1] <- "gene_id"
colnames(str.LT_D12_nt_vs_e13_deseq.norm)[1] <- "gene_id"
#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames

LT_d12_nt_vs_e13.merged_df <- merge(str.LT_D12_nt_vs_e13_deseq.res, str.LT_D12_nt_vs_e13_deseq.norm, by = "gene_id")
LT_d12_nt_vs_e13.merged_df_deseq_count <- merge(LT_d12_nt_vs_e13.merged_df,stringtie_count,by = "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(LT_d12_nt_vs_e13.merged_df_deseq_count)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


strLT_d12_nt_vs_e13.mer.df_deseq_count <- na.omit(LT_d12_nt_vs_e13.merged_df_deseq_count)

##### mean of Deseq normalized for nt
colnames(strLT_d12_nt_vs_e13.mer.df_deseq_count)
strLT_d12_nt_vs_e13.mer.df_deseq_count$mean.nt.deseq <- rowMeans(strLT_d12_nt_vs_e13.mer.df_deseq_count[, 8:9], na.rm = TRUE)



library(ggrepel)
library(dplyr)

##### basemean cutoff
strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean <- strLT_d12_nt_vs_e13.mer.df_deseq_count %>% filter(mean.nt.deseq >= 10)

###########
count_below_001 <- sum(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj <= 0.01)

min(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 10.58
max(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 271195.3
############## labelling 
# Create an empty vector of labels
### add index 
#strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean<- cbind(index = 1:13853, strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean)
#rownames(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean) <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id


##### writing out the files
str.d12.nt_Vs_e13_LT.control.upregulated_genes <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange > 0 & strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]
str.d12.nt_Vs_e13_LT.control.downregulated_genes <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange < 0 & strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]

write.csv(str.d12.nt_Vs_e13_LT.control.upregulated_genes,file = "str.d12.nt_Vs_e13_LT.control.upregulated_genes.csv")
write.csv(str.d12.nt_Vs_e13_LT.control.downregulated_genes,file = "str.d12.nt_Vs_e13_LT.control.downregulated_genes.csv")


library(ggplot2)

# Your data
str.xval <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq
str.yval <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange
str.padj <- strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj

# Create a data frame
str.df <- data.frame(x = str.xval, y = str.yval, padj = str.padj)

# MA plot
# Calculate counts
str.upregulated_count <- sum(str.df$y > 0 & str.df$padj < 0.01)
str.downregulated_count <- sum(str.df$y < 0 & str.df$padj < 0.01)
str.nonsignificant_count <- sum(str.df$padj >= 0.01)

library(ggplot2)
library(ggrepel)
library(scales)

str.gg <- ggplot(str.df, aes(x = x, y = y, color = factor(ifelse(padj < 0.01, 
                                                                 ifelse(y > 0, "Upregulated", "Downregulated"),
                                                                 "Nonsignificant")))) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "Stringtie counts",
       subtitle = "Human iNeurons Day 12 (Wildtype vs. ZSWIM8-KO)") +
  geom_text_repel(data = data.frame(x = str.xval[which(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    y = str.yval[which(strLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    label = "ZSWIM8"),
                  aes(x = x, y = y, label = label),
                  color = "black",
                  nudge_x = 0.5,  # changing ZSWIM8 label position
                  nudge_y = -0.5,
                  size = 4) +
  scale_x_log10(labels = scales::trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm"),  # Adjust tick length here
        legend.position = "top",
        legend.justification = c(1, 1),
        legend.box.just = "right")

str.gg <- str.gg +
  geom_text(data = data.frame(x = 400000, y = 4.5, label = paste("Upregulated:", str.upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 4, label = paste("Downregulated:", str.downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 3.5, label = paste("Nonsignificant:", str.nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 1)

# Print the plot
print(str.gg)




############################ Salmon counts

### salmons counts
salmon_count <- read.csv("./3_rounded_salmon.merged.gene_counts.csv",header = T) 
colnames(salmon_count) <- paste(colnames(salmon_count),"count",sep = "_")
colnames(salmon_count)[1] <- "gene_id"
###


#### deseq normalized 
## 
setwd("./")
sal.LT_D12_nt_vs_e13_deseq.norm <- read.table("./salmonall.normalised_counts.tsv",header = T)
sal.LT_D12_nt_vs_e13_deseq.res <- read.table("./salmoncondition_control_treated.deseq2.results.tsv",header=T)
##### marking colnames with id before merging 
colnames(sal.LT_D12_nt_vs_e13_deseq.norm) <- paste(colnames(sal.LT_D12_nt_vs_e13_deseq.norm),"sal_LT_nt_vs_e13_deseq.norm",sep="-")
colnames(sal.LT_D12_nt_vs_e13_deseq.res) <- paste(colnames(sal.LT_D12_nt_vs_e13_deseq.res),"sal_LT_nt_vs_e13_deseq",sep="-")
colnames(sal.LT_D12_nt_vs_e13_deseq.res)[1] <- "gene_id"
colnames(sal.LT_D12_nt_vs_e13_deseq.norm)[1] <- "gene_id"
#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames
salLT_d12_nt_vs_e13.merged_df <- merge(sal.LT_D12_nt_vs_e13_deseq.res, sal.LT_D12_nt_vs_e13_deseq.norm, by = "gene_id")
salLT_d12_nt_vs_e13.merged_df_deseq_count <- merge(salLT_d12_nt_vs_e13.merged_df,salmon_count,by = "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(salLT_d12_nt_vs_e13.merged_df_deseq_count)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


salLT_d12_nt_vs_e13.merged_df_deseq_count <- na.omit(salLT_d12_nt_vs_e13.merged_df_deseq_count)

##### mean of Deseq normalized for nt
colnames(salLT_d12_nt_vs_e13.merged_df_deseq_count)
salLT_d12_nt_vs_e13.merged_df_deseq_count$mean.nt.deseq <- rowMeans(salLT_d12_nt_vs_e13.merged_df_deseq_count[, 8:9], na.rm = TRUE)



library(ggrepel)
library(dplyr)

##### basemean cutoff
salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean <- salLT_d12_nt_vs_e13.merged_df_deseq_count %>% filter(mean.nt.deseq >= 10)

###########
count_below_001 <- sum(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$`padj-sal_LT_nt_vs_e13_deseq` <= 0.01)

min(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 10.00204
max(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 108043.6
############## labelling 
# Create an empty vector of labels
### add index 
#salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean<- cbind(index = 1:13853, salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean)
#rownames(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean) <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id


##### writing out the files
sal.d12.nt_Vs_e13_LT.control.upregulated_genes <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange > 0 & salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]
sal.str.d12.nt_Vs_e13_LT.downregulated_genes <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange < 0 & salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]

write.csv(sal.d12.nt_Vs_e13_LT.control.upregulated_genes,file = "sal.d12.nt_Vs_e13_LT.control.upregulated_genes")
write.csv(sal.str.d12.nt_Vs_e13_LT.downregulated_genes,file = "sal.str.d12.nt_Vs_e13_LT.downregulated_genes.csv")


library(ggplot2)
library(ggrepel)
library(scales)

# Your data
sal.xval <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq
sal.yval <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange
sal.padj <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj

# Create a data frame
sal.df <- data.frame(x = sal.xval, y = sal.yval, padj = sal.padj)

# MA plot
# Calculate counts
sal.upregulated_count <- sum(sal.df$y > 0 & sal.df$padj < 0.01)
sal.downregulated_count <- sum(sal.df$y < 0 & sal.df$padj < 0.01)
sal.nonsignificant_count <- sum(sal.df$padj >= 0.01)

sal_gg <- ggplot(sal.df, aes(x = x, y = y, color = factor(ifelse(padj < 0.01, 
                                                                 ifelse(y > 0, "Upregulated", "Downregulated"),
                                                                 "Nonsignificant")))) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "salmon counts",
       subtitle = "Human iNeurons Day 12 (Wildtype vs. ZSWIM8-KO)") +
  geom_text_repel(data = data.frame(x = sal.xval[which(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    y = sal.yval[which(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    label = "ZSWIM8"),
                  aes(x = x, y = y, label = label),
                  color = "black",
                  nudge_x = 0.5,  # changing ZSWIM8 label position
                  nudge_y = -0.5,
                  size = 4) +
  scale_x_log10(labels = scales::trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm"),  # Adjust tick length here
        legend.position = "top",
        legend.justification = c(1, 1),
        legend.box.just = "right")

sal_gg <- sal_gg +
  geom_text(data = data.frame(x = 400000, y = 4.5, label = paste("Upregulated:", sal.upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 4, label = paste("Downregulated:", sal.downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 3.5, label = paste("Nonsignificant:", sal.nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 1)

# Print the plot
print(sal_gg)

########### featurecount_FC

### FCs counts
FC_count <- read.csv("./Feature_gene_counts.csv",header = T) 
colnames(FC_count) <- paste(colnames(FC_count),"count",sep = "_")
colnames(FC_count)[1] <- "gene_id"
###


#### deseq normalized 
## 
setwd("./")
FC.LT_D12_nt_vs_e13_deseq.norm <- read.table("./all.normalised_counts.tsv_featurecount.tsv",header = T)
FC.LT_D12_nt_vs_e13_deseq.res <- read.table("./condition_control_treated.deseq2.results.tsv_featurecount.tsv",header=T)
##### marking colnames with id before merging 
colnames(FC.LT_D12_nt_vs_e13_deseq.norm) <- paste(colnames(FC.LT_D12_nt_vs_e13_deseq.norm),"FC_LT_nt_vs_e13_deseq.norm",sep="-")
colnames(FC.LT_D12_nt_vs_e13_deseq.res) <- paste(colnames(FC.LT_D12_nt_vs_e13_deseq.res),"FC_LT_nt_vs_e13_deseq",sep="-")
colnames(FC.LT_D12_nt_vs_e13_deseq.res)[1] <- "gene_id"
colnames(FC.LT_D12_nt_vs_e13_deseq.norm)[1] <- "gene_id"
#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames
FCLT_d12_nt_vs_e13.merged_df.norm.res <- merge(FC.LT_D12_nt_vs_e13_deseq.res, FC.LT_D12_nt_vs_e13_deseq.norm, by = "gene_id")
FCLT_d12_nt_vs_e13.merged_df.norm.res.count <- merge(FCLT_d12_nt_vs_e13.merged_df.norm.res,FC_count,by = "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(FCLT_d12_nt_vs_e13.merged_df.norm.res.count)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


FCLT_d12_nt_vs_e13.merged_df.norm.res.count <- na.omit(FCLT_d12_nt_vs_e13.merged_df.norm.res.count)

##### mean of Deseq normalized for nt
colnames(FCLT_d12_nt_vs_e13.merged_df.norm.res.count)
FCLT_d12_nt_vs_e13.merged_df.norm.res.count$mean.nt.deseq <- rowMeans(FCLT_d12_nt_vs_e13.merged_df.norm.res.count[, 8:9], na.rm = TRUE)



library(ggrepel)
library(dplyr)

##### basemean cutoff
FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean <- FCLT_d12_nt_vs_e13.merged_df.norm.res.count %>% filter(mean.nt.deseq >= 1)

###########
count_below_001 <- sum(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$`padj-FC_LT_nt_vs_e13_deseq` <= 0.01)

min(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 8.084946
max(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 101107.7
############## labelling 
# Create an empty vector of labels
### add index 
#FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean<- cbind(index = 1:13853, FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean)
#rownames(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean) <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id


##### writing out the files
fc.d12.nt_Vs_e13_LT.upregulated_genes <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange > 0 & FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]
fc.d12.nt_Vs_e13_LT.downregulated_genes <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange < 0 & FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj < 0.01, ]

write.csv(fc.d12.nt_Vs_e13_LT.upregulated_genes,file = "fc.d12.nt_Vs_e13_LT.upregulated_genes.csv")
write.csv(fc.d12.nt_Vs_e13_LT.downregulated_genes,file = "fc.d12.nt_Vs_e13_LT.downregulated_genes")


library(ggplot2)
library(ggrepel)
library(scales)

# Your data
fc.xval <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq
fc.yval <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$log2FoldChange
fc.padj <- FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$padj

# Create a data frame
fc.df <- data.frame(x = fc.xval, y = fc.yval, padj = fc.padj)

# MA plot
# Calculate counts
fc.upregulated_count <- sum(fc.df$y > 0 & fc.df$padj < 0.01)
fc.downregulated_count <- sum(fc.df$y < 0 & fc.df$padj < 0.01)
fc.nonsignificant_count <- sum(fc.df$padj >= 0.01)

fc_gg <- ggplot(fc.df, aes(x = x, y = y, color = factor(ifelse(padj < 0.01, 
                                                               ifelse(y > 0, "Upregulated", "Downregulated"),
                                                               "Nonsignificant")))) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "Featurecount",
       subtitle = "Human iNeurons Day 12 (Wildtype vs. ZSWIM8-KO)") +
  geom_text_repel(data = data.frame(x = fc.xval[which(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    y = fc.yval[which(FCLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    label = "ZSWIM8"),
                  aes(x = x, y = y, label = label),
                  color = "black",
                  nudge_x = 0.5,  # changing ZSWIM8 label position
                  nudge_y = -0.5,
                  size = 4) +
  scale_x_log10(labels = scales::trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm"),  # Adjust tick length here
        legend.position = "top",
        legend.justification = c(1, 1),
        legend.box.just = "right")

fc_gg <- fc_gg +
  geom_text(data = data.frame(x = 400000, y = 4.5, label = paste("Upregulated:", fc.upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 4, label = paste("Downregulated:", fc.downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 3.5, label = paste("Nonsignificant:", fc.nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 1)

# Print the plot
print(fc_gg)

maplot_list <- list(str.gg,sal_gg,fc_gg)
ggarrange(plotlist = maplot_list, ncol = 3, nrow = 1)

############################################################# comparing the counts between str,sal and FC ---- part
### stringties counts
stringtie_count <- read.csv("./stringtie_gene_count_matrix.csv",header = T) 
colnames(stringtie_count) <- paste("str",colnames(stringtie_count),sep=".")
colnames(stringtie_count)[1] <- "gene.id"

### salmons counts
salmon_count <- read.csv("./3_rounded_salmon.merged.gene_counts.csv",header = T) 
colnames(salmon_count) <- paste("sal",colnames(salmon_count),sep=".")
colnames(salmon_count)[1] <- "gene.id"

###FC
FC_count <- read.csv("./Feature_gene_counts.csv",header = T) 
colnames(FC_count) <- paste("FC",colnames(FC_count),sep=".")
colnames(FC_count)[1] <- "gene.id"

###

str-62703
sal-62703
FC -63086 
### genes only in 
FC <- FC_count$gene_id
sal <- salmon_count$gene.id
only_in_sal <- setdiff(sal,FC) #### 2614
only_in_FC <- setdiff(FC,sal)#### 2997
common_in_FC_and_salmon <- intersect(stringtie_count$gene_id,FC_count$gene_id)#### 60089

####### coloring the only_in_sal in sal MA plot 

sal.LT_D12_nt_vs_e13_deseq.norm <- read.table("./salmonall.normalised_counts.tsv",header = T)
sal.LT_D12_nt_vs_e13_deseq.res <- read.table("./salmoncondition_control_treated.deseq2.results.tsv",header=T)
##### marking colnames with id before merging 
colnames(sal.LT_D12_nt_vs_e13_deseq.norm) <- paste(colnames(sal.LT_D12_nt_vs_e13_deseq.norm),"sal_LT_nt_vs_e13_deseq.norm",sep="-")
colnames(sal.LT_D12_nt_vs_e13_deseq.res) <- paste(colnames(sal.LT_D12_nt_vs_e13_deseq.res),"sal_LT_nt_vs_e13_deseq",sep="-")
colnames(sal.LT_D12_nt_vs_e13_deseq.res)[1] <- "gene_id"
colnames(sal.LT_D12_nt_vs_e13_deseq.norm)[1] <- "gene_id"
#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames
salLT_d12_nt_vs_e13.merged_df <- merge(sal.LT_D12_nt_vs_e13_deseq.res, sal.LT_D12_nt_vs_e13_deseq.norm, by = "gene_id")
salLT_d12_nt_vs_e13.merged_df_deseq_count <- merge(salLT_d12_nt_vs_e13.merged_df,salmon_count,by = "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(salLT_d12_nt_vs_e13.merged_df_deseq_count)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


salLT_d12_nt_vs_e13.merged_df_deseq_count <- na.omit(salLT_d12_nt_vs_e13.merged_df_deseq_count)

##### mean of Deseq normalized for nt
colnames(salLT_d12_nt_vs_e13.merged_df_deseq_count)
salLT_d12_nt_vs_e13.merged_df_deseq_count$mean.nt.deseq <- rowMeans(salLT_d12_nt_vs_e13.merged_df_deseq_count[, 8:9], na.rm = TRUE)



library(ggrepel)
library(dplyr)

##### basemean cutoff
salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean <- salLT_d12_nt_vs_e13.merged_df_deseq_count %>% filter(mean.nt.deseq >= 10)
write.csv(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean,file="salmon_DEG_filter_basemean_greater_10.csv")

###########
count_below_001 <- sum(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$`padj-sal_LT_nt_vs_e13_deseq` <= 0.01)

min(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 10.00204
max(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$mean.nt.deseq) ## 108043.6
############## labelling 
# Create an empty vector of labels
### add index 
#salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean<- cbind(index = 1:13853, salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean)
#rownames(salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean) <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id

####### intersect only_in_sal with salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean

only_in_sal_in_sal_MA <- intersect(only_in_sal,salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id)
only_in_sal_in_sal_MA_df <- salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean[salLT_d12_nt_vs_e13.mer.df_deseq_count.basemean$gene_id %in% only_in_sal_in_sal_MA, ]

library(ggplot2)
library(ggrepel)
library(scales)

library(ggplot2)
library(ggrepel)
library(scales)

# Example vector of gene IDs to highlight
only <-  # Replace with your actual gene IDs
  
  # Create a new column 'highlight' that marks points to be highlighted
  sal.df$highlight <- ifelse(sal.df$sal.gene_id %in% only_in_sal, "Highlight", 
                             ifelse(sal.df$padj < 0.01, 
                                    ifelse(sal.df$y > 0, "Upregulated", "Downregulated"), 
                                    "Nonsignificant"))



# Calculate the counts for the annotations
sal.upregulated_count <- sum(sal.df$highlight == "Upregulated")
sal.downregulated_count <- sum(sal.df$highlight == "Downregulated")
sal.nonsignificant_count <- sum(sal.df$highlight == "Nonsignificant")
only_in_salmon_counts <- sum(sal.df$highlight == "Highlight")

# Modify the plot
sal_gg <- ggplot(sal.df, aes(x = x, y = y, color = highlight)) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray",
                                "Highlight" = "blue"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "salmon counts",
       subtitle = "Human iNeurons Day 12 (Wildtype vs. ZSWIM8-KO)") +
  scale_x_log10(labels = scales::trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm"),  # Adjust tick length here
        legend.position = "top",
        legend.justification = c(1, 1),
        legend.box.just = "right")

sal_gg <- sal_gg +
  geom_text(data = data.frame(x = 400000, y = 4.5, label = paste("Upregulated:", sal.upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 4, label = paste("Downregulated:", sal.downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 3.5, label = paste("Nonsignificant:", sal.nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = 400000, y = 3, label = paste("Only in salmon counts:", only_in_salmon_counts)),
            aes(x = x, y = y, label = label), color = "blue", hjust = 1, vjust = 1)

# Print the plot
print(sal_gg)

























###
count <- read.csv("./gene_count_matrix.csv",header = T)
colnames(count) <- paste(colnames(count),"count",sep = "_")
colnames(count)[1] <- "gene_id"
###
TPM <- read.csv("./Final.bulk.bc_v5z8_mchyTPM.csv",header = T,row.names = 1)
colnames(TPM)[1] <-"gene_id"


#### deseq normalized -----mcherry treated and V5z8 as control ----------------------------- first 
## 
chy.V5.deseq <- read.table("./1_chy_vs_v5z8_deseq2.with_shrink_true.tsv",header = T)
chy.V5.deseq.norm_shrink <- read.table("./1.all.chy_vs_v5z8_norm_chy.treated_vs_v5z8.controlwith_shrink_counts.tsv",header=T)
colnames(chy.V5.deseq.norm_shrink) <- paste(colnames(chy.V5.deseq.norm_shrink),"Dseq.norm.shrinked",sep="-")
colnames(chy.V5.deseq.norm_shrink)[1] <- "gene_id"

#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames
chy.V5.merged_df <- merge(count, TPM, by = "gene_id")
chy.V5.merged_df <- merge(chy.V5.merged_df, chy.V5.deseq, by = "gene_id")
chy.V5.fi.merge_df.with.na <- merge(chy.V5.merged_df, chy.V5.deseq.norm_shrink, by= "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(chy.V5.fi.merge_df.with.na)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


chy.V5.fi.merge_df <- na.omit(chy.V5.fi.merge_df.with.na)

##### mean of Deseq normalized for V5z8 
chy.V5.fi.merge_df$mean.V5z8.deseq <- rowMeans(chy.V5.fi.merge_df[, 28:30], na.rm = TRUE)


#/ then install the package from Github:
#install.packages("remotes",force=T)
#remotes::install_github("ATpoint/vizzy")
library("vizzy")
library(ggrepel)


##### basemean cutoff
chy.V5.fi.basemean<- chy.V5.fi.merge_df %>%
  filter(mean.V5z8.deseq >= 100)

###########
count_below_005 <- sum(chy.V5.fi.basemean$padj <= 0.01)

min(chy.V5.fi.basemean$mean.V5z8.deseq) ## 100
max(chy.V5.fi.basemean$mean.V5z8.deseq) ## 470534.3
############## labelling 
# Create an empty vector of labels
### add index 
chy.V5.fi.basemean <- cbind(index = 1:11487, chy.V5.fi.basemean)
rownames(chy.V5.fi.basemean) <- chy.V5.fi.basemean$gene_id


library(ggplot2)

# Your data
xval <- chy.V5.fi.basemean$mean.V5z8.deseq
yval <- chy.V5.fi.basemean$log2FoldChange
padj <- chy.V5.fi.basemean$padj

# Create a data frame
df <- data.frame(x = xval, y = yval, padj = padj)

# MA plot
# Calculate counts
upregulated_count <- sum(df$y > 0 & df$padj < 0.01)
downregulated_count <- sum(df$y < 0 & df$padj < 0.01)
nonsignificant_count <- sum(df$padj >= 0.01)


# Create the plot
gg <- ggplot(df, aes(x = x, y = y, color = factor(ifelse(padj < 0.01, 
                                                         ifelse(y > 0, "Upregulated", "Downregulated"),
                                                         "Nonsignificant")))) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "MA Plot",
       subtitle = "v5z8.control_vs__shrink_chy.treated_padj_0.01") +
  geom_text_repel(data = data.frame(x = xval[which(chy.V5.fi.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    y = yval[which(chy.V5.fi.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    label = "V5-ZSWIM8"),
                  aes(x = x, y = y, label = label),
                  color = "black",
                  nudge_x = -2,
                  nudge_y = -0.5,
                  size = 4) +
  scale_x_continuous(limits = c(100.0217, 100000), breaks = c(100, 25000, 50000, 75000, 100000)) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal()

# Add custom legend labels
gg <- gg +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste("Upregulated:", upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = Inf, y = -Inf, label = paste("Downregulated:", downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 0) +
  geom_text(data = data.frame(x = Inf, y = -5, label = paste("Nonsignificant:", nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 0.5) +
  theme(legend.position = "bottom")  # Remove default legend

# Print the plot
print(gg)





#### deseq normalized -----mcherry treated and V5z8 as control-----no shrink ----------------------------- first 
## 

## 
no_shrink_chy.V5.deseq <- read.table("./2_condition_control_treated.deseq2.false.shrink.tsv",header = T)
no_shrink_chy.V5.deseq.norm_shrink <- read.table("./2_all.norm_chy.treated_v5z8.treated_no_shrink__counts.tsv",header=T)
colnames(no_shrink_chy.V5.deseq.norm_shrink) <- paste(colnames(no_shrink_chy.V5.deseq.norm_shrink),"Dseq.norm.shrinked",sep="-")
colnames(no_shrink_chy.V5.deseq.norm_shrink)[1] <- "gene_id"

#### merge dataframe 
# Assuming count, TPM, Deseq.normalized, and deseq are your data frames
# Merge count and TPM data frames
no_shrink_chy.V5.merged_df <- merge(count, TPM, by = "gene_id")
no_shrink_chy.V5.merged_df <- merge(no_shrink_chy.V5.merged_df, no_shrink_chy.V5.deseq, by = "gene_id")
no_shrink_chy.V5.fi.merge_df.with.na <- merge(no_shrink_chy.V5.merged_df, no_shrink_chy.V5.deseq.norm_shrink, by= "gene_id")
## Assuming merged_df is your dataframe

# Check if any NA values exist in the dataframe
if (anyNA(no_shrink_chy.V5.fi.merge_df.with.na)) {
  result <- "yes"
} else {
  result <- "no"
}

# Print the result
print(result)


no_shrink_chy.V5.fi.merge_df <- na.omit(no_shrink_chy.V5.fi.merge_df.with.na)

##### mean of Deseq normalized for V5z8 
no_shrink_chy.V5.fi.merge_df$mean.V5z8.deseq <- rowMeans(no_shrink_chy.V5.fi.merge_df[, 28:30], na.rm = TRUE)


#/ then install the package from Github:
#install.packages("remotes",force=T)
#remotes::install_github("ATpoint/vizzy")
library("vizzy")
library(ggrepel)


##### basemean cutoff
no_shrink_chy.V5.fi.basemean<- no_shrink_chy.V5.fi.merge_df %>%
  filter(mean.V5z8.deseq >= 100)

###########
count_below_001 <- sum(no_shrink_chy.V5.fi.basemean$padj <= 0.01) # 924

min(no_shrink_chy.V5.fi.basemean$mean.V5z8.deseq) ## 100.0305
max(no_shrink_chy.V5.fi.basemean$mean.V5z8.deseq) ## 473799.3
############## labelling 
# Create an empty vector of labels
### add index 
no_shrink_chy.V5.fi.basemean <- cbind(index = 1:11487, no_shrink_chy.V5.fi.basemean)
rownames(no_shrink_chy.V5.fi.basemean) <- no_shrink_chy.V5.fi.basemean$gene_id




library(ggplot2)

# Your data
xval <- no_shrink_chy.V5.fi.basemean$mean.V5z8.deseq
yval <- no_shrink_chy.V5.fi.basemean$log2FoldChange
padj <- no_shrink_chy.V5.fi.basemean$padj

# Create a data frame
df <- data.frame(x = xval, y = yval, padj = padj)

# MA plot
# Calculate counts
upregulated_count <- sum(df$y > 0 & df$padj < 0.01)
downregulated_count <- sum(df$y < 0 & df$padj < 0.01)
nonsignificant_count <- sum(df$padj >= 0.01)

# Create the plot

library(ggplot2)



# Create the plot
gg <- ggplot(df, aes(x = x, y = y, color = factor(ifelse(padj < 0.01, 
                                                         ifelse(y > 0, "Upregulated", "Downregulated"),
                                                         "Nonsignificant")))) +
  geom_point(alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Upregulated" = "cyan", 
                                "Downregulated" = "deeppink2",
                                "Nonsignificant" = "darkgray"),
                     guide = FALSE) +
  labs(x = "WT normalized counts", y = "Fold change KO/WT (log2)",
       title = "MA Plot",
       subtitle = "v5z8.control_vs_mno_shrink_chy.treated_padj_0.01") +
  geom_text_repel(data = data.frame(x = xval[which(no_shrink_chy.V5.fi.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    y = yval[which(no_shrink_chy.V5.fi.basemean$gene_id == "ENSG00000214655.11|ZSWIM8")],
                                    label = "V5-ZSWIM8"),
                  aes(x = x, y = y, label = label),
                  color = "black",
                  nudge_x = -2,
                  nudge_y = -0.5,
                  size = 4) +
  scale_x_continuous(limits = c(100.0217, 100000), breaks = c(100, 25000, 50000, 75000, 100000)) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_minimal()

# Add custom legend labels
gg <- gg +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste("Upregulated:", upregulated_count)),
            aes(x = x, y = y, label = label), color = "cyan", hjust = 1, vjust = 1) +
  geom_text(data = data.frame(x = Inf, y = -Inf, label = paste("Downregulated:", downregulated_count)),
            aes(x = x, y = y, label = label), color = "deeppink2", hjust = 1, vjust = 0) +
  geom_text(data = data.frame(x = Inf, y = -5, label = paste("Nonsignificant:", nonsignificant_count)),
            aes(x = x, y = y, label = label), color = "darkgray", hjust = 1, vjust = 0.5) +
  theme(legend.position = "bottom")  # Remove default legend

# Print the plot
print(gg)
























