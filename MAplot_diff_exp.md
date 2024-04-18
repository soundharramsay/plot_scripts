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
























