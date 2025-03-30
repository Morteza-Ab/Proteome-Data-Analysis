library(DESeq2)
library(readxl)
library(dplyr)
library(openxlsx)

setwd("/Users/mortezaabyadeh/Documents/npc-fty/senescence omics data/proteome data analysis new")


counts_df <- read_excel("Cleaned_Proteomic_count.xlsx")
counts_df <- as.data.frame(counts_df)
head(counts_df)


dim(counts_df)
sample_names <- colnames(counts_df)[2:ncol(counts_df)] # Extract sample names from counts_df
colnames(counts_df)
# Create col_data with sample names and group information
col_data <- data.frame(
  sampleName = sample_names,
  group = factor(rep(c("CTRL", "Eto", "FTY", "FTY_Eto"), each = 3))
)

countData <- counts_df[, -1]
rownames(countData) <- make.unique(counts_df[[1]])
rownames(countData)
head(countData)
dim(countData)

str(countData)
countData <- round(countData)

countData <- as.data.frame(lapply(countData, as.integer))
warning()
any(is.na(countData))
sum(is.na(countData))
which(is.na(countData), arr.ind = TRUE)
countData[] <- lapply(countData, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)  # Replace NAs with the median of the column
  return(x)
})



# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = col_data,
                              design = ~ group)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)


results <- results(dds)


results_Eto_vs_CTRL <- results(dds, contrast = c("group", "Eto", "CTRL"))
results_FTY_vs_CTRL <- results(dds, contrast = c("group", "FTY", "CTRL"))
results_FTY_Eto_vs_CTRL <- results(dds, contrast = c("group", "FTY_Eto", "CTRL"))
results_FTY_vs_Eto <- results(dds, contrast = c("group", "FTY", "Eto"))
results_FTY_Eto_vs_Eto <- results(dds, contrast = c("group", "FTY_Eto", "Eto"))
results_FTY_Eto_vs_FTY <- results(dds, contrast = c("group", "FTY_Eto", "FTY"))

results_Eto_vs_CTRL$Gene <- counts_df$`Gene Symbol`
results_FTY_vs_CTRL$Gene <- counts_df$`Gene Symbol`
results_FTY_Eto_vs_CTRL$Gene <- counts_df$`Gene Symbol`
results_FTY_vs_Eto$Gene <- counts_df$`Gene Symbol`
results_FTY_Eto_vs_Eto$Gene <- counts_df$`Gene Symbol`
results_FTY_Eto_vs_FTY$Gene <- counts_df$`Gene Symbol`

# Convert results to data frames
results_Eto_vs_CTRL_df <- as.data.frame(results_Eto_vs_CTRL)
results_FTY_vs_CTRL_df <- as.data.frame(results_FTY_vs_CTRL)
results_FTY_Eto_vs_CTRL_df <- as.data.frame(results_FTY_Eto_vs_CTRL)
results_FTY_vs_Eto_df <- as.data.frame(results_FTY_vs_Eto)
results_FTY_Eto_vs_Eto_df <- as.data.frame(results_FTY_Eto_vs_Eto)
results_FTY_Eto_vs_FTY_df <- as.data.frame(results_FTY_Eto_vs_FTY)



results_Eto_vs_CTRL_df <- cbind(results_Eto_vs_CTRL_df, normalized_counts[, 1:3], normalized_counts[, 4:6])

# For FTY vs CTRL (first 3 columns for CTRL and next 3 for FTY)
results_FTY_vs_CTRL_df <- cbind(results_FTY_vs_CTRL_df, normalized_counts[, 1:3], normalized_counts[, 7:9])

# For FTY_Eto vs CTRL (first 3 columns for CTRL and next 3 for FTY_Eto)
results_FTY_Eto_vs_CTRL_df <- cbind(results_FTY_Eto_vs_CTRL_df, normalized_counts[, 1:3], normalized_counts[, 10:12])

# For FTY vs Eto (first 3 columns for Eto and next 3 for FTY)
results_FTY_vs_Eto_df <- cbind(results_FTY_vs_Eto_df, normalized_counts[, 4:6], normalized_counts[, 7:9])

# For FTY_Eto vs Eto (first 3 columns for Eto and next 3 for FTY_Eto)
results_FTY_Eto_vs_Eto_df <- cbind(results_FTY_Eto_vs_Eto_df, normalized_counts[, 4:6], normalized_counts[, 10:12])

# For FTY_Eto vs FTY (first 3 columns for FTY and next 3 for FTY_Eto)
results_FTY_Eto_vs_FTY_df <- cbind(results_FTY_Eto_vs_FTY_df, normalized_counts[, 7:9], normalized_counts[, 10:12])



head(results_FTY_Eto_vs_FTY_df)

write.xlsx(results_Eto_vs_CTRL_df, file = "results_Eto_vs_CTRL.xlsx")
write.xlsx(results_FTY_vs_CTRL_df, file = "results_FTY_vs_CTRL.xlsx")
write.xlsx(results_FTY_Eto_vs_CTRL_df, file = "results_FTY_Eto_vs_CTRL.xlsx")
write.xlsx(results_FTY_vs_Eto_df, file = "results_FTY_vs_Eto.xlsx")
write.xlsx(results_FTY_Eto_vs_Eto_df, file = "results_FTY_Eto_vs_Eto.xlsx")
write.xlsx(results_FTY_Eto_vs_FTY_df, file = "results_FTY_Eto_vs_FTY.xlsx")




