library(readxl)
library(limma)
library(openxlsx)  
library(stringi)


setwd("/Users/mortezaabyadeh/Documents/npc-fty/senescence omics data/proteome data analysis new")


data <- read.xlsx("Cleaned_Proteomic_data.xlsx")
abundance_data <- data[, 2:13]  
rownames(abundance_data) <- make.unique(data[[1]])  


#abundance_data_log <- log2(abundance_data + 1)

abundance_data_log <- abundance_data

condition <- factor(rep(c("CTRL", "Eto", "FTY", "FTY+Eto"), each = 3))
design <- model.matrix(~ 0 + condition)  


colnames(design) <- levels(condition)


fit <- lmFit(abundance_data_log, design)
fit <- eBayes(fit)


condition_pairs <- list(
  c("CTRL", "Eto"),
  c("CTRL", "FTY"),
  c("CTRL", "FTY+Eto"),
  c("Eto", "FTY"),
  c("Eto", "FTY+Eto"),
  c("FTY", "FTY+Eto")
)


wb <- createWorkbook()

for (pair in condition_pairs) {
  cond1 <- pair[1]
  cond2 <- pair[2]
  cond1_indices <- which(condition == cond1)
  cond2_indices <- which(condition == cond2)
  
  cond1_data <- abundance_data_log[, cond1_indices]
  cond2_data <- abundance_data_log[, cond2_indices]
  
  fold_change <- rowMeans(cond2_data) / rowMeans(cond1_data)
  log2FC <- log2(fold_change)
  
  p_values <- apply(abundance_data_log, 1, function(x) t.test(x[cond1_indices], x[cond2_indices])$p.value)
  
  adj_p_values <- p.adjust(p_values, method = "BH")
  
  result_df <- data.frame(
    GeneName = rownames(abundance_data_log),
    Cond1_Mean = rowMeans(cond1_data),
    Cond2_Mean = rowMeans(cond2_data),
    log2FC = log2FC,
    P_Value = p_values,
    Adjusted_P_Value = adj_p_values,
    cond1_data,  
    cond2_data   
  )
  

  addWorksheet(wb, paste(cond1, "vs", cond2))
  writeData(wb, sheet = paste(cond1, "vs", cond2), result_df)
}


saveWorkbook(wb, "Differential_Expression_Results.xlsx", overwrite = TRUE)
