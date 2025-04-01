library(readxl)
library(dplyr)
library(data.table)
library(tidyverse)
library(reshape2)
library(writexl)
library(multcomp)
install.packages("multcomp")
library(stats)
setwd("/Users/mortezaabyadeh/Documents/npc-fty/senescence omics data/proteome data analysis new/Anova")
data <- read_excel("Cleaned_Proteomic_data.xlsx")
head(data)
library(dplyr)
library(tidyr)
library(multcomp)
library(stats)
library(openxlsx) 




library(dplyr)
library(tidyr)
library(multcomp)
library(stats)
library(openxlsx) 

group <- factor(rep(c("CTRL", "Eto", "FTY", "FTY_Eto"), each = 3))

comparisons <- list(
  "Eto-CTRL" = c("Eto", "CTRL"),
  "FTY-CTRL" = c("FTY", "CTRL"),
  "FTY-Eto" = c("FTY", "Eto"),
  "FTY_Eto-Eto" = c("FTY_Eto", "Eto"),
  "FTY_Eto-CTRL" = c("FTY_Eto", "CTRL"),
  "FTY_Eto-FTY" = c("FTY_Eto", "FTY")
)

perform_tests <- function(expression_values) {
  df <- data.frame(
    Expression = as.numeric(expression_values),
    Group = group
  )
  
  anova_result <- aov(Expression ~ Group, data = df)
  anova_pvalue <- summary(anova_result)[[1]][["Pr(>F)"]][1]
  
  tukey_result <- TukeyHSD(anova_result)
  tukey_pvalues <- as.data.frame(tukey_result$Group)
  
  tukey_pvalues <- tibble::rownames_to_column(tukey_pvalues, "Comparison")
  
  colnames(tukey_pvalues)[colnames(tukey_pvalues) == "p adj"] <- "p_value"
  
  tukey_pvalues$Comparison <- gsub(" ", "", tukey_pvalues$Comparison) # Remove spaces
  
  tukey_results <- tukey_pvalues %>% 
    filter(Comparison %in% names(comparisons))
  
  for (comp_name in names(comparisons)) {
    if (!comp_name %in% tukey_results$Comparison) {
      tukey_results <- bind_rows(tukey_results, data.frame(Comparison = comp_name, p_value = NA))
    }
  }
  

  tukey_results$Gene <- data$Gene[i]
  
  return(list(anova_pvalue, tukey_results))
}


anova_results <- list()
tukey_results_list <- vector("list", length(comparisons))


for (i in 1:length(comparisons)) {
  tukey_results_list[[i]] <- data.frame(Gene = character(), p_value = numeric(), adj_p_value = numeric())
}


for (i in 1:nrow(data)) {
  result <- perform_tests(data[i, -1])
  

  anova_results[[i]] <- data.frame(
    Gene = data$Gene[i],
    p_value = result[[1]]
  )
  
 
  tukey_df <- result[[2]]
  
  for (j in 1:length(comparisons)) {
    comp_name <- names(comparisons)[j]
    subset_df <- tukey_df %>% filter(Comparison == comp_name)
    
    if (nrow(subset_df) > 0) {
      tukey_results_list[[j]] <- bind_rows(tukey_results_list[[j]], subset_df)
    }
  }
}


anova_results_df <- bind_rows(anova_results)


for (j in 1:length(comparisons)) {
  tukey_results_list[[j]] <- tukey_results_list[[j]] %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
}


file_name <- "anova_tukey_results.xlsx"
wb <- createWorkbook()


addWorksheet(wb, "ANOVA Results")
writeData(wb, "ANOVA Results", anova_results_df)


for (j in 1:length(comparisons)) {
  sheet_name <- names(comparisons)[j]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tukey_results_list[[j]])
}

saveWorkbook(wb, file_name, overwrite = TRUE)

print(paste("Results saved to", file_name))

