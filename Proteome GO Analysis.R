library(openxlsx)

Low_gene <- read.xlsx("/Users/mortezaabyadeh/Desktop/Cleaned_Proteome.xlsx", sheet = "Low Abundance")
print("Here is the head of the data:")
head(Low_gene)
dim(Low_gene)
High_gene <- read.xlsx("/Users/mortezaabyadeh/Desktop/Cleaned_Proteome.xlsx", sheet = "High Abundance")
head(High_gene)
dim(High_gene)


low_abundance_genes <- Low_gene[["Gene.Symbol"]]
high_abundance_genes <- High_gene[["Gene.Symbol"]]
head(low_abundance_genes)
low = as.data.frame(low_abundance_genes)
high = as.data.frame(high_abundance_genes)
dim(low)

gene_names_string <- paste(low, collapse = ", ")
head(gene_names_string,5)
gene_names_string <- gsub("\\\"", "", gene_names_string)
print(gene_names_string)

library(devtools)
install_github("wjawaid/enrichR")
install.packages("enrichR")
library(enrichR)
listEnrichrSites()

setEnrichrSite("WormEnrichr")

dbs <- listEnrichrDbs()
head(dbs, 30)
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019")
enriched <- enrichr(gene_names_string, dbs)
dim(enriched[["GO_Biological_Process_2018"]])
dim(enriched[["KEGG_2019"]])
plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
enriched[["KEGG_2019"]]

kegg_results <- enriched[["KEGG_2019"]]
sorted_kegg_results <- kegg_results[order(kegg_results$P.value), ]
kegg_results