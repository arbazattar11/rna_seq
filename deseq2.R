# Load DESeq2 library
library(DESeq2)

# Load count data
count_data <- read.table("count_matrix.csv", header=TRUE, row.names=1)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=count_data, colData=col_data, design=~condition)

# Perform normalization and differential expression analysis
dds <- DESeq(dds)

# Extract differential expression results
results <- results(dds)

# Filter results based on adjusted p-value and fold change
significant_results <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

# View top differentially expressed genes
top_genes <- head(significant_results[order(significant_results$padj),], n=10)
print(top_genes)
