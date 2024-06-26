# RNA-seq Analysis Pipeline

Welcome to the RNA-seq Analysis Pipeline repository! This repository contains a linux script for performing RNA-seq data analysis. RNA sequencing (RNA-seq) is a powerful technique used to measure the expression of genes in a biological sample. This pipeline takes RNA-seq expression data as input, preprocesses it and visualizes the results.

## Introduction

RNA-seq has become a standard tool for studying gene expression due to its high-throughput nature and ability to quantify transcript levels accurately. However, analyzing RNA-seq data can be complex and requires a series of computational steps to extract meaningful insights.

This pipeline aims to simplify the RNA-seq analysis process by providing a modular and user-friendly linux script. Users can easily adapt this script to their specific datasets and research questions.

## Features

- **Data Preprocessing**: Normalize and log-transform the expression data to prepare it for analysis.
- **Customizable**: Easily modify the script to incorporate additional analysis steps or adapt it to different experimental designs.
- **Visualization**: Generate plots to visualize the results and explore the underlying structure of the data.

## Usage

1. **Clone the Repository**: Clone this repository to your local machine using `git clone https://github.com/arbazattar11/rna_seq.git`.
2. **Prepare Data**: Prepare your RNA-seq expression data in tabular format.
3. **Run the Script**: Update the file path in the script to point to your data file and execute.
4. **Explore Results**: Explore the generated plots to gain insights into the structure of your RNA-seq data.

## Dependencies

- Java
- R
- Pandas
- Fastqc
- Fastp
- Hisat2
- Samtools
- Stringtie
- Deseq2

## Pipeline

### FastQC

**Installation:**
```bash
sudo apt update
sudo apt install fastqc
```

**Usage:**
```bash
fastqc *.fastq.gz
```

### Fastp

**Installation:**
```bash
# Fastp is available via apt or you need to download and install it manually or via bioconda
```
**Installation:**
```bash
sudo apt install fastp
```

**Usage - Single End:**
```bash
fastp -i input.fastq.gz -o output.fastq.gz
```

**Usage - Paired End:**
```bash
fastp -i input_R1.fastq.gz -I input_R2.fastq.gz -o output_R1.fastq.gz -O output_R2.fastq.gz --detect_adapter_for_pe
```

**Some common options and parameters**

- -i: Input fastq file.
- -o: Output fastq file.
- --detect_adapter_for_pe: Automatically detect and remove adapter sequences for paired-end reads.
- --cut_front: Remove a specified number of bases from the 5' end of reads.
- --cut_tail: Remove a specified number of bases from the 3' end of reads.
- --cut_mean_quality: Trim low-quality bases from the ends of reads until the mean quality score meets a specified threshold.
- --length_required: Filter out reads shorter than a specified length.
- --qualified_quality_phred: Specify the minimum Phred quality score required to keep a base.
- --unqualified_percent_limit: Specify the maximum percentage of bases allowed with quality scores below the specified threshold.
- --thread: Specify the number of threads to use for processing (default is 1).

### Hisat2 
**Installation:**
```bash
sudo apt install Hisat2
```
**Usage - build a HISAT2 index:**
```bash
hisat2-build genome.fa genome_index
```
**Usage - alignment and mapping using HISAT2:**
```bash
hisat2 -x index_prefix -1 read1.fastq.gz -2 read2.fastq.gz -S output.sam
```
### Samtools & Stringtie
**Installation:**
```bash
sudo apt install Samtools
```
```bash
sudo apt install Stringtie
```
**Usage - Binary Alignment / Map format using Samtools & Stringtie:**
- Convert sam file to bam format
```bash
samtools view -Sb input.sam > output.bam
```
- Sorting
```bash
samtools sort input.bam -o output.bam
```
- Indexing
```bash
samtools index input.bam
```
- Merged files
```bash
stringtie --merge -G reference.gtf -o merged.gtf input_list.txt
```
- Creating gtf filename
```bash
stringtie input.bam -G reference.gtf -o output.gtf
```
- Transcript assembly

```bash
stringtie [options] -G reference.gtf -o output.gtf input.bam
```
-[options]: Additional options to customize the behavior of StringTie. In your command, -e enables the output of gene-level expression estimates, and -B outputs read coverages for each transcript.

**Usage - For Expreesion Quantification Running the python script:**
```bash
import pandas as pd

def quantify_expression(gtf_file):
    # Load the GTF file into a DataFrame
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#', names=columns)

    # Filter rows to retain only transcripts (feature == 'transcript')
    transcripts = df[df['feature'] == 'transcript']

    # Extract gene_id and transcript_id from attributes column
    transcripts['gene_id'] = transcripts['attributes'].str.extract(r'gene_id "(.*?)";')
    transcripts['transcript_id'] = transcripts['attributes'].str.extract(r'transcript_id "(.*?)";')

    # Group by gene_id and sum up the length of transcripts to get expression
    expression = transcripts.groupby('gene_id')['end'].sum()

    return expression

if __name__ == "__main__":
    # Path to the output GTF file from StringTie
    gtf_file = 'path/to/your/output.gtf'

    # Perform expression quantification
    expression = quantify_expression(gtf_file)

    # Save expression values to a CSV file
    expression.to_csv('expression.csv')
```
### Deseq2
**Usage - Differential gene expression analysis using DESeq2 in R:**
```bash
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
```
**Usage -explore and visualize gene expression patterns and differential expression between conditions.**
- ggplot2
```bash
library(ggplot2)

# Set up multiple plots in a grid
par(mfrow=c(2,3))

# List of genes to plot
genes <- c("NM_024330_1", "NM_002508", "NM_007046", "NM_012249", "NM_014899", "NM_016495")

# Loop through each gene and plot its counts
for (gene in genes) {
    plotCounts(DDS, gene=gene, intgroup="Condition")
}
```
- PlotPCA
```bash
# Perform VST transformation
pcaData <- vst(DDSeq, blind = FALSE)

# Plot PCA
plotPCA(pcaData, intgroup = "Condition")
```
- PlotMA
```bash
# Perform differential expression analysis
Result <- nbinomTest(DDS, "primary", "metastatic")

# Plot MA plot
plotMA(Result)
```
### Further Analysis:
- Identification of Differentially Expressed Genes (DEGs): Based on the plots and any additional statistical tests, identify genes that are significantly differentially expressed between conditions. DESeq2 typically provides statistical tests and adjusted p-values for this purpose.

- Annotation and Functional Analysis: Once DEGs are identified, you might want to annotate them with biological information such as gene ontology (GO) terms, pathways, or functional annotations. Tools like clusterProfiler or enrichR in R can be used for functional enrichment analysis.

- Validation: Validate the differential expression of selected genes using alternative methods such as qPCR or validation datasets if available.

- Biological Interpretation: Interpret the results in the context of your biological question or hypothesis. Consider the biological relevance of identified DEGs and their potential roles in the studied biological process or disease.

- Further Analysis: Depending on your research goals, you may proceed with downstream analyses such as network analysis, pathway analysis, or integration with other omics data.

## Contributing

Contributions to this repository are welcome! If you have suggestions for improvements, bug fixes, or additional features, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
