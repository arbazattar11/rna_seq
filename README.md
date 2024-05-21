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

1. **Clone the Repository**: Clone this repository to your local machine using `git clone https://github.com/username/repository.git`.
2. **Prepare Data**: Prepare your RNA-seq expression data in tabular format.
3. **Run the Script**: Update the file path in the script to point to your data file and execute `python rna_seq_analysis.py`.
4. **Explore Results**: Explore the generated plots to gain insights into the structure of your RNA-seq data.

## Dependencies

- Java
- R
- Fastqc
- Fastp
- Hisat2
- samtools
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
# Fastp is not available via apt, you need to download and install it manually or via bioconda
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

**Usage - build a HISAT2 index:**
```bash
hisat2-build genome.fa genome_index
```
**Usage - alignment and mapping using HISAT2:**
```bash
hisat2 -x index_prefix -1 read1.fastq.gz -2 read2.fastq.gz -S output.sam
```
### Samtools & Stringtie
**Usage - Binary Alignment / Map format using Samtools & Stringtie:**
```bash
samtools view -Sb input.sam > output.bam
```
```bash
samtools sort input.bam -o output.bam
```
```bash
samtools index input.bam
``````
```bash
stringtie --merge -G reference.gtf -o merged.gtf input_list.txt
```
```bash
stringtie input.bam -G reference.gtf -o output.gtf
```
```bash
stringtie input.bam -G reference.gtf -o output.gtf
```
## Contributing

Contributions to this repository are welcome! If you have suggestions for improvements, bug fixes, or additional features, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
