sudo apt update
sudo apt install fastqc
fastqc *.fastq.gz
# Fastp is available via apt or you need to download and install it manually or via bioconda
sudo apt install fastp
fastp -i input.fastq.gz -o output.fastq.gz
sudo apt install Hisat2
hisat2-build genome.fa genome_index
hisat2 -x index_prefix -1 read1.fastq.gz -2 read2.fastq.gz -S output.sam
sudo apt install Samtools
sudo apt install Stringtie
samtools view -Sb input.sam > output.bam
samtools sort input.bam -o output.bam
samtools index input.bam
stringtie -G reference.gtf -o output.gtf input.bam
stringtie --merge -G reference.gtf -o merged.gtf input_list.txt
stringtie input.bam -G reference.gtf -o output.gtf
stringtie [options] -G reference.gtf -o output.gtf input.bam
