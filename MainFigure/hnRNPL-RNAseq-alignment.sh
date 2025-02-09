#!/bin/bash
# hnRNPL RNA-seq alignment pipline

########## step 1 fastqc
fastqc -o ./ -d ./ -f fastq -t 30 $sample_1.clean.fq.gz $$sample_2.clean.fq.gz 


########## step 2 tophat2
genome='/home/yiping/genome/grch38'
tophat2 -p 5 -G /home/genome/grch38/gencode.v35.primary_assembly.annotation.gtf \
-o ./Sample/ \
/home/genome/grch38/GRCh38.bowtie2/GRCh38.primary_assembly.genome \
Sample_R1.fq.gz Sample_R2.fq.gz 
samtools_0.1.18 sort -n accepted_hits.bam Nsort.bam

########## step 3 Count
htseq-count -f bam Nsort.bam.bam $genome/gencode.v35.primary_assembly.annotation.gtf 1>Sample.htseq 2>Sample.log