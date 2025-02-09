#!/bin/bash
# hnRNPL ChIP-seq alignment pipline

########## step 1 fastqc
fastqc -o ./ -d ./ -f fastq -t 30 $sample_1.clean.fq.gz $$sample_2.clean.fq.gz 

########## step 2 cutadapter
cutadapt -f fastq --match-read-wildcards --times 2 -e 0.1 -O 5 -m 35 \
-a CTGTCTCTTATACAC -a TGTCTCTTATACACA -a GTCTCTTATACACAT -a TCTCTTATACACATC -a CTCTTATACACATCT -a TCTTATACACATCTC -a CTTATACACATCTCC -a TTATACACATCTCCG -a TATACACATCTCCGA -a ATACACATCTCCGAG -a TACACATCTCCGAGC -a ACACATCTCCGAGCC -a CACATCTCCGAGCCC -a ACATCTCCGAGCCCA -a CATCTCCGAGCCCAC -a ATCTCCGAGCCCACG -a TCTCCGAGCCCACGA -a CTCCGAGCCCACGAG -a TCCGAGCCCACGAGA -a CCGAGCCCACGAGAC \
-A CTGTCTCTTATACAC -A TGTCTCTTATACACA -A GTCTCTTATACACAT -A TCTCTTATACACATC -A CTCTTATACACATCT -A TCTTATACACATCTG -A CTTATACACATCTGA -A TTATACACATCTGAC -A TATACACATCTGACG -A ATACACATCTGACGC -A TACACATCTGACGCT -A ACACATCTGACGCTG -A CACATCTGACGCTGC -A ACATCTGACGCTGCC -A CATCTGACGCTGCCG -A ATCTGACGCTGCCGA -A TCTGACGCTGCCGAC -A CTGACGCTGCCGACG -A TGACGCTGCCGACGA \
-o ./$sample.cutadapt_1.fq.gz -p ./$sample.cutadapt_2.fq.gz ./$sample_1.clean.fq.gz ./$sample_2.clean.fq.gz >./$sample_cutadaptor.metrics

cutadapt -f fastq --match-read-wildcards --times 2 -e 0.1 -O 5 -m 35 \
-g GTGTATAAGAGACAG -g TGTGTATAAGAGACA -g ATGTGTATAAGAGAC -g GATGTGTATAAGAGA -g AGATGTGTATAAGAG -g CAGATGTGTATAAGA -g TCAGATGTGTATAAG -g GTCAGATGTGTATAA -g CGTCAGATGTGTATA -g GCGTCAGATGTGTAT -g AGCGTCAGATGTGTA -g CAGCGTCAGATGTGT -g GCAGCGTCAGATGTG -g GGCAGCGTCAGATGT -g CGGCAGCGTCAGATG -g TCGGCAGCGTCAGAT -g GTCGGCAGCGTCAGA -g CGTCGGCAGCGTCAG -g TCGTCGGCAGCGTCA \
-G GTGTATAAGAGACAG -G TGTGTATAAGAGACA -G ATGTGTATAAGAGAC -G GATGTGTATAAGAGA -G AGATGTGTATAAGAG -G GAGATGTGTATAAGA -G GGAGATGTGTATAAG -G CGGAGATGTGTATAA -G TCGGAGATGTGTATA -G CTCGGAGATGTGTAT -G GCTCGGAGATGTGTA -G GGCTCGGAGATGTGT -G GGGCTCGGAGATGTG -G TGGGCTCGGAGATGT -G GTGGGCTCGGAGATG -G CGTGGGCTCGGAGAT -G TCGTGGGCTCGGAGA -G CTCGTGGGCTCGGAG -G TCTCGTGGGCTCGGA -G GTCTCGTGGGCTCGG \
-o ./$sample.cutadapt_R2_1.fq.gz -p ./$sample.cutadapt_R2_2.fq.gz ./$sample.cutadapt_1.fq.gz ./$sample.cutadapt_2.fq.gz >./$sample_cutadaptor.metrics
done


########## step 3 bowtie2
genome='/home/genome/grch38/GRCh38.bowtie2'
path="/home/software/samtools-1.11/bin"
core=30

bowtie2 -x $genome/GRCh38.primary_assembly.genome \
--un-conc ./$sample.un.fastq \
-1 $sample.cutadapt_R2_1.fq.gz -2 $sample.cutadapt_R2_2.fq.gz \
-p 20 --sensitive-local -N 1 | $path/samtools sort -@ $core -O bam -o $sample.bam

$path/samtools view -@ $core -hf 0x2 $sample.bam | grep -v "XS:i:" > $sample.uni.sam
$path/samtools view -@ $core -b $sample.uni.sam >$sample.uni.bam
$path/samtools flagstat -@ $core $sample.uni.bam >$sample.uni.txt
$path/samtools sort -@ $core $sample.uni.bam >$sample.uni.srt.bam

########## step 4 picard
java -jar /home/software/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=sample.uni.srt.bam OUTPUT=sample.picard.bam METRICS_FILE=sample.metrics
samtools sort -@ 20 sample.picard.bam >sample.srt.bam
samtools index -@ 15 sample.srt.bam
samtools flagstat -@ 15 sample.picard.bam >sample.picard.txt

########## step 5 macs2
path='/home/software/MACS-2.2.7.1/bin'
$path/macs2 callpeak --call-summits -q 0.1 \
-t IP1.srt.bam IP2.srt.bam \
-c IgG.srt.bam \
--outdir ./ -g hs -B -f BAMPE -n pool

$path/macs2 callpeak --broad \
-t IP1.srt.bam IP2.srt.bam \
--outdir ./ -g hs -B -f BAMPE -n pool


########## step 5 bw
core=30
bamCoverage -p $core --scaleFactor 12.54 -e 300 -b IgG.srt.bam -o IgG.bw
bamCoverage -p $core --scaleFactor 3.25 -e 300 -b IP1.srt.bam -o IP1.srt.bw
bamCoverage -p $core --scaleFactor 4.09 -e 300 -b IP2.srt.bam -o IP2.srt.bw

