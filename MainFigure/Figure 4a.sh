#!/bin/bash
# figure 4a, bw file from 'step 5 bw' output file of /Users/carolwei/Downloads/github-plot/hnRNPL-Chipseq-alignment.sh 
# reigion file calculate by deeptools using gencode_grch35 genome

genome='/home/miniconda3/envs/rna_2.7/share/homer/data/genomes/gencode_grch35'
bedtools='/home/software/bedtools2/bin'
core=30
mkdir ./profile_rna
cd ./profile_rna
computeMatrix scale-regions -p $core \
  -R $genome/gencode_grch35.Rna.for.deeptools \
  -S ../*.bw \
  -b 2000 -a 2000 \
  --regionBodyLength 5000 \
  --skipZeros -o bw_morm.gz \
  --outFileNameMatrix bw_morm.tab \
  --outFileSortedRegions bw_morm.bed
plotProfile --dpi 300 -m bw_morm.gz \
              -out profile.pdf \
              --numPlotsPerRow 2 \
              --plotTitle "ChIP-seq Profile"