#!/bin/bash
#### for plot
ls *.bed |while read id
do
sort -k1,1 -k2,2n $id >${id%%.*}.srt.bed
bedtools merge -i ${id%%.*}.srt.bed -c 4 -d 10 -o collapse -delim "|"  >${id%%.*}.merge.bed
sort -k1,1 -k2,2n ${id%%.*}.merge.bed >${id%%.*}.merge.srt.bed
done


bedtools intersect -a NC.merge.srt.bed -b hnRNPL.merge.srt.bed -wa |wc -l 
bedtools intersect -a NC.merge.srt.bed -b hnRNPL.merge.srt.bed -wb |wc l 
bedtools intersect -a NC.merge.srt.bed -b hnRNPL.merge.srt.bed -wa -wb > intersect.bed

###### Profile
out='/data/yiping/hnRNPL/08_20231219_chipseq/results/08_merge/02_Profile/01_Profile_shNC'
dir='/home/yiping/wql/hnRNPL/04_20210809_chipseq/chip-seq/07_bw_IGV'
ref='/data/yiping/hnRNPL/08_20231219_chipseq/results/08_merge/02_Profile/NC_summits.merge.srt.bed'
computeMatrix reference-point -p 40 \
       --referencePoint center \
       -b 5000 -a 5000 \
       -R $ref \
       -S $dir/*bw \
       --skipZeros \
       -o $out/HnRNPL_itChIPoverNC_PolII.gz \
       --outFileSortedRegions $out/HnRNPL_itChIPoverNC_PolII.bed

plotHeatmap --dpi 300 -m $out/HnRNPL_itChIPoverNC_PolII.gz \
     -out $out/HnRNPL_itChIPoverNC_PolII.pdf \
     --colorMap Blues 


###### Profile
out='/data/yiping/hnRNPL/08_20231219_chipseq/results/08_merge/02_Profile/02_PolII_QC'
dir='/data/yiping/hnRNPL/08_20231219_chipseq/results/04_bw/IGV'
ref='/data/yiping/hnRNPL/08_20231219_chipseq/results/08_merge/02_Profile/merge_summits.merge.srt.bed'
computeMatrix reference-point -p 40 \
       --referencePoint center \
       -b 8000 -a 8000 \
       -R $ref \
       -S $dir/*bw \
       --skipZeros \
       -o $out/PolII_QC.gz \
       --outFileSortedRegions $out/PolII_QC.bed

plotHeatmap --dpi 300 -m $out/PolII_QC.gz \
     -out $out/PolII_QC.pdf \
     --colorMap Blues 
