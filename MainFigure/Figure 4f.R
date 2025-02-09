chip <- read.delim('./05_20210810_itchipseq/06_macs2/broad/anno/chip_peak_list_6287peak_5937gene.txt')
rna <- read.delim('./02_20200825_rnaseq/statistics/DGE_all_ID.txt') %>%
  mutate(targets = ifelse(ENSEMBL %in% chip$ENSEMBL, 'Targets', 'Non_targets')) %>%
  mutate(targets = as.factor(targets))

rna$targets <- Replace(rna$targets, pattern = c('Targets:ChIP_targets(5213)', 'Non_targets:Non_targets(11657)'))

###### config
a <- read.delim('./gencode_grch35.ensg.enst') %>%
  mutate(ID = gene_id) %>%
  separate(ID, into = 'ID')

b <- read.delim('./grch38_v35_deeptools/gencode_grch35.Rna.for.deeptools', header = F)

up <- a %>%
  filter(ID %in% rna$ENSEMBL[rna$change=='UP'])
down <- a %>%
  filter(ID %in% rna$ENSEMBL[rna$change=='DOWN'])
not <- a %>%
  filter(ID %in% rna$ENSEMBL[rna$change=='NOT'])

up_for <- b %>%
  filter(V4 %in% up$transcript_id)

down_for <- b %>%
  filter(V4 %in% down$transcript_id)

not_for <- b %>%
  filter(V4 %in% not$transcript_id)

write.table(up_for, file = 'up_regulated_all.bed', row.names = F, col.names = F, sep = '\t', quote = F)
write.table(down_for, file = 'down_regulated_all.bed', row.names = F, col.names = F, sep = '\t', quote = F)
write.table(not_for, file = 'not_changed_all.bed', row.names = F, col.names = F, sep = '\t', quote = F)

rna1 <- rna %>%
  filter(targets == 'ChIP_targets(5213)')
data <- rna1[, c(1, 10)] %>%
  right_join(chip[, c(8, 14, 15)], by = 'ENSEMBL') %>%
  mutate(part = annotation) %>%
  separate(part, into = 'part', sep = ' \\(ENST') %>%
  group_by(change, part) %>%
  dplyr::summarise(n = n())

s <- rna1 %>%
  right_join(chip, by = 'ENSEMBL') %>%
  filter(!is.na(change))
s1 <- s %>%
  filter(change != 'NOT')
getwd()
write.table(s, file = '5537_rna_chip_peaks.txt', row.names = F, sep = '\t')
write.table(s1, file = '450_rna_chip_diff_peaks.txt', row.names = F, sep = '\t')
m <- s [, 1:10] %>%
  distinct(ENSEMBL, .keep_all = T)

m1 <- m%>%
  filter(change != 'NOT')

up <- a %>%
  filter(transcript_id %in% s$transcriptId[s$change=='UP'])
down <- a %>%
  filter(transcript_id %in% s$transcriptId[s$change=='DOWN'])
not <- a %>%
  filter(transcript_id %in% s$transcriptId[s$change=='NOT'])

up_for <- b %>%
  filter(V4 %in% up$transcript_id)

down_for <- b %>%
  filter(V4 %in% down$transcript_id)

not_for <- b %>%
  filter(V4 %in% not$transcript_id)

write.table(up_for, file = 'up_regulated_chip.bed', row.names = F, col.names = F, sep = '\t', quote = F)
write.table(down_for, file = 'down_regulated_chip.bed', row.names = F, col.names = F, sep = '\t', quote = F)
write.table(not_for, file = 'not_changed_chip.bed', row.names = F, col.names = F, sep = '\t', quote = F)




###### profile for RNA_group
dir='/home/yiping/wql/hnRNPL/04_20210809_chipseq/chip-seq/09_group_profile/01_rna_all'
cd $dir
computeMatrix reference-point \
       --referencePoint TSS \
       -b 5000 -a 5000 \
       -R down_regulated_all.bed  not_changed_all.bed  up_regulated_all.bed \
       -S ../HnRNPL_itChIP.bw \
       --skipZeros \
       -o HnRNPL_itChIP.gz \
       --outFileSortedRegions HnRNPL_itChIP.bed
plotProfile --dpi 300 -m HnRNPL_itChIP.gz \
              -out profile_group_rna.pdf \
              --perGroup \
              --plotTitle "ChIP-seq counts around RNA-seq genes TSS"


##### profile for RNA_ChIP_group
dir='/home/yiping/wql/hnRNPL/04_20210809_chipseq/chip-seq/09_group_profile/02_rna_chip'
cd $dir
computeMatrix reference-point \
       --referencePoint TSS \
       -b 5000 -a 5000 \
       -R down_regulated_chip.bed  not_changed_chip.bed  up_regulated_chip.bed \
       -S ../HnRNPL_itChIP.bw \
       --skipZeros \
       -o HnRNPL_itChIP.gz \
       --outFileSortedRegions HnRNPL_itChIP.bed
plotProfile --dpi 300 -m HnRNPL_itChIP.gz \
              -out profile_group_bind.pdf \
              --perGroup \
              --plotTitle "ChIP-seq peaks around RNA-seq genes TSS"
