
####### pool_peaks.broadPeak file from 'step 5 macs2' of /Users/carolwei/Downloads/github-plot/hnRNPL-Chipseq-alignment.sh
a <- read.delim('pool_peaks.broadPeak', header = F)
names(a) <- c("chr", "start","end",'name', '10*-log10(pvalue)', 'strand', 'fold_enrichment', '-log10(pvalue)', "-log10(qvalue)")

dat1 <- a[, 1:5] %>% 
  mutate(a = substr(chr, 1,2)) %>%
  filter(a=='ch') %>%
  dplyr::select(-'a')
names(dat1) <- paste0('V', 1:5)
dir.create('./anno')
write.table(dat1, file = './anno/narowpeak_for_chipseeker', quote = F, sep = '\t', row.names = F)


#### load genome
setwd('/Users/carolwei/GH/bioinfomatics/Genome/Gencode/')
gtf <- "gencode.v35.primary_assembly.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf, format = "gtf") 

setwd('/Users/carolwei/GH/bioinfomatics/hnRNPL/05_20210810_itchipseq/06_macs2/broad/anno/')
name <- 'narowpeak_for_chipseeker'

### annotation
peakAnno <- annotatePeak(name, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
RNA <- read.delim('/Users/carolwei/GH/bioinfomatics/Genome/gencode_v35.txt') %>%
  separate(gene_id, into = 'ENSEMBL') %>%
  distinct(ENSEMBL, .keep_all = T)
bed1 <- read.delim('../pool_peaks.broadPeak', header = F)
data <- as.data.frame(peakAnno) %>%
  separate(geneId, into = 'ENSEMBL') %>%
  left_join(RNA[, 1:2], by = 'ENSEMBL') %>%
  left_join(bed1[, c(4,7,8,9)], by = 'V4') 

names(data)[21:23] <- c('fold_enrichment', '-log10(pvalue)', "-log10(qvalue)")
write.table(data, file = 'chip_peak_list_all.txt', sep = '\t', row.names = F)

a <- data %>%
  filter(`-log10(qvalue)` >-log10(0.00000001))

#### DGE_all_ID.txt is the differential expressed gene list of hnRNPL RNAseq
rna <- read.delim('DGE_all_ID.txt')
length(intersect(a$ENSEMBL, rna$ENSEMBL))

write.table(a, file = 'chip_peak_list_6287peak_5937gene.txt', sep = '\t', row.names = F)


b <- a[, 1:8] %>%
  mutate(part = annotation) %>%
  separate(part, into = 'part', sep = ' \\(ENST') %>%
  group_by(part) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n)) %>%
  mutate(Feature = part) %>%
  mutate(Feature = ifelse(substr(Feature, 1, 3) =='Dow', 'Downstream', Feature)) %>%
  mutate(Feature = gsub('1-2kb', '1-3kb', Feature)) %>%
  mutate(Feature = gsub('2-3kb', '1-3kb', Feature)) %>%
  group_by(Feature) %>%
  dplyr::summarise(Frequency = sum(n)) %>%
  mutate(sum = sum(Frequency)) %>%
  mutate(Freq = round(Frequency/sum*100, 4)) 
  
b <- as.data.frame(b)
c <- b[, c(1,2, 4)] %>%
  dplyr::arrange(desc(Freq)) %>%
  mutate(label =  as.factor(paste0(Feature, ' (', Freq, " %)"))) %>%
  mutate(sample = 'HnRNPL')
c$label <- factor(c$label, levels = rev(c$label))



ggplot(c, aes(sample, Freq, fill = label)) +
  geom_bar(stat='identity', position = 'stack', width = 1, color = 'white') +
  scale_fill_brewer() 

chip <- read.delim('chip_peak_list_6287peak_5937gene.txt')[, c(1:3, 14, 20:23)] %>%
  distinct(ENSEMBL, .keep_all = T)
write.xlsx(chip, file = 'chip_list_5937.xlsx', row.names = F)