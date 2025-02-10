
###### chip_peak_list_6287peak_5937gene.txt file from 'Figure 4b.R'
before_HL = read.delim('/Users/carolwei/GH/bioinfomatics/03_hnRNPL/05_20210810_itchipseq/06_macs2/broad/anno/chip_peak_list_6287peak_5937gene.txt')


####### all_52729_5991Up_13411Down.txt file from 'Figure S3g.R'
PolII = read.delim('./all_52729_5991Up_13411Down.txt')

up = PolII %>% filter(change == 'UP')
down = PolII %>% filter(change == 'DOWN')

genes = list(pol_uP = unique(up$ENSEMBL),
             pol_down = unique(down$ENSEMBL))


inter = intersect(unique(up$ENSEMBL), unique(down$ENSEMBL))
up_uniqe = up %>% filter(!(ENSEMBL %in% inter))
down_unique = down %>% filter(!(ENSEMBL %in% inter))

###### DGE_all_ID.txt file from 'Figure 4c.R'
RNA = read.delim('/Users/carolwei/GH/bioinfomatics/03_hnRNPL/02_20200825_rnaseq/statistics/DGE_all_ID.txt') %>%
  mutate(type = ifelse(ENSEMBL %in% inter, 'BOTH', '')) %>%
  mutate(type = ifelse(ENSEMBL %in% up_uniqe$ENSEMBL, 'shHL_ChIP_UP', type)) %>%
  mutate(type = ifelse(ENSEMBL %in% down_unique$ENSEMBL, 'shHL_ChIP_DOWN', type)) %>%
  mutate(type = ifelse(is.na(type), 'NOT', type)) %>%
  mutate(type = ifelse(type == '', 'NOT', type))
  
genes = list(polII_DOWN = down_unique$ENSEMBL,
             HL_ChIP = unique(before_HL$ENSEMBL),
             RNA = unique(RNA$ENSEMBL[RNA$change != 'NOT'])
             )

pdf('pol_HL_RNA_39UP_113DOWN.pdf')
vennplot(genes)
dev.off()
