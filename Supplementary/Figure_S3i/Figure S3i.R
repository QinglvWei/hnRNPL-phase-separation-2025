
b1 <- read.delim('./results/all_52729_5991Up_13411Down.txt') %>%
  mutate(change = as.factor(change)) %>%
  mutate(al = ifelse(change == 'NOT', 0.5, 0.6)) %>% 
  filter(change != 'NOT') %>% 
  separate(geneId, into = 'ENSEMBL') %>% 
  dplyr::select(c('ENSEMBL', 'change'))


c = bitr(unique(b1$ENSEMBL), fromType = 'ENSEMBL',
         toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db') %>%
  na.omit() %>%
  full_join(b1, by = 'ENSEMBL') %>%
  na.omit()


######### 
before_HL = read.delim('./chip_peak_list_6287peak_5937gene.txt')
PolII = read.delim('./all_52729_5991Up_13411Down.txt')

up = PolII %>% filter(change == 'UP')
down = PolII %>% filter(change == 'DOWN')

inter = intersect(unique(up$ENSEMBL), unique(down$ENSEMBL))
up_uniqe = up %>% filter(!(ENSEMBL %in% inter))
down_unique = down %>% filter(!(ENSEMBL %in% inter))

RNA = read.delim('./DGE_all_ID.txt') %>%
  mutate(type = ifelse(ENSEMBL %in% inter, 'BOTH', '')) %>%
  mutate(type = ifelse(ENSEMBL %in% up_uniqe$ENSEMBL, 'shHL_ChIP_UP', type)) %>%
  mutate(type = ifelse(ENSEMBL %in% down_unique$ENSEMBL, 'shHL_ChIP_DOWN', type)) %>%
  mutate(type = ifelse(is.na(type), 'NOT', type)) %>%
  mutate(type = ifelse(type == '', 'NOT', type))
  

a = RNA %>% filter(ENSEMBL %in% before_HL$ENSEMBL) %>%
  filter(change != 'NOT') %>%
  filter(ENSEMBL %in% down_unique$ENSEMBL)
write.table(a, file = './Final_152inter_39up_113down.txt', row.names = F, sep = '\t', quote = F)


b = c[, 1:2] %>%
  mutate(w = paste0(ENSEMBL, ENTREZID)) %>%
  distinct(w, .keep_all = T) %>%
  filter(ENSEMBL %in% a$ENSEMBL) %>%
  right_join(a, by = 'ENSEMBL') %>% na.omit()

genes = list(up = unique(b$ENTREZID[b$change == 'UP']),
             down = unique(b$ENTREZID[b$change == 'DOWN']))

compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.5,
                           pAdjustMethod = "BH")
compKEGG = setReadable(compKEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

a = as.data.frame(compKEGG)
write.table(a, file = './Final_enrichKEGG_152.txt', row.names = F, sep = '\t', quote = F)







