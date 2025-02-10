
######## Deseq2.Rda and Deseq_result.Rda calculated from Figure 4c.R
load("Deseq2.Rda")
load("Deseq_result.Rda")

data <- as.data.frame(vsd) %>%
  mutate(ENSEMBL=rownames(.)) %>%
  separate(ENSEMBL, into = 'ENSEMBL') %>%
  inner_join(DGE_diff[, 1:2], by = 'ENSEMBL') %>%
  tibble::column_to_rownames('gene_name') %>%
  dplyr::select(-'ENSEMBL')
  
vsd_DGE_ID <- data
## heatmap
names(vsd_DGE_ID) <- c("shNC_1", "shNC_2", "sh1_1", "sh1_2", 'sh2_1', 'sh2_2')

class <- as.factor(c("shNC","shNC", "shHL-1","shHL-1", "shHL-2", 'shHL-2')) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(vsd_DGE_ID)

pheatmap(vsd_DGE_ID,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col =annotation_col, 
         annotation_legend=TRUE, 
         show_rownames = F,
         scale = "row", 
         color =colorRampPalette(c("blue", "white","red"))(100),
         cellwidth = 50, cellheight = 0.18,
         fontsize = 10,
         cutree_rows = 2, cutree_cols = 3)
## save 7.43*7.43
