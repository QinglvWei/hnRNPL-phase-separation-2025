######### deseq2 for hnRNPL RNAseq
setwd('./03_htseq/')
dir()
directory <- getwd()
sampleFiles <- grep("*", list.files(directory), value=TRUE)
sampleCondition <- as.factor(c(rep("NC", 2), rep("sh", 4)))

sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                       directory = directory, 
                                       design = ~condition)

dds <- DESeq(ddsHTSeq)
counts <- assay(dds) ## counts数
vsd <- vst(dds, blind = FALSE)
vsd <- assay(vsd)

count <- as.data.frame(counts)

d = estimateSizeFactors(dds)
dds2 = counts(d, normalized = TRUE)

save(sampleTable, ddsHTSeq, counts, vsd, dds2, num, file = "../statistics/Deseq2.Rda")
res <- results(dds, contrast = c("condition", "sh", "NC"))

mcols(res)$description
metadata(res)$filterThreshold

resOrdered <- res[order(res$pvalue), ]
resOrdered_p <- as.data.frame(resOrdered)

DGE_all <- na.omit(resOrdered_p)
dim(DGE_all)

logFC_cutoff <- log2(1.5)
DGE_all$change = as.factor(ifelse
                           (DGE_all$pvalue < 0.05&abs(DGE_all$log2FoldChange) > logFC_cutoff,
                             ifelse(DGE_all$log2FoldChange > -logFC_cutoff, "UP", "DOWN"),
                             "NOT"))
table(DGE_all$change)
RNA <- data.table::fread('./Genome/Gencode/gencode_v35_annotation_gtf_RNA_type.txt')

DGE_all1 <- DGE_all %>%
  mutate(gene_id = rownames(.)) %>%
  inner_join(RNA, by = 'gene_id') %>%
  dplyr::select(c(8:10, 1:7)) %>%
  separate(gene_id, into = 'ENSEMBL') 
  
DGE_diff <- DGE_all1 %>%
  filter(change != "NOT")

save(res, resOrdered_p, DGE_all1, DGE_diff, file = "../statistics/Deseq_result.Rda")

write.table(DGE_all1, file='../statistics/DGE_all_ID.txt', sep = '\t', row.names = F)
write.table(DGE_diff, file='../statistics/DGE_diff_ID.txt', sep = '\t', row.names = F)



load(file = "Deseq_result.Rda")

data <- DGE_all1
logFC_cutoff <- log2(1.5)
this_title <- paste0("Cutoff for logFC_cutoff is ", round(logFC_cutoff, 3), 
                     "\nThe number of up gene is ", nrow(data[data$change == "UP", ]), "\nThe number of down gene is ", nrow(data[data$change == "DOWN", ]))

title =  paste0("Cutoff for logFC_cutoff is ", round(logFC_cutoff, 3))
up =  paste0("UP: ", nrow(data[data$change == "UP", ]))
down = paste0("DOWN: ", nrow(data[data$change == "DOWN", ]))

col3 <- c("#00B2EE", "grey","#FF4500")
ggplot(data = data, aes(x = log2FoldChange, y= -log10(pvalue), color = change))+
  scale_color_manual(values = col3)+
  geom_point(size=1)+
  geom_hline(yintercept = -log10(0.05), lty = 2)+
  geom_vline(xintercept = c(-log2(1.2), log2(1.2)), lty = 5)+
  theme_bw()+
  coord_fixed(ratio = 0.4) +
  labs (title = title, x="log2 fold change", y="-log10(Pvalue)")+
  theme(plot.title = element_text(color='black', hjust = 0.5, size = 20)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(colour="black",size=15, hjust=0.5,vjust=0.5),
        axis.text.y = element_text(colour="black",size=15, hjust=0.5,vjust=0.5),
        axis.title = element_text(colour = "black", size=20)) +
  ylim(0, 15) +
  xlim(-4, 4) +
  annotate('text',label = up,  x=3.5, y=15, colour = "#FF4500", size=5) +
  annotate('text',label = down, x= -3.2, y=15, colour = "#00B2EE", size=5)

## save 7.43*5.46