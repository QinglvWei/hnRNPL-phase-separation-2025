chip <- read.delim('chip_peak_list_6287peak_5937gene.txt')
rna <- read.delim('DGE_all_ID.txt') %>%
  mutate(targets = ifelse(ENSEMBL %in% chip$ENSEMBL, 'Targets', 'Non_targets')) %>%
  mutate(targets = as.factor(targets))

rna$targets <- Replace(rna$targets, pattern = c('Targets:ChIP_targets(5213)', 'Non_targets:Non_targets(11657)'))

data <- split(rna, rna$targets)
title = paste0('p=',signif(wilcox.test(data[[1]]$log2FoldChange, data[[2]]$log2FoldChange, alternative = "two.sided")$p.value,3))
ggplot(rna, aes(log2FoldChange, colour = targets))+
  scale_color_manual(values = c("indianred1", 'lightblue'))+
  #scale_color_manual(values = c("#EE7621", '#008B00'))+
  stat_ecdf(lwd=2) +
  labs(title = 'ChIP_RNA', x = 'shhnRNPL log2(FC)', y='Cumulative proportion')  +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  guides(color=guide_legend(title=NULL)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(plot.title=element_text(size=rel(2),hjust=0.5), ##标题
        axis.title=element_text(size=rel(1.5)), ## 坐标轴标题
        axis.text=element_text(size=rel(1.5)), ## 坐标轴横纵坐标
        panel.grid.major=element_line(color="white"), ## 边框
        panel.grid.minor=element_line(color="white"),
        legend.key.size=unit(6,'mm')) + ## 内部框架
  guides(fill=guide_legend(title=NULL))  + ## 图注去标题
  geom_vline(xintercept = c(-log2(1.5), 0, log2(1.5)), lwd=0.5,col="gray",lty=2) +
  geom_hline(yintercept = c(0,0.5, 1), lwd=0.5,col="gray",lty=2) +
  scale_x_continuous(limits = c(-3, 3),breaks = c(-3, -2,-1, 0, 1,2, 3)) +
  xlim(-2, 2) +
  annotate("text", x=1.5,y=0.05, label=title, col="red", size = 4) ### 画在内部，可以调整图片的大小
#### save(5*5)