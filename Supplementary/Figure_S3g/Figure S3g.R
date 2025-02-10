
data1 = read.delim('./results/all_52729_5991Up_13411Down.txt')

b <- data1 %>%
  mutate(change = as.factor(change)) %>%
  mutate(al = ifelse(change == 'NOT', 0.5, 0.6))
table(b$al)
up=length(b$name[(b$change == 'UP')])
down=length(b$name[(b$change == 'DOWN')])
table(b$al)

col3 <- c("#00B2EE","grey","#FF4500")
p <- ggplot(b, aes(log2(shNC_FC), log2(shHL_FC))) +
  geom_point(aes(color=change, alpha = al)) +
  scale_color_manual(values = col3) +
  labs(title = 'FoldChange Threshold: 2', x= "shNC (EF)", y="shHL (EF)") +
  theme(legend.position = 'none') +
  xlim(0, 10)+
  ylim(0, 10) +
  annotate('text', label = paste0('DOWN: ',down), x=9, y=3, size=5) +
  annotate('text', label = paste0('UP: ', up), x=3, y=9, size=5) 

p
pdf('volcano.pdf')
p
dev.off()