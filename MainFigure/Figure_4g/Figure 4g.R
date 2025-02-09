
a <- read.delim('321_RNA_ChIP_down.txt') %>%
  filter(PValue < 0.05) %>%
  arrange(PValue)

go <- a %>%
  separate(Term, into = c('a', 'Term'), sep = '~') %>%
  filter(Category != 'KEGG_PATHWAY')

kegg <- a %>%
  separate(Term, into = c('a', 'Term'), sep = ':') %>%
  filter(Category == 'KEGG_PATHWAY')

a <- rbind(go, kegg) %>%
  arrange(PValue)

s <- split(a, a$Category)
b <- rbind(s[[1]][1:10, ], s[[2]][1:10, ], s[[3]][1:10, ], s[[4]][1:10, ])
b <- na.omit(b)
data <- b %>%
  arrange(Category, PValue) %>%
  mutate(number=factor(rev(1:nrow(.))))

data$Term <- as.factor(data$Term)

shorten_names <- function(x, n_word=40, n_char=400){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 400))
  {
    if (nchar(x) > 400) x <- substr(x, 1, 400)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply(
  levels(data$Term)[as.numeric(data$Term)],shorten_names))
names(labels) = rev(1:nrow(data))
labels


ggplot(data, aes(number, -log10(PValue), fill=Category)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  theme_bw() + 
  scale_x_discrete(labels=labels) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = "", y="-log10(PValue)", title = "hnRNPL 321 RNA down & itChIP targets") +
  theme(plot.title = element_text(color='black', hjust = 0.5, size = 18)) +
  theme(axis.text.x = element_text(colour="black",size=13),
        axis.text.y = element_text(colour="black",size=13),
        legend.key.size=unit(5,'mm'),
        legend.title=element_text(colour = "black", size=14),
        legend.text = element_text(colour = "black", size=14))+
  theme(legend.position = c(0.8, 0.2)) +
  theme(title=element_text(size=16,color="black")) +
  geom_hline(yintercept = -log10(0.05), lty = 2)

### 7.5 *10
