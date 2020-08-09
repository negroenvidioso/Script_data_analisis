ncRNA_human <- read.delim("~/Documentos/tesis/ncRNA_human.txt")

ncRNA_human[order(ncRNA_human$x),]

library(ggplot2)
ggplot(data = ncRNA_human, aes(x = x, y = y, fill = x))+
  theme_classic()+
  geom_bar(stat = "identity")+
  ggtitle("Predicción de lncRNA por Cromosoma")+
  geom_text(aes(label = ncRNA_human$y), vjust = -0.3, size = rel(3.9))+
  labs(y = "Número de lncRNA", x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, color = "#666666", size = rel(3)),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)))
