# Graficos con el numero de mirna en humanos en mirBase

setwd("D:/Documentos/Magister/tesis")
mirna_humanos <- read_excel("mirna_humanos.xlsx")

a <- mirna_humanos$ID
b <- mirna_humanos$Mature1_ID
c <- mirna_humanos$Mature2_ID

a1 <- sum(!is.na(a))
b1 <- sum(!is.na(b))
c1 <- sum(!is.na(c))

mirna <- c("miRNA Total", "miRNA Mature", "miRNA Mature 2")
numero <- as.numeric(c(a1, b1, c1))
todo <- data.frame(mirna, numero)
todo$mirna = factor(todo$mirna, levels=c("miRNA Mature 2", "miRNA Mature", "miRNA Total"))  # ordena los calores como yo quiero

library(ggplot2)

ggplot(data = todo, aes(x = todo$mirna, y = todo$numero, fill = todo$mirna))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label= todo$numero), hjust = 1.1, color = "black", size = rel(5.5))+
  theme_classic()+
  ggtitle("miRNA reportated in miRBase \n for Homo sapiens")+
  theme(plot.title = element_text(hjust = 0.5,vjust = 3.5, color = "#666666", size = 25, face = "italic"),
        axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = "none")+
  labs(x = element_blank(), y = element_blank() )+
  coord_flip()
