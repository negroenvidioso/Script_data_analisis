setwd('F:/respaldo/Documentos/tesis/microarrays')

library(readxl)
datos_utilizados <- read_excel("datos_utilizados.xlsx")
View(datos_utilizados)
names(datos_utilizados) = c('GEO', 'Muestras analizadas', 'Muestras No analizadas')
library(reshape2)
DF1 <- melt(datos_utilizados, id.var="GEO")

library(ggplot2)

ggplot(DF1, aes(x = DF1$GEO, y = DF1$value, fill = DF1$variable)) +
  geom_bar(stat="identity")+
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values=c('#20B2AA', '#808080'))+
  xlab(element_text('GEO accession'))+
  ylab(element_text('Numero de muestras'))+
  guides(fill=FALSE)
