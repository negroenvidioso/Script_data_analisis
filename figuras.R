# hipertension arterial
# GSE112093
# GSE118578
# GSE75670
# hipertension arterial pulmonar
# GSE108707
# GSE21281
# GSE67597
# hipertension pulmonar
# GSE56914

# graficosgrafico d

### grafico D

mirna_32 = intersect(mirna_comb$gene, mirna_vote$gene)

library(dplyr)

todos = rbind(GSE108707,
              GSE112093,
              GSE118578,
              GSE21281,
              GSE56914,
              GSE67597,
              GSE75670)

resumen = todos %>%
  group_by(gene) %>%
  summarise(estudios = length(gene))

by_todos = group_by(todos, gene)
intento = summarise(by_todos, mean(logFC), mean(CI.L), mean(CI.R))

datos = mirna_vote[mirna_vote$gene %in% mirna_32,]

# grasfico d

tabla_grafico_cesar = merge(intento, datos, by = 'gene')

library(ggplot2)
library(ggrepel)

grafico_d = merge(tabla_grafico_cesar, resumen, by = 'gene')

lm = lm(grafico_d$`mean(logFC)` ~ round(grafico_d$`mean(logFC)`))

ggplot(data = grafico_d, aes(x = round(grafico_d$`mean(logFC)`, digits = 0),
                             y = grafico_d$`mean(logFC)`,
                             color = grafico_d$estudios))+
  theme_classic()+
  geom_point(size = 3)+
  coord_flip()+
  scale_color_gradient2("Studies")+
  geom_smooth(method='lm', se = F)+
  geom_errorbar(aes(ymin = grafico_d$`mean(CI.L)`, ymax = grafico_d$`mean(CI.R)`))+
  xlab(element_text('Promedio fold change reportado (log2)'))+
  ylab(element_text('Resumen fold change (log2)'))+
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# grafico e

grafico_e = merge(tabla_grafico_cesar, resumen, by = 'gene')

ggplot(data = grafico_e, aes(x = grafico_e$ndeg,
                             y = grafico_e$estudios,
                             color = grafico_e$`mean(logFC)`,
                             label = grafico_e$gene))+
  theme_classic()+
  geom_point(size = 3)+
  geom_text(aes(label = ifelse(grafico_e$gene %in% mirna_comb$gene,'')),hjust=2, vjust=1,
            size = 0.5)+
  geom_text_repel()+
  scale_color_gradient2("metafc", limits = c(-4, 4), 
                        low = "blue", mid = "grey", high = "red")+
  scale_x_continuous(limit = c(1,4))+
  xlab(element_text('papers citados'))+
  ylab(element_text('Papers tablas suplementarias'))+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20))

### grafico f

grafico_f = merge(todos, resumen, by = 'gene')
grafico_f = grafico_f[!duplicated(grafico_f$gene),]

ggplot(data = grafico_f, aes(x = grafico_f$estudios,
                             y = grafico_f$logFC,
                             color = grafico_f$logFC,
                             #fill = grafico_f$estudios
))+
  theme_classic()+
  geom_point(size = grafico_f$estudios)+
  scale_color_gradient2("metafc", limits = c(-8, 8), 
                        low = "blue", mid = "grey", high = "red")+
  xlab(element_text('Estudios'))+
  ylab(element_text('Señal de consistencia'))+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7))

## arreglar
ggplot(data = grafico_f)+
  theme_classic()+
  geom_point(aes(x = grafico_f$estudios,
                 y = grafico_f$logFC,
                 size = grafico_f$estudios,
                 fill = grafico_f$logFC), shape = 21)+
  guides(fill = guide_legend(override.aes = list(size=8)))+
  scale_color_gradient2("metafc", limits = c(-8, 8), 
                        low = "blue", mid = "grey", high = "red")+
  xlab(element_text('Studies'))+
  ylab(element_text('Sing consistency'))
