x = meta_degs_comb_all@metaresult
x = subset(x, x$metap < 0.05)
y = subset(x, x$metafc > 0.5)
y1 = subset(x, x$metafc < -0.5)

je = rbind(y, y1)

write.table(je,
            file = 'listado_avance_comb.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

deg = meta_degs_vote_all@metaresult
deg1 = subset(deg, deg$degvcount != '1.Unperturbed')

intersect(je$gene, deg1$gene)

library(VennDiagram)
library(scales)
venn.plot <- draw.pairwise.venn(area1      = length(je$gene),
                                area2      = length(deg1$gene),
                                cross.area = length(intersect(je$gene, deg1$gene)),
                                category   = c("Combinig Method", "DEGs Vote"),
                                scaled     = FALSE,
                                fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                                cat.col = c("#440154ff", '#21908dff'),
                                cex = 3,
                                cat.cex = 3,
                                cat.pos = c(-20, 20),
                                cat.dist = c(0.055, 0.055))
require(gridExtra)
grid.arrange(gTree(children = venn.plot),
             top = textGrob('miRNA compartidos por MetaVolcanoR', gp = gpar(fontsize = 30, font = 8)),
             bottom = textGrob(''))

sub_3 = subset(deg, deg$ndeg == 3)          
sub_3$gene
sub_2 = subset(deg, deg$ndeg == 2)

intersect(je$gene, sub_3$gene)
length(intersect(je$gene, sub_2$gene))

datos_combining = je[je$gene %in% intersect(je$gene, deg1$gene),]

write.table(datos,
            file = 'subset_mirna_compartidos_deg_combining.txt',
            quote = F,
            sep = '\t',
            col.names = T,
            row.names = F)

datos = read.delim('~/Documentos/tesis/microarrays/metavolcano_mirna_pclib/subset_mirna_compartidos_deg_combining.txt',
                   header = T,
                   sep = '\t')
library(ggplot2)
library(ggrepel)

ggplot(data = datos, aes(x = datos$gene, y = datos$ndeg, fill = datos$ddeg))+
  theme_classic()+
  coord_flip()+
  geom_bar(stat = 'identity')

datos_combining_2 = merge(datos_combining, datos, by = 'gene')

ggplot(data = datos_combining_2, aes(x = datos_combining_2$metafc, y = datos_combining_2$ndeg, color = datos_combining_2$metafc))+
  theme_classic()+
  geom_point(size = 3)+
  scale_y_continuous(limit = c(1,4))+
  scale_color_gradient2("metafc", limits = c(-4, 4), 
                       low = "blue", mid = "grey", high = "red")+
  xlab(element_text('summary fold-change'))+
  ylab(element_text('Studies'))

### tabla con promedios ###
GSE108707_sub = GSE108707[GSE108707$gene %in% datos$gene,]
GSE112093_sub = GSE112093[GSE112093$gene %in% datos$gene,]
GSE118578_sub = GSE118578[GSE118578$gene %in% datos$gene,]
GSE21281_sub = GSE21281[GSE21281$gene %in% datos$gene,]
GSE56914_sub = GSE56914[GSE56914$gene %in% datos$gene,]
GSE67597_sub = GSE67597[GSE67597$gene %in% datos$gene,]
GSE75670_sub = GSE75670[GSE75670$gene %in% datos$gene,]

library(dplyr)

todos = rbind(GSE108707_sub,
              GSE112093_sub,
              GSE118578_sub,
              GSE21281_sub,
              GSE56914_sub,
              GSE67597_sub,
              GSE75670_sub)
resumen = todos %>%
  group_by(gene) %>%
  summarise(estudios = length(gene))
  
by_todos = group_by(todos, gene)
intento = summarise(by_todos, mean(logFC), mean(CI.L), mean(CI.R))
write.table(intento,
            file = 'promedio_de_mirna_estudios.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

tabla_grafico_cesar = merge(intento, datos, by = 'gene')

ggplot(data = tabla_grafico_cesar, aes(x = tabla_grafico_cesar$`mean(logFC)`, y = datos_combining_2$ndeg, color = datos_combining_2$metafc))+
  theme_classic()+
  geom_point(size = 3)+
  scale_y_continuous(limit = c(1,4))+
  scale_color_gradient2("metafc", limits = c(-4, 4), 
                        low = "blue", mid = "white", high = "red")+
  xlab(element_text('summary fold-change'))+
  ylab(element_text('Studies'))

### grafico D
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

### Grafico E
grafico_e = merge(tabla_grafico_cesar, resumen, by = 'gene')

ggplot(data = grafico_e, aes(x = grafico_e$ndeg,
                             y = grafico_e$estudios,
                             color = grafico_e$`mean(logFC)`,
                             label = grafico_e$gene))+
  theme_classic()+
  geom_point(size = 3)+
  geom_text(aes(label = ifelse(grafico_e$gene %in% datos_combining$gene,'')),hjust=2, vjust=1,
            size = 20)+
  geom_text_repel()+
  scale_color_gradient2("metafc", limits = c(-4, 4), 
                        low = "blue", mid = "grey", high = "red")+
  scale_x_continuous(limit = c(1,4))+
  xlab(element_text('Papers tablas suplementarias'))+
  ylab(element_text('papers citados'))+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20))

### grafico f

grafico_f = merge(todos, resumen, by = 'gene')

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
        axis.text = element_text(size = 20))

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

### grafico d_1

# arterial hipertension
# GSE112093, GSE118578, GSE75670
arterial = rbind(GSE112093, GSE118578, GSE75670)
arteria_0.05 = subset(arterial, arterial$P.Value < 0.05)
arterial_up = subset(arteria_0.05, arteria_0.05$logFC > 0.5)
arterial_down = subset(arteria_0.05, arteria_0.05$logFC < -0.5)

# pulmonar arterial hipertension
# GSE108707, GSE67597, GSE21281
arterial_pulmonar = rbind(GSE108707, GSE67597, GSE21281)
arterial_pulmonar_0.05 = subset(arterial_pulmonar, arterial_pulmonar$P.Value < 0.05)
arterial_pulmonar_up = subset(arterial_pulmonar_0.05, arterial_pulmonar_0.05$logFC > 0.5)
arterial_pulmonar_down = subset(arterial_pulmonar_0.05, arterial_pulmonar_0.05$logFC < -0.5)

# pulmonar hipertension
# GSE56914
library(scales)

pulmonar = GSE56914
pulmonar_0.05 = subset(pulmonar, pulmonar$P.Value < 0.05)
pulmonar_up = subset(pulmonar_0.05, pulmonar_0.05$logFC > 0.5)
pulmonar_down = subset(pulmonar_0.05, pulmonar_0.05$logFC < -0.5)


tipos_hipertension = c('arterial hipertension', 'pulmonar arterial hipertension', 'pulmonar hipertension')
up = c(length(unique(arterial_up$gene)),
       length(unique(arterial_pulmonar_up$gene)),
       length(unique(pulmonar_up$gene)))
down = c(-length(unique(arterial_down$gene)),
         -length(unique(arterial_pulmonar_down$gene)),
         -length(unique(pulmonar_down$gene)))
error_up = c(length(arterial_up$gene),
          length(arterial_pulmonar_up$gene),
          length(pulmonar_up$gene))
error_down = c(length(arterial_down$gene),
               length(arterial_pulmonar_down$gene),
               length(pulmonar_down$gene))
grafico_d1 = data.frame(tipos_hipertension,
                        up,
                        down,
                        error_up,
                        error_down)

ggplot(data = grafico_d1, aes(x = grafico_d1$tipos_hipertension,
                              y = grafico_d1$up))+
  theme_classic()+
  coord_flip()+
  geom_bar(stat = 'identity', position="stack", fill = alpha('#B22222'))+
  geom_bar(aes(y = grafico_d1$down), stat = 'identity', position="stack", fill = alpha('#00008B'))+
  geom_errorbar(aes(ymin = c(270, 90, 0), ymax = c(274, 94, 0)))+
  geom_errorbar(aes(ymin = c(-110, 0, 0), ymax = c(-126, 0, 0)))+
  scale_x_discrete(limits = c('pulmonar hipertension',
                              'pulmonar arterial hipertension',
                              'arterial hipertension'))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size = 20)
        )+
  ylab(element_text('gene'))
        