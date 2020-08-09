### lugar de trabajo ###
setwd("~/Documentos/tesis/plataformas/")

### generar lista para obtener los mirna de las plataformas ###
GPL21572 = read.delim("GPL21572-124634.txt", sep = "\t", header = T)
GPL21572_hsa = as.character(GPL21572$miRNA_ID[grep("hsa", GPL21572$miRNA_ID)])

GPL20712 = read.delim("GPL20712-tbl-1.txt", sep = "\t", header = T)
GPL20712_hsa = as.character(GPL20712$X[grepl("hsa", GPL20712$X)])

GPL18402 = read.delim("GPL18402-tbl-1.txt", sep = "\t", header = F)
GPL18402_hsa = as.character(GPL18402$V1[grepl("hsa", GPL18402$V1)])

GPL17292 = read.delim("GPL17292-tbl-1.txt", sep = "\t", header = T)
GPL17292_hsa = as.character(GPL17292$Hy3[grepl("hsa", GPL17292$Hy3)])

GPL18587 = read.delim("GPL18587-tbl-1.txt", header = T, sep = "\t")
GPL18587_hsa = as.character(GPL18587$Hy3[grepl("hsa", GPL18587$Hy3)])

GPL16601 = read.delim("GPL16601-tbl-1.txt", header = F, sep = "\t")
GPL16601_hsa = as.character(GPL16601$V1[grepl("hsa", GPL16601$V1)])

GPL10305 = read.delim("GPL10305-tbl-1.txt", header = F, sep = "\t")
GPL10305_hsa = as.character(GPL10305$V3[grepl("hsa", GPL10305$V3)])

GPL4133 = read.delim("GPL4133-12599.txt", header = T, sep = "\t")  # es de mRNA

GPL10558 = read.delim("GPL10558-50081.txt", header = T, sep = "\t")
GPL10558_hsa = subset(GPL10558, grepl(GPL10558))  # es de mRNA

### listas ###

lista_hsa_unique =  unique(c(GPL21572_hsa,
                             GPL20712_hsa,
                             GPL18402_hsa,
                             GPL17292_hsa,
                             GPL18587_hsa,
                             GPL16601_hsa,
                             GPL10305_hsa))

### matriz binaria ### 

write.table(unique(lista_hsa_unique),
            file = "lista_hsa_unique.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL21572_hsa),
            file = "GPL21572_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL20712_hsa),
            file = "GPL20712_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL18402_hsa),
            file = "GPL18402_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL17292_hsa), #
            file = "GPL17292_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL18587_hsa), #
            file = "GPL18587_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL16601_hsa),
            file = "GPL16601_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

write.table(unique(GPL10305_hsa),
              file = "GPL10305_hsa_mirna.txt",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

library(UpSetR)
setwd("~/Documentos/tesis/plataformas/matriz_binaria/")
matriz = read.delim("matrix_plataformas_hypertension.txt", header = T, sep = "\t")
names(matriz) = c("miRNA", "GPL10305", "GPL16601", "GPL17292",
                  "GPL18402", "GPL18587", "GPL20712", "GPL21572")
upset(matriz,
      nsets = 20,
      mb.ratio = c(0.55, 0.45),
      order.by = "freq",
         matrix.color="steelblue", # pelotas de la matriz
      sets.bar.color = "steelblue", # barras laterales
      main.bar.color = "#FD0101", # barras superiores
      mainbar.y.label = "miRNAs unicos por plataforma",
      sets.x.label = "Cantidad total de miRNAs presentes",
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))
   