setwd('/home/ignacio/Documentos/tesis_allan/sra/PRJNA329554/')
library(DESeq2)

# Leemos conteos 
diseño = read.table('samples.table', header = T)
PRJNA329554 = read.table("genes.fpkm_table", header = T)
names(PRJNA329554)
names(PRJNA329554) = c('gene_id',
                       'control_2',
                       'control_1',
                       'control_3',
                       'hipertension_1',
                       'hipertension_2',
                       'hipertension_3')
View(PRJNA480270)
dim(PRJNA480270)

#Subseting de los datos. S?lo la columna de los conteos:
countdata = PRJNA329554[, c(2:7)]
countdata = round(countdata, digits = 0)
dim(countdata)

#coldata: tabla con informacion de las muestras
coldata = data.frame(
  row.names = colnames( countdata ),
  tipo = c( rep('control', 3), rep("hipertension", 3)))

# Construimos el objeto DESeqDataSet a partir de la matriz de conteos y la informaci?n de las muestras
dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~ tipo)
dds <- DESeq(dds)

# Si quereremos extraer los conteos normalizados:
conteosnormalizados = counts(dds, normalized = TRUE)
dim(conteosnormalizados)
View(conteosnormalizados) #Conteos normalizados
View(countdata) #Conteos brutos

# Extraemos los resultados co el logFC estimado y p values
resDESeq2 <- results(dds)

#Asignamos los nombres de las filas igual a los nombres de genes
rownames(resDESeq2) = PRJNA329554$gene_id

#Guardamos el archivo de resultados sin filtro
write.table(resDESeq2,
            file="DESeq2_rna-seq_PRJNA329554.txt",
            sep="\t")

# filtramos la tabla
resDESeq2 = as.data.frame(resDESeq2)
resDESeq2 = na.omit(resDESeq2)
DEGsdeseq2 = resDESeq2[resDESeq2$padj< 0.01,]
genenames = rownames(resDESeq2[resDESeq2$padj< 0.05,])
View(DEGsdeseq2)

# anotacion de los nombres
anotacion = read.table('genes.attr_table', header = TRUE)
anotacion = subset(anotacion, selec = c('tracking_id', 'gene_short_name'))
resultado_nombres = anotacion[anotacion$tracking_id %in% genenames,]