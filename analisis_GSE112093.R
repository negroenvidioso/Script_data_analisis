setwd("D:/Documentos/Magister/tesis/pc_lib/pc_lib/Documentos/tesis/microarrays/GSE112093")

library(oligo)
library(ggplot2)
library(limma)

celpath = "D:/Documentos/Magister/tesis/pc_lib/pc_lib/Documentos/tesis/microarrays/GSE112093_RAW"

list = list.files(celpath, full.names = T)
data = oligo::read.celfiles(list)

int = oligo::intensity(data)

ph = data@phenoData

ph@data[,1] = c('control1',
                'control2',
                'control3',
                'control4',
                'control5',
                'hipertension1',
                'hipertension2',
                'hipertension3',
                'hipertension4',
                'hipertension5')

ph@data

### imagenes de la placa ## se pueden eliminar revisar mas tarde
#for (i in 1:10)
#{
#  name = paste("image", i, ".jpg", sep = "")
#  jpeg(name)
#  image(data[,i], main = ph@data$index[i])
#  dev.off()
#}

## pseudo imagenes de la placa ## se puede eliminar revisar mas tarde
#Pset = oligo::fitProbeLevelModel(data)
#for (i in 1:10)
#{
#  name = paste("pseudoimage", i, ".jpg", sep = "")
#  jpeg(name)
#  image(Pset, which = i, type = 'residuals', main = ph@data$index[i])
#  dev.off()
#}

### histogramas de distribucion de los datos

pmexp = oligo::pm(data)

sampleNames = vector()
logs = vector()
for (i in 1:10)
{
  sampleNames = c(sampleNames, rep(ph@data[i,1], dim(pmexp)[1]))
  logs = c(logs, log2(pmexp[,i]))
}

logData = data.frame(logInt = logs,
                     sampleName = sampleNames)

### normalizacion 
data.rma = rma(data)
data.matrix = exprs(data.rma)

### histograma
ggplot(logData, aes(logInt, color = sampleName))+
  geom_density()+
  theme_classic()

### box plot
### sin normalizar
ggplot(logData, aes(sampleName, logInt))+
  geom_boxplot()+
  theme_classic()

### normalizado
names = vector()
normlogs = vector()
for (i in 1:10)
{
  names = c(names, rep(ph@data[i,1], dim (data.matrix)[1]))
  normlogs = c(normlogs, data.matrix[,1])
}
normdata = data.frame(norm_logInt = normlogs, sampleName = names)

ggplot(normdata, aes(sampleName, norm_logInt))+
  geom_boxplot()+
  theme_classic()+
  ylim (2, 16)+
  ggtitle("GSE1129093 normalizada")

### grafico de MA para microarray
## sin normalizar
#for (i in 1:10)
#{
#  name = paste ("MAplot", i, ".jpg", sep = "")
#  jpeg(name)
#  MAplot(data, which = i)
#  dev.off()
#}

### normalizado
#for (i in 1:10)
#{
#  name = paste("MAplot_norm", i, ".jpg", sep = "")
#  jpeg(name)
#  oligo::MAplot(data.rma, which = i)
#  dev.off()
#}

### PCA ###

color = c('Magenta',
          'Magenta',
          'Magenta',
          'Magenta',
          'Magenta',
          'steelblue',
          'steelblue',
          'steelblue',
          'steelblue',
          'steelblue')

data.pc = prcomp(t(data.matrix),scale. = T)
#View(data.pc$x)

plot(data.pc$x[,2], col = color)

ph@data[,2] = c("control",
                "control",
                "control",
                "control",
                "control",
                "hipertension",
                "hipertension",
                "hipertension",
                "hipertension",
                "hipertension")

colnames(ph@data)[2] = 'source'

groups = ph@data$source     
f = factor(groups, levels = c("control", "hipertension"))

desing = model.matrix(~ 0 + f)
colnames(desing) = c('control', 'hipertension')

data.fit = lmFit(data.matrix, desing)

data.fit$coefficients
contrast.matrix = makeContrasts(hipertension - control, levels = desing)
data.fit.con = contrasts.fit(data.fit, contrast.matrix)

data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients

volcanoplot(data.fit.eb, coef = 1, highlight = 100)

### IDs genes 

tab = topTable(data.fit.eb, coef = 1, number = 250, adjust.method = 'BH', confint = T)

tab1 = DataFrame(gene = row.names(tab),
                 logFC = tab$logFC,
                 CI.L = tab$CI.L,
                 CI.R = tab$CI.R,
                 AveExpr = tab$AveExpr,
                 t = tab$t,
                 P.Value = tab$P.Value,
                 adj.P.Val = tab$adj.P.Val,
                 B = tab$B)

write.table(tab1,
            file = 'GSE112093_table.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

topgenes = tab[tab[,'P.Value'] < 0.001,]
dim(topgenes)

topups = topgenes[topgenes[, 'logFC'] > 0.5,]
dim(topups)
topdowns = topgenes[topgenes[, 'logFC'] < -0.5,]
dim(topdowns)

IDs.up = rownames(topups)
IDs.down = rownames(topdowns)

### DEresults

DEresults = decideTests(data.fit.eb, method = 'global', adjust.method = 'BH', p.value = 0.05, lfc = 0.5)

data.matrix.up = data.matrix[(rownames(topups)),]

sampleNames = vector()
featureNames = vector()
heatlogs = vector()
for (i in 1:10)
{
  sampleNames = c(sampleNames, rep(ph@data[i, 1], dim(topups)[1]))
  featureNames = c(featureNames, row.names(data.matrix.up[1:dim(topups)[1],]))
  heatlogs = c(heatlogs, data.matrix.up[1:dim(topups)[1], i])
}

heatdata = data.frame(norm_logint = heatlogs,
                      sampleName = sampleNames,
                      featureName = featureNames)

ggplot(heatdata, aes(sampleName, featureName))+
  geom_tile(aes(fill = norm_logint))+
  scale_fill_gradient(low = 'green', high = 'red')


vennDiagram(DEresults)  
