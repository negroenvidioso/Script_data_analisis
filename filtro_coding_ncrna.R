setwd('/home/ignacio/Documentos/tesis_allan')

bedcoding = read.delim('codings_Human.bed', header = F, sep = '\t')

bedcoding = subset(bedcoding, selec = c('V4'))

bedcoding = gsub('[_PAR_Y]', '', bedcoding$V4)
bedcoding = unique(bedcoding)

write.table(bedcoding,
            file = 'listados_est_coding.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = F)

anotacion = read.delim('coding_anotation.txt', header = T, sep = '\t')
length(unique(anotacion$Transcript.stable.ID.version))

DEGsdeseq2$gene = row.names(DEGsdeseq2)

coding = DEGsdeseq2[DEGsdeseq2$gene %in% anotacion$Gene.stable.ID.version,]
up_coding = subset(coding, coding$log2FoldChange > 0)
down_coding = subset(coding, coding$log2FoldChange < 0)
ncrna = DEGsdeseq2[!DEGsdeseq2$gene %in% anotacion$Gene.stable.ID.version,]
up_ncrna = subset(ncrna, coding$log2FoldChange > 0)
down_ncrna = subset(ncrna, coding$log2FoldChange < 0)
