setwd("D:/Documentos/Magister/tesis/plataformas/")

GPL10558 = read.delim("GPL10558-50081.txt", header = T, sep = "\t")


library(seqinr)
GPL10558_data = data.frame(GPL10558$ID, GPL10558$SEQUENCE)

write.fasta(as.list(GPL10558$SEQUENCE),GPL10558$ID,"fasta_proves_GPL10558.txt")

GPL21572 = read.delim("GPL21572-124634.txt", header = T, sep = "\t")

write.fasta(as.list(GPL21572$Sequence),GPL21572$Accession_ID,"fasta_proves_GPL21572.txt")

GPL4133 = read.delim("anotacion GPL4133-12599.txt", header = T, sep = "\t")
