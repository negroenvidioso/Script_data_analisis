setwd("~/Documentos/tesis/")


lnc_predicted = read.delim("Human_lncRNAs.bed12",
                           header = F)

View(lnc_predicted)

x = lnc_predicted$V4
View(x)

ensemblid = subset(lnc_predicted, !grepl("[.]", lnc_predicted$V4))
otras = subset(lnc_predicted, grepl("[.]", lnc_predicted$V4))

y = unique(ensemblid$V4)
x = subset(y, grepl("ENST", y))
no_x = subset(y, !grepl("ENST", y))

write.table(x,
            file = "ENST_id.txt",
            append = F,
            quote = F,
            sep = " \t ",
            row.names = F,
            col.names = F)
  
z = gsub("[_PAR_Y]", "", otras)

write.table(z,
            file = "gencode_lncrna_predicted.txt",
            append = F,
            quote = F,
            sep = " \t ",
            row.names = F,
            col.names = F)

trans = read.delim("gencode_lncrna_predicted.txt",
                   header = F,
                   sep = ".")
write.table(unique(trans$V1),
            file = "fix_gencode_lncrna_predicted.txt",
            append = F,
            quote = F,
            sep = " \t ",
            row.names = F,
            col.names = F)
