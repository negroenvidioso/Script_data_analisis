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


mirna_vote = meta_degs_vote_all@metaresult
mirna_vote = subset(mirna_vote, !mirna_vote$degvcount == '1.Unperturbed')

mirna_comb = meta_degs_comb_all@metaresult
mirna_comb = subset(mirna_comb, mirna_comb$metap <= 0.05)
mirna_comb = subset(mirna_comb, abs(mirna_comb$metafc) >= 0.5)

write.table(mirna_vote$gene,
            file = 'mirna_vote_all.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

write.table(mirna_comb$gene,
            file = 'mirna_comb_all.txt',
            quote = F,
            row.names = F,
            col.names = F)


write.table(intersect(mirna_comb$gene, mirna_vote$gene),
            file = '32_mirna_compartidos_correccion.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

mirna_total = c(mirna_vote$gene, mirna_comb$gene)
mirna_total = unique(mirna_total)

write.table(mirna_total,
            file = 'mirna_total.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

GSE112093 = subset(GSE112093, GSE112093$P.Value <= 0.05)
GSE112093 = subset(GSE112093, abs(GSE112093$logFC) >= 0.5)
GSE118578 = subset(GSE118578, GSE118578$P.Value <= 0.05)
GSE118578 = subset(GSE118578, abs(GSE118578$logFC) >= 0.5)
GSE75670 = subset(GSE75670, GSE75670$P.Value <= 0.05)
GSE75670 = subset(GSE75670, abs(GSE75670$logFC) >= 0.5)
GSE108707 = subset(GSE108707, GSE108707$P.Value <= 0.05)
GSE108707 = subset(GSE108707, abs(GSE108707$logFC) >= 0.5)
GSE21281 = subset(GSE21281, GSE21281$P.Value <= 0.05)
GSE21281 = subset(GSE21281, abs(GSE21281$logFC) >= 0.5)
GSE67597 = subset(GSE67597, GSE67597$P.Value <= 0.05)
GSE67597 = subset(GSE67597, abs(GSE67597$logFC) >= 0.5)
GSE56914 = subset(GSE56914, GSE56914$P.Value <= 0.05)
GSE56914 = subset(GSE56914, abs(GSE56914$logFC) >= 0.5)

arterial = c(GSE112093$gene, GSE118578$gene, GSE75670$gene)
arterial = unique(arterial)

arterial_pulmonar = c(GSE108707$gene, GSE21281$gene, GSE67597$gene)
arterial_pulmonar = unique(arterial_pulmonar)

pulmonar = c(GSE56914$gene)
pulmonar = unique(pulmonar)

write.table(arterial,
            file = 'mirnas_arterial_correccion.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

write.table(arterial_pulmonar,
            file = 'mirnas_arterial_pulmonar_correccion.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

write.table(pulmonar,
            file = 'mirnas_pulmonar_correccion.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)
    