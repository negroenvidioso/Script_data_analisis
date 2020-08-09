# reanalisis
setwd('~/Documentos/tesis/microarrays/metavolcano_mirna_pclib/')

GSE112093 = read.delim('GSE112093_table_2.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE112093) <- GSE112093 == ''
GSE112093 = GSE112093[!is.na(GSE112093),]
GSE112093 = subset(GSE112093, grepl('hsa', GSE112093$miRNA_ID))
GSE112093 = data.frame(gene = GSE112093$miRNA_ID,
                       logFC = GSE112093$logFC,
                       CI.L = GSE112093$CI.L,
                       CI.R = GSE112093$CI.R,
                       P.Value = GSE112093$P.Value,
                       adj.P.Val = GSE112093$adj.P.Val)
GSE112093$gene = gsub('[*]', '', GSE112093$gene)
GSE112093 = GSE112093[!duplicated(GSE112093$gene),]
names(GSE112093) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE112093$logFC = 2^GSE112093$logFC

GSE118578 = read.delim('GSE118578_table.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE118578) <- GSE118578 == ''
GSE118578 = GSE118578[!is.na(GSE118578),]
GSE118578 = subset(GSE118578, grepl('hsa', GSE118578$gene))
GSE118578 = data.frame(gene = GSE118578$gene,
                       logFC = GSE118578$logFC,
                       CI.L = GSE118578$CI.L,
                       CI.R = GSE118578$CI.R,
                       P.Value = GSE118578$P.Value,
                       adj.P.Val = GSE118578$adj.P.Val)
GSE118578$gene = gsub('[*]', '', GSE118578$gene)
GSE118578 = GSE118578[!duplicated(GSE118578$gene),]
names(GSE118578) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE118578$logFC = 2^GSE118578$logFC

GSE108707	= read.delim('GSE108707_table_2.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE108707) <- GSE108707 == ''
GSE108707 = GSE108707[!is.na(GSE108707),]
GSE108707 = subset(GSE108707, grepl('hsa', GSE108707$miRNA_ID))
GSE108707 = data.frame(gene = GSE108707$miRNA_ID,
                       logFC = GSE108707$logFC,
                       CI.L = GSE108707$CI.L,
                       CI.R = GSE108707$CI.R,
                       P.Value = GSE108707$P.Value,
                       adj.P.Val = GSE108707$adj.P.Val)
GSE108707$gene = gsub('[*]', '', GSE108707$gene)
GSE108707 = GSE108707[!duplicated(GSE108707$gene),]
names(GSE108707) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE108707$logFC = 2^GSE108707$logFC

GSE67597	= read.delim('GSE67597_table_2.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE67597) <- GSE67597 == ''
GSE67597 = GSE67597[!is.na(GSE67597),]
GSE67597 = subset(GSE67597, grepl('hsa', GSE67597$miRNA_ID))
GSE67597 = data.frame(gene = GSE67597$miRNA_ID,
                      logFC = GSE67597$logFC,
                      CI.L = GSE67597$CI.L,
                      CI.R = GSE67597$CI.R,
                      P.Value = GSE67597$P.Value,
                      adj.P.Val = GSE67597$adj.P.Val)
GSE67597$gene = gsub('[*]', '', GSE67597$gene)
GSE67597 = GSE67597[!duplicated(GSE67597$gene),]
names(GSE67597) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE67597$logFC = 2^GSE67597$logFC

GSE75670	= read.delim('GSE75670_table.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE75670) <- GSE75670 == ''
GSE75670 = GSE75670[!is.na(GSE75670),]
GSE75670 = subset(GSE75670, grepl('hsa', GSE75670$gene))
GSE75670 = data.frame(gene = GSE75670$gene,
                      logFC = GSE75670$logFC,
                      CI.L = GSE75670$CI.L,
                      CI.R = GSE75670$CI.R,
                      P.Value = GSE75670$P.Value,
                      adj.P.Val = GSE75670$adj.P.Val)
GSE75670$gene = gsub('[*]', '', GSE75670$gene)
GSE75670 = GSE75670[!duplicated(GSE75670$gene),]
names(GSE75670) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE75670$logFC = 2^GSE75670$logFC

GSE56914	= read.delim('GSE56914_table.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE56914) <- GSE56914 == ''
GSE56914 = GSE56914[!is.na(GSE56914),]
GSE56914 = subset(GSE56914, grepl('hsa', GSE56914$miRNA_ID))
GSE56914 = data.frame(gene = GSE56914$miRNA_ID,
                      logFC = GSE56914$logFC,
                      CI.L = GSE56914$CI.L,
                      CI.R = GSE56914$CI.R,
                      P.Value = GSE56914$P.Value,
                      adj.P.Val = GSE56914$adj.P.Val)
GSE56914$gene = gsub('[*]', '', GSE56914$gene)
GSE56914 = GSE56914[!duplicated(GSE56914$gene),]
names(GSE56914) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE56914$logFC = 2^GSE56914$logFC

GSE21281	= read.delim('GSE21281_table.txt', header = T, sep = '\t', stringsAsFactors = F)
is.na(GSE21281) <- GSE21281 == ''
GSE21281 = GSE21281[!is.na(GSE21281),]
GSE21281 = subset(GSE21281, grepl('hsa', GSE21281$miRNA_ID))
GSE21281 = data.frame(gene = GSE21281$miRNA_ID,
                      logFC = GSE21281$logFC,
                      CI.L = GSE21281$CI.L,
                      CI.R = GSE21281$CI.R,
                      P.Value = GSE21281$P.Value,
                      adj.P.Val = GSE21281$adj.P.Val)
GSE21281$gene = gsub('[*]', '', GSE21281$gene)
GSE21281 = GSE21281[!duplicated(GSE21281$gene),]
names(GSE21281) = c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")
#GSE56914$logFC = 2^GSE56914$logFC

setwd('~/Documentos/tesis/microarrays/metavolcano_mirna_pclib/reanalisis/')

# todas las patologias
library(knitr)
library(MetaVolcanoR)
name = c("GSE112093", "GSE118578", "GSE108707", "GSE67597", "GSE75670", "GSE56914", "GSE21281")
data = list(GSE112093, GSE118578, GSE108707, GSE67597, GSE75670, GSE56914, GSE21281)
names(data) = name

meta_degs_vote_all <- votecount_mv(data,
                                   pcriteria="P.Value",
                                   foldchangecol="logFC",
                                   genenamecol="gene",
                                   geneidcol=NULL,
                                   pvalue=0.05,
                                   foldchange=0,
                                   metathr=0.05,
                                   collaps=FALSE,
                                   jobname="MetaVolcano_all_0.05",
                                   outputfolder=".",
                                   draw='HTML')

meta_degs_comb_all <- combining_mv(data,
                                   pcriteria='P.Value',
                                   foldchangecol='logFC',
                                   genenamecol='gene',
                                   geneidcol=NULL,
                                   metafc='Mean',
                                   metathr=0.05,
                                   collaps=TRUE,
                                   jobname="MetaVolcano_all_0.05_comb",
                                   outputfolder=".",
                                   draw='HTML')

# solo arterial
# hipertension arterial
# GSE112093
# GSE118578
# GSE75670
arterial = list(GSE112093, GSE118578, GSE75670)
nombre = c('GSE112093', 'GSE118578', 'GSE75670')
names(arterial) = nombre

setwd('/home/negro/Documentos/tesis/microarrays/metavolcano_mirna_pclib/reanalisis/arterial')

meta_degs_vote_all <- votecount_mv(arterial,
                                   pcriteria="P.Value",
                                   foldchangecol="logFC",
                                   genenamecol="gene",
                                   geneidcol=NULL,
                                   pvalue=0.05,
                                   foldchange=0,
                                   metathr=0.05,
                                   collaps=FALSE,
                                   jobname="MetaVolcano_arterial_0.05",
                                   outputfolder=".",
                                   draw='HTML')

meta_degs_comb_all <- combining_mv(arterial,
                                   pcriteria='P.Value',
                                   foldchangecol='logFC',
                                   genenamecol='gene',
                                   geneidcol=NULL,
                                   metafc='Mean',
                                   metathr=0.05,
                                   collaps=TRUE,
                                   jobname="MetaVolcano_arterial_0.05_comb",
                                   outputfolder=".",
                                   draw='HTML')

arterial_comb = meta_degs_comb_all@metaresult
arterial_comb = subset(arterial_comb, arterial_comb$metap <= 0.05)
arterial_comb = subset(arterial_comb, abs(arterial_comb$metafc) >= 0.5)

arterial_vote = meta_degs_vote_all@metaresult
arterial_vote = subset(arterial_vote, !arterial_vote$degvcount == '1.Unperturbed')

write.table(arterial_comb$gene,
            file = 'arterial_comb.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

write.table(arterial_vote$gene,
            file = 'arterial_vote.txt',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F)

intersect(arterial_comb$gene, arterial_vote$gene)

# hipertension arterial pulmonar
# GSE108707
# GSE21281
# GSE67597
setwd('/home/negro/Documentos/tesis/microarrays/metavolcano_mirna_pclib/reanalisis/arterial_pulmonar/')
arterial_pulmonar = list(GSE108707, GSE21281, GSE67597)
nombre = c('GSE108707', 'GSE21281', 'GSE67597')
names(arterial_pulmonar) = nombre

meta_degs_vote_all <- votecount_mv(arterial_pulmonar,
                                   pcriteria="P.Value",
                                   foldchangecol="logFC",
                                   genenamecol="gene",
                                   geneidcol=NULL,
                                   pvalue=0.05,
                                   foldchange=0,
                                   metathr=0.05,
                                   collaps=FALSE,
                                   jobname="MetaVolcano_arterial_pulmonar_0.05",
                                   outputfolder=".",
                                   draw='HTML')

meta_degs_comb_all <- combining_mv(arterial_pulmonar,
                                   pcriteria='P.Value',
                                   foldchangecol='logFC',
                                   genenamecol='gene',
                                   geneidcol=NULL,
                                   metafc='Mean',
                                   metathr=0.05,
                                   collaps=TRUE,
                                   jobname="MetaVolcano_arterial_pulmonar_0.05_comb",
                                   outputfolder=".",
                                   draw='HTML')