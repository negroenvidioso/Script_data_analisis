setwd('F:/RNA-seq')

GSE75360 = read.delim('expresion_diferencial_anotado.txt', header = T, sep = '\t')
PRJNA480270 = read.delim('DESeq2_rna-seq_PRJNA480270.txt', header = T, sep = '\t')
write.table(PRJNA480270$gene,
            file = 'gene_PRJNA480270.txt',
            sep = '\t',
            quote = F,
            col.names = F,
            row.names = F)

gene_PRJNA480270 = read.delim('gene_PRJNA480270.txt', sep = '.', header = F)

PRJNA480270_fix = cbind(PRJNA480270, gene_PRJNA480270$V1)
PRJNA480270_fix = subset(PRJNA480270_fix, select = c('gene_PRJNA480270$V1',
                                                     'log2FoldChange',
                                                     'pvalue',
                                                     'padj'))
datos_PRJNA480270 = PRJNA480270_fix[!is.na(PRJNA480270_fix$log2FoldChange),]
names(datos_PRJNA480270) = c('gene', 'logFC', 'P.Value', 'adj.P.Val')
names(GSE75360) = c('gene', 'logFC', 'P.Value', 'adj.P.Val')

# metavolcano
library(knitr)
library(MetaVolcanoR)
name = c('PRJNA480270', 'GSE75360')
data = list(datos_PRJNA480270, GSE75360)
names(data) = name

meta_degs_vote_all <- votecount_mv(data,
                                   pcriteria="P.Value",
                                   foldchangecol="logFC",
                                   genenamecol="gene",
                                   geneidcol=NULL,
                                   pvalue=0.05,
                                   foldchange=0,
                                   metathr=0.05,
                                   collaps=TRUE,
                                   jobname="MetaVolcano_all_0.05_vote_RNA-Seq",
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
                                   jobname="MetaVolcano_all_0.05_comb-RNA-Seq",
                                   outputfolder=".",
                                   draw='HTML')

comb_rna_seq = meta_degs_comb_all@metaresult
comb_rna_seq = subset(comb_rna_seq, comb_rna_seq$metap <= 0.05)
comb_rna_seq = subset(comb_rna_seq, abs(comb_rna_seq$metafc) >= 0.5)

vote_rna_seq = meta_degs_vote_all@metaresult
vote_rna_seq = subset(vote_rna_seq, vote_rna_seq$degvcount != '1.Unperturbed')

write.table(comb_rna_seq,
            file = 'gene_comb_rna_seq.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = F)

write.table(vote_rna_seq,
            file = 'gene_vote_rna_seq.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = F)

intersect = intersect(comb_rna_seq$gene, vote_rna_seq$gene)

hsa_miR_106b_5p = read.delim('hsa-miR-106b-5p_lncrna_interaction.txt', header = F, sep = '\t')

hsa_miR_106b_5p_lncrna = comb_rna_seq[comb_rna_seq$gene %in% hsa_miR_106b_5p$V1,]
hsa_miR_106b_5p_lncrna_gse = vote_rna_seq[vote_rna_seq$gene %in% hsa_miR_106b_5p$V1,]
