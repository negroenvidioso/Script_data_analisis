# GSE75670

library(limma)
library(affy)

setwd('D:/Documentos/Magister/tesis/pc_lib/Documentos/tesis/microarrays/GSE75670')
targets = readTargets('target.txt')

setwd('D:/Documentos/Magister/tesis/pc_lib/Documentos/tesis/microarrays/GSE75670_RAW/')
raw = read.maimages(targets, source = 'agilent')

setwd("D:/Documentos/Magister/tesis/pc_lib/Documentos/tesis/microarrays/GSE75670")

raw_BGcorrected = backgroundCorrect(raw, method = 'normexp')

raw_BGandNormalized = normalizeBetweenArrays(raw_BGcorrected, method = 'quantile')
raw_aver = avereps(raw_BGandNormalized, ID = raw_BGandNormalized$genes$ProbeName)

### box plot
### nto normalizated
boxplot(log(as.matrix(raw_BGcorrected)), las = 2, ylab = 'Log (intensity)')
### normalizated
boxplot(as.matrix(raw_BGandNormalized), las = 2, ylab = 'log2 (intensity)')

affy::mva.pairs(as.matrix(raw))
affy::mva.pairs(as.matrix(raw_BGcorrected))
affy::mva.pairs(as.matrix(raw_BGandNormalized))                

f = factor(targets$Condicion)
design = model.matrix(~f)

desing = cbind(hipertension = c(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),
               control = c(1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0))

fit = lmFit(raw_aver, desing)
contrastMatrix = makeContrasts(hipertension - control, levels = desing)

fit2 = contrasts.fit(fit, contrastMatrix)
fit2 = eBayes(fit2, 0.01)

topTable(fit2, coef = 'hipertension - control', confint = T, number = 250)  # confint = T

signif = topTable(fit2, coef = 'hipertension - control', number = nrow(fit2), confint = T)

table = data.frame(gene = signif$GeneName,
                   logFC = signif$logFC,
                   CI.L = signif$CI.L,
                   CI.R = signif$CI.R,
                   AveExpr = signif$AveExpr,
                   t = signif$t,
                   P.Value = signif$P.Value,
                   adj.P.Val = signif$adj.P.Val,
                   B = signif$B)

write.table(table,
            file = 'GSE75670_table.txt',
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

upregulated = signif[which(signif[,11] > 0),]
downregulated = signif[which(signif[,11] < 0),]
