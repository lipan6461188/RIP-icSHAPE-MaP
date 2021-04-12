
library(edgeR)

############################
#######  Normal
############################

outFn <- '/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/4.limmaDE/limma.txt'
genecountFn <- "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/genecount.txt"

counts <- read.delim(genecountFn, row.names = 1)
head(counts)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

group <- c('INPUT','INPUT','INPUT','RIP','RIP')
plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupRIP - groupINPUT, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = outFn, row.names = F, sep = "\t", quote = F)


"""
总结所有的基因类型：
awk 'NR==1{print $0}$2>1&&$6<0.05{print $0}' limma.txt | awk '{split($1,arr,"_"); print arr[1]}' | sort | uniq -c
查看相关的基因：
awk 'NR==1{print $0}$2>1&&$6<0.05{print $0}' limma.txt | expand -t 50 | les
"""

