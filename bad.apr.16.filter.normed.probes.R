## THIS IS BAD. DON'T USE IT
## ------------------------------
## load REFSEQ filtered, SCAN normalized data.
## Re-construct expression using past results.
## ------------------------------
library(Biobase)
library(celeganscerefseq.db)
load("../GSE2180.BRAINARRAY.REFSEQ.SCAN.RData")
load("../GSE2180.BRAINARRAY.REFSEQ.UPC.RData")
# all(rownames(exprs(GSE2180.SCAN))==rownames(exprs(GSE2180.UPC)))
# TRUE
load("../celegans.apr8.expr.RData")
# E.expr, E.expr.trans, E.gold
# ls("package:celeganscerefseq.db")
gsms.wtms <- colnames(exprs(E.gold))       # previously filtered GSM names in wt/ms set

trans.factors <- unique(as.vector(unlist(read.table("celegans.transcription.factors.txt", stringsAsFactors=F))))
# 606 +3 = 609 transcription factor names
# from http://www.macwormlab.net/ntfdb/index.php?class=1&state=1 (plus mab-21, cwn-1, scrt-1)
gold.genes <- c('pal-1', 'tbx-8', 'tbx-9', 'elt-1', 'hnd-1', 'scrt-1', 'cwn-1', 'unc-120', 'hlh-1', 'nob-1', 'elt-3', 'nhr-25', 'mab-21', 'lin-26', 'vab-7')
# lin-26 probes
lin26 <- c('174037_at', '190309_at') # use 190309_at per
# http://www.wormbook.org/chapters/www_transcriptionalregulation/transcriptionalregulation.html

pal1 <- c('174043_at', '193341_at', '193342_s_at')
elt3 <- c('175801_at') # 193640_s_at ? 

ids <- rownames(exprs(GSE2180.SCAN))
genes <- mget(ids, celeganscerefseqSYMBOL, ifnotfound = NA)
featureData(GSE2180.SCAN)$gene <- genes
featureData(GSE2180.UPC)$gene <- genes

# 0) Inspect gold probe expressions
# ------------------------------
match(gold.genes, featureData(E)$gene) # all gold genes are in platform
gold.qq <- featureData(E)$gene %in% gold.genes
summary(as.factor(unlist(featureData(E[gold.qq,])$gene)))
  cwn-1   elt-1   elt-3   hlh-1   hnd-1  lin-26  mab-21  nhr-25   nob-1   pal-1 
      2       1       1       1       1       1       1       1       1       1 
 scrt-1   tbx-8   tbx-9 unc-120   vab-7 
      1       1       1       1       1
apply(exprs(P[gold.qq,]) >= 0.5, 1, sum)
# WARNING: no elt-3???


# 1) Select samples, add back pheno data
# ------------------------------
colnames(exprs(GSE2180.SCAN)) <- sub(".CEL.gz","",colnames(exprs(GSE2180.SCAN)))
colnames(exprs(GSE2180.UPC)) <- sub(".CEL.gz","",colnames(exprs(GSE2180.UPC)))
all(colnames(exprs(GSE2180.UPC)) == colnames(exprs(GSE2180.SCAN))) # TRUE
gsms <- colnames(exprs(GSE2180.SCAN))
qq.col <- match(gsms.wtms, gsms)
E <- GSE2180.SCAN[,qq.col]
P <- GSE2180.UPC[,qq.col]
all(rownames(pData(phenoData(E.expr))) == colnames(rownames(E))) # TRUE
phenoData(E) <- phenoData(E.expr)
phenoData(P) <- phenoData(E.expr)

# 2) Select expressed probes (this might be too conservative)
# ------------------------------
upc.num.over.5 <- apply(exprs(P) > 0.5,1,sum)
qq.row <- upc.num.over.5>=4
sum(qq.row)
#[1] 5191 of 17623 features
E.exp <- E[qq.row,]
length(unique(featureData(E.exp)$gene))
#[1] 5053 unique gene symbols (-1 for NA)
sum(is.na(featureData(E.exp)$gene))
#[1] 38 without gene symbol

# 3) Select Expressed Gold List
# ------------------------------
sum(featureData(E)$gene %in% gold.genes)

match(gold.genes, featureData(E.exp)$gene)


# 4) get expressed transcription factor list

