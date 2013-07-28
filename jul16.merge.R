load("../jun27.GSE2180.gsplom.RData")
load("../jul1.NAfiltered.tf.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/nec.net.plots.R")

gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")

# k=420
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k420.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z420 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=420)
zz420 <- split(names(z420),z420)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z420.sig <- collapse.cls(BOOL.TF.F, z420, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z420.sig$MIX.SIGN + R.z420.sig$MIX.DIR
CRIT <- (log2(R.z420.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.420 <- clust.names.to.idx(zz420,syms)
plot.rects.coords(clust.coords.420)
map.plot.rects(clust.coords.420, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.420 <- get.coh.M.score(R.z420.sig, min.dcor=0.36)
S.420$tri.wavg # [1] 0.975891
S.420$edge.flaws # 98
S.420$edge.all.n # 87990
S.420$edge.n # 47638
S.420$loser.clusters # 0

# 1) check for "loser clusers": DONE

# 2) pick out gold genes, create BRCA-like network

# 2.1) get clusters that contain a gold gene
qq <- sapply(zz420, function(z) any(z %in% gold.genes))
zz420[qq]
## $`23`
## [1] "chd-7"      "K05F1.5"    "hil-7"      "elt-1"      "ldb-1"     
## [6] "hil-6"      "Y82E9BR.17" "ceh-20"    
## $`57`
## [1] "unc-120" "cdc-14"  "sma-9"  
## $`70`
## [1] "nhr-69" "zip-8"  "nob-1"  "zip-12"
## $`92`
## [1] "cwn-1"
## $`107`
## [1] "hlh-1"  "somi-1" "nhr-25" "ztf-19"
## $`154`
## [1] "zip-2"  "tbx-8"  "hlh-14" "hnd-1" 
## $`166`
## [1] "mab-21"         "CELE_Y41D4B.26" "CELE_M03D4.4"   "nhr-152"       
## $`169`
## [1] "pal-1"
## $`206`
## [1] "tbx-9"
## $`343`
## [1] "elt-3"
## $`374`
## [1] "scrt-1"

# Save relevant files in correct formats for network generation
dirpath <- "/Users/z/Dropbox/biostat/local_c.elegans/jul16.gold.clusts/"
# BOOL.TF.F, DCOR.TF.F, WEAK.TF.F
rownames(BOOL.TF.F) <- rownames(DCOR.TF.F)
colnames(BOOL.TF.F) <- rownames(DCOR.TF.F)
colnames(DCOR.TF.F) <- rownames(DCOR.TF.F)
rownames(WEAK.TF.F) <- rownames(DCOR.TF.F)
colnames(WEAK.TF.F) <- rownames(DCOR.TF.F)

write.table(BOOL.TF.F, file=paste0(dirpath,"BOOL.TF.F.tab"), sep="\t")
write.table(DCOR.TF.F, file=paste0(dirpath,"DCOR.TF.F.tab"), sep="\t")
write.table(WEAK.TF.F, file=paste0(dirpath,"WEAK.TF.F.tab"), sep="\t")

# 3) SUMMARIZE NETWORK
# IDEAS
# Merge branches that don't change any edges per sub-tree
# partition by Mutual Exclusion glyphs
