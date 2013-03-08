library(celegansceentrezg.db)
library(celegansceentrezgprobe) 
library(celegans.db)
library(sva)
library(lumi)

GPL200 <- read.table("../GPL200-2880.txt", sep="\t", quote="", comment="", header=TRUE, row.names=1)

genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
"lin-26"
"unc-120"
load("../GSE2180.clean.RData")
load("../GSE9665.clean.RData")

ids <- rownames(exprs(GSE2180))
GSE2180.syms <- mget(ids, celegansceentrezgSYMBOL, ifnotfound=NA)
GSE2180.syms2 <- mget(ids, celegansSYMBOL, ifnotfound=NA)
match(genes, GSE2180.syms)

x <- celegansceentrezgSYMBOL
xx <- as.list(x[mappedkeys(x)])
match(genes, xx)

190309_at
GPL200[rownames(GPL200)=="190309_at",]$Gene.Symbol
[1] lin-26 /// lir-1 /// WBGene00003012

which(celegansceentrezgprobe$Probe.Set.Name == 187454_at)
Error: unexpected input in "which(celegansceentrezgprobe$Probe.Set.Name == 187454_"
> which(celegansceentrezgprobe$Probe.Set.Name == "187454_at")
[1] 153710 153711 153712 153713 153714 153715 153716 153717
> qq <- celegansceentrezgprobe$Probe.Set.Name == "187454_at"
> celegansceentrezgprobe$sequence[qq]
[1] "GAACTCTCGACTTTAATCCGTGGAT" "GTTAGCCAGCAACTAGCAGCACAGC"
[3] "GCAGCACAGCTATCACAGAACGGTA" "AAAATCGACCTGGTGCTCCAAATCA"
[5] "CTCCAAATCAAGTTCCGCTTCTGAA" "TCCTCGGCTTTCATGGTTGTACCCG"
[7] "TTGTACCCGTACATGAGCAAGTCCC" "GCGTCAGCAAGAAGCCTCGGACAGC"
> 

mget(c("190309_at"), celegansSYMBOL, ifnotfound=NA)
190309_at


[1] "GAACTCTCGACTTTAATCCGTGGAT" "GTTAGCCAGCAACTAGCAGCACAGC"
[3] "GCAGCACAGCTATCACAGAACGGTA" "AAAATCGACCTGGTGCTCCAAATCA"
[5] "CTCCAAATCAAGTTCCGCTTCTGAA" "TCCTCGGCTTTCATGGTTGTACCCG"
[7] "TTGTACCCGTACATGAGCAAGTCCC" "GCGTCAGCAAGAAGCCTCGGACAGC"

atggccaatc agtttatcat tctagctata gtaatagtga ttactttgct aaattttatg
ctctggagtg ggaatagtta cggaagcaca caactaactg agttctacac atcactattg
gacagaaatt caaacaataa cagtgtggac caacaattgg ctcaaatccg gagccaattg
gatgctgaaa tggagatttt acggaatttc gaatccaaat ttggtcaggg agccggagac
agttttgagc aacactttga gaaaatgaaa gcattctcaa tagcccaggc acctcttcgc
aagagaactc tagcccagat cttcaagcca acagagtttt acttgctcta tggagccatg
ggaccagagg tcttttgtcc tgaaaaagta cgcatcggaa cagttggaga cggcggaaaa
tgggtttgca atccgtggaa gtgcccgaac aattctgtaa tgttctcgct aggtctcaac
aattggataa ctttcgagga ggaatggcag aagcttacag ataatcgaaa tattttatat
ggttttgatg cggctgacca aaacgacaga acccgtcaaa cctactcaaa aatccgtgga
acctcgaaaa aagctctaat ctccgtggca accgatccag caactagcaa atacacaatt
gatgatcttg ccaagcactt taatgtttcc aatattgaaa ttctaaaaat cgatattgaa
ggcgccgagt tgacgtgttt aattccattc cttgagaaat acgaagtatg tcaaatctac
ttggagctcc acggcggagc acaggaacac gcgaaattgc tcagggaaat tggtcactta
aattatcgat tgttttcata tgaggttaat ggttttgagt tgaaagcctg tgaatatagc
tttatacatg aaaagtgtgt ggagaagtac ggtggaatga gaattgctaa ttatttggat
tataggaatt aaaatttgta cgtttttgac ttgttgaatt gattgtcggg tgatatttat
ttattgattg attacaataa agttattaaa attt

GPL200$RefSeq.Transcript.ID[grep("lin-26", GPL200$Gene.Symbol)]
[1] NM_001026921 /// NM_001026922 /// NM_001026923 /// NM_001026924 /// NM_001026925 /// NM_001047300
[2] NM_001026921 /// NM_001026922 /// NM_001026923 /// NM_001026924 /// NM_001026925 /// NM_001047300

GPL200$Sequence.Type[grep("lin-26", GPL200$Gene.Symbol)]
[1] Exemplar sequence Exemplar sequence

GPL200$Representative.Public.ID[grep("lin-26", GPL200$Gene.Symbol)]
[1] 551223  F18A1.2

GPL200$Sequence.Source[grep("lin-26", GPL200$Gene.Symbol)]
[1] GenBank                         Affymetrix Proprietary Database

GPL200$Target.Description[grep("lin-26", GPL200$Gene.Symbol)]
[1] g551223 /REP_DB=GenBank Identifier /CHR=2 /FEA=Genomic Cluster
[2] F18A1.2 /REP_DB=WormBase Gene ID /WP=CE27972 /GEN=lin-26 /GB=AAK67251.1 /SUBMIT=ST.LOUIS /CHR=2 /FEA=Sanger Annotation /DEF=transcription factor

GPL200$Gene.Symbol[grep("lin-26", GPL200$Gene.Symbol)]
[1] lin-26 /// lir-1 /// WBGene00003012 lin-26 /// lir-1 /// WBGene00003012

rownames(GPL200)[grep("lin-26", GPL200$Gene.Symbol)]
[1] "174037_at" "190309_at"

GPL200$Gene.Title[grep("lin-26", GPL200$Gene.Symbol)]
[1] abnormal cell LINeage /// LIn-26 Related /// locus:lin-26
[2] abnormal cell LINeage /// LIn-26 Related /// locus:lin-26

> GPL200$Gene.Symbol[grep("lir-1", GPL200$Gene.Symbol)]
[1] lin-26 /// lir-1 /// WBGene00003012 lir-1                              
[3] lin-26 /// lir-1 /// WBGene00003012
