library(diffHic)

library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")


setwd ("/Volumes/Joanna_HD/bams_for_diffHiC/")
bamFiles <- list.files("/Volumes/Joanna_HD/bams_for_diffHiC/")
print(bamFiles)

#Date: 12.04.2017

## bam files for diffHiC must be name sorted - sort files using: samtools sort -n - done on cluster using gi/samtools/1.2 module


### cutGenome with BSgenome object has different chromosomes names, use fasta from hicup as a string factor
#hs.farg <- cutGenome(BSgenome.Hsapiens.UCSC.hg38, "CCATGG", overhang=4L) ## NcoI enzyme
#hs.farg
#hs.param <- pairParam(hs.farg)
#sessionInfo()

hs.frag <- cutGenome("/Volumes/Joanna_HD/Annotations/hg38/fasta_hicup/WholeGenomeFasta/genome.fa", "CCATGG", overhang=4L) ## NcoI enzyme
hs.frag
as.character(runValue(seqnames(hs.frag)))
hs.param <- pairParam(hs.frag)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("diffHic", version = "3.8")



names(Rsamtools::scanBamHeader("TKCC_MCF7C_NcolII_1.asd.bam")[[1]]$targets)
as.character(runValue(seqnames(hs.frag)))

### preparePairs - matching mapped reads to restriction fragments - one sample at the time
min.inward <- 1000
min.outward <- 25000

##MCF7 T0 a1

diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_T0_a1_sorted.asd.bam", hs.param, file="MCF7_T0_a1.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_T0_a1.h5", hs.param, file.out="MCF7_T0_a1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#total   length   inward  outward retained 
#22070709  1439100   117407   670903 19853149 

##   MCF7 T0 a2
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_T0_a2_sorted.asd.bam", hs.param, file="MCF7_T0_a2.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_T0_a2.h5", hs.param, file.out="MCF7_T0_a2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#14369355  1659998    82946   452555 12180888  

## MCF7_p16_a1
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_p16_a1_sorted.asd.bam", hs.param,
                            file="MCF7_p16_a1.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_p16_a1.h5", hs.param, file.out="MCF7_p16_a1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#total   length   inward  outward retained 
#23856659   941040   116639   682828 22125785 


##MCF7 p16_a2
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_p16_a2_sorted.asd.bam", hs.param, file="MCF7_p16_a2.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_p16_a2.h5", hs.param, file.out="MCF7_p16_a2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#total   length   inward  outward retained 
#18453125  1031148    96646   543232 16789755

########## 23/04/2019

##MCF7 p32_a1
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_p32_a1_sorted.asd.bam", hs.param, file="MCF7_p32_a1.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_p32_a1.h5", hs.param, file.out="MCF7_p32_a1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#total   length   inward  outward retained 
#22112316  1257193   113819   649111 20101284 

##MCF7 p32_a2
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/20190418_Hi-C_MCF7Caldon_p32_a2_sorted.asd.bam", hs.param, file="MCF7_p32_a2.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_p32_a2.h5", hs.param, file.out="MCF7_p32_a2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#8227134   886468    52527   235084  7056833 


##Cardiff MCF71 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_MCF7C_NcolII_1_sorted.asd.bam", hs.param, file="MCF7_1.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_1.h5", hs.param, file.out="MCF7_1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#26967758  2244937   153061   793527 23789243 


##Cardiff MCF72 - sort by name
#diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/", hs.param, file="MCF7_2.h5", dedup=TRUE, minq=10)
#prunePairs("MCF7_2.h5", hs.param, file.out="MCF7_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

##Cardiff MCF73 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_MCF7C_NcolII_3_sorted.asd.bam", hs.param, file="MCF7_3.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_3.h5", hs.param, file.out="MCF7_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#86033049 12540783   650493  2825272 70079910 

##Cardiff TAMR1 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_TAMR_pooled_NcolII_1_sorted.asd.bam", hs.param, file="TAMR_1.h5", dedup=TRUE, minq=10)
prunePairs("TAMR_1.h5", hs.param, file.out="TAMR_1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

#   total   length   inward  outward retained 
#54525142  5597850   503209  1862828 46598069 

##Cardiff TAMR2 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_TAMR_pooled_NcolII_2_sorted.asd.bam", hs.param, file="TAMR_2.h5", dedup=TRUE, minq=10)
prunePairs("TAMR_2.h5", hs.param, file.out="TAMR_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

##Cardiff TAMR3 - sort by name
#diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/", hs.param, file="TAMR_3.h5", dedup=TRUE, minq=10)
#prunePairs("TAMR_3.h5", hs.param, file.out="TAMR_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

##Cardiff FASR1 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_FASR_pooled_NcolII_1_sorted.asd.bam", hs.param, file="FASR_1.h5", dedup=TRUE, minq=10)
prunePairs("FASR_1.h5", hs.param, file.out="FASR_1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#37848889  6063765   215322  1118206 30471807 

##Cardiff FASR2 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_FASR_pooled_NcolII_2_sorted.asd.bam", hs.param, file="FASR_2.h5", dedup=TRUE, minq=10)
prunePairs("FASR_2.h5", hs.param, file.out="FASR_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
#61428192  4903361   499345  2248241 53821168 


##Cardiff FASR3 - sort by name
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_FASR_pooled_NcolII_3_sorted.asd.bam", hs.param, file="FASR_3.h5", dedup=TRUE, minq=10)
prunePairs("FASR_3.h5", hs.param, file.out="FASR_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

#   total   length   inward  outward retained 
#57660319  5801323   404106  1874623 49623236 


### 23/04/2019



## counting read pais into interactions - do not use merged biological replicates - merge technical replicates and keep biol separarte
## technical replicates were merged at fastq level

# save datafiles as index files
anchor1.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
anchor2.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
dummy <- data.frame(anchor1.id, anchor2.id, other.data=as.integer(runif(100, 1, 100)))
savePairs(dummy, "example.h5", hs.param)
## for next step use MCF7 1 2 and 3 (3 biological replicates) and TAMR 1 2 and 3 (3 biological replicates)


#### resolution - 1Mb


### need to re-make h5 files from Cardiff data using a new fasta file as chromosomes are different

## this step has to be done at 1Mb bin - to much memory to do 100kb with filter=1 and not enough bin pairs to estamite threshold
input <- c("MCF7_1_trimmed.h5", "MCF7_3_trimmed.h5", "TAMR_1_trimmed.h5", "TAMR_2_trimmed.h5", "FASR_1_trimmed.h5", "FASR_2_trimmed.h5", "FASR_3_trimmed.h5", "MCF7_p16_a2_trimmed.h5", "MCF7_p16_a1_trimmed.h5", "MCF7_T0_a1_trimmed.h5", "MCF7_T0_a2_trimmed.h5", "MCF7_p32_a1_trimmed.h5", "MCF7_p32_a2_trimmed.h5" )

hs.param
bin.size <- 5e5
## filter - not sure what it it, but in LNCaP/PrEC data Aaron L used 10 --> 10
data <- squareCounts(input, hs.param, width=bin.size, filter=10)
data



head(anchors(data, type="first"))

head(anchors(data, type="second"))


## choosing a bin width
head(regions(data))

margin.data <- marginCounts(input, hs.param, width=bin.size)

# prior count??
require(edgeR)
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5) #prior count 5 in user guide, Aaron used 3 in our analysis
#MDS plot mfrow -plot 2 by 2 graphs together, mar - margins

## do-do MDS plots for corrected data - after filtering
par(mfrow=c(2,2), mar=c(5,4,2,2))
adjc <- cpm(asDGEList(data), log=TRUE, prior.count=5)
for (top in c(100, 500, 1000, 5000)) {out <- plotMDS(adjc, main=top, top=top, label=1:ncol(data))}

par(mfrow=c(2,2), mar=c(5,4,2,2))
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5)
for (top in c(100, 500, 1000, 5000)) {out <- plotMDS(adjc, main=top, top=top, label=1:ncol(data))}

smoothScatter(0.5*(adjc[,1]+adjc[,4]), adjc[,1]-adjc[,4], xlab="A", ylab="M", main="TAMR1 vs MCF71")
smoothScatter(0.5*(adjc[,2]+adjc[,5]), adjc[,2]-adjc[,5], xlab="A", ylab="M", main="TAMR2 vs MCF72")
smoothScatter(0.5*(adjc[,3]+adjc[,6]), adjc[,3]-adjc[,6], xlab="A", ylab="M", main="TAMR3 vs MCF73")

## check differences within relicates
smoothScatter(0.5*(adjc[,1]+adjc[,2]), adjc[,1]-adjc[,2], xlab="A", ylab="M", main="TAMR1 vs TAMR2")
smoothScatter(0.5*(adjc[,2]+adjc[,3]), adjc[,2]-adjc[,3], xlab="A", ylab="M", main="TAMR2 vs TAMR3")
smoothScatter(0.5*(adjc[,4]+adjc[,5]), adjc[,4]-adjc[,5], xlab="A", ylab="M", main="MCF71 vs MCF72")


#filtering out uninteresitng interactions
require(edgeR)

ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance")

count.keep <- ave.ab >= aveLogCPM(2, lib.size=mean(data$totals))
summary(count.keep)



dummy <- data[count.keep,] #interesing interactions
dummy


#remove low abundance interactions
direct <- filterDirect(data)

direct$threshold


direct.keep <- direct$abundance > log2(2) + direct$threshold ## threshold in user guide 5 
summary (direct.keep)


direct.keep2 <- direct.keep & count.keep # to ensure that retained bin pairs have large absolute counts
summary(direct.keep2)

## other filer? trended? OMIT - Aaron suggested to use direct filter only



###### peak calling in the interaction space
flank.width <- 5
enrichments <- enrichedPairs(data, flank=flank.width)
summary(enrichments)

#min.count 5 in user guide
peak.keep <- filterPeaks(data, enrichments, min.enrich=0.5, min.count=5, min.diag=2L)
sum(peak.keep)



## peak-calling diagnostics
neighborhood <- (2*flank.width + 1) * metadata(data)$width
boxed <- boxPairs(data[peak.keep], reference=neighborhood)

out <- tabulate(tabulate(boxed$indices[[1]]))
setNames(round(c(out[1:5], sum(out[5:length(out)]))/sum(out)*100, 1), c(1:5, ">5"))


peak.keep2 <-filterPeaks(data, enrichments, min.enrich=0, min.count=5, min.diag=2L) # no enrichment cut off to see if needed
boxed <- boxPairs(data[peak.keep2], reference=neighborhood)
out <- tabulate(tabulate(boxed$indices[[1]]))
setNames(round(c(out[1:5], sum(out[5:length(out)]))/sum(out)*100, 1), c(1:5, ">5"))
# larger proportions suggest that more agressive filtering on the enrichment score is required


## creating the "smaller data" object for testing sig interactions - same resolution here
bin.size <- 4e4
smaller.data <- squareCounts(input, hs.param, width=bin.size, filter=10)
direct <- filterDirect(smaller.data, reference=data)
direct$threshold

small.keep <- direct$abundances > direct$threshold + log2(5)
summary(small.keep)
    

# special cases summary
orginal <- data
data <- data[direct.keep2, ]
data
smaller.data <- data[direct.keep2, ]




## normalization

normfacs <- normOffsets(data)
normfacs

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- cpm(asDGEList(data), log=TRUE)
mval <- adj.counts[,1]-adj.counts[,4] #
smoothScatter(ab, mval, xlab="A", ylab="M", main="TAMR1 vs MCF71")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

mval <- adj.counts[,2]-adj.counts[,5] 
smoothScatter(ab, mval, xlab="A", ylab="M", main="TAMR2 vs MCF72")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

nb.off <- normOffsets(data, type="loess")

head(nb.off)


neardiag <- filterDiag(data, by.dist=1.5e6)
nb.off <- matrix(0, nrow=nrow(data), ncol=ncol(data))
nb.off[neardiag] <- normOffsets(data[neardiag,], type="loess")
nb.off[!neardiag] <- normOffsets(data[!neardiag,], type="loess")

adj.counts <- log2(assay(data) + 0.5) - nb.off/log(2)
mval <- adj.counts[,4]-adj.counts[,1] ## numbers of datasets to compare
smoothScatter(ab, mval, xlab="A", ylab="M", main="MCF71 vs TAMR1")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")


# iterative correction for interaction intensities
corrected <- correctedContact(orginal, winsor.high=0.02, ignore.low=0.02)
head(corrected$truth)

## "some NA will be present due to removal of low-abundance bins that do not exibit stable behavior during correction"
corrected$max

corrected <- correctedContact(orginal, average=FALSE)
anchor1.bias <- corrected$bias[anchors(orginal, type="first", id=TRUE),]
anchor2.bias <- corrected$bias[anchors(orginal, type="second", id=TRUE),]
iter.off <- log(anchor1.bias * anchor2.bias)

## correcting for CNVs - not neccessery here
cnv.offs <- normalizeCNV(data, margin.data)
head(cnv.offs)


#compare between one replicate of each - shouldn't be necceserry for the same cell line in both conditions?
count.files <- c("TKCC_MCF7C_NcolII_1_sorted.bam.bam.h5", "TKCC_TAMR_NcolII_1_sorted.bam.bam.h5")
cnvcheck.data <- squareCounts(count.files, hs.param, width=1e5)
cnvcheck.marg <- marginCounts(count.files, hs.param, width=1e5)
cnvcheck.data <- cnvcheck.data[aveLogCPM(asDGEList(cnvcheck.data)) > 0,]

matched <- matchMargins(cnvcheck.data, cnvcheck.marg)
m.adjc <- cpm(asDGEList(cnvcheck.marg), log=TRUE)
sum.madjc <- m.adjc[matched$anchor1,] + m.adjc[matched$anchor2,]
margin.lr <- sum.madjc[,1] - sum.madjc[,2]

before <- cpm(asDGEList(cnvcheck.data), log=TRUE)
after <- log2(assay(cnvcheck.data)+0.5) - normalizeCNV(cnvcheck.data, cnvcheck.marg, maxk=1000)/log(2)
par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.4)
smoothScatter(margin.lr, before[,1]-before[,2], ylim=c(-4,4), main="Before", xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")
smoothScatter(margin.lr, after[,1]-after[,2], ylim=c(-4,4), main="After", xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")

## 26/09/2016 - start                                                        
## modelling biological variability
design <- model.matrix(~factor(c("TAMR", "TAMR", "TAMR", "MCF7", "MCF7", "MCF7")))
colnames(design) <- c("TAMR", "MCF7")
design



y <- asDGEList(data)
y$offset <- nb.off

## estimate NB dispertion
y <- estimateDisp(y, design)

y$common.dispersion
#[1] 0.02244597
plotBCV(y)

## estimate QL dispersion
fit <- glmQLFit(y, design, robust=TRUE)

plotQLDisp(fit)

summary(fit$df.prior)

### testing for significant interactions

## using the quasi-likelihood F-test
result <- glmQLFTest(fit, coef=2)
topTags(result)

## multiplicity correction and the FDR

# cluster bins that are no more than 1bp apart - not bigger than 100kb for 40kb resolution
clustered.small <- clusterPairs(data, tol=1, upper=1e5)

length(clustered.small$interactions)
head(clustered.small$interactions)



## clustering based on significant bin pairs at FDR 5%

clustered.sig <- diClusters(data, result$table, target=0.05, cluster.args=list(tol=1))
length(clustered.sig$interactions)
##significant at ax 2Mb clustering at 1% FDR: 
#at 5% FDR



clustered.sig$FDR


library(csaw)
tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
head(tabcom)

tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
head(tabbest)

##save the output table - RESULT - DIFF INTERACTIONS CLUSTERS
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
results.d <- data.frame(as.data.frame(clustered.sig$interactions) [,useful.cols], tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
o.d <- order(results.d$PValue)
write.table(results.d[o.d,], file="DIclusteres_40kb_FDR0.05.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## independent clustering of adj bin pairs - to make sparse data more interpretable? - RESULT TABLE 2
y.small <- asDGEList(data)
y.small$offset <- normOffsets(data, type="loess")
y.small <- estimateDisp(y.small, design)
fit.small <- glmQLFit(y.small, design, rebust=TRUE)
result.small <- glmQLFTest(fit.small)

tabcluster <- combineTests(clustered.small$indices[[1]], result.small$table)
head(tabcluster)


sum(tabcluster$FDR <= 0.05)

inter.frame <- as.data.frame(clustered.small$interactions) [,useful.cols]
results.i <- data.frame(inter.frame, tabcluster)
o.i <- order(results.i$PValue)
write.table(results.i[o.i,], file="Independent_40kb_FDR0.05.tsv", sep="\t", quote=FALSE, row.names=FALSE)


useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
results.d <- data.frame(as.data.frame(clustered.sig$interactions) [,useful.cols], tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
o.d <- order(results.d$PValue)


## direct application of BH method - for large bins - easier to interpret each bin pair on its own, without clustering
## to many single bins to interpret?

adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= 0.01)

inter.frame <- as.data.frame(interactions(data)) [,useful.cols]
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
o.r <- order(results.r$PValue)
write.table(results.r[o.r,], file="BinPairs_20kb_FDR0.01.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## identify most significant bin pairs within clustres or other bins - based on p value

#inside <- getBestTest(boxed$indices$smaller, result.small$table)
#best.interactions <- interactions(smaller.data) [inside$best,]

#inter.frame <- as.data.frame(best.interactions) [,useful.cols[-c(1,4)]]
#nested <- data.frame(inter.frame, inside[,c("logFC", "F")])
#head(nested)
#o.r <- order(results.r$PValue)
#write.table(results.r[o.r,], file="BinPairs.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## visualise
chosen <- o.r[1]
chosen.a1 <- anchors(data[chosen], type="first")
chosen.a2 <- anchors(data[chosen], type="second")
expanded1 <- resize(chosen.a1, fix="center", width=bin.size*5)
expanded2 <- resize(chosen.a2, fix="center", width=bin.size*5)
cap.nt <- 30
cap.kd <- cap.nt*data$totals[1]/data$totals[3] ## ko / wt
plotPlaid(input[1], first=expanded1, second=expanded2, max.count=cap.kd, width=5e4, param=hs.param, main="CTCF_KD_merged_1")
rect(start(chosen.a1), start(chosen.a2), end(chosen.a1), end(chosen.a2))

plotPlaid(input[3], first=expanded1, second=expanded2, max.count=cap.nt, width=5e4, param=hs.param, main="CTCF_NT_1")
rect(start(chosen.a1), start(chosen.a2), end(chosen.a1), end(chosen.a2))


## diff plot plaids
chosen <- o.r[5]
chosen.a1 <- anchors(data[chosen], type="first")
chosen.a2 <- anchors(data[chosen], type="second")
expanded1 <- resize(chosen.a1, fix="center", width=5e6)
expanded2 <- resize(chosen.a2, fix="center", width=5e6)
colfun <- plotDI(data, result$table$logFC, expanded1, expanded2, diag=FALSE)


## identify TADs
finder <- domainDirections(input, hs.param, width=4e4, span=10)
prior.c <- 10
directionality <- log2((finder$Up + prior.c)/(finder$Down + prior.c))
onchr1 <- as.logical(seqnames(finder)=="chr1")
chr1.dir <- directionality[onchr1]
hist(chr1.dir, breaks=100)

