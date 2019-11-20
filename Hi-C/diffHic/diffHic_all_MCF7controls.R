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
hs.param <- pairParam(hs.frag)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("diffHic", version = "3.8")


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
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_MCF7C_NcolII_2_sorted.asd.bam", hs.param, file="MCF7_2.h5", dedup=TRUE, minq=10)
prunePairs("MCF7_2.h5", hs.param, file.out="MCF7_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

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
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_for_diffHiC/bams_for_diffhic/TKCC_TAMR_pooled_NcolII_3_sorted.asd.bam", hs.param, file="TAMR_3.h5", dedup=TRUE, minq=10)
prunePairs("TAMR_3.h5", hs.param, file.out="TAMR_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

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


### 



## counting read pais into interactions - do not use merged biological replicates - merge technical replicates and keep biol separarte
## technical replicates were merged at fastq level

# save datafiles as index files
anchor1.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
anchor2.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
dummy <- data.frame(anchor1.id, anchor2.id, other.data=as.integer(runif(100, 1, 100)))
savePairs(dummy, "example.h5", hs.param)
## for next step use MCF7 1 2 and 3 (3 biological replicates) and TAMR 1 2 and 3 (3 biological replicates)


#### resolution - 100kb


## 
input <- c("MCF7_1_trimmed.h5",
           "MCF7_2_trimmed.h5",
           "MCF7_3_trimmed.h5",
           "TAMR_1_trimmed.h5",
           "TAMR_2_trimmed.h5",
           "TAMR_3_trimmed.h5",
           "FASR_1_trimmed.h5",
           "FASR_2_trimmed.h5",
           "FASR_3_trimmed.h5",
           "MCF7_p16_a2_trimmed.h5",
           "MCF7_p16_a1_trimmed.h5",
           "MCF7_T0_a1_trimmed.h5",
           "MCF7_T0_a2_trimmed.h5",
           "MCF7_p32_a1_trimmed.h5",
           "MCF7_p32_a2_trimmed.h5" )

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


require(edgeR)
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5) #conts per bin = library size

#prior count 5 in user guide, Aaron used 3 in our analysis

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

##Aaron suggested to use direct filter only

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


## do MDS plots for corrected data - after filtering
par(mfrow=c(2,2), mar=c(5,4,2,2))
adjc <- cpm(asDGEList(data), log=TRUE, prior.count=5)
for (top in c(100, 500, 1000, 5000)) {out <- plotMDS(adjc, main=top, top=top, label=1:ncol(data))}

par(mfrow=c(2,2), mar=c(5,4,2,2))
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5)
for (top in c(100, 500, 1000, 5000)) {out <- plotMDS(adjc, main=top, top=top, label=1:ncol(data))}









