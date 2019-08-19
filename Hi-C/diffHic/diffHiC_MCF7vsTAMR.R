library(diffHic)

library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

setwd ("/Users/joanna/Desktop/PROJECTS/Endocrine_Resistance_Project/Analysis/Hi-C_Level_3_hg38/diffHic")
bamFiles <- list.files("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/")
print(bamFiles)

#Date: 12.04.2017

## bam files for diffHiC must be name sorted - sort files using: samtools sort -n - done on cluster using gi/samtools/1.2 module

load("/Volumes/Joanna_HD/diffHic/MCF7vsTAMR_20kb/diffHiC_MCF7vsTAMR_20kb_19042017.RData")
load("/Volumes/Joanna_HD/diffHic/MCF7vsFASR_20kb/diffHic_MCF7vsFASR_20kb_24042017.RData")
hs.farg <- cutGenome(BSgenome.Hsapiens.UCSC.hg38, "CCATGG", overhang=4L) ## NcoI enzyme
hs.farg
hs.param <- pairParam(hs.farg)




### preparePairs - matching mapped reads to restriction fragments - one sample at the time

##MCF71

diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_MCF7C_NcolII_1_sorted.bam.bam", hs.param, file="MCF7_1.h5", dedup=TRUE, minq=10)
diagnostics
#$pairs


diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#

## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("MCF7_1.h5", hs.param, file.out="MCF7_1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)


#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("MCF7_1.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")
# save in /figures

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

#     0      1      2      3 
996360 997830 998049 997219 
llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
col=c("darkgreen", "red", "blue"), lwd=2)

## observ: no self-circulation and no dangling ends?

##2   MCF72
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_MCF7C_NcolII_2_sorted.bam.bam", hs.param, file="MCF7_2.h5", dedup=TRUE, minq=10)
diagnostics
#$pairs
total   marked filtered   mapped 
33944496        0        0 33944496 

$same.id
dangling self.circle 
0          68 

$singles
[1] 0

$chimeras
total  mapped   multi invalid 
0       0       0       0 

diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#[1] NaN
## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("MCF7_2.h5", hs.param, file.out="MCF7_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#   total   length   inward  outward retained 
33944428  4780630   203159   965137 28013558  

#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("MCF7_2.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

#      0       1       2       3 
2520820 2523086 2518982 2528421 

llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
                plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
       col=c("darkgreen", "red", "blue"), lwd=2)


##3   MCF73
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_MCF7C_NcolII_3_sorted.bam.bam", hs.param,
                            file="MCF7_3.h5", dedup=TRUE, minq=10)
diagnostics

#$pairs
total   marked filtered   mapped 
48539432        0        0 48539432 

$same.id
dangling self.circle 
0          39 

$singles
[1] 0

$chimeras
total  mapped   multi invalid 
0       0       0       0 

diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#[1] NaN

## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("MCF7_3.h5", hs.param, file.out="MCF7_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

#   total   length   inward  outward retained 
48539393 10996013   188654  1221083 36160590 

#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("MCF7_3.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

#      0       1       2       3 
4377718 4387855 4385067 4393525 

llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
                plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
       col=c("darkgreen", "red", "blue"), lwd=2)


##TAMR1

diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_TAMR_NcolII_1_sorted.bam.bam", hs.param, file="TAMR_1.h5", dedup=TRUE, minq=10)
diagnostics
#$pairs
total   marked filtered   mapped 
26905931        0        0 26905931 

$same.id
dangling self.circle 
0          49 

$singles
[1] 0

$chimeras
total  mapped   multi invalid 
0       0       0       0 


diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#[1] NaN
## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("TAMR_1.h5", hs.param, file.out="TAMR_1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

#   total   length   inward  outward retained 
26905882  2246080   152824   791712 23728321

#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("TAMR_1.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")
# save in /figures

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

#      0       1       2       3 
1702930 1702983 1701451 1702895 
llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
                plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
       col=c("darkgreen", "red", "blue"), lwd=2)

## observ: no self-circulation and no dangling ends?

##2   TAMR2
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_TAMR_NcolII_2_sorted.bam.bam", hs.param, file="TAMR_2.h5", dedup=TRUE, minq=10)
diagnostics
#$pairs
total    marked  filtered    mapped 
101556344         0         0 101556344 

$same.id
dangling self.circle 
0         205 

$singles
[1] 0

$chimeras
total  mapped   multi invalid 
0       0       0       0 

diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#[1] NaN
## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("TAMR_2.h5", hs.param, file.out="TAMR_2_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)
#    total    length    inward   outward  retained 
101556139  14162511    627846   3146431  83687483 

#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("TAMR_2.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])


0       1       2       3 
7308753 7315071 7318594 7323049 

llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
                plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
       col=c("darkgreen", "red", "blue"), lwd=2)


##3   MCF7_3
diagnostics <- preparePairs("/Volumes/Joanna_HD/bams_diffHiC/hiccup_bams_new/TKCC_TAMR_NcolII_3_sorted.bam.bam", hs.param,
                            file="TAMR_3.h5", dedup=TRUE, minq=10)
diagnostics

#$pairs
total   marked filtered   mapped 
85791440        0        0 85791440 

$same.id
dangling self.circle 
0         227 

$singles
[1] 0

$chimeras
total  mapped   multi invalid 
0       0       0       0 

diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]
#[1] NaN

## filtering artifacts
min.inward <- 1000
min.outward <- 25000

# 25000 outward in user guide
prunePairs("TAMR_3.h5", hs.param, file.out="TAMR_3_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

#   total   length   inward  outward retained 
85791213 12524527   649108  2818215 69862992 

#examine the distribution of inferred frag lentghs ( length is 1000 in User Guide)
diags <- getPairData("TAMR_3.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

#      0       1       2       3 
6140969 6153298 6156634 6161256

llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert) 
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L], plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L],
                plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"),
       col=c("darkgreen", "red", "blue"), lwd=2)





## merging technical replicates


### 20/09/2016 



## counting read pais into interactions - do not use merged biological replicates - merge technical replicates and keep biol separarte
## technical replicates were merged at fastq level

# save datafiles as index files
anchor1.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
anchor2.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
dummy <- data.frame(anchor1.id, anchor2.id, other.data=as.integer(runif(100, 1, 100)))
savePairs(dummy, "example.h5", hs.param)
## for next step use MCF7 1 2 and 3 (3 biological replicates) and TAMR 1 2 and 3 (3 biological replicates)


#### resolution - 1Mb 19/09/2016

#1-3 TAMRs
#4-6 MCF7s

## this step has to be done at 1Mb bin - to much memory to do 100kb with filter=1 nand not enough bin pairs to estamite threshold
input <- c("TKCC_TAMR_NcolII_1_sorted.bam.bam.h5", "TKCC_TAMR_NcolII_2_sorted.bam.bam.h5", "TKCC_TAMR_NcolII_3_sorted.bam.bam.h5", "TKCC_MCF7C_NcolII_1_sorted.bam.bam.h5", "TKCC_MCF7C_NcolII_2_sorted.bam.bam.h5", "TKCC_MCF7C_NcolII_3_sorted.bam.bam.h5")
hs.param
bin.size <- 4e4
## filter - not sure what it it, but in LNCaP/PrEC data Aaron L used 10 --> 10
data <- squareCounts(input, hs.param, width=bin.size, filter=10)
data

#class: InteractionSet 
dim: 1709425 6 
metadata(2): param width
assays(1): counts
rownames: NULL
rowData names(0):
  colnames: NULL
colData names(1): totals
type: ReverseStrictGInteractions
regions: 72903


head(anchors(data, type="first"))
GRanges object with 6 ranges and 1 metadata column:
  seqnames            ranges strand |    nfrags
<Rle>         <IRanges>  <Rle> | <integer>
  [1]     chr1 [705584,  792748]      * |        25
[2]     chr1 [792745,  900367]      * |        33
[3]     chr1 [792745,  900367]      * |        33
[4]     chr1 [792745,  900367]      * |        33
[5]     chr1 [900364, 1002627]      * |        34
[6]     chr1 [900364, 1002627]      * |        34
-------
  seqinfo: 93 sequences from an unspecified genome

head(anchors(data, type="second"))

GRanges object with 6 ranges and 1 metadata column:
  seqnames           ranges strand |    nfrags
<Rle>        <IRanges>  <Rle> | <integer>
  [1]     chr1 [705584, 792748]      * |        25
[2]     chr1 [529439, 600190]      * |        27
[3]     chr1 [705584, 792748]      * |        25
[4]     chr1 [792745, 900367]      * |        33
[5]     chr1 [705584, 792748]      * |        25
[6]     chr1 [792745, 900367]      * |        33
-------
  seqinfo: 93 sequences from an unspecified genome

## choosing a bin width
head(regions(data))
GRanges object with 6 ranges and 1 metadata column:
  seqnames           ranges strand |    nfrags
<Rle>        <IRanges>  <Rle> | <integer>
  [1]     chr1 [     1, 102784]      * |        29
[2]     chr1 [102781, 171571]      * |        13
[3]     chr1 [171568, 320626]      * |        12
[4]     chr1 [320623, 400884]      * |        21
[5]     chr1 [400881, 529442]      * |        29
[6]     chr1 [529439, 600190]      * |        27
-------
  seqinfo: 93 sequences from an unspecified genome

margin.data <- marginCounts(input, hs.param, width=bin.size)
margin.data
class: RangedSummarizedExperiment 
dim: 72903 6 
metadata(1): param
assays(1): counts
rownames: NULL
rowData names(1): nfrags
colnames: NULL
colData names(1): totals
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
#   Mode   FALSE    TRUE    NA's 
logical  284274 1425151       0 


dummy <- data[count.keep,] #interesing interactions
dummy
#class: InteractionSet 
class: InteractionSet 
dim: 1425151 6 
metadata(2): param width
assays(1): counts
rownames: NULL
rowData names(0):
  colnames: NULL
colData names(1): totals
type: ReverseStrictGInteractions
regions: 72903

#remove low abundance interactions
direct <- filterDirect(data)

direct$threshold
#[1] -4.714415

direct.keep <- direct$abundance > log2(2) + direct$threshold ## threshold in user guide 5 
summary (direct.keep)
#    Mode   FALSE    TRUE    NA's 
logical  284274 1425151       0

direct.keep2 <- direct.keep & count.keep # to ensure that retained bin pairs have large absolute counts
summary(direct.keep2)
#   Mode   FALSE    TRUE    NA's 
logical  284274 1425151       0 

## other filer? trended? OMIT - Aaron suggested to use direct filter only



###### peak calling in the interaction space
flank.width <- 5
enrichments <- enrichedPairs(data, flank=flank.width)
summary(enrichments)
#       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-8.6430 -1.0900 -0.3107 -0.5137  0.2242  8.2900 
#min.count 5 in user guide
peak.keep <- filterPeaks(data, enrichments, min.enrich=0.5, min.count=5, min.diag=2L)
sum(peak.keep)
#[1] 39816


## peak-calling diagnostics
neighborhood <- (2*flank.width + 1) * metadata(data)$width
boxed <- boxPairs(data[peak.keep], reference=neighborhood)

out <- tabulate(tabulate(boxed$indices[[1]]))
setNames(round(c(out[1:5], sum(out[5:length(out)]))/sum(out)*100, 1), c(1:5, ">5"))
#   1    2    3    4    5   >5 
35.3 23.0 14.3  9.8  6.3 17.6 

peak.keep2 <-filterPeaks(data, enrichments, min.enrich=0, min.count=5, min.diag=2L) # no enrichment cut off to see if needed
boxed <- boxPairs(data[peak.keep2], reference=neighborhood)
out <- tabulate(tabulate(boxed$indices[[1]]))
setNames(round(c(out[1:5], sum(out[5:length(out)]))/sum(out)*100, 1), c(1:5, ">5"))
# larger proportions suggest that more agressive filtering on the enrichment score is required
#   1    2    3    4    5   >5 
23.3 13.0  8.1  6.8  5.6 48.8 

## creating the "smaller data" object for testing sig interactions - same resolution here
bin.size <- 4e4
smaller.data <- squareCounts(input, hs.param, width=bin.size, filter=10)
direct <- filterDirect(smaller.data, reference=data)
direct$threshold

small.keep <- direct$abundances > direct$threshold + log2(5)
summary(small.keep)
#   Mode   FALSE    TRUE    NA's 
logical 1206209  503216      

# special cases summary
orginal <- data
data <- data[direct.keep2, ]
data
smaller.data <- data[direct.keep2, ]

#class: InteractionSet 
dim: 330472 6 
metadata(2): param width
assays(1): counts
rownames: NULL
rowData names(0):
  colnames: NULL
colData names(1): totals
type: ReverseStrictGInteractions
regions: 29273



## normalization

normfacs <- normOffsets(data)
normfacs
[1] 1.2264851 0.9755719 0.9518191 1.0928354 0.9227664 0.8707179

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
[,1]       [,2]       [,3]       [,4]      [,5]      [,6]
[1,] -0.7672489 -0.2691330 0.02380741 -0.3459529 0.8047603 0.5537670
[2,] -0.6801966 -0.3041508 0.04143132 -0.4044173 0.7986340 0.5486993
[3,] -1.0382868 -0.1749132 0.08231419 -0.1893320 0.7777377 0.5424800
[4,] -0.7450431 -0.2755880 0.03485234 -0.4341037 0.8110478 0.6088347
[5,] -0.6692401 -0.3055382 0.02852524 -0.3971231 0.7983636 0.5450124
[6,] -0.7403362 -0.2770442 0.03637798 -0.4320599 0.8092837 0.6037786

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
#[1] [1] 2.2633348        NA 0.8126069 1.6963809 0.2270919 0.8045230

## "some NA will be present due to removal of low-abundance bins that do not exibit stable behavior during correction"
corrected$max
[1] 195.393994   3.003061   2.066425   1.588233   1.346186   1.218860   1.157997   1.137339   1.122417   1.110288   1.100249   1.091812   1.084632
[14]   1.078453   1.073084   1.068380   1.064227   1.060536   1.057235   1.054267   1.051585   1.049149   1.046929   1.044897   1.043031   1.041311
[27]   1.039721   1.038247   1.036877   1.035601   1.034408   1.033292   1.032246   1.031262   1.030336   1.029463   1.028638   1.027857   1.027118
[40]   1.026416   1.025749   1.025115   1.024511   1.023935   1.023386   1.022861   1.022358   1.021878   1.021417   1.020975
corrected <- correctedContact(orginal, average=FALSE)
anchor1.bias <- corrected$bias[anchors(orginal, type="first", id=TRUE),]
anchor2.bias <- corrected$bias[anchors(orginal, type="second", id=TRUE),]
iter.off <- log(anchor1.bias * anchor2.bias)

## correcting for CNVs - not neccessery here
cnv.offs <- normalizeCNV(data, margin.data)
head(cnv.offs)
#          [,1]      [,2]       [,3]       [,4]      [,5]      [,6]
[1,] -1.482656 0.1627746 0.04271005 -0.4730445 0.9060872 0.8441291
[2,] -1.450210 0.1625738 0.04352591 -0.4942120 0.9047285 0.8335942
[3,] -1.262446 0.1489471 0.06122216 -0.6134148 0.8608634 0.8048277
[4,] -1.617903 0.1479011 0.02064781 -0.4133300 0.9833637 0.8793206
[5,] -1.311013 0.1411140 0.03923753 -0.5670814 0.8784256 0.8193170
[6,] -1.118944 0.1094946 0.04177934 -0.6828798 0.8335770 0.8169730

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
#  TAMR MCF7
1    1    1
2    1    1
3    1    1
4    1    0
5    1    0
6    1    0
attr(,"assign")
[1] 0 1
attr(,"contrasts")
attr(,"contrasts")$`factor(c("TAMR", "TAMR", "TAMR", "MCF7", "MCF7", "MCF7"))`
[1] "contr.treatment"


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
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.3695 52.5200 52.5200 52.2000 52.5200 52.5200 

### testing for significant interactions

## using the quasi-likelihood F-test
result <- glmQLFTest(fit, coef=2)
topTags(result)

#Coefficient:  MCF7 
logFC     logCPM        F       PValue          FDR
72499   7.635960  2.1580594 640.3640 3.008953e-13 9.943747e-08
72500   7.108720  1.3757530 405.6611 7.156884e-12 1.182575e-06
205183  7.844171  0.2351482 381.2912 1.096844e-11 1.208254e-06
129980  5.416137  0.5746072 305.8375 4.974969e-11 3.418292e-06
205179  7.926698 -0.5181231 291.6002 6.886441e-11 3.418292e-06
106307 -7.856399  1.5230655 286.1685 7.827278e-11 3.418292e-06
101306 -6.402354  1.0295264 285.9515 7.867812e-11 3.418292e-06
136683  7.415012 -0.2880344 283.8398 8.274933e-11 3.418292e-06
205188  6.944458 -0.7914673 274.8935 1.028849e-10 3.777841e-06
205178  7.804691  0.1642020 269.3064 1.182832e-10 3.908928e-06
## multiplicity correction and the FDR

# cluster bins that are no more than 1bp apart - not bigger than 100kb for 40kb resolution
clustered.small <- clusterPairs(data, tol=1, upper=1e5)

length(clustered.small$interactions)
[1] 373856
head(clustered.small$interactions)
#ReverseStrictGInteractions object with 6 interactions and 0 metadata columns:
seqnames1            ranges1     seqnames2            ranges2
<Rle>          <IRanges>         <Rle>          <IRanges>
  [1]      chr1 [ 705584, 1700952] ---      chr1 [ 705584, 1700952]
[2]      chr1 [1700949, 2800914] ---      chr1 [ 792745, 1700952]
[3]      chr1 [1700949, 2800914] ---      chr1 [1700949, 2800914]
[4]      chr1 [2800911, 3695417] ---      chr1 [ 792745, 1700952]
[5]      chr1 [2800911, 3801020] ---      chr1 [1700949, 2800914]
[6]      chr1 [2800911, 3801020] ---      chr1 [2800911, 3801020]
-------
  regions: 21094 ranges and 0 metadata columns
seqinfo: 93 sequences from an unspecified genome

## clustering based on significant bin pairs at FDR 5%

clustered.sig <- diClusters(data, result$table, target=0.05, cluster.args=list(tol=1))
length(clustered.sig$interactions)
##significant at ax 2Mb clustering at 1% FDR: 
#at 5% FDR
[1] 678




clustered.sig$FDR
#[1] 0.05014749

library(csaw)
tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
head(tabcom)
#  nWindows logFC.up logFC.down       PValue          FDR
1        1        1          0 2.252309e-06 6.973496e-06
2        1        1          0 1.065524e-05 1.784607e-05
3        4        0          4 2.245813e-05 2.425009e-05
4        9        0          9 8.141687e-06 1.613545e-05
5       12        0         12 7.293840e-06 1.530576e-05
6        1        0          1 1.558589e-05 2.040104e-05
tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
head(tabbest)
#   best      logFC   logCPM        F       PValue          FDR
1  4650  1.5838936 1.451067 57.60483 2.252309e-06 6.973496e-06
2 22325  0.6492815 6.199059 43.73494 1.065524e-05 1.902273e-05
3 25713 -1.2383987 2.446763 47.43325 2.720979e-05 2.790303e-05
4 24597 -0.9025305 5.138642 65.95920 9.157574e-06 1.797144e-05
5 24678 -1.0480495 4.185895 64.00501 1.458988e-05 2.060501e-05
6 24739 -1.3709322 1.165626 40.77890 1.558589e-05 2.108679e-05
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
#  nWindows logFC.up logFC.down    PValue       FDR
1       53        0          1 0.7837882 0.9447686
2       87       27          0 0.4502493 0.7641082
3       58       11          1 0.3171903 0.6662096
4       66       11          5 0.9110106 0.9898437
5       95       12          1 0.7643100 0.9357873
6       55        1          0 0.7071970 0.9120371

sum(tabcluster$FDR <= 0.05)
#[1] 1092

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
#[1] 8058

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

