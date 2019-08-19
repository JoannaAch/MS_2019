
module load joaach/homer/4.7
module load qiadu/R/3.2.3
module load gi/samtools/1.2
module load gi/seqlogo/2.8.2

makeTagDirectory TAMR_1 TKCC_TAMR_pooled_NcolII_1_R1_genome.bwt2merged.bam TKCC_TAMR_pooled_NcolII_1_R2_genome.bwt2merged.bam -illuminaPE -tb

qsub -cwd -V -pe smp 16 -N homer_makeDir -b y ./makeTagDirectory.sh


runHiCpca.pl pcaOut_TAMR_1 TAMR_1/ -res 50000 -cpu 16 -genome hg38

qsub -cwd -V -pe smp 16 -N homer_2PCA -b y ./pca.sh

needs R in path

submitted - 50KB bins

The genome was then separated into 25K base pair consecutive bins, with each bin marked as type A (active) or type B (inactive). We observed bin switching following AZA treatment and differential loop formations (e.g. promoter-enhancer). Bins that switched from B to A at both time points were found enriched in genes that are related to acute phase response, interferon pathway, and cancer including PI3KCB, DDX58, CD274 and CDKN2A.


#TAMR
new:
makeTagDirectory TAMR TKCC_TAMR_pooled_NcolII_1_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_1_R2_genome.bwt2merged.bam TKCC_TAMR_pooled_NcolII_2_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_2_R2_genome.bwt2merged.bam TKCC_TAMR_pooled_NcolII_3_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE


qsub -cwd -V -pe smp 16 -N homer_tag2 -b y ./makeTagDirectory_2.sh


submitted - 25KB bins

runHiCpca.pl pcaOut_TAMR TAMR/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca2 -b y ./pca.sh

#MCF7
makeTagDirectory MCF7 TKCC_MCF7C_NcolII_1_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_1_R2_genome.bwt2merged.bam TKCC_MCF7C_NcolII_2_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_2_R2_genome.bwt2merged.bam TKCC_MCF7C_NcolII_3_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE

qsub -cwd -V -pe smp 16 -N homer_makeDir3 -b y ./makeTagDirectory_3.sh


submitted - 25KB bins

runHiCpca.pl pcaOut_MCF7 MCF7/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca3 -b y ./pca_mcf7.sh


#FASR

makeTagDirectory FASR TKCC_FASR_pooled_NcolII_1_R1_genome.bwt2merged.bam,TTKCC_FASR_pooled_NcolII_1_R2_genome.bwt2merged.bam TKCC_FASR_pooled_NcolII_2_R1_genome.bwt2merged.bam,TKCC_FASR_pooled_NcolII_2_R2_genome.bwt2merged.bam TKCC_FASR_pooled_NcolII_3_R1_genome.bwt2merged.bam,TKCC_FASR_pooled_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE

qsub -cwd -V -pe smp 16 -N homer_tag2 -b y ./makeTagDirectory_4.sh
makeTagDirectory_4.sh


runHiCpca.pl pcaOut_FASR FASR/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca4 -b y ./pca_fasr.sh

#differences between PCA

getHiCcorrDiff.pl out_MCF7_TAMR MCF7 TAMR -cpu 16 -res 25000

getHiCcorrDiff.pl out_MCF7_FASR MCF7 FASR -cpu 16 -res 25000


qsub -cwd -V -pe smp 16 -N homer_corr_1 -b y ./getHiCcorrDiff_fasr.sh

qsub -cwd -V -pe smp 16 -N homer_corr_2 -b y ./getHiCcorrDiff_tamr.sh



##find compartments A

findHiCCompartments.pl pcaOut_MCF7.PC1.txt > compartmentsA_MCF7.txt
findHiCCompartments.pl pcaOut_TAMR.PC1.txt > compartmentsA_TAMR.txt
findHiCCompartments.pl pcaOut_FASR.PC1.txt > compartmentsA_FASR.txt

find_comp.sh

qsub -cwd -V -pe smp 16 -N homer_find_comp -b y ./find_comp.sh


## find compartments B

findHiCCompartments.pl pcaOut_MCF7.PC1.txt -opp > compartmentsB_MCF7.txt
findHiCCompartments.pl pcaOut_TAMR.PC1.txt -opp > compartmentsB_TAMR.txt
findHiCCompartments.pl pcaOut_FASR.PC1.txt -opp > compartmentsB_FASR.txt

qsub -cwd -V -pe smp 8 -N homer_find_comp -b y ./find_comp_opp.sh


## make bedfiles of compartments

cat compartmentsA_MCF7.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_MCF7.bed
cat compartmentsA_TAMR.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_TAMR.bed
cat compartmentsA_FASR.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_FASR.bed

cat compartmentsB_MCF7.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_MCF7.bed
cat compartmentsB_TAMR.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_TAMR.bed
cat compartmentsB_FASR.txt | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_FASR.bed

# different compartments
Differential Compartments (i.e. Flipping)
Usage: findHiCCompartments.pl <PC1.txt file> [options]

	Options:
		-opp (return inactive, not active regions)
		-thresh <#> (threshold for active regions, default: 0)
		-bg <2nd PC1.txt file> (for differential domains)
			-diff <#> (difference threshold, default: 50)
		-corr <corrDiff.txt file> (for differential domains)
			-corrDiff <#> (correlation threshold, default: 0.4)
		-peaks (output as peaks/regions, not continuous domains)

##MCF7 vs TAMR



qsub -cwd -V -pe smp 8 -N homer_diff1 -b y ./find_comp_opp.sh


findHiCCompartments.pl pcaOut_TAMR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_TAMR.corrDiff.txt  -corrDiff 0.2 > compartmentsA_MCF7-TAMR.bed
findHiCCompartments.pl pcaOut_FASR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_FASR.corrDiff.txt  -corrDiff 0.2 > compartmentsA_MCF7-FASR.bed

cat compartmentsA_MCF7-TAMR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_MCF7-TAMR_bed.bed
cat compartmentsA_MCF7-FASR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_MCF7-FASR_bed.bed



findHiCCompartments.pl pcaOut_TAMR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_TAMR.corrDiff.txt  -corrDiff 0.2 -opp > compartmentsB_MCF7-TAMR.bed
findHiCCompartments.pl pcaOut_FASR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_FASR.corrDiff.txt  -corrDiff 0.2 -opp > compartmentsB_MCF7-FASR.bed

cat compartmentsB_MCF7-TAMR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_MCF7-TAMR_bed.bed
cat compartmentsB_MCF7-FASR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_MCF7-FASR_bed.bed

qsub -cwd -V -pe smp 8 -N homer_diffB -b y ./diff_B.sh


## histogram around features:
## TSS:

annotatePeaks.pl tss hg38 -size 1000000 -hist 1000 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram.txt

qsub -cwd -V -pe smp 8 -N homer_annots1 -b y ./annots_tss.sh

## ER binding sites

annotatePeaks.pl ESR1_nr.bed hg38 -ann ESRbs_custom.txt -size 1000000 -hist 1000 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_er.txt

annotatePeaks.pl JC_MCF7_ER_pooled.bed hg38 -ann ER_MCF7_custom.txt -size 5000 -hist 10 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_er_mcf7.txt

ER_MCF7_custom.txt


ESRbs_custom.txt
ESR1_nr.bed


assignGenomeAnnotation annotations.txt annotations.txt -prioritize annotations.final.txt > stats.txt

assignGenomeAnnotation ESR1_nr.bed ESR1_nr.bed -prioritize ESRbs_custom.txt > stats.txt


## MCF7 ER

JC_MCF7_ER_pooled.bed

assignGenomeAnnotation JC_MCF7_ER_pooled.bed JC_MCF7_ER_pooled.bed -prioritize ER_MCF7_custom.txt > stats.txt



## MCF7 CTCF

Encode_MCF7_CTCF_1.bed

assignGenomeAnnotation Encode_MCF7_CTCF_1.bed Encode_MCF7_CTCF_1.bed -prioritize CTCF_MCF7_custom.txt > stats.txt

##ReMap CTCF

CTCF_nr.bed

assignGenomeAnnotation CTCF_nr.bed CTCF_nr.bed -prioritize CTCF_ReMap_custom.txt > stats2.txt

make_annot_ctcf.sh 

qsub -cwd -V -pe smp 16 -N make_annot2 -b y ./make_annot_ctcf.sh


annotatePeaks.pl CTCF_nr.bed hg38 -ann CTCF_ReMap_custom.txt -size 1000000 -hist 1000 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_ctcf_remap.txt

annotatePeaks.pl Encode_MCF7_CTCF_1.bed hg38 -ann CTCF_MCF7_custom.txt -size 1000000 -hist 1000 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_ctcf_mcf7.txt





annotatePeaks.pl JC_MCF7_ER_pooled.bed hg38 -ann ER_MCF7_custom.txt -size 1000000 -hist 10 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_er_mcf7.txt

annotatePeaks.pl JC_TAMR_ER_pooled.bed hg38 -ann ER_TAMR_custom.txt -size 50000 -hist 10 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_er_tamr.txt


annotatePeaks.pl JC_MCF7_FOXA1.bed hg38 -ann FOXA1_MCF7_custom.txt -size 50000 -hist 10 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_foxa1_mcf7.txt

annotatePeaks.pl JC_TAMR_FOXA1.bed hg38 -ann FOXA1_TAMR_custom.txt -size 50000 -hist 10 -bedGraph pcaOut_MCF7.PC1.bedGraph pcaOut_TAMR.PC1.bedGraph pcaOut_FASR.PC1.bedGraph > output.histogram_foxa1_tamr.txt


FOXA1_MCF7_custom.txt
FOXA1_TAMR_custom.txt
ER_TAMR_custom.txt
ER_MCF7_custom.txt


JC_MCF7_ER_pooled.bed
JC_MCF7_FOXA1.bed
JC_TAMR_ER_pooled.bed
JC_TAMR_FOXA1.bed














##Other features:
Usage: annotatePeaks.pl <peak file | tss> <genome version>  [additional options...]

	Available Genomes (required argument): (name,org,directory,default promoter set)
		hg19	human	/share/ClusterShare/software/contrib/joaach/homer/4.7/.//data/genomes/hg19/	default
		hg38	human	/share/ClusterShare/software/contrib/joaach/homer/4.7/.//data/genomes/hg38/	default
			-- or --
		Custom: provide the path to genome FASTA files (directory or single file)
		If no genome is available, specify 'none'.
		If using FASTA file or none, may want to specify '-organism <...>'

	User defined annotation files (default is UCSC refGene annotation):
		annotatePeaks.pl accepts GTF (gene transfer formatted) files to annotate positions relative
		to custom annotations, such as those from de novo transcript discovery or Gencode.
		-gtf <gtf format file> (Use -gff and -gff3 if appropriate, but GTF is better)
		-gid (by default the GTF file is processed by transcript_id, use this option for gene_id)
		-ann <custom homer annotation file> (created by assignGenomeAnnotation, see website)

	Peak vs. tss/tts/rna mode (works with custom GTF file):
		If the first argument is "tss" (i.e. annotatePeaks.pl tss hg18 ...) then a TSS centric
		analysis will be carried out.  Tag counts and motifs will be found relative to the TSS.
		(no position file needed) ["tts" now works too - e.g. 3' end of gene]
		["rna" specifies gene bodies, will automaticall set "-size given"]
		NOTE: The default TSS peak size is 4000 bp, i.e. +/- 2kb (change with -size option)
		-list <gene id list> (subset of genes to perform analysis [unigene, gene id, accession,
			 probe, etc.], default = all promoters)
		-cTSS <promoter position file i.e. peak file> (should be centered on TSS)

	Primary Annotation Options:
		-mask (Masked repeats, can also add 'r' to end of genome name)
		-m <motif file 1> [motif file 2] ... (list of motifs to find in peaks)
			-mscore (reports the highest log-odds score within the peak)
			-nmotifs (reports the number of motifs per peak)
			-mdist (reports distance to closest motif)
			-mfasta <filename> (reports sites in a fasta file - for building new motifs)
			-fm <motif file 1> [motif file 2] (list of motifs to filter from above)
			-rmrevopp <#> (only count sites found within <#> on both strands once, i.e. palindromic)
			-matrix <prefix> (outputs a motif co-occurrence files:
				prefix.count.matrix.txt - number of peaks with motif co-occurrence
				prefix.ratio.matrix.txt - ratio of observed vs. expected  co-occurrence
				prefix.logPvalue.matrix.txt - co-occurrence enrichment
				prefix.stats.txt - table of pair-wise motif co-occurrence statistics
				additional options:
				-matrixMinDist <#> (minimum distance between motif pairs - to avoid overlap, default: 4)
				-matrixMaxDist <#> (maximum distance between motif pairs)
			-mbed <filename> (Output motif positions to a BED file to load at UCSC (or -mpeak))
			-mlogic <filename> (will output stats on common motif orientations)
		-d <tag directory 1> [tag directory 2] ... (list of experiment directories to show
			tag counts for) NOTE: -dfile <file> where file is a list of directories in first column
		-bedGraph <bedGraph file 1> [bedGraph file 2] ... (read coverage counts from bedGraph files)
		-wig <wiggle file 1> [wiggle file 2] ... (read coverage counts from wiggle files)
		-p <peak file> [peak file 2] ... (to find nearest peaks)
			-pdist to report only distance (-pdist2 gives directional distance)
			-pcount to report number of peaks within region
		-vcf <VCF file> (annotate peaks with genetic variation infomation, one col per individual)
			-editDistance (Computes the # bp changes relative to reference)
			-individuals <name1> [name2] ... (restrict analysis to these individuals)
			-editDistance (Computes the # bp changes relative to reference)
			-individuals <name1> [name2] ... (restrict analysis to these individuals)
		-gene <data file> ... (Adds additional data to result based on the closest gene.
			This is useful for adding gene expression data.  The file must have a header,
			and the first column must be a GeneID, Accession number, etc.  If the peak
			cannot be mapped to data in the file then the entry will be left empty.
		-go <output directory> (perform GO analysis using genes near peaks)
		-genomeOntology <output directory> (perform genomeOntology analysis on peaks)
			-gsize <#> (Genome size for genomeOntology analysis, default: 2e9)

	Annotation vs. Histogram mode:
		-hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)
		The -hist option can be used to generate histograms of position dependent features relative
		to the center of peaks.  This is primarily meant to be used with -d and -m options to map
		distribution of motifs and ChIP-Seq tags.  For ChIP-Seq peaks for a Transcription factor
		you might want to use the -center option (below) to center peaks on the known motif
		** If using "-size given", histogram will be scaled to each region (i.e. 0-100%), with
		the -hist parameter being the number of bins to divide each region into.
			Histogram Mode specific Options:
			-nuc (calculated mononucleotide frequencies at each position,
				Will report by default if extracting sequence for other purposes like motifs)
			-di (calculated dinucleotide frequencies at each position)
			-histNorm <#> (normalize the total tag count for each region to 1, where <#> is the
				minimum tag total per region - use to avoid tag spikes from low coverage
			-ghist (outputs profiles for each gene, for peak shape clustering)
			-rm <#> (remove occurrences of same motif that occur within # bp)

	Peak Centering: (other options are ignored)
		-center <motif file> (This will re-center peaks on the specified motif, or remove peak
			if there is no motif in the peak.  ONLY recentering will be performed, and all other
			options will be ignored.  This will output a new peak file that can then be reanalyzed
			to reveal fine-grain structure in peaks (It is advised to use -size < 200) with this
			to keep peaks from moving too far (-mirror flips the position)
		-multi (returns genomic positions of all sites instead of just the closest to center)

	Genome comparisons (need genome & liftOver)
		-cmpGenome <genome1> [genome2] (Genomes to compare for sequence/motifs)
		-cmpLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)

	Normalization options:
		-fpkm (normalize read counts to million reads or fragments per kilobase mapped)
		-raw (do not adjust the tag counts based on total tags sequenced, -noadj works too)
		-norm <#> (normalize tags to this tag count, default=1e7, 0=average tag count in all directories)
		-normLength <#> (Fragment length to normlize to for experiments with different lens, def: 100)
		-log (output tag counts as log2(x+1+rand) values - for scatter plots)
		-sqrt (output tag counts as sqrt(x+rand) values - for scatter plots)
		-ratio (process tag values as ratios - i.e. chip-seq, or mCpG/CpG)

	Advanced normalization options: (-rlog and -vst require R and DESeq2 to be installed)
		-rlog (quantile/variance normalization on reported genes using DESeq2 rlog funcition, needs R)
		-vst (quantile/variance normalization on reported genes using DESeq2 vst function, needs R)

	Advanced Options:
		-len <#> / -fragLength <#> (Fragment length, default=auto, might want to set to 1 for 5'RNA)
		-size <#> (Peak size[from center of peak], default=inferred from peak file)
			-size #,# (i.e. -size -10,50 count tags from -10 bp to +50 bp from center)
			-size "given" (count tags etc. using the actual regions - for variable length regions)
		-strand <+|-|both> (Count tags on specific strands relative to peak, default: both)
		-pc <#> (maximum number of tags to count per bp, default=0 [no maximum], -tbp <#> works too)
		-CpG (Calculate CpG/GC content)
		-nfr (report nuclesome free region scores instead of tag counts, also -nfrSize <#>)
		-norevopp (do not search for motifs on the opposite strand [works with -center too])
		-gwasCatalog <gwasCatalog file from UCSC> (list overlapping GWAS risk SNPs)
		-pdist (only report distance to nearest peak using -p, not peak name)
		-map <mapping file> (mapping between peak IDs and promoter IDs, overrides closest assignment)
		-noann, -nogene (skip genome annotation step, skip TSS annotation)
		-homer1/-homer2 (by default, the new version of homer [-homer2] is used for finding motifs)
		-cpu <#> (Number of processors to use when possible - only some parts utilize multiple cores)
		-noblanks (remove peaks/rows with missing data)

