
module load joaach/homer/4.7
module load qiadu/R/3.2.3
module load gi/samtools/1.2
module load gi/seqlogo/2.8.2



#TAMR
##nano makeTagDirectory.sh
makeTagDirectory TAMR TKCC_TAMR_pooled_NcolII_1_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_1_R2_genome.bwt2merged.bam TKCC_TAMR_pooled_NcolII_2_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_2_R2_genome.bwt2merged.bam TKCC_TAMR_pooled_NcolII_3_R1_genome.bwt2merged.bam,TKCC_TAMR_pooled_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE


qsub -cwd -V -pe smp 16 -N homer_tag2 -b y ./makeTagDirectory_2.sh

##runHicpca
runHiCpca.pl pcaOut_TAMR TAMR/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca2 -b y ./pca.sh

#MCF7
##nano makeTagDirectory.sh
makeTagDirectory MCF7 TKCC_MCF7C_NcolII_1_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_1_R2_genome.bwt2merged.bam TKCC_MCF7C_NcolII_2_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_2_R2_genome.bwt2merged.bam TKCC_MCF7C_NcolII_3_R1_genome.bwt2merged.bam,TKCC_MCF7C_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE

qsub -cwd -V -pe smp 16 -N homer_makeDir3 -b y ./makeTagDirectory_3.sh


##runHicpca
runHiCpca.pl pcaOut_MCF7 MCF7/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca3 -b y ./pca_mcf7.sh


#FASR
##nano makeTagDirectory.sh
makeTagDirectory FASR TKCC_FASR_pooled_NcolII_1_R1_genome.bwt2merged.bam,TTKCC_FASR_pooled_NcolII_1_R2_genome.bwt2merged.bam TKCC_FASR_pooled_NcolII_2_R1_genome.bwt2merged.bam,TKCC_FASR_pooled_NcolII_2_R2_genome.bwt2merged.bam TKCC_FASR_pooled_NcolII_3_R1_genome.bwt2merged.bam,TKCC_FASR_pooled_NcolII_3_R2_genome.bwt2merged.bam -illuminaPE

qsub -cwd -V -pe smp 16 -N homer_tag2 -b y ./makeTagDirectory_4.sh


##runHicpca
runHiCpca.pl pcaOut_FASR FASR/ -res 25000 -cpu 16 -genome hg38


qsub -cwd -V -pe smp 16 -N homer_pca4 -b y ./pca_fasr.sh

#differences between PCA

getHiCcorrDiff.pl out_MCF7_TAMR MCF7 TAMR -cpu 16 -res 25000

getHiCcorrDiff.pl out_MCF7_FASR MCF7 FASR -cpu 16 -res 25000


qsub -cwd -V -pe smp 16 -N homer_corr_1 -b y ./getHiCcorrDiff_fasr.sh

qsub -cwd -V -pe smp 16 -N homer_corr_2 -b y ./getHiCcorrDiff_tamr.sh



##find compartments A
##nano find_comp.sh
findHiCCompartments.pl pcaOut_MCF7.PC1.txt > compartmentsA_MCF7.txt
findHiCCompartments.pl pcaOut_TAMR.PC1.txt > compartmentsA_TAMR.txt
findHiCCompartments.pl pcaOut_FASR.PC1.txt > compartmentsA_FASR.txt

qsub -cwd -V -pe smp 16 -N homer_find_comp -b y ./find_comp.sh


## find compartments B
###find_comp_opp.sh
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

##MCF7 vs TAMR

findHiCCompartments.pl pcaOut_TAMR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_TAMR.corrDiff.txt  -corrDiff 0.2 > compartmentsA_MCF7-TAMR.bed
findHiCCompartments.pl pcaOut_FASR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_FASR.corrDiff.txt  -corrDiff 0.2 > compartmentsA_MCF7-FASR.bed

cat compartmentsA_MCF7-TAMR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_MCF7-TAMR_bed.bed
cat compartmentsA_MCF7-FASR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsA_MCF7-FASR_bed.bed



findHiCCompartments.pl pcaOut_TAMR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_TAMR.corrDiff.txt  -corrDiff 0.2 -opp > compartmentsB_MCF7-TAMR.bed
findHiCCompartments.pl pcaOut_FASR.PC1.txt -bg pcaOut_MCF7.PC1.txt -diff 5 -corr out_MCF7_FASR.corrDiff.txt  -corrDiff 0.2 -opp > compartmentsB_MCF7-FASR.bed

cat compartmentsB_MCF7-TAMR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_MCF7-TAMR_bed.bed
cat compartmentsB_MCF7-FASR.bed | awk '{OFS="\t"; print $2,$3,$4}' > compartmentsB_MCF7-FASR_bed.bed

qsub -cwd -V -pe smp 8 -N homer_diffB -b y ./diff_B.sh

qsub -cwd -V -pe smp 8 -N homer_diff1 -b y ./find_comp_opp.sh
##END
