java -jar GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R /share/ClusterShare/biodata/contrib/genomeIndices_garvan/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     -I:tumor /share/ScratchGeneral/joaach/WGS_endocrine_hg38/FASR/bwa/LizCaldon_FASR_WGS_sorted_markedduplicates_addgroup.asd.bam \
     -I:normal /share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7C/bwa/LizCaldon_MCF7_WGS_sorted_markedduplicates_addgroup.asd.bam \
     -o /share/ScratchGeneral/joaach/WGS_endocrine_hg38/output_mutect2_fasr-mcf7_chr13-24.vcf \
     --dbsnp /share/ScratchGeneral/joaach/temp/dbsnp_146.hg38.vcf \
     --cosmic /share/ScratchGeneral/joaach/temp/CosmicCodingMuts.hg38.final2.vcf \
     --cosmic /share/ScratchGeneral/joaach/temp/CosmicNonCodingVariants.hg38.final2.vcf \
     -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY


java -jar GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R /share/ClusterShare/biodata/contrib/genomeIndices_garvan/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     -I:tumor /share/ScratchGeneral/joaach/WGS_endocrine_hg38/TAMR/bwa/LizCaldon_TAMR_WGS_sorted_markedduplicates_addgroup.asd.bam \
     -I:normal /share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7C/bwa/LizCaldon_MCF7_WGS_sorted_markedduplicates_addgroup.asd.bam \
     -o /share/ScratchGeneral/joaach/WGS_endocrine_hg38/output_mutect2_tamr-mcf7_chr.vcf \
     --dbsnp /share/ScratchGeneral/joaach/temp/dbsnp_146.hg38.vcf \
     --cosmic /share/ScratchGeneral/joaach/temp/CosmicCodingMuts.hg38.final2.vcf \
     --cosmic /share/ScratchGeneral/joaach/temp/CosmicNonCodingVariants.hg38.final2.vcf \
     -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -nct 14 \
