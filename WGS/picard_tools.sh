#MCF7
java -jar picard.jar MarkDuplicates \
      I= /share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7/bwa/LizCaldon_MCF7_WGS_sorted.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7/bwa/LizCaldon_MCF7_WGS_sorted_markedduplicates.asd.bam \
      M=MCF7_marked_dup_metrics.txt

java -jar picard.jar AddOrReplaceReadGroups \
      I=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7C/bwa/LizCaldon_MCF7_WGS_sorted_markedduplicates.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/MCF7C/bwa/LizCaldon_MCF7_WGS_sorted_markedduplicates_addgroup.asd.bam \
      RGID=1 \
      RGLB=1 \
      RGPL=illumina \
      RGPU=None \
      RGSM=MCF7


#TAMR
java -jar picard.jar MarkDuplicates \
      I= /share/ScratchGeneral/joaach/WGS_endocrine_hg38/TAMR/bwa/LizCaldon_TAMR_WGS_sorted.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/TAMR/bwa/LizCaldon_TAMR_WGS_sorted_markedduplicates.asd.bam \
      M=TAMR_marked_dup_metrics.txt

java -jar picard.jar AddOrReplaceReadGroups \
      I=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/TAMR/bwa/LizCaldon_TAMR_WGS_sorted_markedduplicates.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/TAMR/bwa/LizCaldon_TAMR_WGS_sorted_markedduplicates_addgroup.asd.bam \
      RGID=1 \
      RGLB=1 \
      RGPL=illumina \
      RGPU=None \
      RGSM=TAMR


#FASR
java -jar picard.jar MarkDuplicates \
      I= /share/ScratchGeneral/joaach/WGS_endocrine_hg38/FASR/bwa/LizCaldon_FASR_WGS_sorted.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/FASR/bwa/LizCaldon_FASR_WGS_sorted_markedduplicates.asd.bam \
      M=FASR_marked_dup_metrics.txt
java -jar picard.jar AddOrReplaceReadGroups \
      I=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/FASR/bwa/LizCaldon_FASR_WGS_sorted_markedduplicates.asd.bam \
      O=/share/ScratchGeneral/joaach/WGS_endocrine_hg38/FASR/bwa/LizCaldon_FASR_WGS_sorted_markedduplicates_addgroup.asd.bam \
      RGID=1 \
      RGLB=1 \
      RGPL=illumina \
      RGPU=None \
      RGSM=FASR
