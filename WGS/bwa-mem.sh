#MCF7
bwa mem -t 8 /share/ClusterShare/biodata/contrib/weejar/Homo_sapiens_assembly38_chromosome.fasta LizCaldon_MCF7C_WGS_R1.fastq.gz \
LizCaldon_MCF7C_WGS_R2.fastq.gz -M \ -R '@RG\tID:LizCaldon_MCF7C_WGS\tPL:ILLUMINA\tPU:None\tLB:1\tSM:1\tCN:hcpcg' | samblaster -M --addMateTags | samtools view -h -S -b - \
> LizCaldon_MCF7C_WGS.bam

#TAMR
bwa mem -t 8 /share/ClusterShare/biodata/contrib/weejar/Homo_sapiens_assembly38_chromosome.fasta LizCaldon_TAMR_WGS_R1.fastq.gz \
LizCaldon_TAMR_WGS_R2.fastq.gz -M \ -R '@RG\tID:LizCaldon_TAMR_WGS\tPL:ILLUMINA\tPU:None\tLB:1\tSM:1\tCN:hcpcg' | samblaster -M --addMateTags | samtools view -h -S -b - \
> LizCaldon_TAMR_WGS.bam

#FASR
bwa mem -t 8 /share/ClusterShare/biodata/contrib/weejar/Homo_sapiens_assembly38_chromosome.fasta LizCaldon_FASR_WGS_R1.fastq.gz \
LizCaldon_FASR_WGS_R2.fastq.gz -M \ -R '@RG\tID:LizCaldon_FASR_WGS\tPL:ILLUMINA\tPU:None\tLB:1\tSM:1\tCN:hcpcg' | samblaster -M --addMateTags | samtools view -h -S -b - \
> LizCaldon_FASR_WGS.bam
