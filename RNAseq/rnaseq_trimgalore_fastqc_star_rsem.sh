#!/bin/bash

#samples=( "LNCaP_KD" )
samples=( "LNCaP" "PrEC" "FASR" "MCF7" "TAMR" )

for sample in ${samples[@]};do
	echo $sample
	for file in `dx ls /fastq/$sample/*_R1.fastq.gz`;do
		echo $file 
		dx run rnaseq_trimgalore_fastqc_star_rsem --dest /$sample -istage-Bb1QJQj0g8Kxv0381844vVf1.reads=/fastq/$sample/$file -istage-Bb1QJQj0g8Kxv0381844vVf1.reads2=/fastq/$sample/${file%_R1.fastq.gz}_R2.fastq.gz --yes;
	done;
done;
