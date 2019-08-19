#!/bin/bash


##### Login to DNAnexus #####

#module load katgou/dx-toolkit/0.141.0
#dxstart
#dx login

##### Upload the raw files #####

raw_file_dir="fastq/LNCaP"

outDir="RNAseq:/fastq/LNCaP/"

dx mkdir -p $outDir
dx cd $outDir
#make raw file directory on the DNAnexus side
for file in `ls -d $raw_file_dir/TKCC20150326_LNCaP_*`; do
	dx upload $file
done;

raw_file_dir="fastq/PrEC"

outDir="RNAseq:/fastq/PrEC/"

dx mkdir -p $outDir
dx cd $outDir
#make raw file directory on the DNAnexus side
for file in `ls -d $raw_file_dir/*`; do
        dx upload $file
done;

raw_file_dir="fastq/FASR"

outDir="RNAseq:/fastq/FASR/"

dx mkdir -p $outDir
dx cd $outDir
#make raw file directory on the DNAnexus side
for file in `ls -d $raw_file_dir/*`; do
        dx upload $file
done;

raw_file_dir="fastq/MCF7"

outDir="RNAseq:/fastq/MCF7/"

dx mkdir -p $outDir
dx cd $outDir
#make raw file directory on the DNAnexus side
for file in `ls -d $raw_file_dir/*`; do
        dx upload $file
done;

raw_file_dir="fastq/TAMR"

outDir="RNAseq:/fastq/TAMR/"

dx mkdir -p $outDir
dx cd $outDir
#make raw file directory on the DNAnexus side
for file in `ls -d $raw_file_dir/*`; do
        dx upload $file
done;




