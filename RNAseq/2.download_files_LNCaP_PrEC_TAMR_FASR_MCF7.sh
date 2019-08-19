#!/bin/bash


##### Login to DNAnexus #####

#source ~/local/lib/dx-toolkit/environment
#dxstart
#dx login

##### download the various files for this RNAseq analysis #####

# FOLDERS_TO_DOWNLOAD=(LNCaP PrEC LNCaP_KD FASR MCF7 TAMR)
FOLDERS_TO_DOWNLOAD=(LNCaP_KD)

ROOT_WORKING_DIR=`pwd`

echo ${ROOT_WORKING_DIR}

for current_folder in ${FOLDERS_TO_DOWNLOAD[*]};
do
  inDir="RNAseq:/"${current_folder}"/"
  outDir="DNAnexus_results/"${current_folder}"/"
  INTERNAL_FOLDERS_LIST=(genome logs rsem)
  for internal_folder in ${INTERNAL_FOLDERS_LIST[*]};
  do
    mkdir -p ${outDir}${internal_folder}
    cd ${outDir}${internal_folder}
    dx download -r ${current_folder}/${internal_folder}/*
    cd ${ROOT_WORKING_DIR}
  done
  # for trimgalore don't download the fastq files
  mkdir -p ${outDir}trimgalore
  cd ${outDir}trimgalore
  dx download -r ${current_folder}/trimgalore/*.txt
  dx download -r ${current_folder}/trimgalore/*.html
  dx download -r ${current_folder}/trimgalore/*.zip
  cd ${ROOT_WORKING_DIR}
done


