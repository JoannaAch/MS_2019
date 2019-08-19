#!/bin/bash -e

source /etc/profile.d/modules.sh
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH

module load gi/bedtools/2.22.0 gi/ucsc_utils/283

INPUT_FOLDER="FASR/"
BAM_FILES=$(find $INPUT_FOLDER -name \*.bam)

for bam_file in ${BAM_FILES[*]};
do

  # parse current file location into array
  IFS='/' read -a array <<< "$bam_file"
  # sample name will be the last index of array length
  BAM_FILE_NAME=${array[${#array[@]}-1]}
  
  # remove "stats" to get bam file name
  SAMPLE_NAME="$(sed 's/.genome.sorted.bam//g' <<< "$BAM_FILE_NAME")"
  
  CMD="genomeCoverageBed -bg -split -ibam ${bam_file} -g hg38.chrom.sizes > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph"
  echo $CMD & eval $CMD
  
  CMD="sed -n '/chr/p' ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt"
  echo $CMD & eval $CMD
  
  CMD="sort -k1,1 -k2,2n ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted"
  echo $CMD & eval $CMD

  CMD="bedGraphToBigWig ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted hg38.chrom.sizes ${INPUT_FOLDER}${SAMPLE_NAME}.bw"
  echo $CMD & eval $CMD
  
  CMD="rm ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph*"
  echo $CMD & eval $CMD
  
done


#!/bin/bash -e

source /etc/profile.d/modules.sh
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH

module load gi/bedtools/2.22.0 gi/ucsc_utils/283

INPUT_FOLDER="MCF7/"
BAM_FILES=$(find $INPUT_FOLDER -name \*.bam)

for bam_file in ${BAM_FILES[*]};
do

  # parse current file location into array
  IFS='/' read -a array <<< "$bam_file"
  # sample name will be the last index of array length
  BAM_FILE_NAME=${array[${#array[@]}-1]}
  
  # remove "stats" to get bam file name
  SAMPLE_NAME="$(sed 's/.genome.sorted.bam//g' <<< "$BAM_FILE_NAME")"
  
  CMD="genomeCoverageBed -bg -split -ibam ${bam_file} -g hg38.chrom.sizes > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph"
  echo $CMD & eval $CMD
  
  CMD="sed -n '/chr/p' ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt"
  echo $CMD & eval $CMD
  
  CMD="sort -k1,1 -k2,2n ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted"
  echo $CMD & eval $CMD

  CMD="bedGraphToBigWig ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted hg38.chrom.sizes ${INPUT_FOLDER}${SAMPLE_NAME}.bw"
  echo $CMD & eval $CMD
  
  CMD="rm ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph*"
  echo $CMD & eval $CMD
  
done




#!/bin/bash -e

source /etc/profile.d/modules.sh
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH

module load gi/bedtools/2.22.0 gi/ucsc_utils/283

INPUT_FOLDER="TAMR/"
BAM_FILES=$(find $INPUT_FOLDER -name \*.bam)

for bam_file in ${BAM_FILES[*]};
do

  # parse current file location into array
  IFS='/' read -a array <<< "$bam_file"
  # sample name will be the last index of array length
  BAM_FILE_NAME=${array[${#array[@]}-1]}
  
  # remove "stats" to get bam file name
  SAMPLE_NAME="$(sed 's/.genome.sorted.bam//g' <<< "$BAM_FILE_NAME")"
  
  CMD="genomeCoverageBed -bg -split -ibam ${bam_file} -g hg38.chrom.sizes > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph"
  echo $CMD & eval $CMD
  
  CMD="sed -n '/chr/p' ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt"
  echo $CMD & eval $CMD
  
  CMD="sort -k1,1 -k2,2n ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt > ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted"
  echo $CMD & eval $CMD

  CMD="bedGraphToBigWig ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph.filt.sorted hg38.chrom.sizes ${INPUT_FOLDER}${SAMPLE_NAME}.bw"
  echo $CMD & eval $CMD
  
  CMD="rm ${INPUT_FOLDER}${SAMPLE_NAME}.bedGraph*"
  echo $CMD & eval $CMD
  
done
