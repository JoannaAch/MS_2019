/Users/joanna/Documents/ChromHMM/ChromHMM.jar

Step 1 (12/12/2016)

Binarize bam files

java -mx1600M -jar ChromHMM.jar BinarizeBam chrom_size_file inputdir cellmarkfiletable outputbinarydir

*cellmarkfiletable

tab delimiated file 

A tab delimited file each row contains the cell type or other identifier for a groups of marks, then the associated mark, then the name of a bam file, and optionally a corresponding control bam file:

cell1 mark1 cell1_mark1.bam cell1_control.bam 
cell1 mark2 cell1_mark2.bam cell1_control.bam
cell2 mark1 cell2_mark1.bam cell2_control.bam 
cell2 mark2 cell2_mark2.bam cell2_control.bam 

MCF7 bam files:

/Volumes/GRIW/Cancer-Epigenetics/joaach/ChIP-seq_hg38/MCF7


TAMR bam files:

/Volumes/GRIW/Cancer-Epigenetics/joaach/ChIP-seq_hg38/TAMR

___ binarizeBam
inputdir: /Volumes/GRIW/Cancer-Epigenetics/joaach/ChIP-seq_hg38/chromHMM_hg38
outputdir: /Users/joanna/Desktop/PROJECTS/Endocrine_Resistance_Project/Analysis/hg38/ChromHMM_h38

chromlegnthfile: /Users/joanna/Desktop/PROJECTS/Endocrine_Resistance_Project/Analysis/hg38/ChromHMM_h38/hg38.chrom.sizes

/Users/joanna/Desktop/PROJECTS/Endocrine_Resistance_Project/Analysis/hg38/ChromHMM_h38/cellmarkfiletable


MCF7	H3K27ac	USC20120622_MCF7_H3K27ac.asd.bam	USC20140314_MCF7_input.asd.bam
MCF7	H3K4me1	USC20120622_MCF7_H3K4me1.asd.bam	USC20140314_MCF7_input.asd.bam
MCF7	H3K4me3	AndrewStone_MCF7_H3K4me3_1.asd.bam	USC20140314_MCF7_input.asd.bam
MCF7	H3K27me3	AndrewStone_MCF7_H3K27me3_1.asd.bam	USC20140314_MCF7_input.asd.bam
MCF7	H2AZac	USC20140228_MCF7_H2AZac.asd.bam	USC20140314_MCF7_input.asd.bam
TAMR	H3K27ac	USC20140306_TAMR_H3K27ac.asd.bam	USC20140228_TAMR_input.asd.bam
TAMR	H3K4me1	USC20140306_TAMR_H3K4me1.asd.bam	USC20140228_TAMR_input.asd.bam
TAMR	H3K4me3	AndrewStone_TAMR_H3K4me3_1.asd.bam	USC20140228_TAMR_input.asd.bam
TAMR	H3K27me3	AndrewStone_TAMR_H3K27me3_1.asd.bam	USC20140228_TAMR_input.asd.bam
TAMR	H2AZac	USC20140228_TAMR_H2AZac.asd.bam	USC20140228_TAMR_input.asd.bam

java -mx8000M -jar ChromHMM.jar BinarizeBam hg38.chrom.sizes /Volumes/GRIW/Cancer-Epigenetics/joaach/ChIP-seq_hg38/chromHMM_hg38 cellmarkfiletable 


__ test with just one file copied into the directory:

java -mx8000M -jar ChromHMM.jar BinarizeBam hg38.chrom.sizes ./input/ cellmarkfiletable_test ./


## works with input bam files stored on external drive
output binary files --> output_binary

java -mx8000M -jar ChromHMM.jar BinarizeBam hg38.chrom.sizes /Volumes/EZ_EXT_HD/chromHMM_hg38/ cellmarkfiletable ./

#### run LearnModel

emissions - check 

## reoroder - to rename states according to emissions

model file: model_10.txt

nano labelMappingFile_10
E1	1_Poised_Enahncer
E2	2_Active_Enhancer
E3	3_Strong_Enhancer
E4	4_Strong_Promoter
E5	5_Promoter
E6	6_Promoter
E7	7_Unmarked
E8	8_Unmarked
E9	9_Polycomb_Repressed
E10	10_Unmarked


java -mx8000M -jar ChromHMM.jar Reorder -m labelMappingFile_10 model_10.txt ./output_reorder

### makeSegmentation - using re-ordered files

model inout file:
./output_reorder/model_10.txt

java -mx8000M -jar ChromHMM.jar MakeSegmentation ./output_reorder/model_10.txt ./output_binary ./output_segmentation


## makeBrowsable files - using segmentation

Writing to file ./output_segmentation/MCF7_10_segments.bed
Writing to file ./output_segmentation/TAMR_10_segments.bed

java -mx8000M -jar ChromHMM.jar MakeBrowserFiles -m labelMappingFile_10 -n 10 ./output_segmentation/MCF7_10_segments.bed MCF7_chromHMM MCF7_10states
java -mx8000M -jar ChromHMM.jar MakeBrowserFiles -m labelMappingFile_10 -n 10 ./output_segmentation/TAMR_10_segments.bed TAMR_chromHMM TAMR_10states

_____ MCF7 vs FASR - 8 sates, 4 marks (27me3 did not work in FASR)

copy bam files into external drive
make cellmarkfiletable_f


java -mx8000M -jar ChromHMM.jar BinarizeBam hg38.chrom.sizes /Volumes/EZ_EXT_HD/chromHMM_hg38_MCF7_FASR/  cellmarkfiletable_f ./output_binary_f/

## how many states? Test 6, 8 and 10
java -mx8000M -jar ChromHMM.jar LearnModel output_binary_f/ ./output_learn_f 10 hg38

8 best

nano labelMappingFile_8
E1	Active_Enahncer
E2	Strong_Enhancer
E3	Poised_Enhancer
E4	Unmarked
E5	Unmarked
E6	Promoter
E7	Active_Promoter
E8	Strong_Promoter

(no enter after last record)

java -mx8000M -jar ChromHMM.jar Reorder -m labelMappingFile_8  ./output_reorder_f


java -mx8000M -jar ChromHMM.jar MakeSegmentation ./output_reorder_f/model_8.txt ./output_binary_f ./output_segmentation_f

Writing to file ./output_segmentation_f/FASR_8_segments.bed
Writing to file ./output_segmentation_f/MCF7_8_segments.bed

java -mx8000M -jar ChromHMM.jar MakeBrowserFiles -m labelMappingFile_8 -n 8 ./output_segmentation_f/FASR_8_segments.bed FASR_chromHMM FASR_8states
java -mx8000M -jar ChromHMM.jar MakeBrowserFiles -m labelMappingFile_8 -n 8 ./output_segmentation_f/MCF7_8_segments.bed MCF7_chromHMM MCF7_8states
