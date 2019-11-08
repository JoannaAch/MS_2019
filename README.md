# MS_2019

Analysis pipeline for manuscript "Epigenetic reprogramming at estrogen-receptor binding sites alters the 3D chromatin landscape in endocrine resistant breast cancer":
WGS
Hi-C
ChIPseq
RNAseq

Publically used software:
Hi-C
- HiC-Pro v2.9.0
- HiC-Pro/hicpro2juicebox.sh
- diffHiC v.1.15
- Bowtie2 v2.3.2
- “domain-caller” (Bing Ren: http://chromosome.sdsc.edu/mouse/hi-c/download.html)
- MATLAB vR2015b
- Homer v4.8

ChIP-seq
- NGSane v0.5.2.0
- Bowtie v1.1.0
- Masc2 v2.1.0
- diffBind v2.4.8
- chromHMM v1.17 (http://compbio.mit.edu/ChromHMM/)

WGS
- bwa-mem v0.7.9
- GATK v3.5
- QualiMap v2.1.3
-  Mutect2 v3.8-0

WGBS
- Meth10X (https://github.com/luuloi/Meth10X)
- Bpipe v0.9.9.2
- Trim Galore v0.2.8
- Bwa-meth v0.20
- bwa v0.7.13
- Picard v2.3.0
- QualiMap v2.1.3
- MethylDackel (https://github.com/dpryan79/MethylDackel)
- Biscuit (https://github.com/zwdzwd/biscuit)
- Samtools v1.2

RNA-seq MCF7 vs. TAMR and MCF7 vs. FASR
- Trim Galore v0.11.2
- STAR v2.4.0j
- edgeR v3.18.1

Motif analysis
- Homer v4.7 findMotifsGenome.pl
