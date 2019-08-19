

cat diffHic_MCF7vsFASR_DIclusteres_20kb_FDR0.05.tsv | awk '{if($1==$4){OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}'> cis_FASR_Diclusters_20kb_sig.txt

cat cis_FASR_Diclusters_20kb_sig.txt | awk '{OFS="\t"; print $1,$2,$3}' > cis_FASR_Diclusters_20kb_sig_anchor1.bed
cat cis_FASR_Diclusters_20kb_sig.txt| awk '{OFS="\t"; print $4,$5,$6}' > cis_FASR_Diclusters_20kb_sig_anchor2.bed

cat cis_FASR_Diclusters_20kb_sig.txt | awk '{if($11>"0"){OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'> cis_MCF7-specific_FASR_Diclusters_20kb_sig.txt
cat cis_FASR_Diclusters_20kb_sig.txt | awk '{if($11<"0"){OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'> cis_FASR-specific_FASR_Diclusters_20kb_sig.txt


cat cis_MCF7-specific_FASR_Diclusters_20kb_sig.txt | awk '{OFS="\t"; print $1,$2,$3}' > cis_MCF7-specific_FASR_Diclusters_20kb_cis_sig_anchor1.bed
cat cis_MCF7-specific_FASR_Diclusters_20kb_sig.txt | awk '{OFS="\t"; print $4,$5,$6}' > cis_MCF7-specific_FASR_Diclusters_20kb_cis_sig_anchor2.bed

cat cis_FASR-specific_FASR_Diclusters_20kb_sig.txt | awk '{OFS="\t"; print $1,$2,$3}' > cis_FASR-specific_FASR_Diclusters_20kb_cis_sig_anchor1.bed
cat cis_FASR-specific_FASR_Diclusters_20kb_sig.txt | awk '{OFS="\t"; print $4,$5,$6}' > cis_FASR-specific_FASR_Diclusters_20kb_cis_sig_anchor2.bed
