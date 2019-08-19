qsub -cwd -N TAMR_bigwigs -pe smp 4 -l h_vmem=12G -b y -j y  ./makewigs_TAMR.sh
qsub -cwd -N MCF7_bigwigs -pe smp 4 -l h_vmem=12G -b y -j y  ./makewigs_MCF7.sh
qsub -cwd -N FASR_bigwigs -pe smp 4 -l h_vmem=12G -b y -j y  ./makewigs_FASR.sh
