
# expression per cell line
# author: Kate Gould
# date: 16/10/2015

library(data.table)

TABLES.FOLDER <- "tables/"

OUTPUT.FOLDER <- "expr_per_cell_line_output/"
if (!(file.exists(OUTPUT.FOLDER))){
  dir.create(OUTPUT.FOLDER)
}

tpm.data <- read.csv(paste0(TABLES.FOLDER, "GENE_TAMR_FASR_MCF7_MCF7x_TPM_table.csv"), row.names=1)

colnames(tpm.data) <- gsub("_", ".", paste0 (colnames(tpm.data), ".TPM"))

# only keep gene names beginning with ENSG ie remove ERCC spike-ins
tpm.data <- tpm.data[grep("ENSG*", rownames(tpm.data)),]

tpm.data$TAMR.mean.TPM <- round(rowMeans(tpm.data[,grep("TAMR", colnames(tpm.data))]), 2)
tpm.data$FASR.mean.TPM <- round(rowMeans(tpm.data[,grep("FASR", colnames(tpm.data))]), 2)
tpm.data$MCF7.mean.TPM <- round(rowMeans(tpm.data[,grep("MCF7\\.", colnames(tpm.data))]), 2)
tpm.data$MCF7x.mean.TPM <- round(rowMeans(tpm.data[,grep("MCF7x", colnames(tpm.data))]), 2)

column.order <- c(grep("TAMR", colnames(tpm.data)), grep("FASR", colnames(tpm.data)),  
                  grep("MCF7\\.", colnames(tpm.data)),  grep("MCF7x", colnames(tpm.data)))
tpm.data <- tpm.data[column.order]

# annotation

# annotation
# FROM GTF FILE

gtf.annotation.file <- "../data/hg38/gtf/gene_annotation.tsv"
gtf.anno <- read.table(gtf.annotation.file, header=TRUE, stringsAsFactors=F)
m <- match(rownames(tpm.data), gtf.anno$gene.id)

# assign the rownames(tpm.data) as the gene_id as more specific with version number at end
# as originates from original gtf
tpm.data$gene.id <- rownames(tpm.data)
tpm.data$gene.symbol <- gtf.anno$gene.name[m]
tpm.data$chr <- gtf.anno$chr[m]
tpm.data$start <- gtf.anno$start[m]
tpm.data$end <- gtf.anno$end[m]
tpm.data$strand <- gtf.anno$strand[m]
tpm.data$gene.type <- gtf.anno$gene.type[m]

# annotation file downloaded from http://dec2013.archive.ensembl.org/biomart/martview/
# has following columns: ensembl.gene.ID, chr, start, end, strand, description, HGNC.symbol, entrez.gene.ID
annotation.file <- "Ensembl_87_GRCh38.p7_biomart_export.tsv"
if (!(file.exists(annotation.file))){
  untar(paste0("input/", annotation.file, ".tgz"))
}

DT <- fread(annotation.file)
setnames(DT, gsub(" ", ".", colnames(DT)))
setkey(DT, ensembl.gene.ID)
m <- match(gsub("\\.[0-9]*", "", rownames(tpm.data)), DT$ensembl.gene.ID)

tpm.data$description <- DT$description[m]
tpm.data$entrez.gene.id <- DT$entrez.gene.ID[m]

col.order <- c(17:25, 4, 1:3, 8, 5:7, 12, 9:11, 16, 13:15)
tpm.data <- tpm.data[col.order]
tpm.data[is.na(tpm.data)] <- ""

write.table(tpm.data, file=paste0(OUTPUT.FOLDER, "GENE_TAMR_FASR_MCF7_MCF7x_expression_data.tsv"), 
            sep="\t", quote=FALSE, row.names=F)
