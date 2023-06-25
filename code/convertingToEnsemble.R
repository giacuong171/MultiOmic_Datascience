gene.exp <- read.csv("/fast/work/users/phgi10_c/Uni/final_project/final_data/htseq_count/TCGA-BLCA.htseq_counts.tsv.gz",sep = "\t", stringsAsFactors=FALSE, header=TRUE)

genecode.id <- gene.exp[,1]
fixed_names <- sapply(strsplit(genecode.id, ".", fixed=T), function(x) x[1])

gene.exp[,1] <- fixed_names

cnv_final <- read.csv("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/cnv_analysis/final_cnv.csv",stringsAsFactors=FALSE, header=TRUE)


require('biomaRt')

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)
head(annotLookup)

filt.annotLookup <- subset(annotLookup, hgnc_symbol != '')

matched.genes <- filt.annotLookup[filt.annotLookup$hgnc_symbol %in% cnv_final$GeneSymbol, ]

matched.genes1 <- matched.genes[matched.genes$ensembl_gene_id %in% gene.exp[,1],]
