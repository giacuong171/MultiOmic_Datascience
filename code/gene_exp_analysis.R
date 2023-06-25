#Preprocess the data
gene.exp <- read.csv("/fast/work/users/phgi10_c/Uni/final_project/final_data/htseq_count/TCGA-BLCA.htseq_counts.tsv.gz",
                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
colnames(gene.exp) <- gsub("\\.","-",colnames(gene.exp))

#ENSG00000000003
gene.exp <- gene.exp[, order(names(gene.exp))]

genecode.id <- gene.exp[,1]
fixed_names <- sapply(strsplit(genecode.id, ".", fixed=T), function(x) x[1])

gene.exp[,1] <- fixed_names

row.names(gene.exp) <- gene.exp[,1]
genes.arr <- row.names(gene.exp)
gene.exp <- gene.exp[-1]
gene.exp <- 2^gene.exp -1
overlapped.samples <- colnames(gene.exp)[-1]

clin.data <- read.csv("/fast/work/users/phgi10_c/Uni/final_project/final_data/TCGA-BLCA.survival.tsv",
                      header = TRUE,sep = "\t",stringsAsFactors = FALSE)

clin.data <- clin.data[order(clin.data$sample), ]

#get only samples in gene.exp
clin.data <- clin.data[clin.data$sample %in% overlapped.samples, ]

gene.exp <- gene.exp[,clin.data$sample]
gene.exp <- as.data.frame(lapply(gene.exp, as.integer))
#gene expresison analysis
require(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=gene.exp,
                           colData=DataFrame(OS=as.factor(clin.data$OS)),
                           design=~OS
                           )

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
normalizedCounts <- as.data.frame(counts(dds, normalized = TRUE))
normalizedCounts$genes <- genes.arr
colnames(normalizedCounts) <- gsub("\\.","-",colnames(normalizedCounts))

sig_indices <-which(res$padj < 0.05 & (res$log2FoldChange < -1 | res$log2FoldChange > 1))

significantGenes <- normalizedCounts[sig_indices,]


#### Intersect genes with result from CNV and from gene expression analysis
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

cnv_genes_exp <- which(significantGenes$genes %in% matched.genes$ensembl_gene_id)

significantGenes_cnv <- significantGenes[cnv_genes_exp,]

write.table(significantGenes_cnv, file = "/fast/work/users/phgi10_c/Uni/final_project/gene_exp_analysis/cnv_genes_exp.combine.csv", sep = "\t", quote = FALSE, row.names = FALSE)
