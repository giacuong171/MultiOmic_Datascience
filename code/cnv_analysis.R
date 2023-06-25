require(data.table)
CN <- read.table("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/final_data/no_FFPE_CNV/cnv.segment.seg.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
clinical <- read.table("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/final_data/survival.txt",sep = "\t", stringsAsFactors=FALSE, header=TRUE)
#gene_exp <- read.table("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/final_data/gene_expression/gene_exp/BLCA.mRNAseq_raw_counts.txt",sep = "\t", stringsAsFactors=FALSE, header=TRUE)
#gene_exp_pp <- read.table("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/data/TCGA-BLCA.htseq_fpkm.txt",sep = "\t", stringsAsFactors=FALSE, header=TRUE)
#   Tumor types range from 01 - 09
#   Normal types from 10 - 19
#   Control samples from 20 - 29
types <- gsub("TCGA-[A-Z0-9]*-[A-Z0-9]*-", "", (gsub("-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN[,1])))
clin_types <- gsub("TCGA-[A-Z0-9]*-[A-Z0-9]*-", "", clinical[,1])
unique(clin_types)
#"01" "11" "06"
unique(types)
#"11A" "01A" "10A" "10B" "06A" "10C" "01B" "10D" "10F" "10G"
# NB - in the next 2 lines, it is important to account for all types that are listed from unique(types) command
matched.normal <- c("10A", "10B","10C","10D","10F","10G", "11A", "11B")
tumor <- c("01A", "01B", "06A")

#Divide up the data into tumor and normal
tumor.indices <- which(types %in% tumor)
matched.normal.indices <- which(types %in% matched.normal)
CN.tumor <- CN[tumor.indices,]
CN.matched.normal <- CN[matched.normal.indices,]

#Change IDs to allow for matching to metadata
CN.tumor[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN.tumor[,1])
CN.matched.normal[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN.matched.normal[,1])
################################################################

require(gaia)

# Retrieve probes meta file from Broad Institute
gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
markersMatrix <- read.delim(paste0(gdac.root,"genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt"), as.is=TRUE, header=FALSE)
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")

# Change sex chr names
unique(markersMatrix$Chromosome)
markersMatrix[which(markersMatrix$Chromosome=="X"),"Chromosome"] <- 23
markersMatrix[which(markersMatrix$Chromosome=="Y"),"Chromosome"] <- 24
markersMatrix$Chromosome <- sapply(markersMatrix$Chromosome, as.integer)

# Create a marker ID
markerID <- apply(markersMatrix, 1, function(x) paste0(x[2], ":", x[3]))

#Remove duplicates
markersMatrix <- markersMatrix[-which(duplicated(markerID)),]

# Filter markersMatrix for common CNV
# The 'common CNV' is derived from a Broad Institute panel of normals
markerID <- apply(markersMatrix, 1, function(x) paste0(x[2], ":", x[3]))
commonCNV <- read.delim(paste0(gdac.root,"CNV.hg19.bypos.111213.txt"), as.is=TRUE)
commonCNV[,2] <- sapply(commonCNV[,2], as.integer)
commonCNV[,3] <- sapply(commonCNV[,3], as.integer)
commonID <- apply(commonCNV, 1, function(x) paste0(x[2], ":", x[3]))
table(commonID %in% markerID)
table(markerID %in% commonID)
markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]

# Create the markers object
markers_obj <- load_markers(markersMatrix_fil)
#############################################################
#Prepare CNV matrix
cnvMatrix <- CN.tumor

#Add label (0 for loss, 1 for gain)
#A segment mean of 0.3 is defined as the cut-off
cnvMatrix <- cbind(cnvMatrix, Label=NA)
cnvMatrix[cnvMatrix$Segment_Mean < -0.3,"Label"] <- 0
cnvMatrix[cnvMatrix$Segment_Mean > 0.3,"Label"] <- 1
cnvMatrix <- cnvMatrix[!is.na(cnvMatrix$Label),]

#Remove segment mean as we now go by the binary classification of gain or loss
cnvMatrix <- cnvMatrix[,-6]
colnames(cnvMatrix) <- c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")

#Substitute Chromosomes "X" and "Y" with "23" and "24"
xidx <- which(cnvMatrix$Chromosome=="X")
yidx <- which(cnvMatrix$Chromosome=="Y")
cnvMatrix[xidx,"Chromosome"] <- 23
cnvMatrix[yidx,"Chromosome"] <- 24
cnvMatrix$Chromosome <- sapply(cnvMatrix$Chromosome, as.integer)

#Run GAIA, which looks for recurrent aberrations in your input file
n <- length(unique(cnvMatrix[,1]))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, n)
results.all <- runGAIA(cnv_obj, markers_obj, output_file_name="/data/gpfs-1/users/phgi10_c/work/Uni/final_project/cnv_analysis/Tumor.All.txt", aberrations=-1, chromosomes=-1, num_iterations=10, threshold=0.05)


######################

results.all <- read.csv("/data/gpfs-1/users/phgi10_c/work/Uni/final_project/cnv_analysis/Tumor.All.txt",sep = "\t", stringsAsFactors=FALSE, header=TRUE)

#Convert the recurrent aberrations to numeric (GAIA saves it as text)
RecCNV <- t(apply(results.all, 1, as.numeric))
colnames(RecCNV) <- colnames(results.all)

#Add a new column for 'score'
RecCNV <- cbind(RecCNV, score=0)

#Determine the minimum Q value that's not equal to 0
minval <- format(min(RecCNV[RecCNV[,"q.value"]!=0, "q.value"]), scientific=FALSE)
minval <- substring(minval,1, nchar(minval)-1)

#Replace Q values of 0 with the minimum, non-zero value
RecCNV[RecCNV[,"q.value"]==0, "q.value"] <- as.numeric(minval)

#Set the score to equal -log base 10 of the Q value
RecCNV[,"score"] <- sapply(RecCNV[,"q.value"], function(x) -log10(as.numeric(x)))

#Create a function for plotting the recurrent copy number variants
gaiaCNVplot <- function (calls, threshold=0.05, main="main") {
  Calls <- calls[order(calls[,"Region.Start..bp."]),]
  Calls <- Calls[order(Calls[,"Chromosome"]),]
  rownames(Calls) <- NULL
  Chromo <- Calls[,"Chromosome"]
  Gains <- apply(Calls,1,function(x) ifelse(x["Aberration.Kind"]==1, x["score"], 0))
  Losses <- apply(Calls, 1,function(x) ifelse(x["Aberration.Kind"]==0, x["score"], 0))
  plot(Gains, ylim=c(-max(Calls [,"score"]+2), max(Calls[,"score"]+2)), type="h", col="red2", xlab="Chromosome", ylab=expression("-log"[10]~italic(Q)~"value"), main=main, cex.main=4, xaxt="n", font=2, font.axis=2, font.lab=2, font.axis=2)
  points(-(Losses), type="h", col="forestgreen")
  abline(h= 0, cex=4)
  abline(h=-log10(threshold), col="black", cex=4, main="test", lty=6, lwd=2)
  abline(h=log10(threshold), col="black", cex=4, main="test", lty=6, lwd=2)
  uni.chr <- unique(Chromo)
  temp <- rep(0, length(uni.chr))
  
  for (i in 1:length(uni.chr)) {
    temp[i] <- max(which(uni.chr[i] == Chromo))
  }
  
  for (i in 1:length(temp)) {
    abline(v = temp[i], col = "black", lty = "dashed", )
  }
  
  nChroms <- length(uni.chr)
  
  begin <- c()
  
  for (d in 1:nChroms) {
    chrom <- sum(Chromo == uni.chr[d])
    begin <- append(begin, chrom)
  }
  
  temp2 <- rep(0, nChroms)
  
  for (i in 1:nChroms) {
    if (i == 1) {
      temp2[1] <- (begin[1] * 0.5)
    }
    else if (i > 1) {
      temp2[i] <- temp[i - 1] + (begin[i] * 0.5)
    }
  }
  
  uni.chr[uni.chr==23] <- "X"
  uni.chr[uni.chr==24] <- "Y"
  
  for (i in 1:length(temp)) {
    axis(1, at = temp2[i], labels = uni.chr[i], cex.axis = 1)
  }
  
  #legend("topright", y.intersp=0.8, c("Amplification"), pch=15, col=c("red2"), text.font=2)
  #legend("bottomright", y.intersp=0.8, c("Deletion"), pch=15, col=c("forestgreen"), text.font=2)
}

gaiaCNVplot(RecCNV, threshold=0.05, "A")
##########################




