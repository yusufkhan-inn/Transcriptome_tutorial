
#
## Transcriptomics analysis 
#
## PhD course Oct-2021
#
## Author - Yusuf Khan
#
# It is suggested to install rtools, R and Rstudio prior to start this lab-class
#
# #   Transcriptomic analysis / RNA-seq analysis
# 
# Either make project or just got he the directory and start analysing
#
# Data belong to part of study where different tissue were compared 
#
# install("Biocmanager")                    ##run if you don't have this pacakge installed
# BiocManager::install("Rsubread")          ##run if you don't have this pacakge installed
# load rsubread package
# 
# for quality filtering, mapping and quantification
# 

library("Rsubread")

#
#cahnge path to work in folder where data and codes are present
#
#setwd("<..yourpath..>/phd2021")
#
# All the sequencing output file 



fastq.files <- list.files(path = "./data1", pattern = ".fastqsanger$", full.names = TRUE)
fastq.files

## upload data information file

metadata <- read.delim("targets1.txt", header = T)

# build index file to convert genome file into binary for fast mapping

buildindex(basename="chr19_hg",reference="data1/chr19.fa.gz")


args(align)

#align(index="chr19_hg",readfile1=fastq.files,input_format="FASTQ",
#      output_format="BAM",nthreads=2)



align(index="chr19_hg", 
      readfile1 =  "data1/S1Adernal_1.fastqsanger", 
      readfile2 = "data1/S1Adernal_2.fastqsanger",
      input_format= "fastq",
      isGTF = "data1/chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S1Adernal.BAM"
      
        )

align(index="chr19_hg", 
      readfile1 =  "data1/S2Adernal_1.fastqsanger", 
      readfile2 = "data1/S2Adernal_2.fastqsanger",
      input_format= "fastq",
      isGTF = "data1/chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S2Adernal.BAM"
      
)

align(index="chr19_hg", 
      readfile1 =  "data1/S3Adernal_1.fastqsanger", 
      readfile2 = "data1/S3Adernal_2.fastqsanger",
      input_format= "fastq",
      isGTF = "data1/chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S3Adernal.BAM"
      
)


align(index="chr19_hg", 
      readfile1 =  "data1/S3Brain_1.fastqsanger", 
      readfile2 = "data1/S3Brain_2.fastqsanger",
      input_format= "fastq",
      isGTF = "data1/chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S3Brain.BAM"
      
)


align(index="chr19_hg", 
      readfile1 =  "data1/S2Brain_1.fastqsanger", 
      readfile2 = "data1/S2Brain_2.fastqsanger",
      input_format= "fastq",
      isGTF = "chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S2Brain.BAM"
      
)


align(index="chr19_hg", 
      readfile1 =  "data1/S1Brain_1.fastqsanger", 
      readfile2 = "data1/S1Brain_2.fastqsanger",
      input_format= "fastq",
      isGTF = "chr19gft.gtf", phredOffset = 20,
      output_format = "BAM",
      output_file = "S1Brain.BAM"
      
)

bam.files <- list.files( pattern = ".BAM$", full.names = TRUE)
bam.files


props <- propmapped(bam.files,countFragments=TRUE,properlyPaired=FALSE)
props

  qs <- qualityScores(filename="data1/S1Adernal_1.fastqsanger",nreads=100)
boxplot(qs)


fc <- featureCounts(bam.files, annot.ext=file.path("data1/chr19gft.gtf"),
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd=TRUE)
head(fc$counts)

##renaming file names

colnames(fc$counts)<-metadata$File


library(edgeR)          ##load edgeR to find differential genes

y <- DGEList(fc$counts)

y

## making groups 

group <- metadata$Tissue
# Take a look
group

##filter out low count or unexpressed reads
selected <- rowSums(fc$count) >= 1
counts <- fc$counts[selected,]

heatmap(counts, col=topo.colors(50), margin=c(10,6))

# transpose the data to have variables (genes) as columns
data_for_PCA <- t(counts)
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA))  
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 


# create design matrix for limma
design <- model.matrix(~0+group)
# substitute "group" from the design column names
colnames(design)<- gsub("group","",colnames(design))
# check your design matrix
design

# calculate normalization factors between libraries
nf <- calcNormFactors(counts)
#
# normalise the read counts with 'voom' function
y <- voom(counts,design,lib.size=colSums(counts)*nf,plot = TRUE)
y1 <- voom(counts,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- y$E


# save normalised expression data into output dir
#write.table(counts.voom,file="~/counts.voom.txt",row.names=T,quote=F,sep="\t")


# fit linear model for each gene given a series of libraries
fit <- lmFit(y,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(Adernal-Brain,levels=design)
cont.matrix 


# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
options(digits=3)

# check the output fit
dim(fit)


# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=1.5 #choose wisely here

# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
mycoef="Adernal - Brain"
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)

plotMD(fit)

# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

# get significant DE genes only (adjusted p-value < mypval)
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
dim(limma.res.pval)


# get significant DE genes with low adjusted p-value high fold change
limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
dim(limma.res.pval.FC)

ids <- rownames(limma.res.pval.FC)


##
library("gprofiler2")

goresult <- gost(ids, organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC")


p<- gostplot(goresult, capped = TRUE, interactive = FALSE)

p

publish_gosttable(
  goresult,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = NULL,
  ggplot = TRUE
)




