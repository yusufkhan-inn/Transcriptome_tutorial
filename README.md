# Transcriptome_tutorial

Before starting transcriptome and assembly lab class. Please install following tools and packages.
Install 
Rtools, R and Rstudio
Install packages 
BiocManager, Rsubread , limma, edgeR, tidyverse, ggplot2, gprofiler2
And
Download example data files 

Open CMD and process quality filtering by trimmomatic : follow command below 
java -jar Trimmomatic-0.39\trimmomatic-0.39.jar PE data\S1Adernal_1.fastqsanger data\S1Adernal_2.fastqsanger data1\S1_Adernal_forward_paired.fastq data1\S1_Adernal_reverse_paired.fastq data1\S1_Adernal_forward_unpaired.fastq data1\S1_Adernal_reverse_unpaired.fastq ILLUMINACLIP:Trimmomatic-0.39\adapters\TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

output of trimmomatic will be used to do further analysis 
