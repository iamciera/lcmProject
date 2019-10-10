#01_normalization.r
#this file is used for import into further analysis, for the knitted
#report see normalization_report.rmd

## Key to Samples

## genotype: either wildtype of *tf2*

## region: A. tip B. early emmerging leaflet C. base

## type: MBR = Marginal Blastozone Region, 
## other = the rachis or midvein region

### libraries

#for RNAseq
library(edgeR)
library(locfit)
library(statmod)

# for GO enrichment
library(goseq)
library(GO.db)


## RNAseq

### Read in Data
## Read in raw count data per gene.

## This path is made for being sourced from a sub directory. 
## May need to be changed depending on the purpose.

counts <- read.delim("../data/Ciera_coveragebed_counts.txt")

colnames(counts)
head(counts)

## Visualize Reads Mapped
reads.mapped <- colSums(counts[2:49])

reads.mapped <- melt(reads.mapped)
reads.mapped['library'] <- row.names(reads.mapped)


ggplot(reads.mapped, aes(library, value)) +
  scale_y_continuous(labels = scales::comma) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 16))

colnames(counts)

## Get rid of low count libraries "wtbother1.4", "wtbmbr8", "tf2ambr3"
counts <- counts[,-c(37,36,3)]

# Normalization 
# make the groups for design table
colnames(counts)

#Categorize each library into sample type

sample <- gsub("[0-9]", "", colnames(counts))
sample <- gsub("[.]", "", sample)
sample <- sample[-1]

sample

# set genotype
designTable <- as.data.frame(sample)

designTable$genotype <- ifelse(
  grepl("wt", designTable$sample, ignore.case = T), "wt",
    ifelse(
      grepl("tf", designTable$sample, ignore.case = T), "tf2", "unknown")
  )

# set type

designTable$tissue <- ifelse(
  grepl("other", designTable$sample, ignore.case = T), "rachis",
    ifelse(
      grepl("mbr", designTable$sample, ignore.case = T), "mbr",
         "unknown")
  )

# Set Region
designTable$region <- ifelse(
  grepl("a", designTable$sample, ignore.case = T), "A", 
    ifelse(
      grepl("c", designTable$sample, ignore.case = T), "C", "B")
  )


genotype <- designTable$genotype 
sample <- designTable$sample
tissue <- designTable$tissue
region <- designTable$region

# write.csv(designTable, "../data/output/designTable_from_01_normalization.R_30Aug2017.csv")
#put into DGE List
dim(counts)
head(counts)
colnames(counts)
y <- DGEList(counts = counts[,2:46], genes = counts[,1], group = sample)

cpm.y <- cpm(y) #counts per million
y <- y[rowSums(cpm.y > 5) >= 3,] # get rid of genes with low counts

y <- estimateCommonDisp(y) #Estimates common negative binomial dispersion by conditional maximum liklihood
y <- calcNormFactors(y)
y <- estimateCommonDisp(y) #Disp = 0.4606678 

plotMDS(cpm(y, log=TRUE), column = 1)

## export normalized read count
normCounts <- as.data.frame(y$pseudo.counts)

ygenes <- y$genes
normCounts <- cbind(ygenes, normCounts)
head(normCounts)
## write.csv(normCounts, file = "../data/output/normalizedReadCount19July2017.csv")

## RPKM

gene_length <- read.table("~/Desktop/seqLengths_.txt")

