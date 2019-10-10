### Cluster Functions 
## Written By: Ciera Martinez
## These are a bunch of functions I created to help me analyze specific 
## clusters when using the Kohonen R package.  These clusters likely will
## only work for this particular project, because of how specific my 
## cluster variables are. Sorry. Maybe someday I will make more generic. 

## Most functions take in the output of SOM or SuperSOM and 1. transform the 
## data and 2. Make a ggplot. 

###########################################
## Prereq annotation files for function
###########################################

clusterVis <- function(clustNum){
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1, 8:13)] # just the sample types
  m.data <- melt(sub_data) 
  p <- ggplot(m.data, aes(x = variable, y = value))
  p + geom_point(alpha = 0.5,position = "jitter", size = 1) + 
    geom_boxplot(alpha = 0.75, outlier.size = 0) + 
    theme_bw() + 
    theme(text = element_text(size = 30),
          axis.text.x = element_text(angle = 90, 
                                     vjust = 1)) +
    xlab("Tissue Region") +
    ylab("Scaled Gene Expression")
}

clusterVis_color <- function(clustNum, palette){
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1, 8:13)] # just the sample types
  m.data <- melt(sub_data) 
  p <- ggplot(m.data, aes(x = variable, y = value, color = variable))
  p + geom_point(alpha = 0.5,position = "jitter", size = 1) + 
    geom_boxplot(alpha = 0.75, outlier.size = 0) + 
    theme_bw(base_size = 28) + 
    scale_color_manual(values = palette) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    scale_y_continuous(breaks = c(0)) +
    xlab("") +
    ylab("")
}

clusterVis_region <- function(clustNum){
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1, 8:13)] # just the sample types
  m.data <- melt(sub_data)
  m.data$region <- ifelse(grepl("wta", m.data$variable, ignore.case = T), "top", 
                          ifelse(grepl("wtb", m.data$variable, ignore.case = T), "mid", "base"))
  
  m.data$tissue <- ifelse(grepl("other", m.data$variable, ignore.case = T), "rachis", 
                          ifelse(grepl("mbr", m.data$variable, ignore.case = T), "margin", "NA"))
  
  p <- ggplot(m.data, aes(y = value, x = tissue, color = tissue))
  p + geom_point(alpha = 0.5,position = "jitter", size = 1) + 
    geom_boxplot(alpha = 0.70, outlier.size = 0) +
    scale_colour_manual(values = c( "#bdbec0","#4da2da")) +
   theme_bw(base_size = 25) + 
   facet_grid(region~.) +
   labs(x = "tissue", y = "scaled value") 
}

## For superSOMs
clusterVis_region_ssom <- function(clustNum){
  sub_cluster <- subset(data.val2, ssom.unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1, 9:14)] # just the sample types
  m.data <- melt(sub_data)
  m.data$region <- ifelse(grepl("A", m.data$variable), "tip", 
                          ifelse(grepl("B", m.data$variable), "middle", "base"))
  
  m.data$tissue <- ifelse(grepl("other", m.data$variable, ignore.case = T), "rachis", 
                          ifelse(grepl("mbr", m.data$variable, ignore.case = T), "margin", "NA"))

  p <- ggplot(m.data, aes(y = value, x = region, color = genotype))
  p + geom_point(alpha = 0.5,
                 position = "jitter", 
                 size = 1) + 
    theme_bw() + 
    scale_colour_manual(values = c( "#d8b365","#5ab4ac")) +
    geom_boxplot(alpha = 0.70, outlier.size = 0) +
    theme(legend.text = element_text(size = 20), 
      text = element_text(size = 20)) + 
    facet_grid(tissue~.)
}


#clusterVis_geno Function
#displays transformed data in a box plot and genotype information 

clusterVis_geno <- function(clustNum){
  
  sub_cluster <- subset(plot.data, som.unit.classif == clustNum)
  
  sub_data <- sub_cluster[,c(1,9:14)] # just the sample types
  names(sub_data)
  m.data <- melt(sub_data) 
  
  m.data$genotype <- as.factor(m.data$genotype)
  
  m.data$region <- ifelse(grepl("A", m.data$variable, ignore.case = T), "A.tip", 
   ifelse(grepl("B", m.data$variable, ignore.case = T), "B.middle", "C.base"))
  #m.data$tissue <- ifelse(grepl("other", m.data$variable, ignore.case = T), "rachis", 
  #ifelse(grepl("mbr", m.data$variable, ignore.case = T), "margin", "base"))

  p <- ggplot(m.data, aes(x = variable, y=value, color = genotype))
  p + geom_point(alpha = 0.5, position = "jitter", size = 1) + 
  geom_boxplot(alpha = 0.75, outlier.size = 0) + 
  theme_bw() + 
  scale_colour_manual(values = c("#ef8a62", "#67a9cf"))  + 
    facet_grid(region ~ .)
}


###clusterGO()

#Prints out how many genes in cluster and performs GO enrichment. 

# must read in genelength file
geneLength <- read.csv("../data/requisite_data/normalized_genes_length.csv")
cate <- read.table("../data/requisite_data/melted.GOTable.txt", header = TRUE)

clusterGO <- function(clustNum){
  ##Sets up plot
  dev.off()
  plot.new()
  
  #sub_cluster
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  
  itag.sc <- as.data.frame(sub_cluster$gene) 
  colnames(itag.sc)[1] <- "itag"
  itag.sc$sc <- 1    
 
  #Since each orthologue between tf2 and wt are represented twice in this set, we have to keep only the unique ITAGs.
  
  itag.sc <- unique(itag.sc) #Check. Should cut the list in half. # dim(itag.sc) before and after
  
  #Merge all by itag
  matrixGO <- merge(itag.sc, geneLength, by = "itag", all = TRUE)
  matrixGO[is.na(matrixGO)] <- 0
  pat <- matrixGO
  
  #Now that we have the data in the right format, we can proceed with GO enrichment.
  
    genes = as.integer(pat[,"sc"])
    names(genes) = pat$itag
    table(genes)
    length(genes)
  
    pwf = nullp(genes,bias.data = pat$length)
  
    GO.wall = goseq(pwf,gene2cat = cate)
    head(GO.wall)
  
  #This is going to correct for multiple testing.  You can specify the p-value cut-off of GO categories you are interested.
  
    enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05]
  
    enriched.GO
  
    my.GO <- as.character(enriched.GO)
    my.GO.table <- Term(my.GO)
    my.GO.table
    t <- as.matrix(my.GO.table)

    print(t) #prints only GO categories that are significant
}


### clusterVis_PCA
## Highlights the cluster in the PCA map 

clusterVis_PCA <- function(clustNum) {
  
  #make dataset for visualization
  data.val3 <- data.val2
  data.val3$cluster[data.val3[,20] == clustNum] <- "subcluster"
  data.val3$cluster[data.val3[,20] != clustNum] <- "other"
  
  #plot

  p <- ggplot(data.val3, aes(PC1, PC2, color = cluster)) 
  p + geom_point(size = I(2), alpha = 0.6) +
    scale_colour_manual(values = c("#cccccc", "#000000")) + 
    theme_bw() + 
    theme(legend.text = element_text(
          size = 30, 
          face = "bold"), 
          text = element_text(size = 30), 
            legend.position = "none")
}

### clusterVis_line

clusterVis_line <- function(clustNum) {
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1,8:13,22)] # just the sample types
  sub_data <- melt(sub_data)
  head(sub_data)
  sub_data <- within(sub_data, lineGroup <- paste(gene,sep = '.'))
  ggplot(sub_data, aes(variable, value, group = lineGroup)) + 
    geom_line(alpha = .1) + 
    geom_point(alpha = .0) +
    theme_bw() 
}

clusterVis_line_ssom <- function(clustNum) {
  sub_cluster <- subset(data.val2, ssom.unit.classif == clustNum)
  sub_data <- sub_cluster[,c(1,2,9:14)] # just the sample types
  sub_data <- melt(sub_data)
  sub_data <- within(sub_data, lineGroup <- paste(gene, genotype, sep = '.'))
  ggplot(sub_data, aes(variable, value, group = lineGroup, color = genotype)) + 
    geom_line(alpha = .1) + 
    geom_point(alpha = .0) +
    theme_bw() +
    scale_color_manual(values = c("#d8b365", "#5ab4ac"
)) +
    theme(text = element_text(size = 20))
}

###genesInClust()
#This looks at how many unique genes are in each cluster. 

genesInClust <- function(clustNum, data.val2, annotation) {
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- as.data.frame(sub_cluster[,1])
  colnames(sub_data) <- "ITAG"
  resultsTable <- merge(sub_data,annotation,by = "ITAG", all.x  = TRUE)
  print(nrow(unique(resultsTable)))
  return(unique(resultsTable))
}

genesInClust_ssom <- function(clustNum) {
  sub_cluster <- subset(data.val2, ssom.unit.classif == clustNum)
  sub_data <- as.data.frame(sub_cluster[,2])
  colnames(sub_data) <- "ITAG"
  resultsTable <- merge(sub_data, annotation, by = "ITAG", all.x = TRUE)
  print(nrow(unique(resultsTable)))
  return(unique(resultsTable))
}


