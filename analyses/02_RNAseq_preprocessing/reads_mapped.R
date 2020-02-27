## Mapped Reads plot
## Date: 2019_12_09

## Set - up
customTheme <- list(theme(axis.title.x = element_text(face="bold", size=20), 
                          axis.title.y = element_text(face="bold", size=20), 
                          axis.text.x  = element_text(size=13),
                          axis.text.y  = element_text(size=13),
                          strip.text.x = element_text(size = 20),
                          strip.text.y = element_text(size = 20)) )

lcmPaletteColors <- c( "#b3a2ce", "#4753a4","#bf9e71", "#956025","#b2d9a6", "#0f7c3e")

## Read in
new_map <- read.table("data/Ciera_coveragebed_counts.txt", header = TRUE, row.names = 1)

#change rowname header
new_map <- cbind(ITAG = rownames(new_map), new_map)
rownames(new_map) <- NULL

new_sum <- data.frame(colSums(new_map[2:49]))
colnames(new_sum)[1] <- "reads"
new_sum <- cbind(library = rownames(new_sum), new_sum)
rownames(new_sum) <- NULL

sample <- gsub("[0-9]", "", new_sum$library)
sample <- gsub("\\.", "", sample)
sample


#set genotype
designTable <- as.data.frame(sample)

dim(designTable)
# assign genotype
designTable$genotype <- ifelse(grepl("wt", designTable$sample, ignore.case = T), "wt", ifelse(grepl("tf", designTable$sample, ignore.case = T), "tf2", "unknown"))

#set type
designTable$tissue <- ifelse(grepl("other", designTable$sample, ignore.case = T), "rachis", ifelse(grepl("mbr", designTable$sample, ignore.case = T), "margin", "unknown"))

#Set Region
designTable$region <- ifelse(grepl("a", designTable$sample, ignore.case = T), "tip", ifelse(grepl("c", designTable$sample, ignore.case = T), "base", "middle"))

designTable

for_vis <- cbind(new_sum, designTable) 
for_vis$sample <- paste(for_vis$tissue, for_vis$region)

for_vis$genotype_2 = factor(for_vis$genotype, levels=c("wt", "tf2"))

for_vis$library_2 <- gsub("mbr", "margin", for_vis$library)
for_vis$library_2 <- gsub("other", "rachis", for_vis$library_2)

## Visualize
ggplot(for_vis, aes(library_2, reads, fill = sample)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = lcmPaletteColors) +
  customTheme +
  facet_wrap(~ genotype_2, nrow = 1,scales = "free_x") +
  ylab("Mapped Reads") +
  xlab("Library Sample") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(labels = scales::comma) 
  

