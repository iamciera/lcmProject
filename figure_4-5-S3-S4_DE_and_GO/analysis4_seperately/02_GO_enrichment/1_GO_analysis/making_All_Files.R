## Making the all files. 


## Read in File and attach filename in a new column
read_table_filename <- function(filename){
  ret <- read.table(paste0("../../requisiteData/data_06Sept2017/", filename),  header = TRUE, fill = TRUE)
  ret$Source <- filename 
  ret
}

# Get list of file names
import.list <- list.files("../../requisiteData/data_06Sept2017/")

all_files <- ldply(import.list, read_table_filename)

head(all_files)

write.csv(all_files, 
          "../../../06diffGeneExp/analysis4_1Sept2017/analysis/R_afterDEAnalysis/allsig_DE_from_1Sept2017_analysis.csv", 
          row.names = FALSE)


