# README.md

All changes in differential expression analysis should be done in the `skeletonDE.Rmd` file. Each analysis should be done in a separate directory.  Include regions that are being analyzed in a `de.yml` file corresponding to each analysis and run `skeletonDE.Rmd` in the directory.  

Requisite files:
1. `ITAG2.3_all_Arabidopsis_annotated.tsv` :  
2. `ITAG2.3_all_Arabidopsis_ITAG_annotations` : 
3. `sam2countsResults.tsv` : Results from mapping and counting reads.

This will output:
1. `sample1_sampe2_DE_all.txt` : all genes in analysis
2. `sample1_sample2_DE_sig.txt` : all significantly differentially expressed genes (FDR < .05)
3. `skeletonDE.pdf` : A knitted full report of all the steps of the analysis.

## How to and workflow

1. This is a the skeleton key script (`skeletonDE.Rmd`) for differential expression analysis.  Place this script in the a appropriatly named folder to do the analysis along with  a file YMAL (`.yml`) file called `de.yml`. This YMAL file will specify which sample types you wish to compare. 

2. **After placing `skeletonDE.Rmd` in the folder where you want to perform your analysis, open in Rstudio and set your working directory to that folder**.  
For instance my working directory before starting analysis:

+-- wtcmbr_wtcother (working directory)
|    +-- de.yml
|    +-- skeletonDE.Rmd


**This is the directory after analysis**:

+-- wtcmbr_wtcother (working directory)
|    +-- de.yml
|    +-- skeletonDE.Rmd
|    +-- wtcmbr_wtother_DE_all.txt
|    +-- wtcmbr_wtother_DE_sig.txt
|    +-- wtcmbr_wtcother_DE.pdf




