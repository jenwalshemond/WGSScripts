##Manhattan Plots
##Combination of bash and R scripts to go from calculating Fst -> Filtering -> Plots
##April 13, 2021

##Start with calculating Fst in sliding windows in VCFTools
#Window size=Size of windows in kb. In this example, Fst will be calculated in 50 kb windows
#Pop files = plain text files with a list of individuals in pop 1 or pop 2. Names as they appear in the VCF
vcftools --vcf FILENAME.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --fst-window-size 50000 --out OUTPUTNAME

#Filter output in R
#Scripts from Sarah Khalil
#Example for how to clean datafiles to make manhattan plot
#15 May 2021

library(scales)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(qqman)
​
options(scipen = 999)
​
#import all files
#all 50kb windows
#In this example, fst_allo = file name throughout
fst_allo<-read.delim("OUTPUTNAME.windowed.weir.fst",header=T)
​
#remove windows were variant<10
fst_allo<- fst_allo %>% filter(N_VARIANTS>9)
​
#change colnames
colnames(fst_allo)[colnames(fst_allo)=="MEAN_FST"] <- "fst_allo"
​
#keep only chrom, bin start, and stat value
fst_allo_sub<-fst_allo[c(1,2,6)]
​
#ended up merging a bunch of files, so if you're only using one file for your stuff, just replace "stat" with "fst_allo"
​
#need this to make sure the BIN_STARTs order correctly
stats <- stats[order(stats$BIN_START),]
​
#remove scaffolds that have less than 3 windows
x.tab <- table(stats$CHROM)
stats_clean<-stats[stats$CHROM %in% names(x.tab[x.tab>=4]),]
​
#first need to add a column to dataframe so manhattan function can work - chrom has to be an integer, not scaffold name
scaffolds<-as.data.frame(unique(fst_allo$CHROM)) #make a dataframe of all scaffolds in order
scaffolds$chromnumber<-1:nrow(scaffolds) #add a column from 1->n
colnames(scaffolds)[colnames(scaffolds)=="unique(fst_allo$CHROM)"] <- "CHROM"
stats_clean_man<-merge(stats_clean,scaffolds,by=c("CHROM"))
​
manhattan(stats_clean_man, chr="chromnumber", bp="BIN_START", p="fst_allo", logp = FALSE, genomewideline = quantile(stats_na_clean_man$fst_allo,.99, na.rm = TRUE), suggestiveline = FALSE, ylab = "Fst", xlab="Scaffold", ylim=c(0,1), col = chr_grey)