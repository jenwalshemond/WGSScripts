library(RColorBrewer)
library(SNPRelate)
library(gdsfmt)
library(scales)

setwd("PATHTOFILES")


#reformats the vcf file to gds file for use in futher analysis
snpgdsVCF2GDS(vcf.fn="PATHTOVCF/input.vcf", out.fn="FileName.gds",
              method = c("biallelic.only"),compress.annotation="ZIP.max", 
              snpfirstdim=FALSE, verbose=TRUE)
snpgdsSummary("FileName.gds")

#Import geno file
genofile <- snpgdsOpen("PATHTOFILE/FileName.gds")

#Calculate % Missing Data
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss
write.csv(miss, file = "MissingData.csv")

#Run PCA
pca <- snpgdsPCA(gdsobj = genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#Write PCA output to table
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
tab

#Plot raw PCA output
#For % explained, see PCA output above and fill in values from your data
plot(tab$EV1, tab$EV2, xlab="PC1, % explained", ylab="PC2, % explained",
     main="PCA TITLE")

#Write output
write.csv(tab, file = "MYPCA.csv")

#To color code points
#Order in vcf: open vcf and check order... then create a list of the species in that same order
#So for example, if you have 6 individuals and the first three are pop A and the last three are pop B your new vec would be:
vec <- c("A","A","A","B","B","B")
vec

##Append sample info to PCA output
tab2 <- cbind(tab,vec)
tab2
write.csv(tab2, file = "oriole_PCA_all_SNPs_rm_incl_species.csv")

colorlist<-c("ColorforpopA","ColorforpopB")

plot(tab2$EV1, tab2$EV2, xlab="PC1, X%", ylab="PC2, X%",
     col=alpha(colorlist[as.integer(tab2$vec)],0.6),pch=16,cex=0.75,
     main="PCA TITLE")
legend("bottomleft", legend=levels(tab2$vec), pch=16,col=colorlist, cex=0.75)

snpgdsClose(genofile)