##GWAS in GEMMA
##March 23, 2021

## prep beagle files - impute missing data
java -Xmx96g -jar /programs/beagle4/beagle4.jar gt=FILENAME.vcf nthreads=20 out=FILENAME_beagle_output impute=true

## create PLINK files
vcftools --gzvcf FILENAME.vcf.gz --plink --out FILENAME_outputPlinkformat

## make .bed files
/programs/plink-1.9-x86_64-beta5/plink --file FILENAME_outputPlinkformat --make-bed --chr-set [number of scaffolds] --allow-extra-chr 0 --out FILENAME_output_bed

## enter phenotypic information into the .fam file
## KEEP TRACK OF COLUMN NUMBERS
# N1 = trait one
# N2 = ...

### RUN GEMMA ###
#https://www.xzlab.org/software/GEMMAmanual.pdf
## GEMMA doesn't seem to like -9 as missing data value. Changing actual missing data to 0. Adding 1 to each plumage score (1-5 instead of 0-4)

## generate relatedness matrix
##
##gk: specify which type of kinship/relatedness matrix to generate (default 1; valid value 1-2; 1: centered matrix; 2: standardized matrix.)
##miss: specify missingness threshold (default 0.05)
##maf: specify minor allele frequency threshold (default 0.01)
##r2: specify r-squared threshold (default 0.9999)
##hwe: specify HWE test p value threshold (default 0; no test)

gemma -bfile PATHTO/output_bed -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o GEMMA_HZ

## run GEMMA: univariate linear models
##
##lmm: linear mixed model (also can run lm or bslmm. read about models and determine the best option for your data)
##lmm (cont): specify frequentist analysis choice (default 1; valid value 1-4; 1: Wald test; 2: likelihood ratio test; 3: score test; 4: all 1-3.)
##n: specify phenotype column in the phenotype file (default 1); or to specify which phenotypes are used in the mvLMM analysis

#Trait one

gemma -bfile PATHTO/output_bed -k PAHTTO/GEMMA_HZ.cXX.txt -lmm 4 -n 1 -o GWAS_HZ_lmm_trait1
