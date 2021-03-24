##WGS analyses - GWAS in GEMMA
##3/23/21

## prep beagle files - impute missing data
java -Xmx96g -jar /programs/beagle4/beagle4.jar gt=FILENAME.vcf nthreads=20 out=NAME_output impute=true

## create PLINK files
vcftools --gzvcf NAME_output.vcf.gz --plink --out FILENAME_outputPlinkformat

## make .bed files
/programs/plink-1.9-x86_64-beta5/plink --file FILENAME_outputPlinkformat --make-bed --chr-set Number_SCAFFOLDS --allow-extra-chr 0 --out FILENAME_output_bed

## enter phenotypic information into the .fam file
## note: append these after the 0 0 0
## note: change any missing data -9 to 0
## KEEP TRACK OF COLUMN NUMBERS
# N1 = Trait 1
# N2 = ...

### RUN GEMMA ###
#https://www.xzlab.org/software/GEMMAmanual.pdf

## GEMMA doesn't seem to like -9 as missing data value. Changing actual missing data to 0. Adding 1 to each plumage score (1-5 instead of 0-4)

## generate relatedness matrix
##
##bfile: specify input plink binary file prefix (requires .fam, .bim, and .bed files)
##gk: specify which type of kinship/relatedness matrix to generate (1 is default). 1 = centered matrix and 2 = standardized matrix
##miss: specify missingness threshold (default 0.05)
##maf: specify minor allele frequency threshold
##r2: specify r2 threashold
##hwe: specify HWE test p value threshold (default 0; no test)

gemma -bfile /PATHtoFILENAME_output_bed -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o GEMMA_HZ

## run GEMMA: univariate linear models
##
##lmm: linear mixed model options (aslo can run lm and bslmm - do some research to identify most appropriate options for your data!)
##lmm (cont): 1 Wald test; 2 likelihood ratio test; 3 score test; 4 all (1-3)
##n: specify phenotype column in the phenotype file

gemma -bfile /PATHtoFILENAME_output_bed -k /PATHtoFILENAME_GEMMA_HZ.cXX.txt -lmm 4 -n 1 -o OUTPUTFILE
