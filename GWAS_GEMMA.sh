## prep beagle files - impute missing data
java -Xmx96g -jar /programs/beagle4/beagle4.jar gt=Orioles_filtered_final_Hybrids.recode.vcf nthreads=20 out=Orioles_filtered_Hybrids_beagle_output impute=true

## create PLINK files
vcftools --gzvcf Orioles_filtered_Hybrids_beagle_output.vcf.gz --plink --out Orioles_filtered_Hybrids_outputPlinkformat

## make .bed files
/programs/plink-1.9-x86_64-beta5/plink --file Orioles_filtered_Hybrids_outputPlinkformat --make-bed --chr-set 28 --allow-extra-chr 0 --out Orioles_filtered_Hybrids_output_bed

## enter phenotypic information into the .fam file 
## KEEP TRACK OF COLUMN NUMBERS	
# N1 = Supercillium
# N2 = Forehead
# N3 = Neck
# N4 = Ears
# N5 = Throat
# N6 = Greater Coverts
# N7 = Lesser Coverts
# N8 = Tail Base
# N9 = Tail Tip

### RUN GEMMA ###

## GEMMA doesn't seem to like -9 as missing data value. Changing actual missing data to 0. Adding 1 to each plumage score (1-5 instead of 0-4)

## generate relatedness matrix
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o GEMMA_HZ
## number of total SNPs/var = 11651297
## number of analyzed SNPs = 11469909

## run GEMMA: univariate linear models
#supercillium
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 1 -o GWAS_HZ_lmm_supercillium
#Forehead
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 2 -o GWAS_HZ_lmm_forehead
#Neck
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 3 -o GWAS_HZ_lmm_neck
#Ears
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 4 -o GWAS_HZ_lmm_ears
#Throat
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 5 -o GWAS_HZ_lmm_throat
#Greater Coverts
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 6 -o GWAS_HZ_lmm_greatercoverts
#Lesser Coverts
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 7 -o GWAS_HZ_lmm_lessercoverts
#Tail Base
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 8 -o GWAS_HZ_lmm_tailbase
#Tail Tip
gemma -bfile /workdir/jlw395/Orioles_filtered_Hybrids_output_bed -k /workdir/jlw395/output/GEMMA_HZ.cXX.txt -lmm 4 -n 9 -o GWAS_HZ_lmm_tailtip