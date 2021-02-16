#### AMPSEQ ANALYSIS pipeline - Run on Oriole Data Set
#### FEB 2021
#### Jen Walsh, STEPFANIE M. AGUILLON, Bronwyn Butcher


# renamed all raw seq data files with the follow command
# sequencing ends up being named SAMPLEID.fastq.gz
mv 11818_3270_121878_J8KBL_PlateB_E04_488431_CGGAGCCT_AAGGAGTA_R1.fastq.gz 488431.fastq.gz


##### Sample QC #####
# run FASTQC (version 0.11.8)
# sample qc was conducted using the following command (for each sample)
fastqc SAMPLEID.fastq.gz

# run MULTIQC (version 1.7)

#activate conda and run multiqc
export LC_ALL=en_US.UTF-8
export PATH=/programs/miniconda3/bin:$PATH
source activate multiqc

multiqc .
# multiqc searches the current directory for results from fastqc analysis and compiles a multiqc report.
# note there were some errors when running this version of multiqc... however it contiued and appears to have created the report.

conda deactivate


##### Remove adaptor contamination #####

# unzip reads using gunzip
# takes only one minute
gunzip SAMPLEID.fastq.gz

# run cutadapt
# used the following command (for each sample) and the perl driver file
cutadapt -a CTGTCTCTTATACACATCT -o SAMPLEID_trimmed.fastq SAMPLEID.fastq


# repeat fastqc and multiqc analyses on the trimmed reads
mkdir trimmed

mv *_trimmed.fastq ./trimmed

fastqc SAMPLEID_trimmed.fastq

export LC_ALL=en_US.UTF-8
export PATH=/programs/miniconda3/bin:$PATH
source activate multiqc

multiqc .

conda deactivate

## I noticed that some of the samples stil have adapter contamination and I'm not sure why this is?
## Also there are samples where the size is smaller than others (especially the negative and empty samples). I think that I should remove all sequences that are shorter than the expected length of the amplicons?
## the smallest amplicon is 101bp. Perhaps remove sequences that are less than 90bp?

## add the -m parameter to the cutadapt script. for now I'll just run cutadapt as follows:
cutadapt -m 90 -o 51068_trim2.fastq 51068_trimmed.fastq #(no adapter listed)

### quality filter of the reads
#100% of the bases in a sequence must have a score of higher than 10 for the sequence to be kept
fastq_quality_filter -q 10 -p 100 -Q33 -i 51068_trim2.fastq -o 51068_tf.fastq &

#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
fastq_quality_filter -q 20 -p 95 -Q33 -i 51068_trim2.fastq -o 51068_tff.fastq &

# I calculated the reads retained by doing a line count in the _tff.fastq files and dividing the #lines by 4 to get #reads
wc -l *_tff.fastq

# Run fastqc on these new files.

##### ALIGN TO REFERENCE #####

# prepare a "reference" - in this case this just means a list of the amplicon sequences in fasta format.
# from geneious export all the regions that Jen sent and create a single fasta file named OrioleAmpSeq.fasta
# notes this inlcudes regions larger than the actual amplicons. perpahs I need to subsample and recreate this reference.

# index reference file
bowtie2-build -f ./OrioleAmpSeq.fasta OrioleAmpSeq2 &

mkdir align2

# run bowtie2 alignment
# alignment was conducted using the following command (for each sample) and piped directly into samtools view (to convert sam file to bam file)
bowtie2 --phred33 --sensitive -p 8 -x OrioleAmpSeq2 --rg-id 51068 --rg SM:51068 -U /workdir/bgb27/trimmed/trimfilter/51068_tff.fastq|samtools view -bS - > /workdir/bgb27/align2/51068.bam

#used grep to get the alignment info from the nohup.out file
grep "aligned exactly 1 time" nohup.out

#sort the bam files using samtools sort
samtools sort SAMPLEID.bam -o SAMPLEID_sorted.bam


# run multibamqc
# need to make sure this command is run from within the directory with the samples
# this runs qualimap on individual samples AND the qc on all of the samples

mkdir qualimap
/programs/qualimap_v2.2.1/qualimap multi-bamqc -d ./multibamqc.txt -outdir /workdir/bgb27/aligned/qualimap -outfile multibam-qc-report -r

# -d directs to file with following structure (one for each sample)
SAMPLEID_sorted SAMPLEID_sorted.bam


##### CALL SNPS #####

# index reference genome
samtools faidx /workdir/bgb27/OrioleAmpSeq.fasta

# run mpileup generate VCF format // run call to generate SNP calls
bcftools mpileup -Ou --threads 7 -a "FORMAT/AD,FORMAT/DP,INFO/AD" -f /workdir/bgb27/OrioleAmpSeq.fasta 57675repeat_sorted.bam 51068_sorted.bam 52153_sorted.bam 52519_sorted.bam 52945_sorted.bam 53095_sorted.bam 53509_sorted.bam 54403_sorted.bam 55009_sorted.bam 55519_sorted.bam 56298_sorted.bam 56769_sorted.bam 56770_sorted.bam 56771_sorted.bam 56772_sorted.bam 56773_sorted.bam 56774_sorted.bam 56775_sorted.bam 56776_sorted.bam 56777_sorted.bam 56778_sorted.bam 56779_sorted.bam 56780_sorted.bam 56781_sorted.bam 56782_sorted.bam 56783_sorted.bam 56784_sorted.bam 56785_sorted.bam 56786_sorted.bam 56787_sorted.bam 56788_sorted.bam 56789_sorted.bam 56790_sorted.bam 56791_sorted.bam 56792_sorted.bam 56793_sorted.bam 56794_sorted.bam 56795_sorted.bam 56796_sorted.bam 56797_sorted.bam 56798_sorted.bam 56799_sorted.bam 56800_sorted.bam 56801_sorted.bam 56802_sorted.bam 56803_sorted.bam 56804_sorted.bam 56805_sorted.bam 56806_sorted.bam 56807_sorted.bam 56808_sorted.bam 56809_sorted.bam 56810_sorted.bam 56811_sorted.bam 56812_sorted.bam 56813_sorted.bam 56814_sorted.bam 56815_sorted.bam 56816_sorted.bam 56817_sorted.bam 56818_sorted.bam 56819_sorted.bam 56820_sorted.bam 56821_sorted.bam 56822_sorted.bam 56823_sorted.bam 56824_sorted.bam 56825_sorted.bam 56826_sorted.bam 56827_sorted.bam 56828_sorted.bam 56829_sorted.bam 56830_sorted.bam 56831_sorted.bam 56832_sorted.bam 56833_sorted.bam 56834_sorted.bam 56835_sorted.bam 56836_sorted.bam 56837_sorted.bam 56838_sorted.bam 56839_sorted.bam 56840_sorted.bam 56841_sorted.bam 56842_sorted.bam 56843_sorted.bam 56844_sorted.bam 56845_sorted.bam 56846_sorted.bam 56847_sorted.bam 56848_sorted.bam 56849_sorted.bam 56850_sorted.bam 56851_sorted.bam 56852_sorted.bam 56853_sorted.bam 56854_sorted.bam 56855_sorted.bam 56856_sorted.bam 56857_sorted.bam 56858_sorted.bam 56859_sorted.bam 56860_sorted.bam 56861_sorted.bam 56862_sorted.bam 56863_sorted.bam 56864_sorted.bam 56865_sorted.bam 56866_sorted.bam 56867_sorted.bam 56868_sorted.bam 56869_sorted.bam 56870_sorted.bam 56871_sorted.bam 56872_sorted.bam 56873_sorted.bam 56874_sorted.bam 56875_sorted.bam 56876_sorted.bam 56877_sorted.bam 56878_sorted.bam 56879_sorted.bam 56880_sorted.bam 56881_sorted.bam 56882_sorted.bam 56883_sorted.bam 56884_sorted.bam 56885_sorted.bam 56886_sorted.bam 56887_sorted.bam 56888_sorted.bam 56889_sorted.bam 56890_sorted.bam 56891_sorted.bam 57611_sorted.bam 57612_sorted.bam 57613_sorted.bam 57614_sorted.bam 57615_sorted.bam 57616_sorted.bam 57617_sorted.bam 57618_sorted.bam 57626_sorted.bam 57627_sorted.bam 57628_sorted.bam 57629_sorted.bam 57646_sorted.bam 57647_sorted.bam 57648_sorted.bam 57649_sorted.bam 57650_sorted.bam 57651_sorted.bam 57652_sorted.bam 57653_sorted.bam 57654_sorted.bam 57655_sorted.bam 57656_sorted.bam 57659_sorted.bam 57660_sorted.bam 57661_sorted.bam 57662_sorted.bam 57663_sorted.bam 57670_sorted.bam 57672_sorted.bam 57673_sorted.bam 57674_sorted.bam 57675_sorted.bam 57681_sorted.bam 57684_sorted.bam 57685_sorted.bam 57687_sorted.bam 57688_sorted.bam 57689_sorted.bam 57691_sorted.bam 57692_sorted.bam 57693_sorted.bam 57698_sorted.bam 57699_sorted.bam 57700_sorted.bam 57710_sorted.bam 57711_sorted.bam 57712_sorted.bam 57713_sorted.bam 57714_sorted.bam 57715_sorted.bam 57716_sorted.bam 57722_sorted.bam 57723_sorted.bam 57724_sorted.bam 57730_sorted.bam 57731_sorted.bam 57732_sorted.bam 57733_sorted.bam 57734_sorted.bam 57736_sorted.bam 57737_sorted.bam 57738_sorted.bam 57739_sorted.bam 57740_sorted.bam 57741_sorted.bam 57747_sorted.bam 57748_sorted.bam 57749_sorted.bam 57750_sorted.bam 57751_sorted.bam 57752_sorted.bam 57753_sorted.bam 57754_sorted.bam 57755_sorted.bam 57756_sorted.bam 57757_sorted.bam 57758_sorted.bam 57759_sorted.bam 57760_sorted.bam 57761_sorted.bam 57762_sorted.bam 57889_sorted.bam 57890_sorted.bam 57891_sorted.bam 57892_sorted.bam 57893_sorted.bam 57895_sorted.bam 57943_sorted.bam 57944_sorted.bam 57945_sorted.bam 57946_sorted.bam 57947_sorted.bam 57948_sorted.bam 57949_sorted.bam 57950_sorted.bam 57951_sorted.bam 57952_sorted.bam 57953_sorted.bam 57954_sorted.bam 57955_sorted.bam 57956_sorted.bam 57957_sorted.bam 57958_sorted.bam 57959_sorted.bam 57960_sorted.bam 57961_sorted.bam 57962_sorted.bam 57963_sorted.bam 57964_sorted.bam 57965_sorted.bam 57966_sorted.bam 57968_sorted.bam 57969_sorted.bam 57970_sorted.bam 57971_sorted.bam 57972_sorted.bam 57973_sorted.bam 57974_sorted.bam 57975_sorted.bam 57976_sorted.bam 57977_sorted.bam 57978_sorted.bam 57979_sorted.bam 57980_sorted.bam 57981_sorted.bam 57982_sorted.bam 57983_sorted.bam 57984_sorted.bam 57989_sorted.bam 57991_sorted.bam 57993_sorted.bam 57994_sorted.bam 57995_sorted.bam 57996_sorted.bam 57997_sorted.bam 57998_sorted.bam 57999_sorted.bam 58000_sorted.bam 58001_sorted.bam 58002_sorted.bam 58003_sorted.bam 58004_sorted.bam 58005_sorted.bam 58006_sorted.bam 58007_sorted.bam 58008_sorted.bam 58009_sorted.bam 58010_sorted.bam 58011_sorted.bam 58012_sorted.bam 58013_sorted.bam 58014_sorted.bam 58015_sorted.bam 58016_sorted.bam 58017_sorted.bam 58018_sorted.bam 58021_sorted.bam 58026_sorted.bam 58046_sorted.bam 58047_sorted.bam 58048_sorted.bam 58049_sorted.bam 58119_sorted.bam 58120_sorted.bam 58121_sorted.bam 58122_sorted.bam 58123_sorted.bam 58124_sorted.bam 58125_sorted.bam 58126_sorted.bam 58127_sorted.bam 58128_sorted.bam 58129_sorted.bam 58137_sorted.bam 58138_sorted.bam 58139_sorted.bam 58140_sorted.bam 58141_sorted.bam DMNH13057p1_1_sorted.bam DMNH13057p1_2_sorted.bam DMNH13057p1_3_sorted.bam DMNH13057p1_45_sorted.bam DMNH13057p1_sorted.bam DMNH13057p2_1_sorted.bam DMNH13057p2_2_sorted.bam DMNH13057p2_3_sorted.bam DMNH13057p2_45_sorted.bam DMNH13057p2_sorted.bam empty1_sorted.bam empty2_sorted.bam empty3_sorted.bam empty4_sorted.bam empty5_sorted.bam empty6_sorted.bam error_sorted.bam neg1_sorted.bam neg2_sorted.bam neg3_sorted.bam neg_p1_45_sorted.bam neg_p1_sorted.bam neg_p1_conc_sorted.bam neg_p2_45_sorted.bam neg_p2_conc_sorted.bam | bcftools call --variants-only --threads 7 -mv -Ov -o orioleampseq2.vcf &

#Getting stats about the sites and coverage: using vcf tools I will start with the full vcf and ask
#--site-depth Generates a file containing the depth per site summed across all individuals. This output file has the suffix ".ldepth".
#--depth Generates a file containing the mean depth per individual. This file has the suffix ".idepth".
vcftools --vcf orioleampseq2.vcf --site-depth --out orioleampseq2 &
vcftools --vcf orioleampseq2.vcf --depth --out orioleampseq2 &

#I also created a file that contains only the 151 SNPs (there are 2 SNPs on Orioles_20) that we designed the amplicons to capture. This file is a tab delimited file with chr and position. I used info from geneious to get these positions.
#I used this file to filter the original vcf to include only these sites.

vcftools --vcf orioleampseq2.vcf --positions all_SNPs --recode --out all_SNPs &
  #After filtering, kept 102 out of a possible 1083 Sites

vcftools --vcf all_SNPs.recode.vcf --site-depth --out all_SNPs &
vcftools --vcf all_SNPs.recode.vcf --depth --out all_SNPs &
vcftools --vcf all_SNPs.recode.vcf --missing-indv --out all_SNPs &
vcftools --vcf all_SNPs.recode.vcf --missing-site --out all_SNPs &

#--missing-indv Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
#--missing-site Generates a file reporting the missingness on a per-site basis. The file has the suffix ".lmiss".

#Open these output files (.ldepth and .idepth) in excel and look at the coverage...

#1. the following 23 amplicons are not present in the outputs from above: Oriole_1, Oriole_2, Oriole_6-7, Oriole_26, Oriole_27, Oriole_22, Oriole_35, Oriole_58, Oriole_59, Oriole_63, Oriole_72, Oriole_91, Oriole_93, Oriole_106, Oriole_124, Oriole_127, Oriole_130, Oriole_132, Oriole_134, Oriole_138, Oriole_141, Oriole_142, Oriole_146
#2. I checked in the output from qualimap (by using grep to create a file containing the coverage information and then searching that for the the amplicons above using grep again) - it looks like these have no coverage in any samples - so were probably not PCRd.
#3. There are about 27 sites (from the vcf that contains only the SNPS we are interested in that have a depth of less than 1000 (sum of the depths across all individuals). If we don't consider the negative controls, etc. we have 312 individuals - so a depth of 1000 would be = to approx a depth of 3 in each individual...I think...
#4. When looking at the depth across indviduals - the negative controls and most of the toepad samples all have a mean depth of less than 1.3. And there are 9 individuals that have a mean depth of less than 5.
#5. Looking at the missing data per individual - the negative controls and toepads have a lot of missing data. I will choose to remove all individuals with >

#Plot PCA in R.
#using the following vcf files:
#1. orioleampseq2.vcf
#2. all_SNPs.recode.vcf

###Filter the vcf:

#I decided that perhaps I should remove the samples with lots of missing data first...

#I looked at the output from the --missing-indv and selected all individuals with a F_MISS >0.5 and created a file with these (one line per individual).

vcftools --vcf all_SNPs.recode.vcf --remove indv_with_missing_data.txt --recode --out all_SNPs_rem_indv &

#check the missing data and depth for sites again: using the all_SNPs_rem_indv.recode.vcf file
vcftools --vcf all_SNPs_rem_indv.recode.vcf --site-depth --out all_SNPs_rem &
vcftools --vcf all_SNPs_rem_indv.recode.vcf --depth --out all_SNPs_rem &
vcftools --vcf all_SNPs_rem_indv.recode.vcf --missing-indv --out all_SNPs_rem &
vcftools --vcf all_SNPs_rem_indv.recode.vcf --missing-site --out all_SNPs_rem &

# >>21 sites have >50% missing data and these correspond to the sites that also have the lowest site depth (<300).
# >>all individuals have a mean depth of >10.

#Now I will filter for mean site depth and missing data.
#Based on the data above, perhaps i need to filter the data in the following ways:
#1. remove samples/sites with too much missing data.
#2. removed samples with a mean depth of less than 5?

#I am hoping that this will remove all the negative controls...

--max-missing <float> Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
--minDP <float> Includes only genotypes greater than or equal to the "--minDP" value. This option requires that the "DP" FORMAT tag is specified for all sites.

vcftools --vcf all_SNPs_rem_indv.recode.vcf --max-missing 0.5 --minDP 5 --recode --out all_SNPs_filtered &
  #After filtering, kept 300 out of 300 Individuals
  #Outputting VCF file...
  #After filtering, kept 68 out of a possible 102 Sites

#Create PCA...

### filter to include only positions/snps from each group (background, fixed, inversion, melanin)

vcftools --vcf all_SNPs_filtered.recode.vcf --positions fixed_SNPs --recode --out fixed_SNPs &
  # After filtering, kept 17 out of a possible 68 Sites (out of 30 amplicons)

vcftools --vcf all_SNPs_filtered.recode.vcf --positions melanin --recode --out melanin_SNPs &
  #After filtering, kept 26 out of a possible 68 Sites (out of 57 sites, 56 amplicons)

vcftools --vcf all_SNPs_filtered.recode.vcf --positions inversion_SNPs --recode --out inversion_SNPs &
  #After filtering, kept 10 out of a possible 68 Sites (out of 34 amplicons)

vcftools --vcf all_SNPs_filtered.recode.vcf --positions background_SNPs --recode --out background_SNPs &
  #After filtering, kept 15 out of a possible 68 Sites (out of 30 amplicons)


# Plot PCAs in R: for background, fixed, inversion and melanin vcf files.
# note: I always lose 1 snp when creating the gds file? - probably becuase I did not filter for biallelic SNPs in vcftools
# note: we have reduced the snps that we have in our analysis to a total of 65 sites (only 43% of the sites we were hoping to get)
