# Script for variant calling in Samtools
# Modified on February 17, 2021

# http://www.htslib.org/doc/bcftools.html#mpileup
# Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files. This is based on the original samtools mpileup command (with the -v or -g options) producing genotype likelihoods in VCF or BCF format, but not the textual pileup output. The mpileup command was transferred to bcftools in order to avoid errors resulting from use of incompatible versions of samtools and bcftools when using in the mpileup+bcftools call pipeline.
# Individuals are identified from the SM tags in the @RG header lines. Multiple individuals can be pooled in one alignment file, also one individual can be separated into multiple files. If sample identifiers are absent, each input file is regarded as one sample.

# -O: output type. u = bcf
# -C: Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50
# -a: annotate. list of common INFO tags for vcf output
# -f: Reference genome
# -multiallelic-caller: alternative model for multialleleic and rare-variant calling designed to overcome known limitations in -c calling mode (previous version)
# -variants-only: output variant sites only
# -O: output type. v=VCF
# -o: output file name

bcftools mpileup -Ou -C50 -threads 12 -a "FORMAT/AD,FORMAT/DP,INFO/AD,FORMAT/DV,FORMAT/DP4,FORMAT/DPR,INFO/DPR" -f REFERENCE.fasta INDIVIDUAL1_sortedRGmark.bam INDIVIDUAL2_sortedRGmark.bam INDIVIDUAL..N_sortedRGmark.bam  | bcftools call --multiallelic-caller --variants-only --threads 12 -vm -Ov -o FILENAME.vcf

# Filter VCF
vcftools --vcf RAWSNPs.vcf --max-missing 0.8 --maf 0.05 --min-meanDP 2 --max-meanDP 50 --minQ 20 --min-alleles 2 --max-alleles 2 --recode --out FilteredVCF.vcf
bgzip FilteredVCF.vcf
tabix Filtered.vcf.gz
