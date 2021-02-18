# Script for prepping mapped reads for variant caller
# Modified on February 17, 2021
# Scripts are an example for one individual. Will have to run for all samples in data set
# A new Bam file will be produced for each step. In the end, you only need to keep samplename_sortedRGmark.bam

# Bam files from mapping script above are indexed and sorted already

# Prepare the genome for GATK: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R= REFERENCE.fa O= REFERENCE.dict
samtools faidx REFERENCE.fa

# add read group information: Assigns all the reads in a file to a single new read-group.
# use samtools view -H sample.bam | grep '@RG' 
# samtools command above will print required fields for the AddOrReplaceReadGroups below
# RGLB: Read-Group library - can be the same for all samples
# RGPL: Read-Group platform
# RGPU: run barcode (will be different for each individual)
# RGSM: Read-Group sample name (use read group info from previous steps)
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar AddOrReplaceReadGroups INPUT=MAPPEDFILE.bam OUTPUT=MAPPED_sortedRG.bam RGID=S1 RGLB=lib1 RGPL=illumina RGPU=ADAPTERSEQUENCE RGSM=READGROUPSAMPLENAME_sorted.bam

# mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000

java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=MAPPED_sortedRG.bam OUTPUT=MAPPED_sortedRGmark.bam METRICS_FILE=SAMPLENAME.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

# index .bam files for HaplotypeCaller
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=MAPPED_sortedRGmark.bam