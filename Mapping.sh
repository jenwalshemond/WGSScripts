#Script for mapping raw reads to a reference genome
#Modified on January 20, 2021

##the-bowtie2-aligner http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner

#First you need to run bowtie2-build to create a bowtie index for your chosen reference.
bowtie2-build -f pathtoREFGENOME NAMEofOUTPUT
# -f: Reference is in FASTA format. FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar.
#note that the output will be a number of files, they will all have the base name that you specify in NAMEofOUTPUT

#Next, align your reads to the reference genome
bowtie2 -p 16 --phred33 --very-sensitive-local -x NAMEofOUTPUT -I 149 -X 900 --rg-id SAMPLEID --rg SM:SAMPLEID -1 /PATH/LIBRARY_NAME.pair1.truncated -2 /PATH/LIBRARY_NAME.pair2.truncated -U /PATH/LIBRARY_NAME.collapsed,/PATTH/LIBRARY_NAME.collapsed.truncated,/PATH/LIBRARY_NAME.singleton.truncated -S ./SAMPLEID.sam
# x = NAMEofOUTPUT from step above (base name of your indexed reference)
# p = number of threads (this can be 16 or 24 on the medium machines at BioHPC)
# phred33: refers to the way the quality information is stored in the files - illumina data is --phred33
# --very-sensitive-local is the same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 (these are options that determine the alignment parameters and effort that the program puts into identifying alignments)
# I: The minimum fragment length for valid paired-end alignments. This will depend on the size of the fragments in your library and the length of the reads. If your reads are 150 and are completely overlapping then the min fragement length would be 150.
# X: The maximum fragment length for valid paired-end alignments. This also depends on the size of the fragments in your library and the length of the reads. eg. If your reads are 150bp and are located 600bp apart then they will be considered valid if X=900
#    Note: The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences between -I and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.
# --rg-id and --rg - these are required for the sam file and adds your SAMPLEID to the read group tag in the sam file. SM means that you are providing the sample name.
# -1 Comma-separated list of files containing first read of mate paired reads (this will be the pair1.truncated file output from AdapterRemoval)
# -2 Comma-separated list of files containing second read of mate paired reads (this will be the pair2.truncated file output from AdapterRemoval)
# -U Comma-separated list of files containing unpaired reads to be aligned (these will be the collapsed, collapsed.truncated and singleton.truncated files output from AdapterRemoval)
# -S File to write SAM alignments to (if you don't specity this your alignments will be written to the standard out file)

# use samtools to convert the human readable sam file into a binary fomat bam file (that takes up less space and is needed for downstream programs)
samtools view -bS SAMPLEID.sam > SAMPLEID.bam
# use samtools to sort the bam file such that the alignments occur in “genome order”. That is, ordered positionally based upon their alignment coordinates on each chromosome.
samtools sort SAMPLEID.bam -o SAMPLEID_sorted.bam
#remove the sam and unsorted bam files, becuase we don't need them and they just take up space.
rm ./SAMPLEID.sam
rm ./SAMPLEID.bam

#note you will need to run bowtie2 and samtools on each of your samples, therefore create a script to do this. (note 2: you only need to index your reference once)

##Qualimap
#http://qualimap.conesalab.org/doc_html/index.html
/programs/qualimap_v2.2.1/qualimap bamqc -bam SAMPLEID_sorted.bam -outfile SAMPLEID.pdf
# bamqc reports information for the evaluation of the quality of the provided alignment data (a BAM file). In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced. This analysis can be performed with any kind of sequencing data, e.g. whole-genome sequencing, exome sequencing, RNA-seq, ChIP-seq, etc.
#-bam: input to bam file
#-outfile: name of output file
