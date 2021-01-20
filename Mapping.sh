Commands used to index the reference genome, -f means fasta input
bowtie2-build -f /workdir/jlw395/final.assembly.fasta SALS_bowtie

bowtie2 -p 16 --phred33 --very-sensitive-local -x SALS_bowtie -I 149 -X 900 --rg-id SAVSOG1 --rg SM:SAVSOG1 -1 /workdir/jlw395/LIBRARY_SAVSOG1.pair1.truncated -2 /workdir/jlw395/LIBRARY_SAVSOG1.pair2.truncated -U /workdir/jlw395/LIBRARY_SAVSOG1.collapsed,/workdir/jlw395/LIBRARY_SAVSOG1.collapsed.truncated,/workdir/jlw395/LIBRARY_SAVSOG1.singleton.truncated -S ./SAVSOG1.sam
samtools view -bS SAVSOG1.sam > SAVSOG1.bam
samtools sort SAVSOG1.bam -o SAVSOG1toSALS_sorted.bam
rm ./SAVSOG1.sam
rm ./SAVSOG1.bam