##Trimming and Adapter Removal: All Samples
##Bullock's Oriole (BUOR)
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 /workdir/jlw395/oriole2_R1.fastq.gz --file2 /workdir/jlw395/oriole2_R2.fastq.gz --trimns --trimqualities --minquality 10 --collapse --threads 16 --adapter-list for_adapter_removal_2.txt --basename LIBRARY_BUOR1
