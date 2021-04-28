
#Phylogenetic Analyses - Trees in RAxML and PAUP
##April 28, 2021


##Filtering data and making input files
##Make tab file
cat INPUTFILE.vcf | vcf-to-tab > OUTPUTFILE.tab

##Swap out indels for missing data
##Scripts are frorm L. Campagna. I am not including these on GitHub, as they aren't mine. But email me if you need perl scripts (jlw395@cornell.edu)
##Before running script below, open fix_tab_file.pl in a text editor and change the name/path for your input and output file names
perl fix_tab_file.pl

##Write Fasta File
perl vcf_tab_to_fasta_alignment.pl -i INPUTFILE.tab > OUTPUTFILE.fixed.fasta

##convert fasta file to phylip format. 
perl fasta2relaxedPhylip.pl -f INPUT.fixed.fasta -o OUTPUT.fixed.phy

#Run RAxML. It will fail, but it will print conflicting sites to nohup. Record how many invariable sites are removed
nohup /programs/RAxML-8.2.4/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -x 12345 -p 12345 -# 10 -m ASC_GTRGAMMA --asc-corr=lewis -s INPUT.fixed.phy -n make_exclude_file &

##Take make_exclude_file.pl and make a file that captures those sites
perl exclude_file_maker.pl

#this second run will make a new phy file without the SNPs in exclude.txt
/programs/RAxML-8.2.4/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -x 12345 -p 12345 -# 10 -m ASC_GTRGAMMA --asc-corr=lewis -E exclude.txt -s INPUT.fixed.phy -n excluded

#So now you have a phy with only the important positions and you're ready to run raxml-ng or svdquartets. USE THIS INPUT FILE FOR BOTH PROGRAMS

#To run RAxML-ng. Medium memory machine.
/programs/raxml-ng_v0.9.0/raxml-ng --bootstrap --threads 40 --bs-trees 100 --msa INPUT.fixed.phy.exclude.txt --tree pars{3},rand{3} --model GTR+G+ASC_LEWIS --seed 1575495751

#Map support values from existing set of replicate trees
/programs/raxml-ng_v0.9.0/raxml-ng --support --tree bestML.tree --bs-trees bootstraps.tree


##Run SVDquartets
##Use the same input file from RAxML, but import into geneious and convert to an interleaved NEXUS file
/programs/paup4a166/paup
execute INPUT.fixed.nex;
SVDQuartets evalQuartets=all bootstrap nreps=100 nthreads=64;
savetrees file=OUTPUTSVD.tre;
savetrees from=1 to =1 savebootp=nodelabels file=fulldatasetbootstrap.tre;