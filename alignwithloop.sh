
#load your sample names into an array. 
#in this example I will get the sample names from the fastq files I generated after a previous step.

INDS=($(for i in /workdir/bgb27/*.fastq; do echo $(basename -s .fastq $i); done))
#basename - will remove the directory path and returns the file name. -s tell it which suffix to remove from the end of the file name (in this case .fastq)
#note: the variable INDS will now contain an array of the sample names extracted from the files names.

#If you want to see what is stored in this variable you can type:
#echo ${INDS[@]}
#@=all the elements in the array

mkdir /workdir/bgb27/align4

for SAMPLEID in ${INDS[@]};
do
	#declare variables. This makes it easier and neater to write your command line and you just have to change these for future projects.
	REFERENCE=OrioleAmpSeq2
	TRIMMEDSEQ=/workdir/bgb27/trimmed/trimfilter/${SAMPLEID}_tff.fastq
	OUTPUT=/workdir/bgb27/align4/${SAMPLEID}.bam
	OUTPUTSORTED=/workdir/bgb27/align4/${SAMPLEID}_sorted.bam


	# align with bowtie - the output is piped directly into samtools to avoid having the intermediate .sam file.
	echo "Aligning $SAMPLEID with bowtie" 
	#this just writes a line telling you which sample is being worked on. 
	bowtie2 --phred33 --very-sensitive-local -p 8 -x $REFERENCE --rg-id $SAMPLEID --rg SM:$SAMPLEID -U $TRIMMEDSEQ|samtools view -bS - > $OUTPUT
	
	#sort the bam files using samtools sort
	samtools sort $OUTPUT -o $OUTPUTSORTED

done

#this loop above will take each element in the INDS variable and use that to replace the file names in the script. 