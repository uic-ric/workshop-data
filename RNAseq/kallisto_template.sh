#!/bin/sh
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ]; then
	echo "Usage: input_folder kallisto_index threads output_prefix"
	exit
fi
# define variables from inputs
folder=$1
index=$2
threads=$3
out=$4
# load kallisto dependency
module load kallisto
echo "Starting RNA-seq analysis" > $out.log
# loop over all fastq files in folder
fastqlist=$folder/*_R1_*.fastq.gz
for read1 in $fastqlist
do
	# find names for read2 (skip if single-end) and sample 
	read2=$(echo $read1 | sed 's/_R1_/_R2_/g')
	name=$(basename $read1 "_L001_R1_001.fastq.gz")
	echo "Analyzing sample $name" >> $out.log
	# run kallisto for paired-end data
	kallisto quant -t $threads -i $index -o $out.$name $read1 $read2
	#command for single-end data: the fragment length (-l) and standard deviation (-s) are required: kallisto quant -t $threads -i $index -o $out.$name --single -l <fragment_length> -s <std_dev> $read1
	# reformat our outputs to combine later
	cut -f1 $out.$name/abundance.tsv > $out.genes.txt
	echo "$name" > $out.$name/counts.txt
	cut -f4 $out.$name/abundance.tsv | sed '1d' >> $out.$name/counts.txt
done
# paste files together
paste $out.genes.txt $out.*/counts.txt > $out.txt