#!/bin/sh
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ] || [ "$5" == "" ]; then
	echo "Usage: input_folder STAR_index threads output_prefix GTF_file"
	exit
fi
# define variables from inputs
folder=$1
index=$2
threads=$3
out=$4
gtf=$5
# load STAR dependency
module load STAR
# load featureCounts dependency
module load Subread
echo "Starting RNA-seq analysis" > $out.log
# loop over all fastq files in folder
fastqlist=$folder/*_R1_*.fastq.gz
for read1 in $fastqlist
do
	# find names for read2 (skip if single-end) and sample
	read2=$(echo $read1 | sed 's/_R1_/_R2_/g')
	name=$(basename $read1 "_L001_R1_001.fastq.gz")
	echo "Analyzing sample $name" >> $out.log
	# run STAR (if single-end, skip $read2)
	STAR --genomeDir $index --runThreadN $threads --readFilesCommand zcat --readFilesIn $read1 $read2 --outFileNamePrefix $out.$name --outSAMtype BAM Unsorted
	featureCounts -T 8 -p -t exon -g gene_id -s 0 -a $gtf -o $out.$name.counts $out.$name\Aligned.out.bam
	# reformat our outputs to combine later
	grep -v "#" $out.$name.counts | cut -f1 > $out.genes.txt
	echo "$name" > $out.$name.counts.txt
	grep -v "#" $out.$name.counts | cut -f7 | sed '1d' >> $out.$name.counts.txt
done
# paste files together
paste $out.genes.txt $out.*.counts.txt > $out.txt