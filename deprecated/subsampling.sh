#!/bin/bash
#this little bash script can be used to subsample fastq or fasta files
#samples_max is the maximum reads to be subsampled and must be smaller than all reads
#samples_min is the minimum read number and must be greater zero
#sample_rate is the ratio with which to subsample
input_file=$1
output_path=$2
samples_max=$3
sample_rate=$4
samples_min=$5

while [ $samples_max -ge $samples_min ]
do
seqtk sample $input_file $samples_max > $output_path/subsample_"$samples_max".fastq
samples_max=$(((samples_max * sample_rate) / 100))
done
