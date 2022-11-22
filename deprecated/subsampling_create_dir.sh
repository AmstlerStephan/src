#!/bin/bash
#this little bash script can be used to to create a directory per fastq file in the output directory 

input_path=$1
output_path=$2

for fastq in $input_path/*.fastq;
do
mkdir $output_path/$(basename "${fastq##*/}" .fastq)
done

