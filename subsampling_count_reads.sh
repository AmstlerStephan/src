#!/bin/bash
#this little bash script can be used to subsample fastq or fasta files
#samples_max is the maximum reads to be subsampled and must be smaller than all reads
#samples_min is the minimum read number and must be greater zero
#sample_rate is the ratio with which to subsample

### set parameters
subsampling_dir=subsampling

count_reads() {
	sample=$1
	subsample=$2
	sample_name=$(basename $sample)
	subsample_name=$(basename $subsample)
	count=$(cat $subsample/*fastq | wc -l)
	echo "$sample_name	$subsample_name	$count" >> subsampling_read_count_2.txt
}

add_header(){
	echo "sample_name	subsample	count" >> subsampling_read_count_2.txt
}

loop_dirs(){
	for sample in $subsampling_dir/A*; do
		for subsample in $sample/*; do
			count_reads $sample $subsample
		done
	done
}

add_header
loop_dirs
