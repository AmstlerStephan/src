#!/bin/bash

### script to convert all user defined BAM-files into SAM files in given input directory

### parse command line arguments
while getopts ":i:" OPTION; do

case "$OPTION" in
        i) input_dir=$OPTARG;;
        ?)
                echo "script usage: $(basename "$0") [-i INPUT_DIR]" >&2
                exit 1
                ;;
esac
done
shift $(($OPTIND -1))

out_dir="$input_dir"/SAM
BAM_dir="$input_dir"/BAM

make_SAM(){
	run=$1
	sample=$2
	BAM=$3
	sample_type=$(basename $BAM .bam)

  	samtools view -o $out_dir/"$sample"_"$sample_type"_"$run".sam $BAM
}

loop_BAM_files(){
  for run in $BAM_dir/run*; do
    for sample in $run/*/; do
	for BAM in $sample*.bam; do
		make_SAM $(basename $run) $(basename $sample) $BAM
	done
    done
  done
}

mkdir -p $out_dir
loop_BAM_files
