#!/bin/bash

### script to loop over all directories of the subsampling dir
### and analyze them with the umi pipeline
### and do variant calling with mutserve afterwards


### Set variables needed in the script (might change this to command line input)
subsampling_dir=sub_test
threads=25
umi_pipeline_dir=pipeline-umi-amplicon
ref_dir_mutserve=ref
mutserve=src/mutserve_LPA_adapted.jar
start_dir=$PWD
fastq_dir=fastq_pass

### UMI PIPELINE

loop_subsampling_dirs() {
	fragment=$1
	for dir in $subsampling_dir/$fastq_dir/*$fragment/; do
		out_dir=$(basename $dir)
		input_dir=$dir
		umi_pipeline_subsampling $input_dir $out_dir $fragment
	done
}

umi_pipeline_subsampling(){
	input_dir=$1
	out_dir=$2
	fragment=$3

	copy_data $input_dir
	run_pipeline $input_dir $fragment
	create_output $out_dir
	clean_up
}

copy_data() {
	input_dir=$1

	echo $input_dir

	mkdir -p $umi_pipeline_dir/data
	mkdir -p $umi_pipeline_dir/output

        cp -r $input_dir $umi_pipeline_dir/data/
        cp -r ref/ $umi_pipeline_dir
}

run_pipeline() {
        ### declare relative paths to use within the container
	input_dir=$1
        fragment=$2

	pipe_work_dir=data/$(basename $input_dir)
        pipe_ref=ref/lpa-ref$fragment.fasta
        pipe_bed=ref/lpa-ref$fragment.bed

        cd $umi_pipeline_dir

        for sample in $pipe_work_dir/*; do
		#### if pipeline is not running anymore switch back to subsample_$(basename $sample).fastq and rename the raw folder to subsample_raw instead of raw!!!
		input_fastq=$sample/*.fastq
		#### The pipeline is not able to use wildcards, echo does..
		input_fastq=$(echo $sample/*.fastq)
 		filename=$(basename $input_fastq .fastq)
		echo $input_fastq
		echo $filename

                snakemake -j $threads -pr --configfile config.yml \
                --config input_fastq=$input_fastq reference_fasta=$pipe_ref sample_name=output/$filename targets_bed=$pipe_bed --use-singularity

        done

        cd $start_dir
}

create_output() {
        output_dir=$subsampling_dir/ont_pl/$1

        mkdir -p $output_dir
        cp -r $umi_pipeline_dir/output/* $output_dir
}

clean_up(){
        rm -r $umi_pipeline_dir/data
	rm -r $umi_pipeline_dir/output
	rm -r $umi_pipeline_dir/ref
}

### MUTSERVE

loop_ont_pl_result_dirs() {
	for sample in $subsampling_dir/ont_pl/*; do
		for subsample in $sample/*/; do
			mutserve $(basename $sample) $(basename $subsample) $subsample
		done
	done
}

mutserve(){
	sample_name=$1
	subsample_name=$2
	subsample=$3
	output_dir=$(get_output_dir $sample_name $subsample_name)
        echo $output_dir
	mkdir -p $output_dir
	output=$output_dir/"$subsample_name"_mutserve.txt
        java -jar $mutserve call --reference $(get_ref $subsample) --contig-name $(get_contig $subsample) --write-raw  --output $output $(get_bam $subsample)
}

get_bam(){
	subsample=$1
        echo $subsample/align/*final*.bam
}

get_ref(){
	subsample=$1
        fragment=$(awk '{ print $3 }' $subsample/targets.bed)
        echo $ref_dir_mutserve/lpa-ref$fragment.fasta
}

get_contig(){
	subsample=$1
        echo $(awk '{ print $1 }' $subsample/targets.bed)
}

get_output_dir(){
	sample_name=$1
	subsample_name=$2
        output_dir=$subsampling_dir/mutserve/$sample_name/$subsample_name
	echo $output_dir
}


loop_subsampling_dirs "2645"
loop_subsampling_dirs "5104"
loop_ont_pl_result_dirs
