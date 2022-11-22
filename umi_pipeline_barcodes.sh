#!/bin/bash
#this little bash script can be used to to create a directory per fastq file in the output directory 


### parse command line arguments
while getopts ":i:o:u:t:r:b:h:" OPTION; do

case "$OPTION" in
	i) input_path=$OPTARG;;
	o) output_path=$OPTARG;;
	u) umi_pipeline_path=$OPTARG;;
	t) threads=$OPTARG;;
	r) ref=$OPTARG;;
	b) bed=$OPTARG;;

	?)
		echo "script usage: $(basename "$0") [-i INPUT_PATH] [-o OUTPUT_PATH] [-u UMI_PIPELINE_PATH] [-r REFERENCE_PATH] [-b BED_PATH] [-t THREADS]" >&2
		exit 1
		;;
	h)
		echo "INPUT_PATH	Path to a directory containing the barcodes and fastq or fastq.gz files"
		echo "OUTPUT_PATH	path to the dire"
		exit 1
		;;
esac
done
shift $(($OPTIND -1))

### initialize needed variables
work_dir="$umi_pipeline_path"$(basename "$input_path")/
start_dir=$PWD

copy_data() {
	cp -r "$input_path" "$umi_pipeline_path"
	mkdir "$work_dir"ref
	cp "$ref" "$bed" "$work_dir"ref/
}

create_dir() {
	for barcode in "$work_dir"barcode*/; do
		mkdir -p "$work_dir"output/$(basename "${barcode}")
	done
}

run_pipeline() {
	### declare relative paths to use within the container 
	pipe_work_dir=$(basename "$input_path")/
	pipe_ref="$pipe_work_dir"ref/$(basename "$ref")
	pipe_bed="$pipe_work_dir"ref/$(basename "$bed")
	
	cd "$umi_pipeline_path"
	
	for barcode in "$pipe_work_dir"barcode*/; do
		filename=$(basename "$barcode")
	
		mkdir "$barcode"data
		cat "$barcode"*fastq* > "$barcode"data/"$filename".fastq
		fastq_file="$barcode"data/"$filename".fastq

		snakemake -j "$threads" -pr --configfile config.yml \
		--config input_fastq="$fastq_file" reference_fasta="$pipe_ref" sample_name="$filename" targets_bed="$pipe_bed" --use-singularity

		cp -r "$filename"/* "$pipe_work_dir"output/"$filename"

	done

	cd "$start_dir"
}

create_output() {

	output_dir="$output_path"ont_pl

	mkdir "$output_dir"
	cp -r "$work_dir"output/* "$output_dir"

}


clean_up(){

	for barcode in "$work_dir"barcode*/; do
		rm -r "$umi_pipeline_path"$(basename "${barcode}")
	done

	rm -r "$work_dir"

}

copy_data
create_dir
run_pipeline
create_output
clean_up
