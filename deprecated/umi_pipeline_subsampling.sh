#!/bin/bash
#this little bash script can be used to to create a directory per fastq file in the output directory 


### parse command line arguments i, o or u
while getopts ":i:o:u:t:r:b:" OPTION; do

case "$OPTION" in
i) input_path=$OPTARG;;
o) output_path=$OPTARG;;
u) umi_pipeline_path=$OPTARG;;
t) threads=$OPTARG;;
r) ref=$OPTARG;;
b) bed=$OPTARG;;

?)
echo "script usage: $(basename "$0") [-i INPUT_PATH] [-o OUTPUT_PATH] [u- UMI_PIPELINE_PATH] [r- REFERENCE_PATH] [b- BED_PATH] [-t THREADS]" >&2
exit 1
;;
:)
echo "script usage: $(basename "$0") [-i INPUT_PATH] [-o OUTPUT_PATH] [u- UMI_PIPELINE_PATH] [r- REFERENCE_PATH] [b- BED_PATH] [-t THREADS]" >&2
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
for fastq in "$work_dir"/*.fastq; do

mkdir -p "$work_dir"output/$(basename "${fastq}" .fastq)

done
}

run_pipeline() {

mkdir "$work_dir"output/

pipe_work_dir=$(basename "$input_path")/
pipe_ref="$pipe_work_dir"ref/$(basename "$ref")
pipe_bed="$pipe_work_dir"ref/$(basename "$bed")

cd "$umi_pipeline_path"

for dir in "$pipe_work_dir"output/*; do
filename=$(basename "$dir")

#echo "$filename"
snakemake -j "$threads" -pr --configfile config.yml \
--config input_fastq="$pipe_work_dir""$filename".fastq reference_fasta="$pipe_ref" sample_name="$filename" targets_bed="$pipe_bed" --use-singularity

cp -r "$filename" "$pipe_work_dir"output/

done

cd "$start_dir"
}

create_output() {

output_dir="$output_path"ont_pl

#echo "$output_dir"

mkdir "$output_dir"
cp -r "$work_dir" "$output_dir"

}


clean_up(){

for fastq in "$work_dir"/*.fastq; do

rm -r "$umi_pipeline_path"$(basename "${fastq}" .fastq)

done 

rm -r "$work_dir"

}

echo "$work_dir"
echo "$input_path"
echo "$output_path"
echo "$threads"


#conda init bash
#conda activate pipeline-umi-amplicon
#copy_data
#create_dir
run_pipeline
create_output
clean_up
#conda deactivate
