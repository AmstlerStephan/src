#!/bin/bash
#this little bash script can be used to to create a directory per fastq file in the output directory 


### parse command line arguments i, o or u
while getopts ":i:o:u:t:r:b:" OPTION; do

case "$OPTION" in
i) input_path=$(realpath $OPTARG);;
o) output_path=$(realpath $OPTARG);;
u) umi_pipeline_path=$(realpath $OPTARG);;
t) threads=$OPTARG;;
r) ref=$(realpath $OPTARG);;
b) bed=$(realpath $OPTARG);;

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
start_dir=$PWD

create_dirs() {
for fastq in $input_path/*.fastq; do

mkdir $output_path/$(basename "${fastq}" .fastq)

done
}

run_pipeline() {
for dir in $output_path/*/; do

filename=$(basename $dir)

echo $input_path
echo $output_path
echo $umi_pipeline_path
echo $bed
echo $ref
echo $filename
echo $dir
echo "end"

cd $dir

snakemake --snakefile $umi_pipeline_path/Snakefile -j $threads reads -pr --configfile $umi_pipeline_path/config.yml \
--config input_fastq=$input_path/$filename.fastq reference_fasta="$ref" sample_name=$filename targets_bed=$bed \
--verbose

cd $start_dir

done
}


#conda init bash
#conda activate pipeline-umi-amplicon
create_dirs
run_pipeline
#conda deactivate











