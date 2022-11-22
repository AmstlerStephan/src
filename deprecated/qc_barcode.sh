#!/bin/bash
# script to do the whole qc for one run and if flag -a (attach to summary) is set it will be attached to the overview of all runs

### parse command line arguments

while getopts ":i:a:t:h:" OPTION; do

case "$OPTION" in
        i) input_path=$OPTARG;;
        a) attach_to_summary=$OPTARG;;
        t) threads=$OPTARG;;
        ?)
                echo "script usage: $(basename "$0") [-i PATH_TO_RUN] [-a ATTACH_TO_OVERALL_SUMMARY] [-t THREADS]" >&2
                exit 1
                ;;
        h)
                echo "INPUT_PATH        Path to a directory containing the barcodes and fastq or fastq.gz files"
                echo "OUTPUT_PATH       path to the dire"
                exit 1
                ;;
esac
done
shift $(($OPTIND -1))

output_path="$input_path"QC
sequencing_summary="$input_path"sequencing_summary_*.txt

nanoplot() {

	mkdir $output_path/Nanoplot
	NanoPlot -t $threads --tsv_stats --info_in_report --barcoded --summary $sequencing_summary -o $output_path/Nanoplot
}

summarize_qc() {

	Rscript src/qc_barcode.R $input_path $output_path
}

add_to_overall_summary() {

	Rscript src/qc_add_to_summary.R

}

mkdir $output_path
nanoplot
summarize_qc

if [ $attach_to_summary ]
then
	echo "attach to summary"
fi
