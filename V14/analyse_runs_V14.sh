#!/bin/bash

### script is used to analyse all sequencing runs including using the QC data and run info to summarize all data into one dataframe 

### parse command line arguments
while getopts ":o:" OPTION; do

case "$OPTION" in
        o) outdir=$OPTARG;;
        ?)
                echo "script usage: $(basename "$0") [-o output folder]" >&2
                exit 1
                ;;
esac
done
shift $(($OPTIND -1))


analyse_runs(){
	for run in run*V14/; do
		Rscript src/V14/Analyse_Run_v1.R $run
	done
}

merge_result_files(){
	bash src/V14/merge_result_files.sh -o $outdir
}

summarize_data(){
	Rscript src/V14/summarize_data.R $outdir
}


analyse_runs
merge_result_files
summarize_data
