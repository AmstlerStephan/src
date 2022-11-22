#!/bin/bash
# script to do the whole qc for all runs


output_path=QC/Nanostat_summary
threads=5

nanostat() {
	echo $(basename $1).txt
	NanoStat -t $threads --tsv --summary $2 -o $output_path -n $(basename $1).txt
}

summarize_qc() {
	Rscript src/qc_summary.R $1 $output_path
	rm $output_path/$(basename $1).txt
}

add_to_overall_summary() {
	Rscript src/qc_add_to_summary.R $output_path
}

loop_runs(){
        for run in run[^_]*/; do
		echo $run
		sequencing_summary="$run"sequencing_summary_*.txt
		nanostat $run $sequencing_summary
		summarize_qc $run
        done
}


mkdir -p $output_path
loop_runs
add_to_overall_summary
