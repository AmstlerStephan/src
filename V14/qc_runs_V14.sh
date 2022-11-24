#!/bin/bash

## script to filter and compare original plus filtered read quality for all runs
## plus directly summarize all QC data with an R script
min_length=500

combine_filter_reads(){
	mkdir -p $outdir
	for barcode in $1/fastq_pass/barcode*; do
		base_barcode=$(basename $barcode)
		mkdir $outdir/$base_barcode
		catfishq $barcode/* -l $min_length | gzip > $outdir/$base_barcode/"$base_barcode"_raw.fastq.gz
	done
}

filter_reads(){
	zcat $1/*raw* | NanoFilt -l $min_length | gzip > $1/$2_filtered.fastq.gz
}

stat_reads(){
	for fastq in $1/*; do
		echo $fastq
		NanoStat --fastq $fastq --tsv --outdir $3/stats/$2 -n $(basename $fastq .fastq.gz) -t 20
	done
}

analyse_barcodes(){
	mkdir $1/stats
	for barcode in $1/barcode*; do
		basename_barcode=$(basename $barcode)
		#filter_reads $barcode $basename_barcode
		stat_reads $barcode $basename_barcode $1
	done
}

loop_runs(){
	for run in ./run*; do
		echo $run
		outdir="QC/Nanostats"/$(basename $run)
		combine_filter_reads $run $outdir
		analyse_barcodes $outdir
	done
}

qc_summary_R(){
	Rscript src/V14/qc_runs_V14.R
}

loop_runs
qc_summary_R
