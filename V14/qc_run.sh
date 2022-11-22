#!/bin/bash

## script to filter and compare original plus filtered read quality

min_length=500
run_dir="run8"
outdir="QC/Nanostats/"$run_dir

combine_reads(){
	mkdir -p $outdir
	for barcode in $run_dir/fastq_pass/barcode*; do
		base_barcode=$(basename $barcode)
		mkdir $outdir/$base_barcode
		cat $barcode/* > $outdir/$base_barcode/"$base_barcode"_raw.fastq.gz
	done
}

filter_reads(){
	zcat $1/*raw* | NanoFilt -l $min_length | gzip > $1/$2_filtered.fastq.gz
}

stat_reads(){
	for fastq in $1/*; do
		echo $fastq
		NanoStat --fastq $fastq --outdir $outdir/stats/$2 -n $(basename $fastq .fastq.gz) -t 6 --tsv
	done
}

analyse_barcodes(){
	mkdir -p $outdir/stats
	for barcode in $outdir/barcode*; do
		#basename_barcode=$(basename $barcode)
		#filter_reads $barcode $basename_barcode
		stat_reads $barcode $basename_barcode
	done
}

#combine_reads
analyse_barcodes
