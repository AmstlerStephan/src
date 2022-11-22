#!/bin/bash
### script to do rebasecalling for all runs

config_file="dna_r9.4.1_450bps_sup.cfg"
reads_per_fastq=4000
min_qscore=9
gpu_device=auto
output_dir_root=rebasecalling_guppy_6_1_5

basecall_runs() {
	mkdir $output_dir_root
	for run in run*; do
		output_dir=$output_dir_root/$run
		input_dir=$run/fast5_pass

		mkdir -p $output_dir

		echo $config_file
		echo $reads_per_fastq
		echo $min_qscore
		echo $gpu_device
		echo $output_dir_root
		echo $output_dir
		echo $input_dir

		guppy_basecaller -i $input_dir -s $output_dir -c $config_file -q $reads_per_fastq_file --min_qscore $min_qscore --device $gpu_device --compress_fastq

		done
}

basecall_runs
