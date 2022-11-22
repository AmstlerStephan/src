#!/bin/bash
# count lines of the consensus_final.bam files

count_consensus_sequences(){

	printf "Barcode \t Consenus_sequences"
	for bam_file in run*/ont_pl/b*/align/*final.bam; do
		echo $bam_file
		n_cons=$(samtools view $bam_file | wc -l)
		echo $(echo $bam_file | awk "/" '{print $3}')
		echo $n_cons
	done

}

count_consensus_sequences
