#!/bin/bash
#this little bash script can be used to to create a directory per fastq file in the output directory 


### parse command line arguments
while getopts ":i:t:r:m:" OPTION; do

case "$OPTION" in
	i) run_path=$OPTARG;;
	t) threads=$OPTARG;;
	r) ref_path=$OPTARG;;
	m) mutserve=$OPTARG;;
	?)
		echo "script usage: $(basename "$0") [-m MUTSERVE.jar] [-i RUN_PATH] [-r REFERENCE_PATH] [-t THREADS]" >&2
		exit 1
		;;
esac
done
shift $(($OPTIND -1))

ont_path="$run_path"/ont_pl

get_bam(){
	bam="$1"/align/final/final.bam
}

get_ref(){
	fragment=$(awk '{ print $3 }' "$1"/*.bed)
	ref="$ref_path"*"$fragment"*.fasta
}

get_contig(){
	contig_name=$(awk '{ print $1 }' "$1"/*.bed)
}

make_output_dir(){
	output_dir="$ont_path"/mutserve
	mkdir -p $output_dir/$1
}

mutserve(){
	echo $output
	echo $ref
	echo $contig_name
	output=$output_dir/$1/$1_mutserve.txt
	java -jar $mutserve call --reference $ref --contig-name $contig_name --insertions --deletions --write-raw --output $output $bam
}

add_barcode_to_output() {
	sed -i -e "s/final.bam/$1/g" $output_dir/$1/*.txt
}

merge_files(){
	rm $output_dir/$(basename $run_path)_summary_mutserve.txt
	cat $output_dir/barcode*/*_raw.txt | head -n1 >> $output_dir/$(basename $run_path)_summary_mutserve.txt
	for barcode in $output_dir/*/; do
		cat $barcode/*_raw.txt | tail -n+2 >> $output_dir/$(basename $run_path)_summary_mutserve.txt
	done
}

main(){
	get_ref $ont_path
	get_contig $ont_path
	for barcode in "$ont_path"/barcode*; do
		get_bam $barcode
		make_output_dir $(basename $barcode)
		mutserve $(basename $barcode)
		add_barcode_to_output $(basename $barcode)
	done
	merge_files

}

main
