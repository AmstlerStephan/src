#!/bin/bash
#this little bash script can be used to subsample fastq or fasta files
#samples_max is the maximum reads to be subsampled and must be smaller than all reads
#samples_min is the minimum read number and must be greater zero
#sample_rate is the ratio with which to subsample

### parse command line arguments
while getopts ":i:s:u:f:l:" OPTION; do

case "$OPTION" in
        i) input_dir=$OPTARG;;
        s) output_path=$OPTARG;;
        u) upper_lim_reads=$OPTARG;;
        f) sampling_factor=$OPTARG;;
        l) lower_lim_reads=$OPTARG;;
        ?)
                echo "script usage: $(basename "$0") [-i input directory (fastq_pass/barcodeXX/)] [-s name of subsampling directory] [-u upper limit reads] [-l lower limit reads] [-f sampling factor (%)]" >&2
                exit 1
                ;;
esac
done
shift $(($OPTIND -1))

### set parameters
fastq_dir=fastq_pass
barcode=$(basename $input_dir)
raw=raw
raw_fastq=$output_path/$fastq_dir/$barcode/$raw/$raw.fastq

do_subsampling() {
	while [ $upper_lim_reads -ge $lower_lim_reads ]; do
		create_output_dir $barcode $upper_lim_reads
		seqtk sample $raw_fastq $upper_lim_reads > $output_path/$fastq_dir/$barcode/$upper_lim_reads/subsample_"$upper_lim_reads".fastq
		echo $upper_lim_reads
		upper_lim_reads=$(((upper_lim_reads * sampling_factor) / 100))
	done
}

create_output_dir(){
	mkdir -p $output_path/$fastq_dir/$1/$2
}

merge_input_into_subsampling_dir(){
	zcat $input_dir* > $raw_fastq
}

echo $barcode
echo $upper_lim_reads
### merge and copy raw fastq file
create_output_dir $barcode $raw
merge_input_into_subsampling_dir
### subsampling
do_subsampling

