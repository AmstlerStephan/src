#!/bin/bash
# merge files with one header line in different directories

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


merge_result_files(){
	for file_type in run5/results/*csv; do
		file_type=$(basename $file_type)
		echo $file_type
		cat run5/results/$file_type | head -n1 >> $outdir/temp_all_$file_type
		for files in run*/results/$file_type; do
			echo $files
			cat $files | tail -n+2 >> $outdir/temp_all_$file_type
        	done

	done

}

echo $outdir
#filtered="filtered"
#raw="raw"
#figures="figures"
#NGS="NGS"

rm -r $outdir
mkdir -p $outdir

#mkdir -p $outdir/$figures/reads_vs_cons_consensus_sequences
#mkdir -p $outdir/$figures/$raw/$NGS/Bland_Altman
#mkdir -p $outdir/$figures/$raw/$NGS/comp_per_pos
#mkdir -p $outdir/$figures/$raw/$NGS/comp_variant_level
#mkdir -p $outdir/$figures/$raw/Plasmid
#mkdir -p $outdir/$figures/$filtered/$NGS/Bland_Altman
#mkdir -p $outdir/$figures/$filtered/$NGS/comp_per_pos
#mkdir -p $outdir/$figures/$filtered/$NGS/comp_variant_level
#mkdir -p $outdir/$figures/$filtered/Plasmid
#mkdir -p $outdir/Plasmid_tables

merge_result_files
