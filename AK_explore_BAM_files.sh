#!/bin/bash

### This is a unifying script to collect BAM files of the AK samples
### convert them into SAM files
### and extract a table with information about the consensus sequences

### parse command line arguments
while getopts ":o:" OPTION; do

case "$OPTION" in
        o) out_dir=$OPTARG;;
        ?)
                echo "script usage: $(basename "$0") [-o OUTPUT_DIR]" >&2
                exit 1
                ;;
esac
done
shift $(($OPTIND -1))

collect_BAM_files() {
  Rscript src/AK_collect_BAM_files.R $out_dir
}

make_SAM_files() {
  bash src/AK_make_SAM_files.sh -i $out_dir
}

analyse_SAM() {
  Rscript src/AK_analyse_SAM_files.R $out_dir
}

mkdir -p $out_dir/
collect_BAM_files
make_SAM_files
analyse_SAM
