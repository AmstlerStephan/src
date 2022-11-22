#!/bin/bash

while getopts "i:o:u:" OPTION; do

case "$OPTION" in
i)
input_path=$OPTARG
echo "input: " $input_path
;;

o)
output_path=$OPTARG
echo "output: " $output_path
;;

u)
umi_pipeline_path=$OPTARG
echo "pipeline: " $umi_pipeline_path
;;

?)
echo "script usage: $(basename "$0") [-i PATH]" >&2
exit 1
;;
esac
done

shift $(($OPTIND -1))
