#!/bin/bash

### Set variables
outdir=info
filename=Barcode_Sample_overview_summary.js
counter=1
DIRS=$(find run*/ -maxdepth 0 | wc -l)


### Loop through runs and attach Barcode_Sample_overview_files

echo '{ "Runs":{' >> $outdir/$filename

for run in run*/; do

	echo "\"$(basename $run)\":" >> $outdir/$filename
	cat $run/lib/Barcode* >> $outdir/$filename

	if(($counter<$DIRS)); then
		echo , >> $outdir/$filename
	fi

	counter=$(($counter + 1))

done

echo '}}' >> $outdir/$filename
