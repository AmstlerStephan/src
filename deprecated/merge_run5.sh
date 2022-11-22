#!/bin/bash
# merge files with one header line in different directories



file_type=run5/run5_minion/sequencing_summary*
echo $file_type

cat $file_type | head -n1 >> run5/sequencing_summary_run5.txt
for summary in run5/run5*/sequencing_summary*; do
        echo $summary
        cat $summary | tail -n+2 >> run5/sequencing_summary_run5.txt

done


