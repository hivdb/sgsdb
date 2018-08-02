#! /bin/bash

for accs in $(tail -n +2 data/SGS.sequences.fact.csv | cut -f 2 -d ',' | xargs -n 500 echo | sed 's/ /,/g')
do
    echo -n '.'
    curl -sSL "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&sendto=on&id=$accs" >> local/SGS.sequences.fas 
done
echo
