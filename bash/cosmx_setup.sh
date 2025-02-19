#!/bin/bash

# REMINDER: make script executable by running `chmod +x ~/devel/spatialhub/bash/cosmx_setup.sh
# Usage: cosmx_setup.sh --directory < projDir > --slide < input_slide >


### Get options

echo "Parsing options"

vars=$(getopt -o d:s: --long directory:,slide: -- "$@")
eval set -- "$vars"

# extract options and their arguments into variables.
for opt; do
    case "$opt" in
      -d|--directory)
        projDir=$2
        shift 2
        ;;
      -s|--slide)
        slideName=$2
        shift 2
        ;;
    esac
done

echo "Project directory: ${projDir}"
echo "Slide unique ID: ${slideName}"

cd ${projDir}/${slideName}/flatFiles/


### Expression matrix

mkdir -p split_exprMat
cd split_exprMat

# split by FOV
gzip -dk ../*_exprMat_file.csv.gz 
awk -F\, '{print>$1}' ../*_exprMat_file.csv
  # print (and save) first field (= first column, or FOV variable) of a CSV (\,) delimited file 

# add back .csv extension to generated files
for i in $(find .); do mv $i "$i.csv"; done

# clarify file names
for i in *.csv; do mv "$i" "${slideName}_exprMat_FOV${i}"; done

# retrieve header from first file
for i in $(find . -type f -name "*_exprMat_*.csv" -not -name "*FOVfov.csv");              
    do echo -e "$(head -1 *FOVfov.csv)\n$(cat $i)" > $i;
done

# remove header only file
rm *FOVfov.csv

# check all FOVs are present
ls -1 | wc -l

echo "Done splitting expression matrix file."
rm ../*_exprMat_file.csv


### Transcript file

cd ..
mkdir -p split_txFile
cd split_txFile

# split by FOV
gzip -dk ../*_tx_file.csv.gz 
awk -F\, '{print>$1}' ../*_tx_file.csv

# add back .csv extension to generated files
for i in $(find .); do mv $i "$i.csv"; done

# clarify file names
for i in *.csv; do mv "$i" "${slideName}_tx_FOV${i}"; done

# retrieve header from first file
for i in $(find . -type f -name "*_tx_*.csv" -not -name "*FOVfov.csv");
    do echo -e "$(head -1 *FOVfov.csv)\n$(cat $i)" > $i;
done

# remove header only file
rm *FOVfov.csv

# check all FOVs are present
ls -1 | wc -l

echo "Done splitting transcripts file."
rm ../*_tx_file.csv 
