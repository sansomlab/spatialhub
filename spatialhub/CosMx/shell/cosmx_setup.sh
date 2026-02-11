#!/bin/bash

# REMINDER: make script executable by running `chmod +x ~/devel/spatialhub/bash/cosmx_setup.sh

set -euo pipefail

### Get options

usage() { 
	echo ""
	echo -e "Usage: $0 <Required arguments>\nRequired arguments:"
	echo -e "\t-s\tCosMx slide name."
	echo -e "\t-o\tOutput directory."
	echo -e "\t-m\tExpression matrix file (e.g., *_exprMat_file.csv.gz)."
	echo -e "\t-t\tTranscript file (e.g., *_tx_file.csv.gz)." 1>&2
	echo ""
	exit
}

[[ $# -eq 0 ]] && usage

while getopts "s:o:m:t:" opt; do
	case $opt in
		s)	slide=$OPTARG
			;;
		o)	outdir=$OPTARG
			;;
		m)	expmat=$OPTARG
			;;
		t)	txfile=$OPTARG
			;;
		*)	echo -e ${RED}"ERROR: invalid option provided"${NOCOLOR}
			usage
			;;
	esac
done

echo ""
echo "Output directory: ${outdir}"
echo "Slide unique ID: ${slide}"
echo "Expression matrix path: ${expmat}"
echo "Transcript file path: ${txfile}"
echo ""

[[ ! -f "${expmat}" ]] && { echo "ERROR: Expression matrix not found: ${expmat}"; exit 1; }
[[ ! -f "${txfile}" ]] && { echo "ERROR: Transcript file not found: ${txfile}"; exit 1; }

### Expression matrix
mkdir -p ${outdir}/split_exprMat
cd ${outdir}/split_exprMat

# create files named by the first field (FOV variable)
header=$(gzip -dc "${expmat}" | head -1 || true)
gzip -dc "${expmat}" | 
	awk -v slide="${slide}" -F',' 'NR > 1 {print > slide"_exprMat_FOV"$1".csv"}'

for f in $(find . -type f -name "*_exprMat_FOV*.csv"); do
	tmpf="${f}.tmp"
    echo "$header" > "$tmpf"
    cat "$f" >> "$tmpf"
    mv -f "$tmpf" "$f"
done

echo "Done splitting expression matrix file."


### Transcript file
mkdir -p ${outdir}/split_txFile
cd ${outdir}/split_txFile

# split by FOV
header=$(gzip -dc "${txfile}" | head -1 || true)
gzip -dc "${txfile}" |
	awk -v slide="${slide}" -F',' 'NR > 1 {print > slide"_tx_FOV"$1".csv"}'

for f in $(find . -type f -name "*_tx_FOV*.csv"); do
	tmpf="${f}.tmp"
    echo "$header" > "$tmpf"
    cat "$f" >> "$tmpf"
    mv -f "$tmpf" "$f"
done

echo "Done splitting transcripts file."
