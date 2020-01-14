#!/usr/bin/env bash

function pinfo {
    echo "$*"  > /dev/stderr
}
function usage {
    echo "Usage: check_adapaters.sh Adapter1  Adapter2 fastq_file "
}

if [ "$3-" == "-" ]; then
    usage    
    exit 1
fi


# replace I by N
A1=$(echo $1|tr "I" "N")
A2=$(echo $2|tr "I" "N")
fastq_file=$3

#A1="GGWACWRGWTGRACWNTNTAYCCYCC"
#A2="TANACYTCNGGRTGNCCRAARAAYCA"
#fastq_file=/home/nf/Research/Projects/WIP/CIBIO/MetaEnv/Minion_data/Mussels/2019_September_001/demux/barcode85/fastq_runid_ff5266d30aff9e0873f8cd2f296a56bc4a0f86a7_0.fastq.gz

set -eux
set -o pipefail
set -eu

R1=$(echo $A1| revseq -sequence stdin -outseq stdout|grep -v "^>") 
R2=$(echo $A2| revseq -sequence stdin -outseq stdout|grep -v "^>")

pinfo "A1=$A1"
pinfo "RA1=$R1"
pinfo "A2=$A2"
pinfo "RA2=$R2"

function run_cut {
    cutadapt -g "$1; max_error_rate=$2;" $fastq_file -o /dev/null| grep "Reads with adapters:"| cut -f 2 -d: | cut -f 1 -d\( | tr -d " ," 
}
function print_csv2tsv {
    echo $@ | tr "," "\t"
}


print_csv2tsv "err,A1-A2,R2-R1,A1-R2,A2-R1"

for error_rate in 0.01 0.05 0.1 0.15 0.2 0.25 0.3; do
    X1=$(run_cut "A1=${A1}...$A2" "$error_rate")
    X2=$(run_cut "A1=${R2}...$R1" "$error_rate")
    X3=$(run_cut "A1=${A1}...$R2" "$error_rate")
    X4=$(run_cut "A1=${A2}...$R1" "$error_rate")
    print_csv2tsv "$error_rate,$X1,$X2,$X3,$X4"
done
exit 0
