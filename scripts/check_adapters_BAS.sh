#!/usr/bin/env bash
# depends on cutadapt and revseq

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

RC1=$(echo $A1| revseq -sequence stdin -outseq stdout|grep -v "^>") 
RC2=$(echo $A2| revseq -sequence stdin -outseq stdout|grep -v "^>")

C1=$(echo -e ">seq\n$A1" | seqkit seq -p -v -t dna | grep -v "^>")
C2=$(echo -e ">seq\n$A2" | seqkit seq -p -v -t dna | grep -v "^>")

R1=$(echo -e ">seq\n$A1" | seqkit seq -r -v -t dna | grep -v "^>")
R2=$(echo -e ">seq\n$A2" | seqkit seq -r -v -t dna | grep -v "^>")

pinfo "A1=$A1"
pinfo "RC_A1=$RC1"
pinfo "A2=$A2"
pinfo "RC_A2=$RC2"
pinfo "C_A1=$C1"
pinfo "C_A2=$C2"
pinfo "R_A1=$R1"
pinfo "R_A2=$R2"

function run_cut {
    cutadapt -g "$1; max_error_rate=$2;" $fastq_file -o /dev/null| grep "Reads with adapters:"| cut -f 2 -d: | cut -f 1 -d\( | tr -d " ," 
}
function print_csv2tsv {
    echo $@ | tr "," "\t"
}


#print_csv2tsv "err,A1-RC2,C2-C1,C1-R2,A2-RC1,C2-R1,"

for error_rate in 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4; do
    X1=$(run_cut "A1=${A1}...$A2" "$error_rate")
    X2=$(run_cut "A1=${A1}...$R2" "$error_rate")
    X3=$(run_cut "A1=${A1}...$C2" "$error_rate")
    X4=$(run_cut "A1=${A1}...$RC2" "$error_rate")

    X5=$(run_cut "A1=${R1}...$A2" "$error_rate")
    X6=$(run_cut "A1=${R1}...$R2" "$error_rate")
    X7=$(run_cut "A1=${R1}...$C2" "$error_rate")
    X8=$(run_cut "A1=${R1}...$RC2" "$error_rate")

    X9=$(run_cut "A1=${C1}...$A2" "$error_rate")
    X10=$(run_cut "A1=${C1}...$R2" "$error_rate")
    X11=$(run_cut "A1=${C1}...$C2" "$error_rate")
    X12=$(run_cut "A1=${C1}...$RC2" "$error_rate")

    X13=$(run_cut "A1=${RC1}...$A2" "$error_rate")
    X14=$(run_cut "A1=${RC1}...$R2" "$error_rate")
    X15=$(run_cut "A1=${RC1}...$C2" "$error_rate")
    X16=$(run_cut "A1=${RC1}...$RC2" "$error_rate")

    X17=$(run_cut "A1=${A2}...$A1" "$error_rate")
    X18=$(run_cut "A1=${A2}...$R1" "$error_rate")
    X19=$(run_cut "A1=${A2}...$C1" "$error_rate")
    X20=$(run_cut "A1=${A2}...$RC1" "$error_rate")

    X21=$(run_cut "A1=${R2}...$A1" "$error_rate")
    X22=$(run_cut "A1=${R2}...$R1" "$error_rate")
    X23=$(run_cut "A1=${R2}...$C1" "$error_rate")
    X24=$(run_cut "A1=${R2}...$RC1" "$error_rate")

    X25=$(run_cut "A1=${C2}...$A1" "$error_rate")
    X26=$(run_cut "A1=${C2}...$R1" "$error_rate")
    X27=$(run_cut "A1=${C2}...$C1" "$error_rate")
    X28=$(run_cut "A1=${C2}...$RC1" "$error_rate")

    X29=$(run_cut "A1=${RC2}...$A1" "$error_rate")
    X30=$(run_cut "A1=${RC2}...$R1" "$error_rate")
    X31=$(run_cut "A1=${RC2}...$C1" "$error_rate")
    X32=$(run_cut "A1=${RC2}...$RC1" "$error_rate")

    print_csv2tsv "$error_rate,$X1,$X2,$X3,$X4,$X5,$X6,$X7,$X8,$X9,$X10,$X11,$X12,$X13,$X14,$X15,$X16,$X17,$X18,$X19,$X20,$X21,$X22,$X23,$X24,$X25,$X26,$X27,$X28,$X29,$X30,$X31,$X32"
done
exit 0
