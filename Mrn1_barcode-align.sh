#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/NIKS022/fastq-gz'
GZDIR=${HOME}${DIR}
WORKDIR="${GZDIR}/work"

MHALO='GAGATCCAGTCACTCGGtcagaCTAT'
MFL='GAGATCCAGTCACTCGGcagctCTAT'
MFRAG='GAGATCCAGTCACTCGGccttgCTAT'

mkdir -p "${WORKDIR}"

SAMPLEFILE="${GZDIR}/samples.txt"

rm -f "${SAMPLEFILE}"
for GZ in `ls "${GZDIR}"/*.fastq.gz`
do
    SAMPLE=`basename ${GZ}.fastq.gz`
    echo "${SAMPLE}" >> "${SAMPLEFILE}"
done

for i in $(cat "${SAMPLEFILE}")
do
    if [[ ! -e "${WORKDIR}/${i}mhalo.fastq.gz" ]]
    then
	cutadapt -a mhalo="${MHALO}" -a mfl="${MFL}" -a mfrag="${MFRAG}" --minimum-length 10 \
		 -o ${WORKDIR}/${i}{name}.fastq.gz ${GZDIR}/${i}.fastq.gz &
#	sleep 1
    fi
done

wait

PROM='mhalo mfl mfrag'

for i in $(cat "${SAMPLEFILE}")
do
    for j in ${PROM}
    do
        if [[ ! -e "${i}${j}-count.txt" ]]
        then
	    zcat ${WORKDIR}/$i$j.fastq.gz \
		 | /mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-count --fastq - --output $i${j}-count.txt --neighborhood $i${j}
        fi
    done
done

MHALOCOUNTS=$(echo ${GZDIR}*mhalo-count.txt)
MFLCOUNTS=$(echo ${GZDIR}*mfl-count.txt)
MFRAGCOUNTS=$(echo ${GZDIR}*mfrag-count.txt)

if [[ ! -e "${GZDIR}halocount.txt" ]]
then
    /mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-tabulate -o ${GZDIR}mhalo_seq1_counts.txt ${MHALOCOUNTS}
    /mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-tabulate -o ${GZDIR}mfl_seq1_counts.txt ${MFLCOUNTS}
    /mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-tabulate -o ${GZDIR}mfrag_seq1_counts.txt ${MFRAGCOUNTS}
fi


