#!/bin/bash
#Author: Ioannis Kampouris
#Purpose: Use cutadapt to remove primers in the R Server (UBUNTU OS)
#Dependencies a path for fastqc and 
# Primer Sequences According to NOVOGENE
# FW: CCTAYGGGRBGCASCAG RE:GGACTACNNGGGTATCTAAT


cd 2020
find -name  "*raw_1*"|sed 's/^..//'  > s_sample_list.txt #Change accordingly to your path 

filename='s_sample_list.txt'
echo "Start removing primers from the samples" 
mkdir -p trimmed_cutadapt
mkdir reports

while read i;
	do 
   SAMPLE=$(echo ${i} | sed "s/.raw_1\.fq\.gz//")
   echo "$SAMPLE"
   cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACNNGGGTATCTAAT   --discard-untrimmed \
   -o ${SAMPLE}.trimmed_cutadapt_1.fq -p   ${SAMPLE}.{SAMPLE}.raw_1.fq.gz   ${SAMPLE}.raw_2.fq.gz \
	 > ${SAMPLE}.log
    
	
	
done <"$filename"

find . -name '*.fq' -exec mv {} trimmed_cutadapt \;
find . -name '*.log' -exec mv {} reports \;
