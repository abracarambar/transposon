#!/bin/bash
BAM=$1
TRANSP_TYPE=$2
FILE_LOC=/g/data1a/jp48/scripts/MELTv2.1.4/${TRANSP_TYPE}_with_priors.txt
GROUPOUTPUTDIR=/g/data1a/jp48/MELT/GroupAnalysis/${TRANSP_TYPE}/

java -Xmx4G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar Genotype \
     -bamfile $BAM \
     -t $FILE_LOC \
     -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
     -w $OUTPUTDIR \
     -p $GROUPOUTPUTDIR;
