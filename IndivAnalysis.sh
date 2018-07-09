#!/bin/bash

##Get the bam name and coverage from the job submission command-line like: ./IndivAnalysis.sh <Preprocessed_Sample.bam>
BAM=$1
COV=$2
TRANSP_TYPE=$3
TRANSP_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/${TRANSP_TYPE}_MELT.zip
#/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/LINE1_MELT.zip
#/g/data1a/jp48/scripts/MELTv2.1.4/prior_files/LINE1.1KGP.sites.vcf
PRIORS_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/prior_files/${TRANSP_TYPE}.1KGP.sites.vcf
OUTPUTDIR=/g/data1a/jp48/MELT/IndivAnalysis/${TRANSP_TYPE}

java -Xmx6G  -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar IndivAnalysis \
     -t $TRANSP_FILE $PRIORS_FILE \
     -w $OUTPUTDIR \
     -c $COV \
     -bamfile $BAM \
     -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
     -b NC_007605/hs37d5
