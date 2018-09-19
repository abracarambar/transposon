#!/bin/bash

######################################################################
#
# Maely Gauthier 19/07/2018
# This is an adaptation of scripts developed by Eugene gardner
# see http://melt.igs.umaryland.edu/manual.php#_Other_MELT_Tools
#
######################################################################


##Get the bam name and coverage from the job submission command-line like: ./IndivAnalysis.sh <Preprocessed_Sample.bam>
BAM=$1
COV=$2
TRANSP_TYPE=$3
BASEDIR=$4
OUTPUTDIR=${BASEDIR}/${TRANSP_TYPE}
TRANSP_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/${TRANSP_TYPE}_MELT.zip
PRIORS_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/prior_files/${TRANSP_TYPE}.1KGP.sites.vcf

java -Xmx6G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar IndivAnalysis \
     -t $TRANSP_FILE $PRIORS_FILE \
     -w $OUTPUTDIR \
     -c $COV \
     -bamfile $BAM \
     -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
     -b NC_007605/hs37d5 \
     -r 150
