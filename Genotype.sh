#!/bin/bash

######################################################################
#
# Maely Gauthier 19/07/2018
# This is an adaptation of scripts developed by Eugene gardner
# see http://melt.igs.umaryland.edu/manual.php#_Other_MELT_Tools
#
######################################################################

BAM=$1
TRANSP_TYPE=$2
OUTPUTDIR=$3
FILE_LOC=/g/data1a/jp48/scripts/MELTv2.1.4/${TRANSP_TYPE}_with_priors.txt
OUTPUTDIR=/g/data1a/jp48/MELT/${TRANSP_TYPE}

java -Xmx4G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar Genotype \
     -bamfile $BAM \
     -t $FILE_LOC \
     -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
     -w $OUTPUTDIR \
     -p $OUTPUTDIR;
