#!/bin/bash

######################################################################
#
# Maely Gauthier 19/07/2018
# This is an adaptation of scripts developed by Eugene gardner
# see http://melt.igs.umaryland.edu/manual.php#_Other_MELT_Tools
#
######################################################################

TRANSP_TYPE=$1
BASEDIR=$2
OUTPUTDIR=${BASEDIR}/${TRANSP_TYPE}
TRANSP_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/${TRANSP_TYPE}_MELT.zip

java -Xmx6G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar MakeVCF \
    -genotypingdir $OUTPUTDIR \
    -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
    -t $TRANSP_FILE \
    -w $OUTPUTDIR \
    -p $OUTPUTDIR \
    -o $OUTPUTDIR;
