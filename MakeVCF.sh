#!/bin/bash
TRANSP_TYPE=$1
#OUTPUTDIR=/g/data1a/jp48/MELT/GroupAnalysis/Genotype_${TRANSP_TYPE}
TRANSP_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/${TRANSP_TYPE}_MELT.zip
#WORKDIR=/g/data1a/jp48/MELT/IndivAnalysis/${TRANSP_TYPE}
OUTPUTDIR=/g/data1a/jp48/MELT/${TRANSP_TYPE}

java -Xmx6G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar MakeVCF \
    -genotypingdir $OUTPUTDIR \
    -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
    -t $TRANSP_FILE \
    -w $OUTPUTDIR \
    -p $OUTPUTDIR \
    -o $OUTPUTDIR;
