#!/bin/bash
TRANSP_TYPE=$1
TRANSP_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/me_refs/1KGP_Hg19/${TRANSP_TYPE}_MELT.zip
PRIORS_FILE=/g/data1a/jp48/scripts/MELTv2.1.4/prior_files/${TRANSP_TYPE}.1KGP.sites.vcf
#provide location of files
FILE_LOC=/g/data1a/jp48/scripts/MELTv2.1.4/${TRANSP_TYPE}_with_priors.txt
#INDOUTPUTDIR=/g/data1a/jp48/MELT/IndivAnalysis/${TRANSP_TYPE}
#GROUPOUTPUTDIR=/g/data1a/jp48/MELT/GroupAnalysis/${TRANSP_TYPE}
OUTPUTDIR=/g/data1a/jp48/MELT/${TRANSP_TYPE}

java -Xmx6G -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar GroupAnalysis\
     -discoverydir $OUTPUTDIR \
     -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta \
     -n /g/data1a/jp48/scripts/MELTv2.1.4/add_bed_files/1KGP_Hg19/hg19.genes.bed \
     -r 150 \
     -t $FILE_LOC \
     -w $OUTPUTDIR
