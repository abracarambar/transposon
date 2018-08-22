#!/bin/bash

SAMPLE=$1
NORMAL=$2
TUMOUR=$3
ELEMENT=$4
OUTDIR=$5
TRANS_DIR=${OUTDIR}/${ELEMENT}
VCF=${ELEMENT}.final_comp.vcf

module load bcftools
module load htslib

##This step filters:
##Chr hs37d5
##FIILTER ac0; lc
##calls that have LP+RP<4
##SR<4
##ASSESS<4
##calls that are 0/0 in both the reference and the tumour sample

if [ ! -d "${TRANS_DIR}/$VCF" ]; then
    echo "Need to compress and index VCF file first"
    bgzip ${TRANS_DIR}/$VCF
    tabix ${TRANS_DIR}/$VCF.gz
fi

#VCF_base=`basename $VCF .vcf`

bcftools view -t ^hs37d5  -s $NORMAL,$TUMOUR -f '%LINE\n' -e '%FILTER=="ac0" || %FILTER=="lc" || LP+RP<4 || SR<4 || ASSESS<4' -U ${TRANS_DIR}/$VCF.gz  | sed 's/0\/0\:/0\/0\ /g' | sed 's/0\/1\:/0\/1\ /g' | sed 's/1\/1\:/1\/1\ /g' | sed 's/1\/0\:/1\/0\ /g' | sed 's/.\/.\:/.\/.\ /g' | tr  ' ' '\t' | awk '{if (($10=="0/0" && $12=="0/1") || ($10=="0/1" && $12=="0/0") || ($10=="0/1" && $12=="./.") || ($10=="0/1" && $12=="1/1") || ($10=="1/1" && $12=="0/1")) print}' > ${TRANS_DIR}/${VCF}_${SAMPLE}_somatic.vcf 

#sleep 5

echo 'Total calls recovered:' 
wc -l ${TRANS_DIR}/${VCF}_${SAMPLE}_somatic.vcf

