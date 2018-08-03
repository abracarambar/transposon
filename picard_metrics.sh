#!/bin/bash

######################################################################
#
# Maely Gauthier 19/07/2018
#by default base and mapping quality is 20.
######################################################################

module load java

cd $PBS_O_WORKDIR
BAM=$1

java -jar /g/data1a/jp48/scripts/picard/build/libs/picard.jar CollectWgsMetrics \
      I=$BAM \
      O=$BAM.metrics \
      REFERENCE_SEQUENCE="/g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta" \
      VALIDATION_STRINGENCY=LENIENT
