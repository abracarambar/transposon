#!/bin/bash

######################################################################
#
# Maely Gauthier 19/07/2018
# This is an adaptation of scripts developed by Eugene gardner
# see http://melt.igs.umaryland.edu/manual.php#_Other_MELT_Tools
#
######################################################################


# Preprocess generates three files, all as suffixes to the current bam file:
#sorted.bam.disc discordant pairs from the current BAM file
#sorted.bam.disc.bai index of these discordant pairs
#sorted.bam.disc.fq fastq version of all discordant pairs for MEI alignment

#These will be added to the folder where the bam files reside



##Get the bam name and coverage from the job submission command-line like: ./IndivAnalysis.sh <Preprocessed_Sample.bam>

BAM=$1

java -Xmx4G  -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar Preprocess -bamfile $BAM -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta
