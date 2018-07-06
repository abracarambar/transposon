#!/bin/bash
#PBS -V
#PBS -N MELT
#PBS -j oe
#PBS -m ae
#PBS -M m.gauthier@garvan.org.au
#PBS -q express
#PBS -l ncpus=8
#PBS -l mem=4gb
#PBS -l walltime=24:00:00
 
module load java
module load bowtie2
 
##Get the bam name and coverage from the job submission command-line like: ./IndivAnalysis.sh <Preprocessed_Sample.bam>
BAM=$1
 
java -Xmx4G  -jar /g/data1a/jp48/scripts/MELTv2.1.4/MELT.jar Preprocess -bamfile $BAM -h /g/data1a/jp48/scripts/human_g1k_v37_decoy.fasta
