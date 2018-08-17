#!/bin/env python
from Cheetah.Template import Template
import argparse
import os
import stat
import sys
import fileinput
import collections
import pandas as pd

######################################################################
#
# Maely Gauthier 19/07/2018
# This script runs the MELT tool on a cohort of samples on torque system
# see http://melt.igs.umaryland.edu/manual.php#_Other_MELT_Tools
#
#EXAMPLE OF COMMAND
#source activate localpython
#first derive coverage
#python pbs_make_all.py -r 'lisa' -u 'm.gauthier@garvan.org.au' -l '/g/data1a/jp48/MELT/bam_files4.txt' -c '/g/data1a/jp48/MELT/Coverage.txt' -o '/g/data1a/jp48/MELT/Lisa' -p '1'
#run the rest of the pipeline including somatic step
#python pbs_make_all.py -r '/g/data1a/jp48/scripts/transposon/ocscc' -u 'm.gauthier@garvan.org.au' -l '/g/data1a/jp48/scripts/transposon/TCGA_BAMs3.txt' -c '/g/data1a/jp48/scripts/transposon/TCGABAMcoverage.txt' -o '/g/data1a/jp48/scripts/transposon/ocscc' -s 'true'

######################################################################
###
#TO DO return error if cannot find coverage
#specify different reference fasta files
#specify different bam file suffixes
###

def main():
    parser = argparse.ArgumentParser(description='Create the parser')
    pbs_path, textfile_path, user, cov_path, outdir_path, picard, somatic = local_parser(parser, require_user=True)
    #print 'Folder where commands will be stored is ' + str(run) + '\nText file specifying the path of bam files is ' + str(textfile_path) + '\nUser email is ' + str(user) + '\nText file specifying coverage is' + str(cov_path) + '\nOutput folder is ' + str(outdir_path)
    email = user
    pbs_path = pbs_path + '/'


    if os.path.exists(textfile_path):
        SampleSheet_df = pd.read_csv(textfile_path, index_col=None, skip_blank_lines=True, header=0, skipinitialspace=True, dtype=str, delimiter="\t")
        #print SampleSheet_df
        samples_path = SampleSheet_df["bam_paths"].tolist()
        df_dict = SampleSheet_df.set_index('sample').T.to_dict('list')

    if not os.path.isdir(pbs_path):
        os.makedirs(pbs_path)
        #might alos need to set permission so everyone can acces files
    if picard:
        #picard_metrics
        for index, row in SampleSheet_df.iterrows():
            sample = row['sample']
            sample_path = row['bam_paths']
            generate_template_per_sample('8', #RAMs used <1 for bam with coverage 100X
                                '4', #CPUs
                                'picard',
                                'picard_metrics.sh',
                                sample,
                                sample_path,
                                pbs_path,
                                email,
                                '8', # took 2h for a bam with coverage 100X
                                'express',
                                cov_path,
                                outdir_path)

    else:

        #creates pbs scripts per sample and per TE
        generate_all_separate_pbs_scripts(SampleSheet_df, samples_path, pbs_path, email, cov_path, outdir_path)

        #create a single script to run all pbs jobs
        create_pbs_launch_script(SampleSheet_df, pbs_path)

    if somatic:
        generate_template_per_patient('4', #might have to bring this up to 8Gb
                    '1', #confirmation from developer that MELT will only use 1 core
                    'somatic_calls',
                    'filter_vcf_for_somatic_calls.sh',
                    pbs_path,
                    email,
                    '6', #optimise, will this be enough for ALU??
                    'normal',
                    outdir_path,
                    SampleSheet_df)
    # print patient
        # TE = ('LINE1', 'ALU', 'SVA')
        # for element in TE:
        #     cmd = "bash '" + full_script_name + "' '" + patient + "' '" + sample_path + "' '" + str(element) + "' '" + str(outdir_path) + "'"
        #     print 'Creating script for '  + step + ' & ' + str(element)
        #     name = '_'.join([str(sample), step, element])
        #     script_name = pbs_path + name + '_pbs.sh'
        #     out = os.path.join(pbs_output_dir, name)
        #     write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)
        #


def names_from_txt_file(textfile_path):
    #sample_no = 0
    sample_name_df = pd.read_table(textfile_path, sep='\t', header=None)
    #print sample_name_df
    samples = sample_name_df.iloc[:, 0].tolist()
    #print samples
    return samples

def bam_file_paths_from_txt_file(textfile_path):
    #bam_path = 0
    sample_path_df = pd.read_table(textfile_path, sep='\t', header=None)
    print sample_path_df
    sample_paths = sample_path_df.iloc[:, 0].tolist()
    #print sample_paths
    return sample_paths

def generate_all_separate_pbs_scripts(SampleSheet_df, samples_path, pbs_path, email, cov_path, outdir_path):

    #for sample_path in samples_path:
    for index, row in SampleSheet_df.iterrows():
        sample = row['sample']
        sample_path = row['bam_paths']

        #Preprocess step
        generate_template_per_sample('10', #RAMs used <1 for bam with coverage 100X #one TCGA bam file of Laveniya failed, required a lot of vmem ~10Gb....
                            '1', #CPUs
                            'preprocess',
                            'preprocess.sh',
                            sample,
                            sample_path,
                            pbs_path,
                            email,
                            '3', # took 2h for a bam with coverage 100X
                            'express',
                            cov_path,
                            outdir_path)

        #Individual analysis step
        generate_template_per_sample('8', #might have to bring this up to 8Gb; failed for both ALu and Line jobs at 6gb
                            '1', #confirmation from developer that MELT will only use 1 core
                            'ind_analysis',
                            'IndivAnalysis.sh',
                            sample,
                            sample_path,
                            pbs_path,
                            email,
                            '24', #optimise, will this be enough for ALU??
                            'normal',
                            cov_path,
                            outdir_path)

        #Group analysis
        generate_template('8', #failed for ALU with 6 Gb and 1 macthed pair sample
                         '1', #confirmation from developer that MELT will only use 1 core
                         'group_analysis',
                         'GroupAnalysis.sh',
                         pbs_path,
                         email,
                         '12', #optimise, will this be enough for ALU??
                         'normal',
                         cov_path,
                         outdir_path)

        #Genotype step
        generate_template_per_sample('4', #might have to bring this up to 8Gb
                            '1', #confirmation from developer that MELT will only use 1 core
                            'genotype',
                            'Genotype.sh',
                            sample,
                            sample_path,
                            pbs_path,
                            email,
                            '6', #optimise, will this be enough for ALU??
                            'normal',
                            cov_path,
                            outdir_path)

        #Make VCFs
        generate_template('6', #might have to bring this up to 8Gb
                         '1', #confirmation from developer that MELT will only use 1 core
                         'MakeVCF',
                         'MakeVCF.sh',
                         pbs_path,
                         email,
                         '4', #optimise, will this be enough for ALU??
                         'normal',
                         cov_path,
                         outdir_path)



def generate_template_per_sample(mem,
                                cpus,
                                step,
                                script_name,
                                sample,
                                sample_path,
                                pbs_path,
                                email,
                                time,
                                batchqueue,
                                cov_path,
                                outdir_path):
    '''
    writes the variable details of the pbs scripts to call all sample specific pbs scripts
    '''
    pbs_output_dir = pbs_path
    scripts_dir='/g/data1a/jp48/scripts/transposon/'

    full_script_name = os.path.join(scripts_dir, script_name)

    #give bam file here, not sample name
    wd = "echo Working directory is ${PBS_O_WORKDIR}\n" + \
         "cd ${PBS_O_WORKDIR}\n"

    if step == 'picard':
        cmd = "bash '" + full_script_name + "' '" + sample_path + "' '" + str(outdir_path) + "'"
        print 'Creating script for '  + step
        name = '_'.join([str(sample), step])
        script_name = pbs_path + name + '_pbs.sh'
        out = os.path.join(pbs_output_dir, name)
        write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)

    elif step == 'ind_analysis':
        #derive coverage specific to sample
        print 'Looking for sample matching ' + str(sample) + ' in coverage file\n'
        cov_for_sample = 0
        with open(cov_path, 'r') as fin:
            for line in fin:
                bits = line.strip().split('\t')
                sample_with_cov = bits[0]
                cov = bits[1]
                #you need to detect when sample has no coverage specified!!!!
                if sample_with_cov == sample:
                    print 'Found a match ' + str(sample_with_cov) + ' ' + str(cov)
                    cov_for_sample = cov
                    #write a script for each transcription element
                    #this is what is different for this step
                    TE = ('LINE1', 'ALU', 'SVA')
                    for element in TE:
                        cmd = "bash '" + full_script_name + "' '" + sample_path + "' '" + str(cov_for_sample) + "' '" + str(element) + "' '" + str(outdir_path) + "'"
                        print 'Creating script for '  + step + ' & ' + str(element)
                        name = '_'.join([str(sample), step, element])
                        script_name = pbs_path + name + '_pbs.sh'
                        out = os.path.join(pbs_output_dir, name)
                        #only for ALU analysis and
                        if element == 'ALU' and cov > 75 and "_T" in sample:
                            time = 36
                        write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)

    elif step == 'genotype':
        TE = ('LINE1', 'ALU', 'SVA')
        for element in TE:
            cmd = "bash '" + full_script_name + "' '" + sample_path + "' '" + str(element) + "' '" + str(outdir_path) + "'"
            print 'Creating script for '  + step + ' & ' + str(element)
            name = '_'.join([str(sample), step, element])
            script_name = pbs_path + name + '_pbs.sh'
            out = os.path.join(pbs_output_dir, name)
            write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)



    #preprocessing step only
    else:
        cmd = "bash '" + full_script_name + "' '" + sample_path + "'"
        name = '_'.join([str(sample), step])
        script_name = pbs_path + name + '_pbs.sh'
        out = os.path.join(pbs_output_dir, name)
        write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)

def generate_template(mem,
                    cpus,
                    step,
                    script_name,
                    pbs_path,
                    email,
                    time,
                    batchqueue,
                    cov_path,
                    outdir_path):
    '''
    writes the variable details of the pbs scripts to call all sample specific pbs scripts
    '''
    pbs_output_dir = pbs_path
    scripts_dir='/g/data1a/jp48/scripts/transposon/'
    full_script_name = os.path.join(scripts_dir, script_name)
    #give bam file here, not sample name
    wd = "echo Working directory is ${PBS_O_WORKDIR}\n" + \
         "cd ${PBS_O_WORKDIR}\n"

    #create a single script
    if step == 'group_analysis' or step == 'MakeVCF':
        TE = ('LINE1', 'ALU', 'SVA')
        for element in TE:
            cmd = "bash '" + full_script_name + "' '" + str(element) + "' '" + str(outdir_path) + "'"
            print 'Creating script for ' + step + ' & ' + str(element)
            name = '_'.join([step, element])
            script_name = pbs_path + name + '_pbs.sh'
            out = os.path.join(pbs_output_dir, name)
            write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)

def generate_template_per_patient(mem,
                    cpus,
                    step,
                    script_name,
                    pbs_path,
                    email,
                    time,
                    batchqueue,
                    cov_path,
                    SampleSheet_df):
    '''
    writes the variable details of the pbs scripts to call all sample specific pbs scripts
    '''
    #pbs_output_dir = pbs_path
    #scripts_dir='/g/data1a/jp48/scripts/transposon/'
    full_script_name = os.path.join(pbs_path, script_name)
    #give bam file here, not sample name
    wd = "echo Working directory is ${PBS_O_WORKDIR}\n" + \
         "cd ${PBS_O_WORKDIR}\n"

    patient_list = []
    for index, row in SampleSheet_df.iterrows():
        patient = row['patient']
        patient_list.append(patient)

    for patient in patient_list:
        tumour_sample = SampleSheet_df.loc[(SampleSheet_df['patient'] == patient) & (SampleSheet_df['test_tissue'] == "Y"), 'sample'].item()
        normal_sample = SampleSheet_df.loc[(SampleSheet_df['patient'] == patient) & (SampleSheet_df['test_tissue'] == "N"), 'sample'].item()
        TE = ('LINE1', 'ALU', 'SVA')
        for element in TE:
            cmd = "bash '" + full_script_name + "' '" + patient + "' '" + normal_sample + "' '" + tumour_sample + "' '" + element + ".final_comp.vcf'"
            #print cmd
            print 'Creating script for '  + step + ' & ' + str(element)
            name = '_'.join([str(patient), step, element])
            script_name = pbs_path + name + '_pbs.sh'
            out = os.path.join(pbs_path, name)
            write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)

def write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name):
    #provide the location of the template and the information to be filled in
    template_file = "/g/data1a/jp48/scripts/transposon/pbs_script.template"
    name_space = {'name' : name,
                  'email' : email,
                  'queue' : batchqueue,
                  'out' : out,
                  'cores' : cpus,
                  'mem' : mem,
                  'time' : time,
                  'command' : cmd,
                  'setwd' : wd,}

    template = Template(file=template_file, searchList=[name_space])
    with open(script_name, "w") as f:
        f.write(str(template))

def create_pbs_launch_script(SampleSheet_df, pbs_path):
    '''
    Manages dependencies in the PBS queue - what waits for what
    '''
    TE = ('LINE1', 'ALU', 'SVA')
    call_all_script = pbs_path + 'call_all_pbs.sh'
    with open(call_all_script, 'w') as fout:
        fout.write('#!/bin/bash' + '\n')
        sample_list = []
        patient_list = []
        for index, row in SampleSheet_df.iterrows():
            patient = row['patient']
            sample = row['sample']
            sample_path = row['bam_paths']

            sample_list.append(sample)
            patient_list.append(patient)

            fout.write('MELT_preprocess_' + str(sample) + '=`qsub ' + pbs_path + str(sample) + '_preprocess_pbs.sh`' + '\n')
            fout.write('echo $MELT_preprocess_' + str(sample) + '\n')
            #individual analysis steps
            for element in TE:
                fout.write('MELT_IndAn_' + str(sample) + '_' + str(element) + '=`qsub -W depend=afterok:$MELT_preprocess_' + str(sample) + ' ' + pbs_path + str(sample) + '_ind_analysis_' + str(element) + '_pbs.sh`' + '\n')
                fout.write('echo $MELT_IndAn_' + str(sample) + '_' + str(element) + '\n')

        #group analysis step
        for element in TE:
            all_samples_dependency = ':'.join([('${MELT_IndAn_' + str(sample) + '_' + str(element) + '}') for sample in sample_list])
            fout.write('MELT_GroupAn_' + str(element) + '=`qsub -W depend=afterok:' + all_samples_dependency + ' ' + pbs_path + 'group_analysis_' + str(element) + '_pbs.sh`' '\n')
            fout.write('echo $MELT_GroupAn_' + str(element) + '\n')
            group_analysis_dependency = ':'.join([('${MELT_preprocess_' + str(element) + '}')])

        #genotype step
        for index, row in SampleSheet_df.iterrows():
            patient = row['patient']
            sample = row['sample']
            sample_path = row['bam_paths']
            #sample_red = sample_path.split('/')[-1].replace(".dedup.realigned.bam", "").replace("_Illumina_sorted_markedDupl.bam","").replace(".bam","").replace(".","_")
            for element in TE:
                fout.write('MELT_Genotype_' + str(sample) + '_' + str(element) + '=`qsub -W depend=afterok:' + '$MELT_GroupAn_' +  str(element) + ' ' + pbs_path + str(sample) + '_genotype_' + str(element) + '_pbs.sh`' '\n')
                fout.write('echo $MELT_Genotype_' + str(sample) + '_' + str(element) + '\n')

        #makeVCF step
        for element in TE:
            genotype_dependency = ':'.join([('${MELT_Genotype_' + str(sample) + '_' + str(element) + '}') for sample in sample_list])
            fout.write('MELT_MakeVCF_' + str(element) + '=`qsub -W depend=afterok:' + genotype_dependency + ' ' + pbs_path + 'MakeVCF_' + str(element) + '_pbs.sh`' '\n')
            fout.write('echo $MELT_MakeVCF_' + str(element) + '\n')
            make_vcf_dependency = ':'.join([('${MELT_MakeVCF_' + str(element) + '}')])
            fout.write('MELT_somatic_' + str(element) + '=`qsub -W depend=afterok:' + make_vcf_dependency +  ' ' + pbs_path + str(patient) + '_somatic_calls_' + str(element) + '_pbs.sh`' '\n')
            fout.write('echo $MELT_somatic_' + str(element) + '\n')

        fout.write('exit 0')

    # Make executable
    st = os.stat(call_all_script)
    os.chmod(call_all_script, st.st_mode | stat.S_IXUSR | stat.S_IXGRP)

def local_parser(parser, require_user=False):
    '''
    parse arguments for all steps of pipe
    '''
    parser.add_argument('-r','--run', help= 'Enter run', required=True)
    parser.add_argument('-l','--sample_list', help= 'Enter sample list', required=True)
    parser.add_argument('-u','--user',
    help='PBS will email the garvan email address of this user for execution start and end of pbs job',
    required=True)
    parser.add_argument('-c','--cov')
    parser.add_argument('-o','--output',
    help='provide the path of file containing mean coverage for each individual samples',
    required=True)
    parser.add_argument('-p','--picard', help='specify if the user needs to derive coverage for bam files first before running MELT', default=False)
    parser.add_argument('-s','--somatic', help='specify if user wants to filter putative somatic calls out of final VCF', default=False)
    args = parser.parse_args()
    return args.run, args.sample_list, args.user, args.cov, args.output, args.picard, args.somatic

def extract_sample_name(sample):
    sample_red = sample.replace(".bam", "")
    return sample_red

if __name__ == "__main__":
    main()
