#!/bin/env python
from Cheetah.Template import Template
import argparse
import os
import stat
import sys
import fileinput
import collections
import pandas as pd


def main():
    #EXAMPLE OF COMMAND
    #python pbs_make_all.py -r 'cscc' -u 'm.gauthier' -l '/g/data1a/jp48/MELT/bam_files.txt'
    #add coverage
    #python pbs_make_all.py -r 'cscc' -u 'm.gauthier' -l '/g/data1a/jp48/MELT/bam_files.txt' -c '/g/data1a/jp48/MELT/Coverage.txt
    #Got error 'function' object has no attribute 'ArgumentParser' as function for arsing called like module

    #TO DO:
    #specify full path of bam files, at the moment hardcoded on line 127

    parser = argparse.ArgumentParser(description='Create the parser')
    run, textfile_path, user, cov_path = local_parser(parser, require_user=True)
    email = user + '@garvan.org.au'
    print email
    print cov_path
    #list of bam files to analyse
    #textfile_path = '/g/data1a/jp48/MELT/bam_files.txt'
    if os.path.exists(textfile_path):
        samples = names_from_txt_file(textfile_path)
    pbs_path = '/g/data1a/jp48/scripts/transposon/' + run + '/'
    if not os.path.isdir(pbs_path):
        os.makedirs(pbs_path)
        #might alos need to set permission so everyone can acces files


    generate_all_separate_pbs_scripts(samples, pbs_path, email, cov_path)

def names_from_txt_file(textfile_path):
    sample_no = 0
    sample_name_df = pd.read_table(textfile_path, sep='\t', header=None)
    print sample_name_df
    samples = sample_name_df.iloc[:, 0].tolist()
    print samples
    return samples

def generate_all_separate_pbs_scripts(samples, pbs_path, email, cov_path):
    RAM = '4'
    CPUs = '1' #try 4 #8 took 1h50 to do the preprocessing step
    for sample in samples:
        #Preprocess step
        generate_template_per_sample('4',
                            '1',
                            'preprocess',
                            'preprocess.sh',
                            sample,
                            pbs_path,
                            email,
                            '3',
                            'express',
                            cov_path)

        #Individual analysis step

        generate_template_per_sample('6', #might have to bring this up to 8Gb
                            '1', #confirmation from developer that MELT will only use 1 core
                            'ind_analysis',
                            'IndivAnalysis.sh',
                            sample,
                            pbs_path,
                            email,
                            '24', #optimise, will this be enough for ALU??
                            'normal',
                            cov_path)

    create_pbs_launch_script(samples, pbs_path)

def generate_template_per_sample(mem,
                                cpus,
                                step,
                                script_name,
                                sample,
                                pbs_path,
                                email,
                                time,
                                batchqueue,
                                cov_path):
    '''
    writes the variable details of the pbs scripts to call all sample specific pbs scripts
    '''
    sample_red = sample.replace(".dedup.realigned.bam", "")
    name = '_'.join([str(sample_red), step])
    #pbs_output_dir = '/g/data1a/jp48/MELT/pbs_logs/'
    pbs_output_dir = pbs_path
    scripts_dir='/g/data1a/jp48/scripts/transposon/'
    out = os.path.join(pbs_output_dir, name)
    full_script_name = os.path.join(scripts_dir, script_name)
    #give bam file here, not sample name

    wd = "echo Working directory is ${PBS_O_WORKDIR}\n" + \
         "cd ${PBS_O_WORKDIR}\n"

    #command for all steps except Individual analysis:
    if step == 'ind_analysis':
        #derive coverage specific to sample
        print 'Looking for sample matching ' + str(sample_red) + ' in coverage file\n'
        cov_for_sample = 0
        with open(cov_path, 'r') as fin:
            for line in fin:
                bits = line.strip().split('\t')
                sample_with_cov = bits[0]
                cov = bits[1]
                if sample_with_cov == sample_red:
                    print 'Found a match ' + str(sample_with_cov) + ' ' + str(cov)
                    cov_for_sample = cov
                #write a script for each transcription element
                #this is what is different for this step
                TE = ('LINE1', 'ALU', 'SVA')
                for element in TE:
                    #at this stage just perform this step for one TE element, but eventually loop through the 3 TEs
                    #cmd += " '" + str(cov_for_sample) + "' 'ALU'"
                    cmd = "bash '" + full_script_name + "' '/g/data3/ba08/softlinks/" + str(sample) + "' '" + str(cov_for_sample) + "' '" + str(element) + "'"
                    print 'Creating script for ' + str(element)
                    script_name = pbs_path + name + '_' + str(element) + '_pbs.sh'
                    write_template(name, email, batchqueue, out, cpus, mem, time, cmd, wd, script_name)
    else:
        cmd = "bash '" + full_script_name + "' '/g/data3/ba08/softlinks/" + str(sample) + "'"
        script_name = pbs_path + name + '_pbs.sh'
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

def create_pbs_launch_script(samples, pbs_path):
    '''
    Manages dependencies in the PBS queue - what waits for what
    '''
    call_all_script = pbs_path + 'call_all_pbs.sh'
    with open(call_all_script, 'w') as fout:
        fout.write('#!/bin/bash' + '\n')
        for sample in samples:
            sample_red = sample.replace(".dedup.realigned.bam", "")
            fout.write('MELT_preprocess_' + str(sample_red) + '=`qsub ' + pbs_path + str(sample_red) + '_preprocess_pbs.sh`' + '\n')
            fout.write('echo $MELT_preprocess_' + str(sample_red) + '\n')
            sample_recal_dependency = ':'.join([('${MELT_preprocess_' + str(sample) + '}')])
            TE = ('LINE1', 'ALU', 'SVA')
            for element in TE:
                fout.write('MELT_IndAn_' + str(sample_red) + '=`qsub -W depend=afterok:$MELT_preprocess_' + str(sample_red) + ' ' + pbs_path +  str(sample_red) + '_ind_analysis_' + str(element) + '_pbs.sh`' + '\n')
                fout.write('echo $MELT_IndAn_' + str(sample_red) + '_' + str(element) + '\n')
        fout.write('exit 0')
    # Make executable.
    st = os.stat(call_all_script)
    os.chmod(call_all_script, st.st_mode | stat.S_IXUSR | stat.S_IXGRP)

def local_parser(parser, require_user=False):
    '''
    parse arguments for all steps of pipe/QC
    '''
    parser.add_argument('-r','--run', help= 'Enter run', required=True)
    parser.add_argument('-l','--sample_list', help= 'Enter sample list')
    parser.add_argument('-u','--user',
    help='PBS will email the garvan email address of this user for execution start and end of pbs job',
    required=require_user)
    parser.add_argument('-c','--cov',
    help='provide the path of file containing mean coverage for each individual samples')

    args = parser.parse_args()
    return args.run, args.sample_list, args.user, args.cov


if __name__ == "__main__":
    main()
         
