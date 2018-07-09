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
    parser = argparse.ArgumentParser()
    run, textfile_path, user = tau_config.tau_argparse(parser, require_email=True)
    email = user + '@garvan.org.au'
    print email
    #list of bam files to analyse
    #textfile_path = '/g/data1a/jp48/MELT/bam_files.txt'
    if os.path.exists(textfile_path):
        samples = names_from_txt_file(textfile_path)
    pbs_path = '/g/data1a/jp48/MELT/' + run
    if not os.path.isdir(pbs_path):
        os.makedirs(pbs_path)
        #might alos need to set permission so everyone can acces files

    generate_all_separate_pbs_scripts(samples, pbs_path, email)

def names_from_txt_file(textfile_path):
    sample_no = 0
    sample_name_df = pd.read_table(textfile_path, sep='\t', header=None)
    print sample_name_df
    samples = sample_name_df.iloc[:, 0].tolist()
    print samples
    return samples

def generate_all_separate_pbs_scripts(samples, pbs_path, email):
    RAM = '4'
    CPUs = '8'
    for sample in samples:
        #Preprocess step
        generate_template_per_sample(RAM,
                            CPUs,
                            'preprocess',
                            'preprocess.sh',
                            sample,
                            pbs_path,
                            email,
                            '2',
                            'express')

        # generate_template_per_sample(RAM,
        #                     CPUs,
        #                     'print_reads',
        #                     'print_reads.py',
        #                     sample,
        #                     pbs_path,
        #                     email)

    create_pbs_launch_script(samples, pbs_path)

def generate_template_per_sample(mem,
                                cpus,
                                what,
                                script_name,
                                sample,
                                pbs_path,
                                email,
                                time,
                                batchqueue):
    '''
    writes the variable details of the pbs scripts to call all sample specific pbs scripts
    '''
    sample_red = sample.replace(".dedup.realigned.bam", "")
    name = '_'.join([str(sample_red), what])
    #pbs_output_dir = '/g/data1a/jp48/MELT/pbs_logs/'
    pbs_output_dir = pbs_path
    scripts_dir='/g/data1a/jp48/scripts/transposon/'
    out = os.path.join(pbs_output_dir, name)
    full_script_name = os.path.join(scripts_dir, script_name)
    #give bam file here, not sample name
    #cmd = "bash '" + full_script_name + "' -s '" + str(sample_red) + "'"
    cmd = "bash '" + full_script_name + "' '/g/data3/ba08/softlinks/" + str(sample) + "'"
    wd = "echo Working directory is ${PBS_O_WORKDIR}\n" + \
         "cd ${PBS_O_WORKDIR}\n"
    #batchqueue = 'normal' #normal or express?
    script_name = pbs_path + name + '_pbs.sh'
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
                  #do we need time too?
    template = Template(file=template_file, searchList=[name_space])
    with open(script_name, "w") as f:
        f.write(str(template))

def create_pbs_launch_script(samples, pbs_path):
    '''
    Manages dependencies in the PBS queue - what waits what
    '''
    call_all_script = pbs_path + 'call_all_pbs.sh'
    with open(call_all_script, 'w') as fout:
        fout.write('#!/bin/bash' + '\n')
        for sample in samples:
            sample_red = sample.replace(".dedup.realigned.bam", "")
            fout.write('MELT_preprocess_' + str(sample_red) + '=`qsub ' + pbs_path + str(sample_red) + '_preprocess_pbs.sh`' + '\n')
            fout.write('echo $MELT_preprocess_' + str(sample_red) + '\n')
            sample_recal_dependency = ':'.join([('${MELT_preprocess_' + str(sample) + '}')])
            # fout.write('GATK_PRINTREADS_' + str(sample_red) + '=`qsub -W depend=afterok:$GATK_RECAL_' + str(sample_red) + ' ' + pbs_path +  str(sample_red)+ '_print_reads_pbs.sh`' + '\n')
            # fout.write('echo $GATK_PRINTREADS_' + str(sample_red) + '\n')
        fout.write('exit 0')
    # Make executable.
    st = os.stat(call_all_script)
    os.chmod(call_all_script, st.st_mode | stat.S_IXUSR | stat.S_IXGRP)

def argparse(parser, require_email=False):
    '''
    parse arguments for all steps of pipe/QC
    '''
    parser.add_argument('-r','--run', help= 'Enter run', required=True)
    parser.add_argument('-list','--sample-list', help= 'Enter sample list')
    parser.add_argument('-u','--user',
    help='PBS will email this address for execution start and end of pbs job',
    required=require_email)
    args = parser.parse_args()
    return args.run, args.list, args.user


if __name__ == "__main__":
    main()
         
