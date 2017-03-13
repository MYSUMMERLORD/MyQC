import os
import commands

### tool function
from Utility import (calculate_total_reads,
                     createDIR,
                     wlog,
                     ewlog,
                     rplog,
                     sp)

# ------------------------------------
# main
# ------------------------------------

def mapping(conf_dict,logfile):
    '''
    main data processing step, including mapping, generate mapping rate
    for fastq or fa format:
        bowtie2 mapping
    for sam format
        do nothing
    :param conf_dict:
    :param logfile:
    :return:
    '''

    ### create mapping dir
    mapping_dir = conf_dict['General']['outputdirectory'] + 'mapping/'
    createDIR(mapping_dir)
    ### check reads file format, start mapping step id format is fastq
    if conf_dict['General']['format'] == "sam":
         wlog('reads file format is sam, skip mapping step',logfile)
         cmd = 'ls' + '\t' + conf_dict['General']['samples_file'] + "*.sam"
         a = commands.getstatusoutput(cmd)
         Temp_Raw = a[1].split('\n')
         Check_Raw_List = [i.split('/')[-1].split('.sam')[0] for i in Temp_Raw]
         conf_dict['General']['sam'] = conf_dict['General']['samples_file']
    else:
        wlog('Now start mapping is %s , all mapping result will be here ' %(mapping_dir),logfile)
        os.chdir(mapping_dir)
        # check input file format
        if conf_dict['General']['format'] == "fastq":
            cmd = 'ls' + '\t' + conf_dict['General']['samples_file'] + "*.fastq"
        elif conf_dict['General']['format'] == "fq":
             cmd = 'ls' + '\t' + conf_dict['General']['samples_file'] + "*.fq"
        elif conf_dict['General']['format'] == "fasta":
             cmd = 'ls' + '\t' + conf_dict['General']['samples_file'] + "*.fasta"
        elif conf_dict['General']['format'] == "fa":
             cmd = 'ls' + '\t' + conf_dict['General']['samples_file'] + "*.fa"
        else:
             ewlog('Input file format Error: should be either FASTQ or FASTA or SAM',logfile)
        a = commands.getstatusoutput(cmd)
        Temp_Raw = a[1].split('\n')
        Check_Raw_List = [i.split('/')[-1].split('.fastq')[0] for i in Temp_Raw]
        Check_Raw_List.sort()
        flag = 0
        sample_name_list = []
        if conf_dict['General']['type'] == "pe":
            for k in Check_Raw_List:
                sample_name = k.split('_')[0]
                if flag == 0:
                    sample_name_list.append(sample_name)
                    flag = 1
                else:
                    flag = 0
            Check_Raw_List = sample_name_list
        ## choose mapping tool from  bowtie2 according to config file
        if conf_dict['General']['mapping_software_main'] == "bowtie":
            wlog('user choose bowtie as alignment software', logfile)
            if sp('which bowtie')[0].strip() == "":
                ewlog('bowtie is not detected in default PATH, make sure you installed bowtie and export it into default PATH',logfile)
            wlog('build the reference', logfile)

            if conf_dict['General']['format'] == "fa" or conf_dict['General']['format'] == "fasta":
                if conf_dict['General']['type'] == "pe":
                    for k in Check_Raw_List:
                        mapping_cmd = 'bowtie -a -m 200 -v 2 %s -f -1 %s%s_1.%s -2 %s%s_2.%s -S %s.sam' % (
                        conf_dict['General']['mapindex'], conf_dict['General']['samples_file'], k,
                        conf_dict['General']['format'], conf_dict['General']['samples_file'], k,
                        conf_dict['General']['format'], k)
                        rplog(mapping_cmd, logfile)
                else:
                    for k in Check_Raw_List:
                        mapping_cmd = 'bowtie -a -m 200 -v 2 %s -f %s%s.%s -S %s.sam' % (
                        conf_dict['General']['mapindex'], conf_dict['General']['samples_file'], k,
                        conf_dict['General']['format'], k)
                        rplog(mapping_cmd, logfile)
            else:
                if conf_dict['General']['type'] == "pe":
                    for k in Check_Raw_List:
                        mapping_cmd = 'bowtie -a -m 200 -v 2 %s -1 %s%s_1.%s -2 %s%s_2.%s -S %s.sam' %(conf_dict['General']['mapindex'],conf_dict['General']['samples_file'],k,conf_dict['General']['format'],conf_dict['General']['samples_file'],k,conf_dict['General']['format'],k)
                        rplog(mapping_cmd,logfile)
                else:
                    for k in Check_Raw_List:
                        mapping_cmd = 'bowtie -a -m 200 -v 2 %s %s%s.%s -S %s.sam' %(conf_dict['General']['mapindex'],conf_dict['General']['samples_file'],k,conf_dict['General']['format'],k)
                        rplog(mapping_cmd,logfile)
            conf_dict['General']['sam'] = mapping_dir
        else:
            ewlog('alignment tools can only be bowtie',logfile)

    # calculate mapping rate
    '''
    Check_Raw_List.sort()
    os.chdir(conf_dict['General']['sam'])
    samples_info = {}
    samples_info['samples_totalN'] = []
    samples_info['samples_mappableN'] = []
    samples_info['samples_mapping_rate'] = []
    samples_info['samples_unique_reads'] = []
    for k in Check_Raw_List:
        totalN,mappableN,mapping_rate, unique_reads = calculate_total_reads('%s.sam'%(k))
        samples_info['samples_totalN'].append(totalN)
        samples_info['samples_mappableN'].append(mappableN)
        samples_info['samples_mapping_rate'].append(mapping_rate)
        samples_info['samples_unique_reads'].append(unique_reads)
    return samples_info
    '''




