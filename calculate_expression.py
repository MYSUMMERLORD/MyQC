import os
import commands
import subprocess
### tool function
from Utility import (createDIR,
                     wlog,
                     ewlog,
                     rwlog,
                     rplog,
                     sp)


# ------------------------------------
# main
# ------------------------------------

def calculate_expression(conf_dict, logfile):
    '''
    generate expression matrix file
    :param conf_dict:
    :param logfile:
    :return:
    '''

    cmd = 'ls' + '\t' + conf_dict['General']['sam'] + "*.sam"
    a = commands.getstatusoutput(cmd)
    Temp_Raw = a[1].split('\n')
    Temp_Raw.sort()
    sample_name_list = [i.split('/')[-1].split('.sam')[0] for i in Temp_Raw]
    sample_name_list.sort()
    ### create mapping dir
    if conf_dict['General']['expression'].rstrip() == "":
        conf_dict['General']['expression'] = conf_dict['General']['outputdirectory'] + 'expression/'
        createDIR(conf_dict['General']['expression'])
        os.chdir(conf_dict['General']['expression'])
        if conf_dict['General']['quantative_software_main'] == "rsem":
            wlog('user choose rsem as quantative software', logfile)
            if sp('which rsem-calculate-expression')[0].strip() == "":
                ewlog('rsem is not detected in default PATH, make sure you installed rsem and export it into default PATH',
                      logfile)
            else:
                reference_dir = conf_dict['General']['expression'] + 'ref/'
                createDIR(reference_dir)
                ## check reference's exist
                if not os.path.exists(conf_dict['General']['reference']):
                    print('rsem reference doesnt exist')
                    return
                cmd = "rsem-prepare-reference %s %srsem_ref" % (conf_dict['General']['reference'], reference_dir)
                rwlog(cmd, logfile)
                processes = set()
                if conf_dict['General']['type'] == "pe":
                    for i in range(len(Temp_Raw)):
                        cmd = 'rsem-calculate-expression --paired-end --alignments %s %srsem_ref %s%s' %(Temp_Raw[i],reference_dir,conf_dict['General']['expression'],sample_name_list[i])
                        # rplog(cmd,logfile)
                        wlog(cmd, logfile)
                        processes.add(subprocess.Popen(cmd, shell=True))
                else:
                    for i in range(len(Temp_Raw)):
                        cmd = 'rsem-calculate-expression --alignments %s %srsem_ref %s%s' %(Temp_Raw[i],reference_dir,conf_dict['General']['expression'],sample_name_list[i])
                        #rplog(cmd,logfile)
                        wlog(cmd, logfile)
                        processes.add(subprocess.Popen(cmd,shell=True))
                count = 0
                for p in processes:
                    count += 1
                    if p.poll() is None:
                        p.wait()
                        print('**************** %s ******************' % count)

        else:
            ewlog('quantative tools can only be rsem', logfile)
    wlog('Start compiling TPM matrix of all samples ...', logfile)
    #print("it's me ")
     # collect gene result
    if conf_dict['General']['quantative_software_main'] == "rsem":
        w_all_TPM = {}
        infile = open('%s%s.genes.results'%(conf_dict['General']['expression'],sample_name_list[0]),'r')
        heeader_temp = infile.readline()
        line = infile.readline()
        while(line):
            # merge different isoform
            if conf_dict['General']['expression_type'].rstrip() == 'Gene':
                Gene_inf = line.split('\t')[0].split('|')[5]
            else:
                Gene_inf = line.split('\t')[0]
            #Gene_inf = line.split('\t')[0].split('|')[5]
            # consider each isoform as a gene
            #Gene_inf = line.split('\t')[0]+'\t'+line.split('\t')[1]
            #Gene_inf = line.split('\t')[0]
            w_all_TPM[Gene_inf] = []
            line = infile.readline()
        infile.close()
        for i in range(len(sample_name_list)):
            infile = open('%s%s.genes.results' % (conf_dict['General']['expression'],sample_name_list[i]), 'r')
            header_temp = infile.readline()
            line = infile.readline()
            while(line):
                if conf_dict['General']['expression_type'].rstrip() == 'Gene':
                    Gene_inf = line.split('\t')[0].split('|')[5]
                else:
                    Gene_inf = line.split('\t')[0]
                #Gene_inf = line.split('\t')[0] + '\t' + line.split('\t')[1]

                TPM =  line.split('\t')[-2]
                if len(w_all_TPM[Gene_inf])==i+1:
                    new_value = float(w_all_TPM[Gene_inf][i]) + float(TPM)
                    w_all_TPM[Gene_inf][i] = str(new_value)
                else:
                    w_all_TPM[Gene_inf].append(TPM)
                line = infile.readline()
            infile.close()
        outfile = open(conf_dict['General']['outputdirectory']+'expression_matrix.txt','w')

        header = 'GeneName'+'\t'+'\t'.join(sample_name_list)+'\n'
        outfile.write(header)
        for k in w_all_TPM.keys():
            outfile.write(k+'\t'+'\t'.join(w_all_TPM[k])+'\n')
        outfile.close()

        # get gene detected in each sample
        sample_detected_genes = []
        for i in range(len(sample_name_list)):
            detected_gene_count = 0
            for k in w_all_TPM.keys():
                if float(w_all_TPM[k][i]) > float(conf_dict['General']['tpm_cutoff']):
                    detected_gene_count += 1
            sample_detected_genes.append(detected_gene_count)
        return sample_detected_genes
    else:
        pass



