#coding=utf-8
'''
Created on 2016年12月19日

@author: summer
'''

import sys
import time
import os
import commands

### tool function
from Utility import (calculate_total_reads,
                     wlog,
                     CMD,
                     read_conf)

# ------------------------------
# main step
# ------------------------------
from mapping import mapping
from calculate_expression import calculate_expression
from Plot_FPR import Plot_FPR
from Estimate_Cutoff_and_Annotate_Artifact import Estimate_Cutoff_and_Annotate_Artifact
from Append_PASS import Append_PASS



# -------------------------------
# Main function
# -------------------------------

def main():
    script_path = os.path.abspath(os.path.dirname(__file__))
    conf_dict = read_conf('%s/template.conf'%script_path)
    #conf_dict['General']['startdir'] = os.getcwd() + '/'
    #conf_dict['General']['startdir'] = '/homeb/liulijuan/single_cell/QC/data/pythonTest/'
    # output directory alse can be designed by user
    #conf_dict['General']['outputdirectory'] = conf_dict['General']['startdir'] + 'output/'
    # samples dir
    #conf_dict['General']['samples_file'] = conf_dict['General']['startdir'] + "samples/"
    # input data format
    #conf_dict['General']['format'] = "sam"
    Parameter_list = sys.argv[1:]
    for i in range(0,len(Parameter_list)-1,2):
        if(Parameter_list[i].split('-')[1] in conf_dict['General']):
            conf_dict['General'][Parameter_list[i].split('-')[1]] = Parameter_list[i+1]
        else:
            pass

    # create output dir
    if(os.path.isfile(conf_dict['General']['outputdirectory'].rstrip("/"))):
        print('ERROR: name of your output dir is exist as a file, cannot create a dir, MyQC exit')
        sys.exit(1)
    elif os.path.isdir(conf_dict['General']['outputdirectory']):
        print('name of your output dir is exist as a dir, overwrite is turned on, write output result in existing dir')
    else:
        CMD("mkdir %s " %(conf_dict['General']['outputdirectory'].rstrip("/")))

    ### move to output dir
    os.chdir(conf_dict['General']['outputdirectory'])
    ### specify the main progress log file
    logfile = conf_dict['General']['outputdirectory'] + 'progress_log.txt'
    ### remove existing log file
    if(os.path.isfile(logfile)):
        # linux delete command is rm
        CMD('rm %s' %logfile)
        # in windows delete command is del
        #CMD('del %s' %logfile)

    ### main step for MyQC
    wlog("Start MyQC", logfile)
    t = time.time()
    start_time = t
    ############## Step1：mapping #####################
    wlog("Step1：mapping",logfile)
    #samples_info = mapping(conf_dict,logfile)
    mapping(conf_dict, logfile)
    s1time = time.time() - t
    wlog("time for alignment : %s"%(s1time),logfile)
    wlog("Step1：alignment DONE",logfile)

    ############## Step 2：calculate expression value##################
    t = time.time()
    wlog("Step2：calculate samples expression value",logfile)
    if conf_dict['General']['expression_matrix'].rstrip() == "":
        sample_detected_genes = calculate_expression(conf_dict,logfile)   
        wlog("time for calculate expression value : %s" % (s2time), logfile)
    wlog("Step2.1：calculate expression DONE", logfile)
    #write sample sequence information
    outfile = open('%sseq_information.txt'% conf_dict['General']['outputdirectory'],'w')
    header = 'Sample.Name'+'\t'+'Total.Reads'+'\t'+'Mapped.Reads'+'\t'+'Mapped.Rate'+'\t'+'Reads.Complexity'+'\t'+'Gene.Detected'+'\n'
    outfile.write(header)

    cmd = 'ls' + '\t' + conf_dict['General']['expression'] + "*.genes.results"
    a = commands.getstatusoutput(cmd)
    Temp_Raw = a[1].split('\n')
    sample_name_list = [i.split('/')[-1].split('.genes.results')[0] for i in Temp_Raw]
    for i in range(len(sample_name_list)):
        totalN, mappedN, mapping_rate, reads_complexity = calculate_total_reads('%s%s.sam' % (conf_dict['General']['sam'],sample_name_list[i]))
        outfile.write(sample_name_list[i]+'\t')
        detected_genes = sample_detected_genes[i]
        outfile.write(str(totalN)+'\t'+str(mappedN)+'\t'+str(mapping_rate)+'\t'+str(reads_complexity)+'\t'+str(detected_genes)+'\n')
        print('The Current sample is %s'%sample_name_list[i])
    outfile.close()
    s2time = time.time() - t
    wlog("Step2.2：calculate sequence information DONE", logfile)
    #################### Step3：calculate distinct p_value ########################
    t = time.time()
    wlog('Step3：calculate distinct p_value',logfile)
    if conf_dict['General']['pvalue_file'].rstrip() == "":
        expression_matirx_file = conf_dict['General']['outputdirectory']+'expression_matrix.txt'

        cmd = 'Rscript'+'\t'+script_path+'/MI.R'+'\t'+expression_matirx_file+'\t'+conf_dict['General']['outputdirectory']
        CMD(cmd)
        s3time = time.time() - t
        wlog("time for calculate distinct p_value : %s" % (s3time), logfile)
    wlog("Step3：calculate calculate distinct p_value DONE", logfile)

    ###################Step4: Calculate SeqQC quantile in all samples###############
    t = time.time()
    wlog('Step4：calculate SeqQC quantile in all samples',logfile)
    print(conf_dict['General']['pvalue_file']+'seq_information.txt')
    cmd = 'Rscript'+'\t'+script_path+'/Calculate_Quantile_in_all_samples_update.R'+'\t'+conf_dict['General']['pvalue_file']+'Distinct_PValue.txt'+'\t'+conf_dict['General']['pvalue_mi_cutoff']+'\t'+conf_dict['General']['pvalue_pearson_cutoff']+'\t'+conf_dict['General']['pvalue_spearman_cutoff']+'\t'+conf_dict['General']['logical_tag']+'\t'+conf_dict['General']['pvalue_file']+'seq_information.txt'+'\t'+conf_dict['General']['outputdirectory']
    CMD(cmd)
    s4time = time.time() - t
    wlog('time for calculate SeqQC quantile in all samples : %s' % (s4time),logfile)
    wlog('Step4：calculate SeqQC quantile in all samples DONE',logfile)

    ####################Step5: calculate Weighted Combined Quality Score and MQS#############
    t = time.time()
    wlog('Step5：calculate Weighted Combined Quality Score and Minimal Quantile Score',logfile)
    cmd = 'Rscript'+'\t'+script_path+'/Weight_Combined_Score_and_MQS_update.R'+'\t'+conf_dict['General']['outputdirectory']+'Samples_Seq_Quality_in_all_Quantile.txt'+'\t'+conf_dict['General']['outputdirectory']
    CMD(cmd)
    s5time = time.time() - t
    wlog('time for calculate Weighted Combined Quality Score and Minimal Quantile Score : %s'  % (s5time),logfile)
    wlog('Step5：calculate Weighted Combined Quality Score and Minimal Quantile Score DONE',logfile)

    ####################Step6: FPR#########################################################
    t = time.time()
    wlog('Step6：FPR',logfile)
    # Check How Many Gene Expression Outliors and MainPopulationCell
    wlog('Start check how many expression outliers and mainpopulationcell',logfile)
    infile = open(conf_dict['General']['outputdirectory']+'Samples_Seq_Quality_in_all_Quantile.txt','r')
    header = infile.readline()
    MainPopulationCell_Count = int(0)
    GeneExpressionOutlier_Count = int(0)
    line = infile.readline()
    while(line):
        type = line.split()[1]
        if(type == 'GeneExpressionOutlier'):
            GeneExpressionOutlier_Count += 1
        elif(type == 'MainPopulationCell'):
            MainPopulationCell_Count += 1
        else:
            print('Error ......')
        line = infile.readline()
    infile.close()
    print('GeneExpressionOutlier_Count = %s, MainPopulationCell_Count = %s'%(GeneExpressionOutlier_Count,MainPopulationCell_Count))
    # Start calculating False Positive Table ...
    wlog('Start calculating false positive table',logfile)
    if(MainPopulationCell_Count > 0 and GeneExpressionOutlier_Count > 0):
        cmd = 'Rscript' + '\t' + script_path + '/False_Positive_Rate_Table_Raw.R' + '\t' + conf_dict['General'][
            'outputdirectory'] + 'Samples_MQS_and_WeightedCombinedQualityScore.txt' + '\t' + conf_dict['General'][
                  'outputdirectory']
        CMD(cmd)
        Plot_FPR(script_path,conf_dict,logfile)
    else:
        pass

    s6time = time.time()
    wlog('time for FPR : %s' %s6time,logfile)
    wlog('Step6：FPR DONE',logfile)

    ###################Step7: Annotate artifact#####################
    t = time.time()
    wlog('Step7：Annotate artifact',logfile)
    if(MainPopulationCell_Count>0 and GeneExpressionOutlier_Count>0):
        Estimate_Cutoff_and_Annotate_Artifact(conf_dict,logfile)
    else:
        Append_PASS(conf_dict,logfile)
    s7time = time.time() - t
    wlog('time for annotate artifact : %s' %s7time,logfile)
    wlog('Step7：Annotate artifact DONE',logfile)

    ############# Ploting PCA###################################
    #cmd = 'Rscript'+'\t'+script_path+'/PCA.R'+'\t'+conf_dict['General']['outputdirectory'] + 'expression_matrix.txt'+'\t'+conf_dict['General']['outputdirectory']+'MyQC_All_Samples_QC_information.txt'+'\t'+conf_dict['General']['outputdirectory']
    CMD(cmd)

    ################ Summary ##########################
    outfile_summary = open(conf_dict['General']['outputdirectory']+'Summary.txt','w')
    outfile_summary.write('Summary of MyQC Run'+'\n')
    outfile_summary.write('\n')
    if(MainPopulationCell_Count > 0 and GeneExpressionOutlier_Count > 0):
        infile = open(conf_dict['General']['outputdirectory']+'Real_cutoff.txt','r')
        header = infile.readline()
        line = infile.readline()
        MQS_Cutoff = line.split()[0].split('(')[1].split(')')[0].split(',')[0]
        WCQS_Cutoff = line.split()[0].split('(')[1].split(')')[0].split(',')[1]
        N_Artifacts = line.split()[1]
        Fraction_Artifacts = line.split()[2]
        Real_FPR = line.split()[-1]
        infile.close()

        total_cells = int(MainPopulationCell_Count)+int(GeneExpressionOutlier_Count)
        w_artifact = {}
        infile = open(conf_dict['General']['outputdirectory']+'MyQC_All_Samples_QC_information.txt','r')
        header = infile.readline()
        line = infile.readline()
        while(line):
            sample_name = line.split()[0]
            QC = line.split()[-1]
            if(QC == 'Artifact'):
                w_artifact[sample_name]=''
            else:
                pass
            line = infile.readline()
        infile.close()

        outfile_summary.write('Total Samples: '+'\t'+str(total_cells)+'\n')
        outfile_summary.write('Total GeneExpressionOutliers: '+'\t'+str(GeneExpressionOutlier_Count)+'\n')
        outfile_summary.write('Maximal False Positive Rate (FPR) Allowed: '+'\t'+str(conf_dict['General']['max_fpr'])+'\n')
        outfile_summary.write('MQS_Cutoff[Estimated]='+str(MQS_Cutoff)+'\t'+'WCQS_Cutoff[Estimated]='+str(WCQS_Cutoff)+'\n')
        outfile_summary.write('Real False Positive Rate (FPR) in the Dataset:'+'\t'+Real_FPR+'\n')
        outfile_summary.write("Total Artifacts:"+'\t'+str(len(w_artifact.keys()))+'\n')
        outfile_summary.write('Artifact Samples:'+'\t'+','.join(w_artifact.keys())+'\n')
        outfile_summary.write('\n')
    else:
        total_cells = int(MainPopulationCell_Count) + int(GeneExpressionOutlier_Count)
        outfile_summary.write('Total Samples: ' + '\t' + str(total_cells) + '\n')
        outfile_summary.write('Total GeneExpressionOutliers: ' + '\t' + str(GeneExpressionOutlier_Count) + '\n')
        outfile_summary.write('Total MainPopulationCell:'+'\t'+str(MainPopulationCell_Count)+'\n')
        if(GeneExpressionOutlier_Count == 0):
            outfile_summary.write('Skip downstream analysis due to No GeneExpressionOutliers found!'+'\n')
        else:
            outfile_summary.write('Skip downstream analysis due to No MainPopulationCells found!'+'\n')


    running_time = time.time() - start_time
    outfile_summary.write('Total Running Time: %.2d:%.2d:%.2d' % (running_time/3600,(running_time%3600)/60,running_time%60)+'\n')
    outfile_summary.close()


if __name__ == '__main__':
    try:
        main()
        print('this is the end')
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt MyQC\n")
        sys.exit(1)