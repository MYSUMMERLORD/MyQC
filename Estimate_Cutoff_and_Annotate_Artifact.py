import commands
import sys

from Utility import ( wlog)


def Estimate_Cutoff_and_Annotate_Artifact(conf_dict,logfile):
    format_FPR_file = conf_dict['General']['outputdirectory']+'Format_False_Positive_Table.txt'
    user_Max_FPR_cutoff = conf_dict['General']['max_fpr']
    all_samples_MQS_WCQS_file = conf_dict['General']['outputdirectory']+'Samples_MQS_and_WeightedCombinedQualityScore.txt'

    meet_criteria_list = []
    infile = open(format_FPR_file,'r')
    header = infile.readline()
    line = infile.readline()
    while(line):
        FPR = float(line.split()[-1])
        if(FPR <= float(user_Max_FPR_cutoff)):
            fraction_artifacts = float(line.split()[2])
            meet_criteria_list.append([fraction_artifacts,line])
        else:
            pass
        line = infile.readline()
    infile.close()

    if(len(meet_criteria_list)>0):
        meet_criteria_list.sort()
        meet_criteria_list.reverse()
        outfile = open(conf_dict['General']['outputdirectory']+'Real_cutoff.txt','w')
        outfile.write(header)
        outfile.write(meet_criteria_list[0][1])
        outfile.close()
        MQS_Cutoff = float(meet_criteria_list[0][1].split()[0].split('(')[1].split(')')[0].split(',')[0])
        WCQS_Cutoff = float(meet_criteria_list[0][1].split()[0].split('(')[1].split(')')[0].split(',')[1])
        infile = open(all_samples_MQS_WCQS_file,'r')
        header = infile.readline()
        header = header.rstrip()+'\t'+'QC'+'\n'
        outfile = open(conf_dict['General']['outputdirectory']+'MyQC_All_Samples_QC_information.txt','w')
        outfile.write(header)
        line = infile.readline()
        while(line):
            outfile.write(line.rstrip()+'\t')
            type = line.split()[1]
            QC = ''
            if(type=='MainPopulationCell'):
                QC = 'PASS'
            elif(type=='GeneExpressionOutlier'):
                MQS = float(line.split()[-2])
                WCQS = float(line.split()[-1])
                if(MQS <= MQS_Cutoff or WCQS <= WCQS_Cutoff):
                    QC = 'Artifact'
                else:
                    QC = 'PASS'
            else:
                wlog('Error',logfile)
            outfile.write(QC+'\n')
            line = infile.readline()
        infile.close()
        outfile.close()
    else:
        outfile = open(conf_dict['General']['outputdirectory']+'Real_cutoff.txt','w')
        outfile.write('None of MQS_WCQS_meet_criteria'+'\n')
        outfile.close()
