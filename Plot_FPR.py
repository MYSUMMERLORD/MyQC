from Utility import (CMD)

def Plot_FPR(script_path,conf_dict,logfile):
    w_FPR = {}
    infile = open(conf_dict['General']['outputdirectory']+'False_Positive_Rate_Table_With_all_combinations.txt','r')
    header = infile.readline()
    line = infile.readline()
    while(line):
        MQS_cufoff = line.split()[0]
        WCQS_cutoff = line.split()[1]
        Artifacts = line.split()[2]
        Fraction_Artifacts = line.split()[3]
        False_Positives = line.split()[4]
        False_Positive_Rate = line.split()[5]
        if(float(False_Positive_Rate) in w_FPR):
            w_FPR[float(False_Positive_Rate)].append([float(Fraction_Artifacts),'('+MQS_cufoff+','+WCQS_cutoff+')'+'\t'+Artifacts+'\t'+Fraction_Artifacts+'\t'+False_Positives+'\t'+False_Positive_Rate])
        else:
            w_FPR[float(False_Positive_Rate)]=[]
            w_FPR[float(False_Positive_Rate)].append([float(Fraction_Artifacts),'(' + MQS_cufoff + ',' + WCQS_cutoff + ')' + '\t' + Artifacts + '\t' + Fraction_Artifacts + '\t' + False_Positives + '\t' + False_Positive_Rate])
        line = infile.readline()
    infile.close()

    FPR_list = [float(i) for i in w_FPR.keys()]
    FPR_list.sort()

    outfile = open(conf_dict['General']['outputdirectory']+'Format_False_Positive_Table.txt','w')
    header = 'Cutoff.Pair.OR[MQS][WCQS]'+'\t'+'Artifacts'+'\t'+'Fraction_Artifacts'+'\t'+'False_Positives'+'\t'+'False_Positive_Rate'+'\n'
    outfile.write(header)
    for k in FPR_list:
        this_FPR_all = w_FPR[k]
        this_FPR_all.sort()
        this_FPR_all.reverse()
        outfile.write(this_FPR_all[0][1]+'\n')
    outfile.close()

    cmd = 'Rscript' + '\t' + script_path + '/Plot_FPR.R' + '\t' + conf_dict['General'][
        'outputdirectory'] + 'Format_False_Positive_Table.txt' + '\t' + conf_dict['General'][
              'outputdirectory']
    CMD(cmd)


