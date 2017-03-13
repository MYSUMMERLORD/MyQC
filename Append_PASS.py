def Append_PASS(conf_dict,logfile):
    all_samples_MQS_WCQS_file = conf_dict['General']['outputdirectory']+'Samples_MQS_and_WeightedCombinedQualityScore.txt'

    infile = open(all_samples_MQS_WCQS_file,'r')
    header = infile.readline()
    header = header.rstrip()+'\t'+'QC'+'\n'
    outfile = open(conf_dict['General']['outputdirectory']+'MyQC_All_Samples_QC_information.txt','w')
    outfile.write(header)
    line = infile.readline()
    while(line):
        outfile.write(line.rstrip()+'\t'+'PASS'+'\n')
        line = infile.readline()
    infile.close()
    outfile.close()