import sys
import time
import os
import commands

def main():
	file_ref = open('exp_name.txt','r')
	all_names = {}
	for line in file_ref:
		all_names[line.strip().split("\t")[0].split(".fastq")[0].split("Teichmann_")[1]] = line.strip().split("\t")[1]

	cmd = 'ls *.fastq'
    	a = commands.getstatusoutput(cmd)
    	Temp_Raw = a[1].split('\n')
   	sample_name_list = [i.split('/')[-1].split('.fastq')[0] for i in Temp_Raw]
	for i in sample_name_list:
        	del all_names[i]
    	cmd = 'ls *.gz'
        a = commands.getstatusoutput(cmd)
        Temp_Raw = a[1].split('\n')
    	sample_name_list = [i.split('/')[-1].split('.fastq.gz')[0] for i in Temp_Raw]
    	for i in sample_name_list:
        	del all_names[i]

	print("the number of missing files:%s"%len(all_names))
    	outfile = open('missing_name.txt','w')
    	for i in all_names:
    		outfile.write(all_names[i]+"\n")


if __name__ == '__main__':
	try:
    		main()
    		print('this is the end')
	except KeyboardInterrupt:
    		sys.stderr.write("User interrupt MyQC\n")
    		sys.exit(1)
