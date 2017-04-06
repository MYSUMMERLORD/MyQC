import sys
import time
import os
import commands

def main():
	file_ref = open('E-MTAB-2600.sdrf.txt','r')
	all_names = {}
    file_ref.readline()
	for line in file_ref:
        id = line.strip().split("\t")[26].split('/')[-1].split(".fastq")[0]
        name = line.strip().split("\t")[24].split("Teichmann_")[1].split(".fastq")[0]
		all_names[id] = name

	#cmd = 'ls *.fastq'
	#a = commands.getstatusoutput(cmd)
	#Temp_Raw = a[1].split('\n')
    #for i in Temp_Raw:
    #    commands.getstatusoutput("mv i %s.fastq" %all_names[i.split(".fastq")[0]])
    cmd = 'ls *.gz'
    a = commands.getstatusoutput(cmd)
    Temp_Raw = a[1].split('\n')
    for i in Temp_Raw:
        commands.getstatusoutput("mv i %s.fastq.gz" %all_names[i.split(".fastq")[0]])

if __name__ == '__main__':
	try:
    		main()
    		print('this is the end')
	except KeyboardInterrupt:
    		sys.stderr.write("User interrupt MyQC\n")
    		sys.exit(1)
