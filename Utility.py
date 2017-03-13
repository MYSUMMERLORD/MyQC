import os
import sys
import ConfigParser
import subprocess


def CMD(cmd):
    os.system(cmd)

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    :param cmd:
    :return:
    '''
    a = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac

def raise_error():
    '''
    raise an error message and exit
    :return:
    '''
    print('error occurs, check log file~!')
    sys.exit(1)

def createDIR(dirname):
    '''
    check dir name and create new dir
    :param dirname:
    :return:
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s' %(dirname))

def wlog(message, logfile):
    '''
    print a message and write the message to logfile
    :param message:
    :param logfile:
    :return:
    '''
    print(message)
    os.system('echo "%s " >> %s'%(message, logfile))

def ewlog(message, logfile):
    '''
    print an error message and write the error message to logfile
    error messages start with [ERROR]
    :param message:
    :param logfile:
    :return:
    '''
    print('ERROR %s' % (message))
    os.system('echo "[ERROR] %s " >> %s' %(message, logfile))
    raise_error()

def rwlog(cmd, logfile):
    '''
    prints an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    :param cmd:
    :param logfile:
    :return:
    '''
    print('[CMD] %s' % (cmd))
    os.system('echo "[CMD] %s" >> %s' %(cmd, logfile))
    CMD(cmd)

def rplog(cmd, logfile):
    '''
    use os.popen to get shell
    :param cmd:
    :param logfile:
    :return:
    '''
    print('[CMD] %s' % (cmd))
    os.system('echo "[CMD] %s" >> %s' % (cmd, logfile))
    os.popen(cmd).read()

# read config
def read_conf(conf_file):
    '''
    Read config file and return a dict containing all information
    :param conf_file:
    :return:
    '''
    conf_dict = {}
    cf = ConfigParser.SafeConfigParser()
    cf.read(conf_file)
    conf_dict['General']={}
    for item in cf.items('General'):
        conf_dict['General'][item[0]] = item[1]
    return conf_dict

def calculate_total_reads(samplefile):
    inf = open(samplefile)
    totalN = 0
    last_seq = ""
    mappableN = 0
    unique_reads = {}
    for line in inf:
        if line.startswith("@"):
            continue
        ll = line.strip().split("\t")
        seq_name = ll[0]
        chrom = ll[2]
        seq = ll[9]
        unique_reads[seq] = []
        if not seq_name == last_seq:
            totalN += 1
            if not chrom == '*':
                mappableN += 1
        last_seq = seq_name
    return totalN,mappableN,(int(mappableN)*1.0/totalN),len(unique_reads)




