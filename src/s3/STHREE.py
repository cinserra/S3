#!/usr/bin/env python

import os,shutil,subprocess,sys
from numpy import *
import pyfits
from optparse import OptionParser
from pyraf import iraf
import s3
import string,os,re,glob
usage= "%prog  -v"

s3_dir = s3.__path__[0]

if __name__ == "__main__":
    parser = OptionParser(usage=usage,version="%prog "+str(s3.__version__))
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    option,args = parser.parse_args()
   
def aiuto(verbose):   ###################                ####################
    import re,string,os
    folder=re.sub('SNAKE\n','',os.popen('which SNAKE').readlines()[0])
    rows=glob.glob(folder+'S*')
    print "#"*30+"   S3   "+"#"*30
    for i in rows:
            prog=string.split(i,'/')[-1]
            print prog,'-h  (show the help message)'
    print "#"*10+" c.inserra@qub.ac.uk "+"#"*10+" github: cinserra "+"#"*10+" Twitter: @COSMO_83"+"#"*10

if __name__ == "__main__":
    aiuto(option.verbose)