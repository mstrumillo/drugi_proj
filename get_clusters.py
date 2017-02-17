#!/usr/bin/env python

import os
import subprocess
from sys import argv
import Bio
from Bio.PDB import *

#/nfs/research2/beltrao/strumill/programs/maxcluster64bit -l list -C 1 >output







def make_list_to_align(domain):
    domains_to_align="/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain
    os.chdir(domains_to_align)
    os.system("ls *pdb >list")
    os.chdir("/nfs/research2/beltrao/strumill/2PROJECT/proteomes/")
    return()






def clustering(domain):
    max_cluster="/nfs/research2/beltrao/strumill/programs/maxcluster64bit"
    domains_to_align="/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain
    os.chdir(domains_to_align)
    os.system("%s -l list -C 1 >output"%(max_cluster))
    os.chdir("/nfs/research2/beltrao/strumill/2PROJECT/proteomes/")
   






#domain=argv[1]
#make_list_to_align(domain)
#clustering(domain)
#
#
#
