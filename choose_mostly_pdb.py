#!/usr/bin/env python

import os
import subprocess
from sys import argv
#creates folder with domain name
#with folders inside including full seq fastas and domain fastas
#creates file needed_fasta_domain

domains=open("most_popular_domains","r").readlines()
output=open("phosphorylated_domains_with_PDB","w")




for i in range(0,len(domains)):
    l=domains[i].split(",")
    id=l[0].split(".")[0]
    znacznik=0
    with open("../pdb_mapping/pdbmap") as f:
        for line in f:
            if id in line:
                znacznik=1
    if znacznik==1:
        output.write(domains[i])


