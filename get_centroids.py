#!/usr/bin/env python
import numpy as np
import glob
import os
import subprocess
from sys import argv




clusters=open("output","r").readlines()
ali="domain_seq_for_mafft_PF00118.22.fasta"
dom_fasta=ali.split("_")[-1]
dom=dom_fasta.split(".")[:-1]
domain=".".join(dom)
how_many_clusters=int(clusters[len(clusters)-2].split()[2])
#print(how_many_clusters)


lista=[]

for i in range(len(clusters)-how_many_clusters-1,len(clusters)-1):
    l=clusters[i].split()
    size=l[5]
    spread=l[6]
    pdb_long=l[7][:-3]
#    pdb=pdb_long.split("_")[1:3]
    krotka=(size,pdb_long)
#    print(size,pdb_long)
    lista.append(krotka)


def check_if_mafft_done():
    os.system("bjobs|grep  \"mafft\" >jobs_mafft")
    i=1
    while i==1:
        if os.stat("jobs_mafft").st_size != 0:
            os.system("bjobs|grep  \"mafft\" >jobs_mafft")
            print("still counting")
            time.sleep(5)
        else:
            i=0
            return(True)



print(lista)

for item in lista:
#PF00118.22_PF00118_1q3q_B_35_526.fasta
    fas=item[1]+"fasta"
    phos=item[1]+"phosp"
    print(fas)
    fasta_file_to_copy=domain+"_"+fas
    phosp_file_to_make=domain+"_"+phos
    nowy="%s_%s"%(ali,fas)
    os.system("cp %s %s"%(ali, nowy))
    os.system("cat %s >> %s"%(fas, nowy))
    os.system("bsub -e e/err%s -o o/out%s mafft %s >%s.ali"%(item, item, nowy, nowy))
#    os.system("mafft %s > %s.ali"%(nowy,nowy))
    os.system("cp %s ../%s"%(fas,fasta_file_to_copy))
    os.system("touch ../%s"%phosp_file_to_make)
    if check_if_mafft_done():
        os.system("cp *ali ../../")

