#!/usr/bin/env python
import numpy as np
import glob
import os
import subprocess
import statistics
from sys import argv




#domain = "PF00069.23"
def get_domain():
    cwd = os.getcwd()
#/nfs/research2/beltrao/strumill/2PROJECT/proteomes/done_proteins/proteins_PF00069.23
    p=cwd.split("_")
    domain=p[-1]
    return(domain)

domain=get_domain()


def get_main_pdb():
    main_cluster=open("domains/pdb/main_cluster","r").readlines()
    name=main_cluster[2].split("_")
    pdb=name[1].strip()
    chain=name[2].strip()
    pdb_start=name[3].strip()
    krotka=(pdb,chain,pdb_start)
    return(krotka)

(pdb,chain,pdb_start)=get_main_pdb()

print(pdb_start)
############maybe run the find the interesting positionsos.system() ##############Or maybe add in rerun

interesting_positions=open("interesting_positions","r").readlines()
in_pos=[]

pdb_name=pdb+"_"+chain+"_"+pdb_start

for i in range(0,len(interesting_positions)):
    if pdb_name in interesting_positions[i]:
        line=interesting_positions[i].split(",")
        pos=line[-1].strip()
        pval=line[-2].strip()
        krotka=(pos,pval)
        in_pos.append(krotka)


print(in_pos)
#[('6', '6.99137072188e-09'), ('7', '3.30947313033e-05'), ('8', '3.63030816559e-09'), ('9', '5.47841997387e-09'), ('10', '0.000259521158416'), ('24', '0.00161598744909'), ('140', '0.00116014461342'), ('141', '0.00156999470887'), ('142', '0.00251804434577'), ('143', '0.00183974878876'), ('144', '2.70416288331e-06'), ('146', '0.00252910378904'), ('147', '8.74518270502e-22'), ('148', '7.29431847303e-163'), ('149', '0.0'), ('150', '0.0'), ('151', '0.0'), ('152', '0.0'), ('153', '0.0'), ('154', '6.50198673715e-122'), ('155', '5.13549392991e-95'), ('156', '5.72002336465e-87')]
##############THIS ARE INTERESTING POSITIONS IN THE CHART#AND I NEED TO PASS ON THE PVALS


ref="done_%s_%s_%s/reference.four"%(pdb,chain,pdb_start)
os.system("sed 1d %s|grep -v \"-\" >%s_2"%(ref,ref))


interesting_alignment_line=[]
ref2=open("done_%s_%s_%s/reference.four_2"%(pdb,chain,pdb_start),"r").readlines()

for ali_pos in in_pos:
    line=int(ali_pos[0])
    pval=float(ali_pos[1])
    l_ref2=ref2[line].split(",")
#    print(l_ref2)i
    pos=int(l_ref2[0].strip())
    krotka=(pos,pval,line)
    interesting_alignment_line.append(krotka)
print(interesting_alignment_line)
print(len(in_pos),len(interesting_alignment_line))
#['16,', '24,', '25,', '26,', '27,', '96,', '759,', '763,', '775,', '782,', '783,', '831,', '856,', '857,', '858,', '859,', '922,', '924,', '952,', '953,', '954,', '955,']
#6#['16,', '50,', 'G,', '0']#7#['24,', '51,', 'T,', '0']#8#['25,', '52,', 'G,', '0']#9#['26,', '53,', 'S,', '0']#10#['27,', '54,', 'F,', '0']#24#['96,', '68,', 'H,', '0']#140#['759,', '184,', 'D,', '0']
######BUT 16,24,25,26,27,96 ARE THE INTERESTING POSITIONS IN THE ALIGNMENT###########
    


protein_list=[]
o_from_cere=open("../Saccharomyces_cerevisiae/output","r").readlines()
for i in range(0,len(o_from_cere)):
    if domain in o_from_cere[i]:
        line=o_from_cere[i].split(",")
        protein=line[1].strip()
        start=line[2].strip()
        end=line[3].strip()
        how_many_phosp=int(line[-1])
        if how_many_phosp>=1:#########################################################UNLESS I WANNA TRACK SOMETHING THATS NOT PHOSPHORYLATED YET
            krotka=(protein,start,end)
            protein_list.append(krotka)

print(protein_list)
#[(' YDL126C', ' 251', ' 380'), (' YDL126C', ' 524', ' 657'), (' YPL074W', ' 507', ' 647'), (' YOR259C', ' 218', ' 351'), (' YER017C', ' 325', ' 458')]



##############now for every protein in the protein_list, I need to find what is the protein that has the same position in the alignment 


output=open("yeast_protein_to_mutate","w")




for protein_krotka in protein_list:##################THIS WILL HAVE TO INCLUDE MORE IDS FOR WHEN I FIX THE NAMING
    protein=protein_krotka[0]
    start=protein_krotka[1]
    end=protein_krotka[2]
#I need to get start and end.    
#PF00102.25_ENSP00000334941.5_1203_1438.ali PF00102.25_ENSP00000396842_2_134.four PF00102.25_FBpp0072601_61_295.four reference.four_2
#FileNotFoundError: [Errno 2] No such file or directory: 'done_4zrt_A_40/PF00102.25_YOR208W.four'
    print("done_%s_%s_%s/%s_%s_%s_%s.four"%(pdb,chain,pdb_start,domain,protein,start,end))
    file_protein=open("done_%s_%s_%s/%s_%s_%s_%s.four"%(pdb,chain,pdb_start,domain,protein,start,end),"r").readlines()
    for j in interesting_alignment_line:#########################IM USING INTERESTIN_ALIGNMENT
        ali_pos=int(j[0])
        pval=float(j[1])
        chart_pos=int(j[2])
        line=file_protein[ali_pos].split(",")
        protein_ali_no=int(line[0].strip())############THIS IS FROM start-end so should be fine, theyre all the sam elength
        protein_no=line[1].strip()#####################THIS CAN BE "-"
        aa=line[2]#####################################THIS CAN BE "-"
        phosp=int(line[3])##################################THIS CAN BE 0/1
        if phosp==1:
#            print("%s, %s, chart_pos: %s, %.12f, %s, protein_no: %s, %s, %s"%(domain, protein, chart_pos, pval,ali_pos, protein_no, aa, phosp))
            output.write("%s, %s, chart_pos: %s,\t%.12f, ali_pos: %s,\t protein_no: %s,\t %s, %s\n"%(domain, protein, chart_pos, pval,ali_pos, protein_no, aa, phosp))
        
    
