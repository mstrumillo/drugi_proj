#!/usr/bin/env python

import os
import subprocess
from sys import argv
import Bio
import glob
from Bio.PDB import *
#creates file needed_pdb_domain
#copies pdb into 2PROJECT/pdb database
#cuts out the domain and puts into dir inside proteins_domain dir



pdb_database="/nfs/ftp/pub/databases/pdb/data/structures/all/pdb"
used_pdb_database="/nfs/research2/beltrao/strumill/2PROJECT/pdb"


def create_needed_pdb_file(domain2):
    domain=domain2.split(".")[0]
    os.system("grep %s /nfs/research2/beltrao/strumill/2PROJECT/pdb_mapping/pdbmap |sort|uniq > needed_pdb_%s"%(domain,domain))



def copy_structures(domain2):
    domain=domain2.split(".")[0]
#    needed_struc=open("proteins_%s/needed_pdb_%s"%(domain2,domain),"r").readlines()
    needed_struc=open("needed_pdb_%s"%(domain),"r").readlines()
    for i in range(0,len(needed_struc)):
        line=needed_struc[i].split(";")
        pdb_id=(line[0].strip()).lower()
        identifier=pdb_id[1:3]
        if not os.path.isdir("%s/%s"%(used_pdb_database,identifier)):
            os.system("mkdir %s/%s"%(used_pdb_database,identifier))
        if not os.path.isfile("%s/%s/pdb%s.ent"%(used_pdb_database,identifier,pdb_id)):
            os.system("cp %s/pdb%s.ent.gz %s/%s"%(pdb_database,pdb_id, used_pdb_database,identifier))
            os.system("gunzip %s/%s/pdb%s.ent.gz"%(used_pdb_database,identifier,pdb_id))
    return()
###########
###########
###########I guess Id parse the file again and pull out the chains and copy them inside the protein_domain/domains_pdb folder
###########
###########

def pull_out_domain_structure(structure_id, chain, start, end, output):
    print(structure_id, chain, start, end, output)
    p = PDBParser(PERMISSIVE=1)
    structure=p.get_structure(structure_id,structure_id)
    print("Im extracting %s, %s, %s, %s"%(structure_id,chain, start, end))
    try:
        Dice.extract(structure,chain,start,end,output)
    except:
        print("I didnt extract  %s, %s, %s, %s"%(structure_id,chain, start, end))
    return()

#pull_out_domain_structure("1a06","A",20,276,"outputcik")

def pull_out_all_domains_from_file(domain2):
#    domain_dir="proteins_%s/domains/pdb"%domain2
    domain_dir="domains/pdb"
    domain=domain2.split(".")[0]
#    domains=open("needed_pdb_%s"%(domain),"r").readlines()
    domains=open("needed_pdb_%s"%(domain),"r").readlines()
    for i in range(0,len(domains)):
        used_pdb="/nfs/research2/beltrao/strumill/2PROJECT/pdb"
        line=domains[i].split(";")
        pdb_id=line[0].strip().lower()
        identifier=pdb_id[1:3]
###############if I actually have a structure for this protein
        if os.path.isfile("%s/%s/pdb%s.ent"%(used_pdb,identifier,pdb_id)):
            chain=line[1].strip()
            start_end=line[5].strip()
            start=int(start_end.split("-")[0])
            end=int(start_end.split("-")[1])
            structure_id="%s/%s/pdb%s.ent"%(used_pdb,identifier,pdb_id)
            output="%s/%s_%s_%s_%s_%s.pdb"%(domain_dir,domain,pdb_id,chain,start,end)
            if not os.path.isfile(output):
                pull_out_domain_structure(structure_id,chain,start,end,output)
        else:
            print("theres no pdb file for %s"%pdb_id)
    return()

def prepare_dir_for_naccess(domain):
    os.system("cp /nfs/research2/beltrao/strumill/2PROJECT/proteomes/backup/standard.data /nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain)
    os.system("cp /nfs/research2/beltrao/strumill/2PROJECT/proteomes/backup/vdw.radii     /nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain)
    os.system("cp /nfs/research2/beltrao/strumill/2PROJECT/proteomes/backup/dssp     /nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain)
    return()


def run_naccess(protein):
    os.system("/nfs/research2/beltrao/software/foreign/naccess-2.1.1/bin/naccess %s"%protein)
    return()

def run_dssp(protein):
    prot=protein[:-4]
    os.system("./dssp -i %s.pdb  -o %s.dssp"%(prot,prot))
    return()

def get_fasta_seq(protein):#######################sequence from dssp file if PDB allows - good pdb filter
    prot=protein[:-4]################################SEQUENCE NEEDS TO COME FROM BIOPYTHON - DSSP IS SHIT, AND DETECTS BREAKS AS !
#    PF00118_3j1e_A_42_532.pdb########################ACKTUALLY I can, just get rid of ! and make a stamp that there's a breakage.
    p=prot.split("_")
    pdb=p[1]
    chain=p[2]
    start=p[3]
    end=p[4]
    if os.path.isfile("%s.dssp"%prot):##############################THIS CANNOT BE LIKE THIS. LETS FIND HOW TO EXTRACT WITH BIPYTHON
        dssp=open("%s.dssp"%prot,"r").readlines()
        output=open("%s.fasta"%prot,"w")
        seq=""
        for i in range(25,len(dssp)):
            l=dssp[i].split()
            aa=l[3]
            if aa.isalpha():
                seq+=aa
        output.write(">%s %s; %s; %s\n"%(prot, chain, start, end))
        output.write(seq)
        output.write("\n")
        output.close()
    return()

def make_list_to_align(domain):
    domains_to_align="/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain
    os.chdir(domains_to_align)
    os.system("ls *fasta |sed \'s/.fasta/.pdb/g\'>list")###########################this has to be done using FASTA, cause pdb might not exist and then what.
    os.chdir("/nfs/research2/beltrao/strumill/2PROJECT/proteomes/")
    return()

def clustering(domain):
    max_cluster="/nfs/research2/beltrao/strumill/programs/maxcluster64bit"
    domains_to_align="/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s/domains/pdb"%domain
    os.chdir(domains_to_align)
    os.system("%s -l list -C 1 >output"%(max_cluster))
#    os.chdir("/nfs/research2/beltrao/strumill/2PROJECT/proteomes/")
    return()   




cwd = os.getcwd()
#cwd
#'/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_PF00118.22'
domenka=cwd.split("/")[-1].strip("\'")
domain= domenka.split("_")[1]
print(domain)
#########################THISSSSSSS@@@@@@@@@@##########################################################################################
#########################THISSSSSSS@@@@@@@@@@##########################################################################################

create_needed_pdb_file(domain)  #########output in proteomes folder
copy_structures(domain)         #########output in ../pdb folder
pull_out_all_domains_from_file(domain) #########output in proteomes/protein
##
##
prepare_dir_for_naccess(domain)
os.chdir("domains/pdb")
        
os.system("find . -name \"*.pdb\" -size -26k -delete") ############################this removes empty PDB files
###
proteins=glob.glob("*pdb")
###
###
for i in proteins:
    run_naccess(i)
    run_dssp(i)
    get_fasta_seq(i)
###
###
###
make_list_to_align(domain)
print("done list, now clustering")
clustering(domain)
############################################GET CENTROIDS STARTS##############################################


clusters=open("output","r").readlines()

ali="domain_seq_for_mafft_%s.fasta"%(domain)
os.system("cp ../../%s ."%ali)
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

print(lista)

for item in lista:
#PF00118.22_1a6d.phosp
#PF00118.22_PF00118_1q3q_B_35_526.fasta
    fas=item[1]+"fasta"
    phos=item[1]+"phosp"
    print(fas)
    fasta_file_to_copy=domain+"_"+fas
    phosp_file_to_make=domain+"_"+phos
    nowy="%s_%s"%(ali,fas)
    os.system("cp %s %s"%(ali, nowy))
    os.system("cat %s >> %s"%(fas, nowy))
    os.system("mafft %s > %s.ali"%(nowy,nowy))
    os.system("cp %s ../%s"%(fas,fasta_file_to_copy))
    os.system("touch ../%s"%phosp_file_to_make)


os.system("cp *ali ../../")


