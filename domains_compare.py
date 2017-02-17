#!/usr/bin/env python
import glob
import os
import subprocess


#creates phosp files in domains folder, so all phosps per protein are known
#runs through alignments and create *four files, that have phosps mapped on them


def get_fasta_letters(fasta_file):
    fasta_file=open(fasta_file+".fasta","r").readlines()
    header=fasta_file[0]
    fasta=fasta_file[1]
    start=fasta_file[0].split(";")[-2]
    end=fasta_file[0].split(";")[-1]
    return(header,start,end,fasta)  #returns header and the sequence


#def create_phosp_files(all_proteins):
#    for i in all_proteins:
##        print(i)
#        domain=i.split("_")[0].lstrip("domains/")
##        print(domain)
##        ens=i.split("_")[1][:-6]#####if theres multiple _ in the file splitting didnt work
#        l_d=int(len(domain))+9#############niektore domeny maja 382742.12 a nie ktore 129381.1 
#        ens=i[l_d:-6]    ##################wiec bylo triky ale dziala
#        print(domain,ens)
#        os.system("grep "+ens+" ../*txt > domains/"+domain+"_"+ens+".phosp")


def make_4_columns(file):
    ali_file=open(file,"r").readlines()
    domain=file.split("_")[0]
    print(file)
###########I dont need the protein name, luckily
    domain_protein=file[:-4]
    print(domain_protein)
    if "ENS" in domain_protein: #########################This is only cause of that ridiculous human 
        domain_protein='.'.join(domain_protein.split(".")[:-1])
        print(domain_protein)
    phosp_file=open("domains/%s.phosp"%domain_protein,"r").readlines()

    (info,start,end,seq)=get_fasta_letters("domains/%s"%domain_protein) 
    output=open(domain_protein+".four","w")
    output.write(info) #start files with fasta info
    j=1
    k=1
    #print(len(seq))
    while j<len(ali_file):
    # print(ali_file[j])
        linijeczka=ali_file[j].split()
        ali_position=linijeczka[1].strip("\"")
        domain_position=linijeczka[2].strip("\"")
        while k<len(seq):
#            print(ali_file[j])
            linijeczka=ali_file[j].split()
            ali_position=linijeczka[1].strip("\"")
            domain_position=linijeczka[2].strip("\"")
            if str(domain_position)==str(k):
                protein_position=int(domain_position)-1+int(start)
                znacznik=0
                for l in range(0,len(phosp_file)):
                    linijka=phosp_file[l].split()
                    pp=linijka[-2].strip(",")
                    aa2=linijka[-1].strip(",")
#                    print(pp,aa2)
                    if int(protein_position)==int(pp) and str(seq[k-1])==str(aa2):
                        znacznik=1
                if znacznik==1:
                    output.write("%s, %s, %s, 1\n"%(ali_position, protein_position, seq[k-1])) #make file with three columns, add the letters to the ali position and protein position
                if znacznik==0:
                    output.write("%s, %s, %s, 0\n"%(ali_position, protein_position, seq[k-1]))
                k+=1
                j+=1
            elif str(domain_position)=="-":
                output.write("%s, %s, %s, 0\n"%(ali_position, domain_position, "-")) #if theres a gap, write -
                j+=1
        if str(domain_position)=="-":
            output.write("%s, %s, %s, 0\n"%(ali_position, domain_position, "-")) #finish file with gaps, after writing all protein
            j+=1



#all_proteins=glob.glob("domains/*fasta")
#create_phosp_files(all_proteins)



alignments=glob.glob("PF*ali")
for i in alignments:
    print(i)
    make_4_columns(i)

#make_4_columns("PF00069.23_ABN64346.ali")

#make_4_columns("PF00069.23_ZC581.1.ali")
