#!/usr/bin/env python


import os
import glob
import subprocess
import time
from sys import argv
from map_charts import *
#creates folder with domain name



def get_fasta_letters(fasta_file):
    fasta_file=open(fasta_file+".fasta","r").readlines()
    header=fasta_file[0]
    fasta=fasta_file[1]
    start=fasta_file[0].split(";")[-2]
    end=fasta_file[0].split(";")[-1]
    return(header,start,end,fasta)  #returns header and the sequence


def run_alignment_mapR(alignment,domain):
    os.system("sed -i \'s/COS/%s/g\' alignment_map.R"%(alignment))
    os.system("sed -i \'s/PF00012.18_/%s_/g\' alignment_map.R"%(domain))
    command='/nfs/research2/beltrao/software-rh7/bin/Rscript'
    path2script='alignment_map.R'
    cmd=[command,path2script]
    subprocess.check_output(cmd,universal_newlines=True)
    os.system("sed -i \'s/%s/COS/g\' alignment_map.R"%(alignment))
    os.system("sed -i \'s/%s_/PF00012.18_/g\' alignment_map.R"%(domain))



def make_4_columns(file):
    ali_file=open(file,"r").readlines()
    domain=file.split("_")[0]
    #######I dont need the protein name, luckily
    domain_protein=file[:-4]#########################DOMAIN_PROTEIN_START_END
    rzeczy=domain_protein.split("_")
    if rzeczy[-4]==rzeczy[-2] and rzeczy[-3]==rzeczy[-1]:###########THEN WE'RE DEALING WITH PDB 
###########PF00004.27_PF00004_3j3u_F_204_339.fasta  PF00004.27_PF00004_3j98_B_539_669.phosp ############THIS IS THE 3rd CASE
###########FileNotFoundError: [Errno 2] No such file or directory: 'domains/PF00004.27_PF00004_5dyi_B_241_371_241_371.fasta'
        domain_protein='_'.join(rzeczy[:-2])
#        print("rzeczy wygladaja tak:"+domain_protein)
        phosp_file=open("domains/%s.phosp"%domain_protein,"r").readlines()
        (info,start,end,seq)=get_fasta_letters("domains/%s"%domain_protein)
    elif "ENS" in domain_protein: #########################This is only cause of that ridiculous human########### 
        ens='.'.join(domain_protein.split(".")[:-1])
#        print(ens)
        phosp_file=open("domains/%s.phosp"%(ens),"r").readlines()
        start=domain_protein.split("_")[-2]
        end=domain_protein.split("_")[-1]
        domain_protein=ens+"_"+start+"_"+end
#        print(domain_protein)
        (info,start,end,seq)=get_fasta_letters("domains/%s"%domain_protein)
#FileNotFoundError: [Errno 2] No such file or directory: 'domains/PF00004.27_ENSP00000157812.1_202_335.phosp'
#PF00004.27_ENSP00000157812_202_335.fasta  PF00004.27_ENSP00000157812.phosp          
    else:
#        print("domains/%s"%domain_protein)
        ens='_'.join(domain_protein.split("_")[:-2])
#        print("to ens: "+ens)
#        print("to domain_protein: "+domain_protein)
        phosp_file=open("domains/%s.phosp"%ens,"r").readlines()
        (info,start,end,seq)=get_fasta_letters("domains/%s"%domain_protein)
#PF00004.27_BRADI1G32310.1_248_379.fasta PF00004.27_BRADI1G32310.1.phosp
    print(domain_protein)
    output=open(domain_protein+".four","w")
    output.write(info) #start files with fasta info
    j=1
    k=1
#    print(len(seq))
    while j<len(ali_file):
#        print(seq)
    # print(ali_file[j])
        linijeczka=ali_file[j].split()
        ali_position=linijeczka[1].strip("\"")
        domain_position=linijeczka[2].strip("\"")
        while k<len(seq):
#            print(ali_file[j])
            linijeczka=ali_file[j].split()
#            print(linijeczka)
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


def check_if_jobs_done():######################3this should be checking number of *four files, comparing to number of *ali files
    os.system("bjobs|grep -v \"/bin/bash\"|grep -v JOBID >jobs_going")
    i=1
    while i==1:
        if os.stat("jobs_going").st_size != 0:
            os.system("bjobs|grep -v \"/bin/bash\"|grep -v JOBID >jobs_going")
            print("still counting")
            time.sleep(5)
        else:
            i=0
            return(True)

def check_if_jobs_done2():######################3this should be checking number of *four files, comparing to number of *ali files
    time.sleep(100)
    alis=glob.glob("PF*ali")
    fours=glob.glob("PF*four")
    if len(alis)==len(fours):
        return(True)
    else:
        print("still counting")
        time.sleep(5)

def clean_up_folder(pdb_chain):
    cleaner="done_%s"%pdb_chain
    if not os.path.isdir("%s"%cleaner):
        os.system("mkdir %s"%cleaner)
    os.system("mv PF*ali %s"%cleaner)
    os.system("mv PF*four %s"%cleaner)
    os.system("mv reference.four %s"%cleaner)
    os.system("mv %s.tar.gz %s"%(pdb_chain,cleaner))
#   'domain_seq_for_mafft*2lck_A_14_111.fasta.ali*.ali'
#    domain_seq_for_mafft_PF00153.25.fasta_PF00153_4c9j_A_220_305.fasta.ali
    os.system("mv domain_seq_for_mafft*%s*fasta.ali %s"%(pdb_chain,cleaner))
    os.system("mv chart                   %s"%cleaner)                          
    os.system("mv permuts                 %s"%cleaner)                         
    os.system("mv e                       %s"%cleaner)                         
    os.system("mv o                       %s"%cleaner)                         
    os.system("mv medians                 %s"%cleaner)                          
    os.system("mv stdev                   %s"%cleaner)                          
    os.system("mv means                   %s"%cleaner)                          
    os.system("mv chart2                  %s"%cleaner)                          
    os.system("mv mapped_chart            %s"%cleaner)                          
    os.system("mv mapped_background       %s"%cleaner)                          
    os.system("mv permutated_background   %s"%cleaner)                          
    os.system("mv chart_window            %s"%cleaner)                          
    return()                                                         

def prepare_a_package(ali):#ali="domain_seq_for_mafft_PF00118.22.fasta_PF00118_1q3q_B_35_526.fasta.ali"
    a=ali.split("_")
    pdb=a[6]
    chain=a[7]
    domain=a[4][:-6]
    start=a[8]
    end=a[9]
    name=pdb+"_"+chain+"_"+start
    b=ali.split("fasta")
    pdb_name=b[1]
    print(name,pdb_name)
#PF00118.22_PF00118_1q3q_B_35_526.four
#PF00118.22_1q3q_B_35_526.fasta.ali.four
    reference_fasta=domain+pdb_name+"four"
###############creates small ali files
    run_alignment_mapR(ali,domain)
#####################creates small four files
    alignments=glob.glob("PF*ali")###################this is why cleaning is SOOOOOOOOOOOO IMPORTANT
    for i in alignments:
#        print(i)
        make_4_columns(i) ###########################this is why cleaning is SOOOOOOOOOOOO IMPORTANT
#####################creates folders chart permuts e o and counts all the permuts
    files=glob.glob("PF*.four")
    if not os.path.isdir("e"):#######################################it has to be here cause of bsub 
        os.system("mkdir e")
        os.system("mkdir o")
    for i in files:
#        print(i)
        os.system("bsub -e e/err%s -o o/out%s ./permutations.py %s"%(i,i,i))

###########waits for the bsub to finish counting permuts
    if check_if_jobs_done():
        print("now done")
    print(reference_fasta)
############################map_charts##########################################################
    create_big_chart_file() #########from map_charts, creates chart2 and permutated background
    map_charts(reference_fasta)######from map_charts maps chart2 and permuatetd bg, into mapped chart and mapped_background
    window_chart()###################from map_charts runs window over mapped_chart and returns chart_window
    write_down_medians()#############from map_charts runs window on permutated backgrounds and returns medians, means and stdev
############
    os.system("tar -zcvf %s.tar.gz medians means stdev chart_window"%name)
    return(name)




def run_for_one_ali(ali):
    n=prepare_a_package(ali)
    clean_up_folder(n)
    return()



ali=glob.glob("domain*ali")
for alignment in ali:
    run_for_one_ali(alignment)


#make_4_columns("PF00118.22_PF00118_1gn1_H_215_376.ali")
##run_for_one_ali("domain_seq_for_mafft_PF00118.22.fasta_PF00118_1a6d_A_34_519.fasta.ali")
##run_for_one_ali("domain_seq_for_mafft_PF00118.22.fasta_PF00118_1a6e_B_33_521.fasta.ali")
