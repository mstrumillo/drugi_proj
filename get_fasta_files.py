#!/usr/bin/env python

import os
import subprocess
from sys import argv
#creates folder with domain name
#with folders inside including full seq fastas and domain fastas
#creates file needed_fasta_domain





def create_needed_fasta_file(domain):
    os.system("grep %s */output|sort|uniq| awk '{if ($7>=1) print $0}' > needed_fasta_%s"%(domain,domain))


def run_alignment_mapR(domain):
    os.system("sed -i \'s/PF00012.18/%s/g\' alignment_map.R"%(domain))
    command='Rscript'
    path2script='alignment_map.R'
    cmd=[command,path2script]
    subprocess.check_output(cmd,universal_newlines=True)







def pull_out_domain(ens,file,start,end):
    fastas=open(file,"r").readlines()
    znacznik=0
    j=0
    sequence=[]
    header=[]
    while j<=len(fastas)-1:
        linijeczka=fastas[j].split()
        id=linijeczka[0].lstrip(">")  #################################I think THIS is where HUMAN fails cause HUMAN have .1 .2. 3. 4. 5.
        if "ENS" in id:                 ######this is added for human
            id=id.split(".")[0]         ###########
        if ens==id:                            #there cant be in here, cause some of the numbers might be .11 instead of .1
            header.append(fastas[j].strip())
            header.append(";%s;%s\n"%(start,end))
            znacznik=1
            j+=1
        if ">" in fastas[j]:
            znacznik=0
        if znacznik==1:
            sequence.append(fastas[j].strip())
        j+=1
    head="".join(header)
    seq="".join(sequence)
    domain_fasta=seq[int(start)-1:int(end)]
    return(head,seq,domain_fasta)



def cut_out_all_domains(domain): ###############I need phosphorylation info at this level, so I wont be doing the alignment 100000 times.
    needed_fasta=open("needed_fasta_%s"%domain,"r").readlines()############and I was so smart I already gathered the phosp info in needed_fasta in LAST COLUMN, coming from process_hmm_file, but thats protein info
    if not os.path.isdir("proteins_%s"%domain):
        os.mkdir("proteins_%s"%domain)
    dir="proteins_%s"%domain
    if not os.path.isdir("%s/full_proteins"%dir):
        os.mkdir("%s/full_proteins"%dir)
        os.mkdir("%s/domains"%dir)
        os.mkdir("%s/domains/pdb"%dir)
    full="%s/full_proteins"%dir
    dom="%s/domains"%dir
    if not os.path.isfile("%s/domain_seq_for_mafft_%s.fasta"%(dir,domain)):
        domains_seqs=open("%s/domain_seq_for_mafft_%s.fasta"%(dir,domain),"w")                                                       #open a file with domains, so you can run alignment on it
        for i in range(0,len(needed_fasta)):
            linijeczka=needed_fasta[i].split(",")
            protein=linijeczka[1].strip()
            start=linijeczka[2].strip()
            end=linijeczka[3].strip()
            domain=linijeczka[4].strip()
            indicator=int(linijeczka[-1].strip())
            print(indicator)
            l=needed_fasta[i].split("/")
            species=l[0].strip()
            if indicator>0:
                os.system("grep "+protein+" /nfs/research2/beltrao/strumill/2PROJECT/proteomes/*txt > "+dir+"/domains/"+domain+"_"+protein+".phosp")
                (header,seq,domain_fasta)=pull_out_domain(protein,species+"/ensembl",start,end)                     #gives you the header of the protein, full sequence, and the domain sequence
##################################### I only want to get the seq if I have the phosphorylation in it:
                phosp_file_path=dir+"/domains/"+domain+"_"+protein+".phosp"
                phosp_file=open(phosp_file_path,"r").readlines()
                is_phosp=0
                for n in range(0,len(phosp_file)):
                    lin2=phosp_file[n].split()
                    nr=int(lin2[-2])
                    if nr>=int(start) and nr<=int(end):
                        is_phosp+=1
                if is_phosp>0:
                    domains_seqs.write(header)                                                                          #add domain sequence to one big file that you run mafft on
                    domains_seqs.write(domain_fasta+"\n")                                                               #
                    
                    domain_fasta_file=open("%s/%s_%s.fasta"%(dom,domain,protein),"w")                                   #this makes tiny fasta file just with the domain (so you can map phosps on them later
                    domain_fasta_file.write(header)                                                                     #
                    domain_fasta_file.write(domain_fasta+"\n")                                                          #
                    
                    output=open("%s/full_protein_%s.fasta"%(full,protein),"w")                                           #this writes down whole fasta to a protein fasta file
                    output.write(header)                                                                                #
                    output.write(seq+"\n")                                                                              #
    
def run_mafft_fastest(domain):
    os.system("mafft domain_seq_for_mafft_%s.fasta > domain_seq_for_mafft_%s.ali"%(domain,domain))

def run_mafft_FFT_NS_i(domain):
    os.system("mafft --retree 2 -maxiterate 1000 domain_seq_for_mafft_%s.fasta > domain_seq_for_mafft_%s_FFT_NS_i.ali"%(domain,domain))

def run_mafft_NW_NS_2(domain):
    os.system("mafft --retree 2 --maxiterate 2 --nofft domain_seq_for_mafft_%s.fasta > domain_seq_for_mafft_%s_no_fft.ali"%(domain,domain))
###############################################
#################provide domain name###########
#######makes folder############################
##########runs alignment#######################
###############################################


#domain=argv[1]
#
#
#create_needed_fasta_file(domain)
#print("1")
#cut_out_all_domains(domain)
#print("2")
#os.system("cp domains_compare.py proteins_%s"%domain)
#print("3")
#os.system("cp alignment_map.R proteins_%s"%domain)
#print("4")
#os.system("cp how_to_make_a_chart proteins_%s"%domain)
#print("5")
#os.system("mv needed_fasta_%s proteins_%s"%(domain, domain))
#print("6")
#os.chdir("proteins_%s"%domain)
#print("7")
#
#print("8")
#run_mafft(domain)
#print("9")
#run_alignment_mapR(domain)
#print("10")
#os.system("python domains_compare.py")
#print("11")
#os.system("./how_to_make_a_chart")
#print("12")
#
#
#
#
#
