#!/usr/bin/env python

import glob

files=glob.glob("ENSP00000*.fasta")

domainID="PF00069.23"


def get_fasta_letters(fasta_file):
    fasta_file=open(fasta_file+".fasta","r").readlines()
    fasta=[]
    for i in range(1,len(fasta_file)):
        for j in range(0,len(fasta_file[i].strip())):
            fasta.append(fasta_file[i][j])
    return(fasta_file[0],fasta)  #returns header and the sequence
        
def start_end_domain(file,domain_ID):
    hmm_file=open(file+".hmm",'r').readlines()
    for i in range(28,len(hmm_file)):
        linijeczka=hmm_file[i].split()
        domID=linijeczka[5]
        if domID==domainID:
            start=int(linijeczka[1])
            end=int(linijeczka[2])
            header="%s, %s, %s, %s"%(file, domainID, start, end)
            return(start,end,header)



def cut_out_domain(file,domain_ID):
    (start,end,head)=start_end_domain(file,domainID)
    (first_line,fasta)=get_fasta_letters(file)
    domain_fasta=fasta[start-1:end]
    output=open(domain_ID+"_"+file,"w")
    output.write(">%s\n"%(head))
    output.write("%s\n"%''.join(domain_fasta))
    return(domain_fasta)




for i in files:
    file=i.strip(".fasta")
    print(i)
    cut_out_domain(file,domainID)





#print(head)
#print(''.join(domain_fasta)) #############theres a difference in fasta[end] and fasta[1:end]










