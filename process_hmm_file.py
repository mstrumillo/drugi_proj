#!/usr/bin/env python
import os.path




def process_output(ensembl,phosp,output):
    hmm_file=open(ensembl,"r").readlines()
    phosp_file=open(phosp,"r").readlines()
    output=open(output,"w")
    for i in range(28,len(hmm_file)):
        linijeczka=hmm_file[i].split()
        protein_id=linijeczka[0].strip().split(".")[0] ###########this works only for humans, for the rest it should be without it
        start=int(linijeczka[1].strip())
        end=int(linijeczka[2].strip())
        domain_id=linijeczka[5].strip()
        domain_name=linijeczka[6].strip()
        indicator=0
        for j in range(0,len(phosp_file)):
            linijka=phosp_file[j].split()
            species=linijka[0].strip()
            protein_id2=linijka[-3].strip()
            number=int(linijka[-2].strip())
            if protein_id==protein_id2:
                if number>=start and number<=end:
                    indicator+=1
        output.write("%s, %s, %s, %s, %s, %s, %s\n"%(species,protein_id,start,end,domain_id,domain_name,indicator))






species=open("species_h","r").readlines()

for i in range(0,len(species)):
    spec=species[i].strip()
    if os.path.isfile("%s/phosp"%spec)  and os.path.isfile("%s/ensembl.hmm"%spec) and os.path.getsize("%s/phosp"%spec):
        print("tu jest %s"%spec)
        process_output("%s/ensembl.hmm"%spec,"%s/phosp"%spec,"%s/output"%spec)
