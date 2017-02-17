#!/usr/bin/env python
import numpy as np
import glob
import os
import subprocess
from decimal import *
from sys import argv

########PERMUTATIONS###################################################################################################################################
#for each sequence count how many pS,pT,pY and then sample the phosphorylation across the STYs that are not pS,pT,pY
#######################################################################################################################################################
########the initial sequence needs to be tuples########################################################################################################
########permute the original sequence for 1 different letter###########################################################################################
def permute_one_letter(letter,phosp_list_S,sequence): ###########################################################################return permutated list
    if len(phosp_list_S)>1:
        x=0
        permutated_list=[]
        while x<len(phosp_list_S):
            for item in sequence:
                if item[0]==letter:
                    krotka=(item[0],phosp_list_S[x])
                    permutated_list.append(krotka)
                    x+=1
                else:
                    permutated_list.append(item)
        return(permutated_list)
    else:
        return(sequence)
#permutated list is still tuples######################################################################################################################
def tuple_list_to_zeros(lista):
    lista2=[]
    for item in lista:
        lista2.append(item[1])
    return(lista2)
#this returns just 0001100000111#######################################################################################################################

##############purely informative function
def extract_sequence_from_file(file):
    file1=open(file,"r").readlines()
    sequence=[]
    for i in range(1,len(file1)):
        l=file1[i].split()
        phosp=l[3].strip(",")
        sequence.append(phosp)
    return(','.join(sequence))

#######################################################################################################################################################
#####from big phos list extract just a one letter phosphorylation list
def extract_one_letter_list(letter,big_phosp_list):
    small_phosp_list=[]
    for item in big_phosp_list:
        if item[0]==letter:
            small_phosp_list.append(item[1])
#    np.random.shuffle(small_phosp_list)#######doesnt make sense to randomize earlier
    return(small_phosp_list)
###################################returns 01010101010 shuffled#########################################################################################
########################################################################################################################################################
########################################################################################################################################################
def create_phosp_lists(file):
    file1=open(file,"r").readlines()
    phosp_list=[]
    sequence=[]
    for i in range(1,len(file1)):
        l=file1[i].split()
        aa=l[2].strip(",")
        phosp=l[3].strip(",")
        krotka=(aa,phosp)
        sequence.append(krotka)#################creates a tuples list of sequence +phosphorylation###################
        if aa=="S" or aa=="T" or aa=="Y":
            phosp_list.append(krotka)###########creates a tuples list of STY no matter the phosp state##############
#    np.random.shuffle(phosp_list)#############doesnt need to be randomized earlier
    phosp_list_S=extract_one_letter_list("S",phosp_list)#################This is not randomized in the function 
    phosp_list_T=extract_one_letter_list("T",phosp_list)#################This is not randomized in the function
    phosp_list_Y=extract_one_letter_list("Y",phosp_list)#################This is not randomized in the function
    return(sequence,phosp_list_S,phosp_list_T,phosp_list_Y)



def return_permutated_lists(sequence,phosp_list_S,phosp_list_T,phosp_list_Y):
    np.random.shuffle(phosp_list_S) 
    np.random.shuffle(phosp_list_T)
    np.random.shuffle(phosp_list_Y)
    l=permute_one_letter("S",phosp_list_S,sequence)
#    l2=tuple_list_to_zeros(l)
    m=permute_one_letter("T",phosp_list_T,l)
#    m2=tuple_list_to_zeros(m)
    n=permute_one_letter("Y",phosp_list_Y,m)
    n2=tuple_list_to_zeros(n)
    p=','.join(n2)
    print(p)
    return(p)



######## print original sequence#############################################


def create_chart2_permuts(file):
    if not os.path.isdir("chart"):
        os.system("mkdir chart")
        os.system("mkdir permuts")
    chart=open("chart/chart2"+file,"a")
    permuts=open("permuts/"+file+"permuts","w")
    seq=extract_sequence_from_file(file)
    chart.write(seq+"\n")
    (sequence,phosp_list_S,phosp_list_T,phosp_list_Y)=create_phosp_lists(file)################what is Y is empty
    print(sequence, phosp_list_S,phosp_list_T,phosp_list_Y)
    for j in range(0,100):
        j=return_permutated_lists(sequence,phosp_list_S,phosp_list_T,phosp_list_Y)
        permuts.write(j+"\n")








create_chart2_permuts(argv[1])#this calls for every *FOUR and makes chart and chartpermut files


