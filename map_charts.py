#!/usr/bin/env python
import numpy as np
import glob
import os
import subprocess
import statistics
from sys import argv


def list_of_ref_pos(plik):
    p=open(plik,"r").readlines()
    ref_seq=[]
    for i in range(1,len(p)):
        l=p[i].split(",")
        ali_pos=int(l[0].strip())
        prot_pos=l[1].strip()
        aa=l[2].strip()
        phosp=l[3].strip()
        if aa!="-":
            ref_seq.append(ali_pos)
    return(ref_seq)


def cut_out_mapped_positions(r,plik,output):
    chart=open(plik,"r").readlines()
    mapped_chart=open(output,"w")
#   print(len(chart))
    for i in range(0,len(chart)):
#       print(i)
        mapped_seq=[]
        l=chart[i].split(",")
        for pos in r:
#            print(pos)
            mapped_seq.append(l[pos-1].strip())
        ready_line=",".join(mapped_seq)+"\n"
        mapped_chart.write(ready_line)
    return()



def create_big_chart_file():
    os.chdir("chart")
    os.system("cat chart*four >../chart2")
    os.chdir("../permuts")
    os.system("cat PF*uts >../permutated_background")
    os.chdir("../")

def map_charts(reference_file):
    os.system("cp %s reference.four"%reference_file)
    r=list_of_ref_pos("reference.four")
    cut_out_mapped_positions(r,"chart2","mapped_chart")
    cut_out_mapped_positions(r,"permutated_background","mapped_background")


#create_big_chart_file()
#map_charts(something.four)

def count_window_for_list(list):
    window_values=[]
    for i in range(2,len(list)-2):
        a=int(list[i-2])
        b=int(list[i-1])
        c=int(list[i])
        d=int(list[i+1])
        e=int(list[i+2])
        how_many=a+b+c+d+e
        bg=float(how_many/5.)
        window_values.append(bg)
    return(window_values)


def window_chart():
    chart=open("mapped_chart","r").readlines()
    output=open("chart_window","w")
    seq_length=len(chart[0].split(","))
    all_values=[]
    for j in range(0,len(chart)):
        line_values=[]
        l=chart[j].split(",")
        for i in range(0,seq_length):
            value=int(l[i].strip())
            line_values.append(value)
        all_values.append(line_values)
    sums=[sum(item) for item in zip(*all_values)]
    w=count_window_for_list(sums)
    for item in w:
        output.write("%.3f "%(item))
    output.write("\n")
    return() 
    

#window_chart()



#######################################################################################################################################################
#################################so sum up every 100* sequence############################################################################
#################################get the window value #################################################################################################
#################################and take a median and average for that 5* (nr protiens) and st dev for that value#####################################
#################################gaps are treated like 0 ##############################################################################################
#######################################################################################################################################################
#permutated list is still tuples######################################################################################################################
def tuple_list_to_zeros(lista):
    lista2=[]
    for item in lista:
        lista2.append(item[1])
    return(lista2)
#this returns just 0001100000111#######################################################################################################################
#######################################################################################################################################################
def count_window_for_list(list):
    window_values=[]
    for i in range(2,len(list)-2):
        a=int(list[i-2])
        b=int(list[i-1])
        c=int(list[i])
        d=int(list[i+1])
        e=int(list[i+2])
        how_many=a+b+c+d+e
        bg=float(how_many/5.)
        window_values.append(bg)
    return(window_values)



def list_of_window_sums(plik2,permutations):
    plik=open(plik2,"r").readlines()
    seq_length=len(plik[0].split(","))
    count_medians_from_this=[]###########this will be returned
    for k in range(0,permutations):###########this will iterate every first, for every permut
        from_all_the_first_lines=[]########for every permutation we need new list
        for j in range(k,len(plik),permutations):###########from 0 to end of file, but jump every 100
            line_values=[]#########get one permutation in one list
            l=plik[j].split(",")
            for i in range(0,seq_length):
                value=int(l[i].strip())
                line_values.append(value)######List of all values from single permutations
            from_all_the_first_lines.append(line_values)#####list of lists of single permut values, for every 100 line
        sums=[sum(item) for item in zip(*from_all_the_first_lines)]###sum up all the values for all 1st permutations
        w=count_window_for_list(sums)######now get the window values for the sums
        count_medians_from_this.append(w)############this is the list of all window values per 100 lines permuts
    return(count_medians_from_this)


def count_medians(list):
    medians_all=[statistics.median(i) for i in zip(*list)]
    averages_all=[statistics.mean(i) for i in zip(*list)]
    stdevs_all=[statistics.stdev(i) for i in zip(*list)]
    krotka=(medians_all,averages_all,stdevs_all)
    return(krotka)


def write_down_medians():
    c=list_of_window_sums("mapped_background",100)
    (m,a,s)=count_medians(c)
#    print(m)
#    print(a)
#    print(s)
    medi=open("medians","w")
    aver=open("means","w")
    st=open("stdev","w")   
    for item in m:
        medi.write("%.3f "%(item))
    medi.write("\n")
    for item in a:
        aver.write("%.3f "%(item))
    aver.write("\n")
    for item in s:
        st.write("%.3f "%(item))
    st.write("\n")
    return()

