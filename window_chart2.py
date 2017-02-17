#!/usr/bin/env python
import numpy as np
import os

###################################################################################################
###################################################################################################
#####################################################################################################################################################################################################
#

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

print(sums)

w=count_window_for_list(sums)


for item in w:
    output.write("%.3f "%(item))
output.write("\n")








