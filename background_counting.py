#!/usr/bin/env python
import numpy as np
import glob
import os
import subprocess
from decimal import *

#                                                                                                                                        ############################FIXED#####################################################
#                                                                                                                                        ##first, make the file linear, not vertical
#                                                                                                                                        #########rewrite ali file#########################
#                                                                                                                                        #def rewrite_file():
#                                                                                                                                        #    chart=open("chart","r").readlines()
#                                                                                                                                        #    output=open("chart_flat","w")
#                                                                                                                                        #    sequences=len(chart[0].split()) #######how many sequences
#                                                                                                                                        #    for i in range(0,sequences):
#                                                                                                                                        #        seq=[]
#                                                                                                                                        #        for j in range(0,len(chart)):
#                                                                                                                                        #            l=chart[j].split()
#                                                                                                                                        #            value=l[i]
#                                                                                                                                        #            seq.append(value)
#                                                                                                                                        #        for item in seq:
#                                                                                                                                        #          output.write("%s," % item)
#                                                                                                                                        #        output.write("\n")
#                                                                                                                                        #    return()
#                                                                                                                                        #####################################################################################
#                                                                                                                                        #rewrite_file()
#                                                                                                                                        ####################################################################################
#



#counts average background for each of the sequences in alignment
#returns list of backgrounds and their average
####################################################################################
def normal_average():
    chart=open("chart_flat","r").readlines()
    sequences=len(chart) #######how many sequences
    list_of_bg=[]
    for j in range(0,len(chart)):
        phosps=0
        gaps=0
        l=chart[j].split(",")
        for i in range(0,len(l)-1):
            item=l[i]
            if item=="1":
                phosps+=1
            elif item=="g":
                gaps+=1
#        print(phosps,gaps,len(l)-1)
        bg=phosps/((len(l)-1)-gaps)
        list_of_bg.append(bg)
    average=np.mean(list_of_bg)
    return(list_of_bg,average)
####################################################################################
####################calling#########################################################
#(l,a)=normal_average()
#print(l,a)
####################################################################################




#define window easily by the argument
###############This takes each sequence separately#######################
#produce new file
#second method takes a window of 5 and assign to position producing new file
def window_singular_sequence(window):
    window=int(window)
    chart=open("chart_flat","r").readlines()
    output=open("chart_window_average","w")
    for j in range(0,len(chart)):
        background=[]
        l=chart[j].split(",")
        for i in range(window,len(l)-2-window):
            item=l[i]
#####################################add iteration for window automatic definition#######################################
            getcontext().prec=4
            a=Decimal(l[i-2])
            b=Decimal(l[i-1])
            c=Decimal(l[i])
            d=Decimal(l[i+1])
            e=Decimal(l[i+2])
            how_many_gaps=a+b+c+d+e
            gaps=((str(float(how_many_gaps)))[-1])
            value=int(how_many_gaps)
    #        print(value,gaps)
    #        bg=(value-int(gaps))/(window*2+1)####################################substracts the gaps so everything went minus################################
            bg=(value)/(window*2+1)
            background.append(bg)
#        print(background)
        for item in background:
            output.write("%s," % item)
        output.write("\n")
    return()
##################################################################################################
##################################################################################################
window_singular_sequence("2")
####################################################################################################################################################################################################

