#!/usr/bin/env python

import os
import subprocess
from sys import argv
import Bio
import glob
import numpy as np
from scipy.stats import norm
import scipy


#done_3wxm_G_4/
def get_chart_values(done_dir):
    start=int(done_dir.split("_")[-1])
#    pdb="_".join(done_dir.split("_")[1:])
    pdb=done_dir[5:]
    print(start)
    means=open('%s/means'%done_dir,'r').readlines() 
    stdev=open('%s/stdev'%done_dir,'r').readlines()
    chart=open('%s/chart_window'%done_dir,'r').readlines()
    m=means[0].split()
    s=stdev[0].split()
    c=chart[0].split()
    return(start,m,s,c,pdb)
#
def count_zscore(value,mean,stdev):
    zscore=(value-mean)/stdev
    return(zscore)
#
def return_in_pos(pdb,start,i,m,s,c,output):
    value=float(c[i])
    mean=float(m[i])
    st=float(s[i])
    if value>mean:
        z=count_zscore(value,mean,st)
        p_val = scipy.stats.norm.sf(abs(z))
        if p_val<=0.005:
            output.write("%s, %s, %s, %s\n"%(pdb, z, p_val, i+start))
    return()



done=glob.glob("done_*")
print(done)
output=open("interesting_positions","w")



for i in done:
#    output.write(i+"\n")
    (start,m,s,c,pdb)=get_chart_values(i)
    ali_length=len(c)
    for i in range(0,ali_length):
        return_in_pos(pdb,start,i,m,s,c,output)

