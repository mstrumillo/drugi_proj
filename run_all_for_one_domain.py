#!/usr/bin/env python

import glob

from get_fasta_files import *
from get_clusters import *


domain=argv[1]


#create_needed_fasta_file(domain)
#print("1")
#cut_out_all_domains(domain)
#print("2")
#os.system("cp domains_compare.py proteins_%s"%domain)
#print("3")
#os.system("cp alignment_map.R proteins_%s"%domain)
#print("4")
##os.system("cp how_to_make_a_chart proteins_%s"%domain)
#os.system("cp permutations.py proteins_%s"%domain)
#print("5")
#os.system("mv needed_fasta_%s proteins_%s"%(domain, domain))
#print("fastest")
#os.chdir("proteins_%s"%domain)
#run_mafft_fastest(domain)
#print("fft")
#run_mafft_FFT_NS_i(domain)
#print("nofft")
#run_mafft_NW_NS_2(domain)
#print("6")

os.chdir("proteins_%s"%domain)
#print("9")
#run_alignment_mapR(domain)
#print("10")
#os.system("python domains_compare.py")
#print("11")


files=glob.glob("*.four")
for i in files:
    print(i)
    os.system("bsub -e e/err%s -o o/out%s ./permutations.py %s"%(i,i,i))

print("12")

