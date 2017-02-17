#!/usr/bin/env python


from get_pdb_files import *



domain=argv[1]
################pdb stuff starts
os.chdir("/nfs/research2/beltrao/strumill/2PROJECT/proteomes/proteins_%s"%domain)
create_needed_pdb_file(domain)  #########output in proteomes folder
print("13")
os.chdir("../")
copy_structures(domain)         #########output in ../pdb folder
print("14")
pull_out_all_domains_from_file(domain) #########output in proteomes/protein
print("15")

make_list_to_align(domain)
print("16")
clustering(domain)




