from ete3 import EvolTree
import os 

directory = './'

for filename in os.listdir(directory):
    if filename.endswith('out'):
        t= EvolTree(newick = os.path.join(directory,"{}.treefile".format(filename)),binpath='/home/users/adriayumi/miniconda3/bin/ete3_apps/bin/')
        t.link_to_alignment(alignment= os.path.join(directory,filename))
        t.workdir='./workdir/'
        OGname = filename.split('_')[0]
        for starting_omega in [0.2, 0.7, 1.2]:
            t.run_model('M0.'+str(starting_omega)+str("_"+OGname), omega=(starting_omega))
