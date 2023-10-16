#Run bsA + bsA1 for all treefiles
from ete3 import EvolTree
import os 

directory = './'

#Test on branch-sites: mark all biotipo-C isolate leaves and run bsA model once

biotipoC = ['MDS6', 'MDS7', 'MDS8', 'MDS11', 'MDS12', 'MDS13', 'MDSP6', 'MDS30', 'MDSB1']

#Função para pegar o nome do isolado no node
def get_species_name(node_name_string):
    spcode = node_name_string.split("_")[0]
    return spcode

#Marcando todas as árvores

for filename in os.listdir(directory):
    if filename.endswith('out'):
        print('Now working with:' +str(filename))
        t= EvolTree(newick = os.path.join(directory,"{}.treefile".format(filename)),binpath='/home/users/adriayumi/miniconda3/bin/ete3_apps/bin/')
        t.link_to_alignment(alignment= os.path.join(directory,filename))
        t.workdir='./workdir/'
        OGname = filename.split('_')[0]
        t.set_species_naming_function(get_species_name)
        matches = []
        for node in t.traverse():
            if node.species in biotipoC:
                matches.append(node.node_id)
                t.mark_tree(matches, ['#1'])
        print(t.write())
        
#Test different starting omega for bsA

        for starting_omega in [0.2, 0.7, 1.2]:
            t.run_model('bsA.'+str(starting_omega)+str("_"+OGname), omega=(starting_omega))

#Test different starting omega for bsA1

        for starting_omega in [0.2, 0.7, 1.2]:
            t.run_model('bsA1.'+str(starting_omega)+str("_"+OGname), omega=(starting_omega))
