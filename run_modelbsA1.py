#Run bsA for all treefiles
from ete3 import EvolTree

treefile = snakemake.input.TREEFILE
alignment = snakemake.input.ALIGNMENT

#Test on branch-sites: mark all biotipo-C isolate leaves and run bsA model once

t = EvolTree(newick = treefile, binpath='/home/users/adriayumi/miniconda3/envs/codeml/bin/ete3_apps/bin')
t.link_to_alignment(alignment = alignment)
t.workdir= '/opt/adri/PR1/pipeline_codeml/workdir/'
name = treefile.split('/')[6].split('.')[0]

print("Now working with " + name + " bsA1")

matches = []
for leaf in t:
    if leaf.name.split('_')[1].split('-')[0] == 'C':
        matches.append(leaf.node_id)
        t.mark_tree(matches, ['#1'])
print(t.write())

#Test different starting omega for branch-site
for starting_omega in [0.2, 0.7, 1.2]:
    t.run_model('bsA1.'+str(starting_omega)+str("_"+name), omega=(starting_omega))
