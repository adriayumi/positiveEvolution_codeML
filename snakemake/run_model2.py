from ete3 import EvolTree

treefile = snakemake.input.TREEFILE
alignment = snakemake.input.ALIGNMENT

tree = EvolTree(newick = treefile, binpath='/home/users/adriayumi/miniconda3/envs/codeml/bin/ete3_apps/bin')
tree.link_to_alignment(alignment = alignment)
tree.workdir='/opt/adri/PR1/pipeline_codeml/workdir/'
name = treefile.split('/')[6].split('.')[0]

print("Now working with " + name + " M2")

for starting_omega in [0.2, 0.7, 1.2]:
    tree.run_model('M2.'+str(starting_omega)+str("_"+name), omega=(starting_omega))
