import pandas as pd
import os
from ete3 import EvolTree

def get_table_codeml(directory, workdir):
    family=[]
    hp=[]
    best_M2=[]
    beb=[]
    omega=[]

    for filename in os.listdir(directory):
        if filename.endswith('out'):
            name = filename.split('.')[0]
            family.append(name)
    
            tree = EvolTree(newick = os.path.join(directory,"{}.treefile".format(filename)))

            best_model_list = []
        
            for model in ['M1', 'M2']:
                best_model = None
                best_lnl = float('-inf')
                for starting_omega in [0.2, 0.7, 1.2]:
                    modelo = model+'.'+str(starting_omega)
                    current_model = tree.link_to_evol_model(workdir+modelo+'_'+name+'/out', modelo+'_'+name)
                    current_model = tree.get_evol_model(modelo+'_'+name)
            
                    if current_model.lnL > best_lnl:
                        best_lnl = current_model.lnL
                        best_lnl_name = modelo +': '+ str(best_lnl)
                        best_model = current_model
                
                best_model_list.append(best_lnl_name.split(':')[0]+'_'+name)
    
            model2= tree.get_evol_model(best_model_list[1])
            best_M2.append(best_model_list[1])
            omega.append(model2.classes['w'][2])
            pval = tree.get_most_likely(best_model_list[1], best_model_list[0])
            if pval < 0.05:
                hp.append('M2 model wins.')
                sitios = []
                for s in range(len(model2.sites['BEB']['aa'])):
                    if model2.sites['BEB']['p2'][s] > 0.95:
                        sitios.append('site %s, position: %s' % (model2.sites['BEB']['aa'][s], s+1))
                        beb_str = '/ '.join(sitios)
                beb.append(beb_str)
            else:
                hp.append('M1 model is not rejected')
                beb.append('NaN')
                    
    data={'family':pd.Series(family), 'best_M2':pd.Series(best_M2), 'omega':pd.Series(omega), 'positive_selection_sites':pd.Series(beb), 'hypothesis':pd.Series(hp)}

    table=pd.DataFrame(data, columns=['family', 'best_M2', 'omega', 'positive_selection_sites','hypothesis'])
    table.to_csv('codeml_sites_output_all.tsv', sep='\t', index=None)
    print(table)
    

get_table_codeml(directory='/opt/adri/PR1/pipeline_codeml/mafft', workdir='./workdir/')
