import os
import pandas as pd
from ete3 import EvolTree


def get_table_codeml(directory, workdir):
    family=[]
    psel=[]
    rel=[]
    hp=[]
    wfrg2a=[]
    bckg2a=[]
    best_bsa=[]
    beb=[]

    for filename in os.listdir(directory):
        if filename.endswith('out'):
            name = filename.split('.')[0]
            family.append(name)
    
            tree = EvolTree(newick = os.path.join(directory,"{}.treefile".format(filename)))

            best_model_list = []
        
            for model in ['M0', 'bsA', 'bsA1']:
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
            #print(best_model_list)
             bsa= tree.get_evol_model(best_model_list[1])
            best_bsa.append(best_model_list[1])
            ps = tree.get_most_likely(best_model_list[1], best_model_list[2])
            rx = tree.get_most_likely(best_model_list[2], best_model_list[0])
            psel.append(str(ps))
            rel.append(str(rx))
            bckg2a.append(str(bsa.classes['background w'][2]))
            wfrg2a.append(str(bsa.classes['foreground w'][2])) 
            bsabeb=pd.DataFrame(bsa.sites['BEB'])
            beb_sites=pd.Series(bsabeb.loc[bsabeb['w'] > 1].loc[bsabeb['pv'] > 0.95].index)+1
            beb.append(beb_sites.tolist())

            if ps<0.05 and float(bsa.classes['foreground w'][2])>1:
                hp.append('positive selection ')
            elif rx<0.05 and ps>=0.05:
                hp.append('relaxation')
            else:
                hp.append('best fit for M0')
    
    data={'family':pd.Series(family), 'best_bsa':pd.Series(best_bsa), 'ps':pd.Series(psel), 'rx':pd.Series(rel),
          'wfrg2a': pd.Series(wfrg2a), 'bckg2a':pd.Series(bckg2a),'beb_sites':pd.Series(beb), 'hypothesis':pd.Series(hp)}

    table=pd.DataFrame(data, columns=['family', 'best_bsa', 'ps', 'rx', 'wfrg2a','bckg2a','beb_sites','hypothesis'])
    table.to_csv('branchsite_output.tsv', sep='\t', index=None)
    print(table)
    
get_table_codeml(directory='./mafft/', workdir='./workdir/')
