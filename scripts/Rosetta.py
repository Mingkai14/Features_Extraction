import os
from scripts.Error import error_obj



def Clean_PDB_by_Rosetta(pdb_path,wt_pdb_path,clean_path,wt_pdb_name):
    output_path = wt_pdb_path + wt_pdb_name + '.pdb'
    files = os.listdir(wt_pdb_path)
    pdbs_names = []
    for file in files:
        pdbs_names.append(file.split('.')[0])
    if wt_pdb_name in pdbs_names:
        return
    os.system(f'{clean_path}clean_pdb.py {pdb_path} ignorechain')
    files=os.listdir('./')
    for file in files:
        if file.split('.')[len(file.split('.'))-1]=='pdb':
            with open(file,'r') as r:
                with open(output_path,'w') as w:
                    w.write(r.read())
            break
    from scripts.Utils import Clean_Main_Directory
    Clean_Main_Directory()





