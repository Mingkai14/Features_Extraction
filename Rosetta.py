import os
from Error import error_obj


def Run_Score_JD2(rostta_bin_path,rosetta_db_path,pdb_path,rosetta_terms_dict:dict):
    from Utils import Clean_Main_Directory
    Clean_Main_Directory()
    os.system(f'{rostta_bin_path}score_jd2.static.linuxgccrelease -database {rosetta_db_path} -in:file:s {pdb_path} -ignore_unrecognized_res')
    score_dict={}
    with open('./score.sc','r') as r:
        lines=r.readlines()
        for line in lines:
            line=line.replace('\n','')
            div=line.split()
            if div[0]=='SCORE:' and line.find('total_score')!=-1 and line.find('fa_atr')!=-1 and line.find('fa_dun')!=-1:
                for name in div:
                    score_dict[name]=''
            if div[0]=='SCORE:' and line.find('total_score')==-1 and line.find('fa_atr')==-1 and line.find('fa_dun')==-1:
                for i in range(len(div)):
                    score_dict[list(score_dict.keys())[i]]=div[i]
    try:
        for key in rosetta_terms_dict.keys():
            rosetta_terms_dict[key]=float(score_dict[key])
    except:
        error_obj.Something_Wrong(Run_Score_JD2.__name__)
        Clean_Main_Directory()
        return False
    Clean_Main_Directory()
    return True




def Clean_PDB_by_Rosetta(pdb_path,output_path,clean_path):
    os.system(f'{clean_path}clean_pdb.py {pdb_path} ignorechain')
    files=os.listdir('./')
    for file in files:
        if file.split('.')[len(file.split('.'))-1]=='pdb':
            with open(file,'r') as r:
                with open(output_path,'w') as w:
                    w.write(r.read())
            break
    from Utils import Clean_Main_Directory
    Clean_Main_Directory()





