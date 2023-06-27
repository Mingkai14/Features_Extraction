import os
from scripts.Error import error_obj
from scripts.Utils import Fetch_Single_Chain_Loc
from scripts.Docker import *
from scripts.Global_Value import Docker_Container_Name,Singularity_Container_Path,Home_Location
import scripts.Global_Value
from scripts.Singularity import *

def Temp_Fasta(seq_dict:dict,chain_id,name,temp_path):
    with open(f'{temp_path}test.fasta','w') as temp:
        temp.write(f'>{name}_{chain_id}\n')
        temp.write(seq_dict[chain_id])


def Trans_Blast_2_Saf(name,blast_path,fasta_path,output_path):
    if scripts.Global_Value.D_or_S=='D':
        Docker_Make_Dir(Docker_Container_Name,'/home/prof_temp/')
        Docker_Import_File(Docker_Container_Name,blast_path,'/home/prof_temp/blast.output')
        Docker_Import_File(Docker_Container_Name, fasta_path, '/home/prof_temp/temp.fasta')
        Docker_Run_Cmd(Docker_Container_Name,f'/usr/share/librg-utils-perl/blastpgp_to_saf.pl fileInBlast=/home/prof_temp/blast.output fileInQuery=/home/prof_temp/temp.fasta fileOutRdb=/home/prof_temp/{name}.useless.rdb fileOutSaf=/home/prof_temp/{name}.saf red=100 maxAli=3000 tile=0')
        Docker_Export_File(Docker_Container_Name,f'/home/prof_temp/{name}.saf',f'{output_path}{name}.saf')
        Docker_Delete_Dir(Docker_Container_Name,'/home/prof_temp/')
    elif scripts.Global_Value.D_or_S=='S':
        Singularity_Make_Dir_on_Home(Singularity_Container_Path,'prof_temp')
        Singularity_Copy_File(blast_path,'~/prof_temp/blast.output')
        Singularity_Copy_File(fasta_path, '~/prof_temp/temp.fasta')
        Singularity_Run_Cmd(Singularity_Container_Path,f'/usr/share/librg-utils-perl/blastpgp_to_saf.pl fileInBlast={Home_Location}/prof_temp/blast.output fileInQuery={Home_Location}/prof_temp/temp.fasta fileOutRdb={Home_Location}/prof_temp/{name}.useless.rdb fileOutSaf={Home_Location}/prof_temp/{name}.saf red=100 maxAli=3000 tile=0')
        Singularity_Copy_File(f'~/prof_temp/{name}.saf',f'{output_path}{name}.saf')
        Singularity_Delete_Dir_on_Home(Singularity_Container_Path,'prof_temp')
    else:
        return False
    #os.system(f'{app_path}blastpgp_to_saf.pl fileInBlast={blast_path} fileInQuery={fasta_path} fileOutRdb={output_path}{name}.useless.rdb fileOutSaf={output_path}{name}.saf red=100 maxAli=3000 tile=0')
    if not os.path.exists(f'{output_path}{name}.saf'):
        return False
    with open(f'{output_path}{name}.saf','r') as saf:
        if saf.read()=='':
            return False
    return True
#blastpgp_to_saf.pl fileInBlast= fileInQuery= fileOutRdb= fileOutSaf= red=100 maxAli=3000 tile=0


#os.system(f'{app_path}copf.pl {saf_path} formatIn=saf formatOut=hssp fileOut={output_path}{name}.hssp exeConvertSeq=convert_seq')
def Trans_Saf_2_Hssp_and_filter(name,saf_path,output_path):
    if scripts.Global_Value.D_or_S=='D':
        Docker_Make_Dir(Docker_Container_Name,'/home/prof_temp/')
        Docker_Import_File(Docker_Container_Name,saf_path,'/home/prof_temp/temp.saf')
        Docker_Run_Cmd(Docker_Container_Name,f'/usr/share/librg-utils-perl/copf.pl /home/prof_temp/temp.saf formatIn=saf formatOut=hssp fileOut=/home/prof_temp/{name}.hssp exeConvertSeq=convert_seq')
        Docker_Export_File(Docker_Container_Name,f'/home/prof_temp/{name}.hssp',f'{output_path}{name}.hssp')
        Docker_Delete_Dir(Docker_Container_Name,'/home/prof_temp/')
    elif scripts.Global_Value.D_or_S=='S':
        Singularity_Make_Dir_on_Home(Singularity_Container_Path,'prof_temp')
        Singularity_Copy_File(saf_path,'~/prof_temp/temp.saf')
        Singularity_Run_Cmd(Singularity_Container_Path,f'/usr/share/librg-utils-perl/copf.pl {Home_Location}/prof_temp/temp.saf formatIn=saf formatOut=hssp fileOut={Home_Location}/prof_temp/{name}.hssp exeConvertSeq=convert_seq')
        Singularity_Copy_File(f'~/prof_temp/{name}.hssp',f'{output_path}{name}.hssp')
        Singularity_Delete_Dir_on_Home(Singularity_Container_Path,'prof_temp')
    else:
        return False
    #os.system(f'{app_path}copf.pl {saf_path} formatIn=saf formatOut=hssp fileOut={output_path}{name}.hssp exeConvertSeq=convert_seq')
    if not os.path.exists(f'{output_path}{name}.hssp'):
        return False
    with open(f'{output_path}{name}.hssp','r') as hssp:
        if hssp.read()=='':
            return False
    if scripts.Global_Value.D_or_S=='D':
        Docker_Make_Dir(Docker_Container_Name,'/home/prof_temp/')
        Docker_Import_File(Docker_Container_Name,f'{output_path}{name}.hssp','/home/prof_temp/temp.hssp')
        Docker_Run_Cmd(Docker_Container_Name,f'/usr/share/librg-utils-perl/hssp_filter.pl red=80 /home/prof_temp/temp.hssp fileOut=/home/prof_temp/{name}_filtered.hssp')
        Docker_Export_File(Docker_Container_Name,f'/home/prof_temp/{name}_filtered.hssp',f'{output_path}{name}_filtered.hssp')
        Docker_Delete_Dir(Docker_Container_Name,'/home/prof_temp/')
    elif scripts.Global_Value.D_or_S=='S':
        Singularity_Make_Dir_on_Home(Singularity_Container_Path,'prof_temp')
        Singularity_Copy_File(f'{output_path}{name}.hssp','~/prof_temp/temp.hssp')
        Singularity_Run_Cmd(Singularity_Container_Path,f'/usr/share/librg-utils-perl/hssp_filter.pl red=80 {Home_Location}/prof_temp/temp.hssp fileOut={Home_Location}/prof_temp/{name}_filtered.hssp')
        Singularity_Copy_File(f'~/prof_temp/{name}_filtered.hssp',f'{output_path}{name}_filtered.hssp')
        Singularity_Delete_Dir_on_Home(Singularity_Container_Path,'prof_temp')
    else:
        return False
    #os.system(f'{app_path}hssp_filter.pl red=80 {output_path}{name}.hssp fileOut={output_path}{name}_filtered.hssp')
    if not os.path.exists(f'{output_path}{name}_filtered.hssp'):
        return False
    with open(f'{output_path}{name}_filtered.hssp','r') as hssp:
        if hssp.read()=='':
            return False
    return True
#copf.pl <saf_formatted_file> formatIn=saf formatOut=hssp fileOut= exeConvertSeq=convert_seq
#hssp_filter.pl red=80 <hssp_formatted_file> fileOut=<filtered_hssp_formatted_file>


def Generate_RDB(name,hssp_path,output_path):
    if scripts.Global_Value.D_or_S=='D':
        Docker_Make_Dir(Docker_Container_Name,'/home/prof_temp/')
        Docker_Import_File(Docker_Container_Name,hssp_path,'/home/prof_temp/test_filtered.hssp')
        Docker_Run_Cmd(Docker_Container_Name,f'prof /home/prof_temp/test_filtered.hssp fileRdb=/home/prof_temp/{name}.hssp.prof')
        Docker_Export_File(Docker_Container_Name,f'/home/prof_temp/{name}.hssp.prof',f'{output_path}{name}.hssp.prof')
        Docker_Delete_Dir(Docker_Container_Name,'/home/prof_temp/')
    elif scripts.Global_Value.D_or_S=='S':
        Singularity_Make_Dir_on_Home(Singularity_Container_Path,'prof_temp')
        Singularity_Copy_File(hssp_path,'~/prof_temp/test_filtered.hssp')
        Singularity_Run_Cmd(Singularity_Container_Path,f'prof {Home_Location}/prof_temp/test_filtered.hssp fileRdb={Home_Location}/prof_temp/{name}.hssp.prof')
        Singularity_Copy_File(f'~/prof_temp/{name}.hssp.prof',f'{output_path}{name}.hssp.prof')
        Singularity_Delete_Dir_on_Home(Singularity_Container_Path,'prof_temp')
    else:
        return False
    #os.system(f'prof {hssp_path} fileRdb={output_path}{name}.hssp.prof')
    if not os.path.exists(f'{output_path}{name}.hssp.prof'):
        return False
    with open(f'{output_path}{name}.hssp.prof','r') as hssp:
        if hssp.read()=='':
            return False
    return True
#prof .hssp fileRdb=.hssp.prof

def Run_Profbval(name,fasta_path,prof_path,hssp_path,output_path):
    if scripts.Global_Value.D_or_S=='D':
        Docker_Make_Dir(Docker_Container_Name, '/home/prof_temp/')
        Docker_Import_File(Docker_Container_Name, fasta_path, '/home/prof_temp/test.fasta')
        Docker_Import_File(Docker_Container_Name,prof_path,'/home/prof_temp/test.hssp.prof')
        Docker_Import_File(Docker_Container_Name, hssp_path, '/home/prof_temp/test_filtered.hssp')
        Docker_Run_Cmd(Docker_Container_Name,f'profbval /home/prof_temp/test.fasta /home/prof_temp/test.hssp.prof /home/prof_temp/test_filtered.hssp /home/prof_temp/{name}.profbval 9 6')
        Docker_Export_File(Docker_Container_Name,f'/home/prof_temp/{name}.profbval',f'{output_path}{name}.profbval')
        Docker_Delete_Dir(Docker_Container_Name,'/home/prof_temp/')
    elif scripts.Global_Value.D_or_S=='S':
        Singularity_Make_Dir_on_Home(Singularity_Container_Path,'prof_temp')
        Singularity_Copy_File(fasta_path,'~/prof_temp/test.fasta')
        Singularity_Copy_File(prof_path, '~/prof_temp/test.hssp.prof')
        Singularity_Copy_File(hssp_path, '~/prof_temp/test_filtered.hssp')
        Singularity_Run_Cmd(Singularity_Container_Path, f'profbval {Home_Location}/prof_temp/test.fasta {Home_Location}/prof_temp/test.hssp.prof {Home_Location}/prof_temp/test_filtered.hssp {Home_Location}/prof_temp/{name}.profbval 9 6')
        Singularity_Copy_File(f'~/prof_temp/{name}.profbval', f'{output_path}{name}.profbval')
        Singularity_Delete_Dir_on_Home(Singularity_Container_Path,'prof_temp')
    #os.system(f'profbval {fasta_path} {prof_path} {hssp_path} {output_path}{name}.profbval 9 6')
    if not os.path.exists(f'{output_path}{name}.profbval'):
        return False
    with open(f'{output_path}{name}.profbval','r') as hssp:
        if hssp.read()=='':
            return False
    return True
#profbval .fasta .hssp.prof .hssp .profbval 9 6


def Compute_B_Factor(seq_dict:dict,chain_id,temp_path,o_folder_name,psi_blast_path,main_path,loc:int,aa_for_test):
    '''
    :purpose: Compute b-factor value of a corresponding AA
    :param seq_dict: Input sequence dict
    :param chain_id: Input chain id of this coresponding AA
    :param temp_path: Input TMP_Path
    :param o_folder_name: Input folder name of outpath
    :param psi_blast_path: Input a path of a psiblast result file
    :param main_path: Input main location
    :param loc: Input a true location to match this AA
    :param aa_for_test: Input an AA type for testing
    :return: B-factor value/False
    :outpath: temp_path/o_folder_name, like TMP_Path/prof_res_ID_WT/
    :process: 1. By sequence dict, chain id to generate fasta file in outpath
              2. By docker/singularity to transfer psiblast result file to saf format in outpath
              3. By docker/singularity to transfer saf file to hssp and filter it in outpath
              4. By docker/singularity to generate rdb file in outpath by hssp file
              5. By docker/singularity to run profbval to get b-factor file in outpath by fasta, rdb and hssp file
              6. Compute position of this AA in a single chain
              7. Read result file and return
    '''
    prof_temp_path=temp_path+o_folder_name+'/'
    if os.path.exists(prof_temp_path):
        import shutil
        shutil.rmtree(prof_temp_path)
    os.mkdir(prof_temp_path)
    Temp_Fasta(seq_dict,chain_id,'test',prof_temp_path)#test.fasta
    if not Trans_Blast_2_Saf('test',psi_blast_path,prof_temp_path+'test.fasta',prof_temp_path):#test.saf
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    if not Trans_Saf_2_Hssp_and_filter('test',prof_temp_path+'test.saf',prof_temp_path):#test.hssp test_filtered.hssp
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    if not Generate_RDB('test',main_path+str(prof_temp_path).replace('./','')+'test_filtered.hssp',main_path+str(prof_temp_path).replace('./','')):#test.hssp.prof
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    if not Run_Profbval('test',main_path+str(prof_temp_path).replace('./','')+'test.fasta',main_path+str(prof_temp_path).replace('./','')+'test.hssp.prof',main_path+str(prof_temp_path).replace('./','')+'test_filtered.hssp',main_path+str(prof_temp_path).replace('./','')):#test.profbval
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    bnorm=''
    count_single_chain=Fetch_Single_Chain_Loc(loc,seq_dict,chain_id)
    try:
        with open(prof_temp_path+'test.profbval','r') as res:
            lines=res.readlines()
            for line in lines[1:]:
                num=int(line.split('\t')[0])
                aa=line.split('\t')[1]
                if num==count_single_chain and aa==aa_for_test:
                    bnorm=line.split('\t')[3]
    except:
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    if bnorm=='':
        error_obj.Something_Wrong(Compute_B_Factor.__name__)
        return False
    # files=os.listdir(prof_temp_path)
    # for file in files:
    #     os.remove(prof_temp_path+file)
    f_bnorm=float(bnorm)
    return f_bnorm






