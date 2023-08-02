import shutil
import xlrd
import xlwt
import requests
import os
from Bio import SeqIO
from Bio.PDB import *

import scripts.Global_Value
from scripts.Error import error_obj
from scripts.Classes import *
from bin.rdkit_2023_3_1.rdkit_compute import Compute_Pharmacophore_with_Rdkit
from bin.Protlego.Hydrophobic_cluster import *
from math import sqrt,pow
from scripts.Rosetta import Clean_PDB_by_Rosetta
from scripts.Docker import Docker_Remove_Container



amino_acid_map={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
                'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',}


uncharged_polar_aa=['G','S','T','C','Y','N','Q']

positively_charged_polar_aa=['K','R','H']

negatively_charged_polar_aa=['D','E']

nonpolar_aa=['A','V','L','I','P','F','W','M']

aromatic_aa=['F','Y']

aliphatic_aa=['A','V','L','I','M','N','Q','K','R','G','S','T','C','D','E']

heterocyclic_aa=['H','W']

sulfur_containing_aa=['M','C']

# secondary_structure_map={'H':0,'E':1,'C':2}
secondary_structure_map={'-':0,'H':1,'B':2,'E':3,'G':4,'I':5,'T':6,'S':7}

secondary_structure_encode={}

amino_acid_num_map={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

amino_acid_encode={}

def Amino_Acid_Encode():
    aa_map=[]
    for i in amino_acid_num_map.keys():
        for j in amino_acid_num_map.keys():
            temp=f'{i}{j}'
            aa_map.append(temp)
    for i in range(len(aa_map)):
        amino_acid_encode[aa_map[i]]=i

def SS_Encode():
    ss_map=[]
    for i in secondary_structure_map.keys():
        for j in secondary_structure_map.keys():
            temp=f'{i}{j}'
            ss_map.append(temp)
    for i in range(len(ss_map)):
        secondary_structure_encode[ss_map[i]]=i



def Get_Mutation_Description(wt_aa:Researched_Amino_Acid,mut_aa:Researched_Amino_Acid,wt_ss,mut_ss):
    wt_encode=amino_acid_num_map[wt_aa.Type_short]
    mut_encode=amino_acid_num_map[mut_aa.Type_short]
    temp=f'{wt_aa.Type_short}{mut_aa.Type_short}'
    mutation_des=amino_acid_encode[temp]

    temp=f'{wt_ss}{mut_ss}'
    mutation_des_by_ss=secondary_structure_encode[temp]
    return [wt_encode,mut_encode,mutation_des,mutation_des_by_ss]

def Prepare_Table(Raw_Data_List,Table_Path,Res_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,Pred_PDB_Path=''):
    count=0
    wrong_count=0
    right_count=0
    for Raw_Data in Raw_Data_List:
        count+=1
        if scripts.Global_Value.Is_Pred==0:
            Raw_PDB_Num=Raw_Data[0]
            Mut_Info=Raw_Data[1]
            if isinstance(Raw_Data[2], str):
                Chain_ID = Raw_Data[2]
            else:
                try:
                    Chain_ID = str(int(Raw_Data[2]))
                except:
                    error_obj.Something_Wrong(Prepare, f'Check xls file {Raw_PDB_Num}')
                    wrong_count += 1
                    continue
            pH=Raw_Data[4]
            Temperature=Raw_Data[5]
            DDG=Raw_Data[3]
            if str(Raw_PDB_Num).find('_')!=-1:
                error_obj.Something_Wrong(Prepare,f'Check xls file {Raw_PDB_Num}')
                wrong_count+=1
                continue
            if not Prepare(Table_Path,Clean_Path,Res_Table_Name,Raw_PDB_Num,Mut_Info,Chain_ID,pH,Temperature,DDG,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path):
                error_obj.Something_Wrong(__name__,Raw_Data[0]+'_'+Raw_Data[1])
                wrong_count+=1
                continue
            right_count+=1
        else:
            Raw_PDB_Name=Raw_Data[0]
            Mut_Info=Raw_Data[1]
            if isinstance(Raw_Data[2], str):
                Chain_ID = Raw_Data[2]
            else:
                try:
                    Chain_ID = str(int(Raw_Data[2]))
                except:
                    error_obj.Something_Wrong(Prepare, f'Check xls file {Raw_PDB_Num}')
                    wrong_count += 1
                    continue
            pH=Raw_Data[3]
            T=Raw_Data[4]
            if str(Raw_PDB_Name).find('_')!=-1:
                error_obj.Something_Wrong(Prepare,f'Check xls file {Raw_PDB_Name}')
                wrong_count += 1
                continue
            from scripts.Global_Value import Pred_Table_Path,Pred_Table_Name
            if not Prepare_for_Pred(Pred_Table_Path,Clean_Path,Pred_Table_Name,Raw_PDB_Name,Pred_PDB_Path,Mut_Info,Chain_ID,pH,T,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path):
                error_obj.Something_Wrong(__name__, Raw_Data[0] + '_' + Raw_Data[1])
                wrong_count += 1
                continue
            right_count += 1

    print(f'Preparing task table: There are whole {count} data in raw dataset, {wrong_count} data has been filter away and {right_count} data has been recorded')
    return True

def Prepare(table_path,clean_path,res_table_name,raw_pdb_num,mut_info,chain_id,pH,temperature,ddg,raw_pdb_path,w_pdb_path,m_pdb_path,raw_fasta_path,m_fasta_path,wt_pssm_data_path,mut_pssm_data_path,wt_psi_blast_data_path,mut_psi_blast_data_path):
    try:
        id=raw_pdb_num+'_'+chain_id+'_'+mut_info
        wt_aa_short=mut_info[0]
        mut_aa_short=mut_info[len(mut_info)-1]
        foo=str(mut_info).replace(wt_aa_short,'').replace(mut_aa_short,'')
        loc=int(foo)
        if wt_aa_short not in amino_acid_map.values() or mut_aa_short not in amino_acid_map.values():
            error_obj.Something_Wrong(Prepare.__name__)
            return False
    except:
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    if not Fetch_PDB(raw_pdb_num,raw_pdb_path):
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    true_loc = Get_True_Loc(loc, wt_aa_short, raw_pdb_path+raw_pdb_num+'.pdb',chain_id)
    if true_loc is False:
        error_obj.Something_Wrong(Prepare.__name__,f'Something wrong in {id} PDB file or Variation info does not match in PDB file')
        return False
    wt_pdb_name = raw_pdb_num
    wt_pdb_path = w_pdb_path + wt_pdb_name + '.pdb'
    Clean_PDB_by_Rosetta(raw_pdb_path+raw_pdb_num+'.pdb',w_pdb_path,clean_path,wt_pdb_name)
    if not os.path.exists(wt_pdb_path):
        error_obj.Something_Wrong(Prepare.__name__, 'PDB can not be cleaned, may only have CA')
        return False
    if not Fetch_Fasta_from_PDB(wt_pdb_path,wt_pdb_name,raw_fasta_path):
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    wt_fasta_path=raw_fasta_path+wt_pdb_name+'.fasta'
    raw_seq_dict=Read_Seq_from_Fasta(wt_fasta_path)
    if raw_seq_dict==False:
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    # if len(list(raw_seq_dict.keys()))>4:
    #     error_obj.Something_Wrong(Prepare.__name__, 'PDB has so many chains')
    #     return False
    # for key in raw_seq_dict.keys():
    #     if len(raw_seq_dict[key])>800:
    #         error_obj.Something_Wrong(Prepare.__name__, 'some chains in PDB have so many aa')
    #         return False
    mut_seq_dict = {}
    count = 0
    is_match = True
    is_find = False
    for key in raw_seq_dict.keys():
        seq_list = list(raw_seq_dict[key])
        for i in range(len(seq_list)):
            count += 1
            if count == true_loc:
                if not seq_list[i] == wt_aa_short:
                    is_match = False
                seq_list[i] = mut_aa_short
                is_find=True
        mut_seq_dict[key] = ''.join(seq_list)
    if not is_match:
        error_obj.Something_Wrong(Prepare.__name__,'There may be HETATM during chain')
        return False

    if not is_find:
        error_obj.Something_Wrong(Prepare.__name__,f'Something wrong in {id} PDB file or Variation info does not match in PDB file')
        return False

    if not Make_Fasta_from_Seq(mut_seq_dict, id, m_fasta_path):
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    mut_fasta_path=m_fasta_path+id+'.fasta'
    mut_pdb_name=id
    mut_pdb_path=''
    wt_pssm_path=''
    mut_pssm_path=''
    wt_psi_blast_path=''
    mut_psi_blast_path=''
    wt_blastp_path=''
    mut_blastp_path=''
    if not os.path.exists(table_path):
        os.mkdir(table_path)
    is_file_existed=os.path.exists(table_path+res_table_name)
    with open(table_path+res_table_name,'a') as table:
        if not is_file_existed:
            table.write('id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,pH,temperature,ddg\n')
        if is_file_existed:
            table.write('\n')
        table.write(id+','+wt_aa_short+','+mut_aa_short+','+str(loc)+','+str(true_loc)+','+wt_pdb_name+','+wt_pdb_path+','+mut_pdb_name+','+mut_pdb_path+','+wt_fasta_path+','+mut_fasta_path+','+wt_pssm_path+','+mut_pssm_path+','+wt_psi_blast_path+','+mut_psi_blast_path+','+wt_blastp_path+','+mut_blastp_path+','+str(pH)+','+str(temperature)+','+str(ddg))
    return True


def Prepare_for_Pred(table_path,clean_path,res_table_name,pdb_name,pdb_path,mut_info,chain_id,pH,temperature,w_pdb_path,m_pdb_path,raw_fasta_path,m_fasta_path,wt_pssm_data_path,mut_pssm_data_path,wt_psi_blast_data_path,mut_psi_blast_data_path):
    try:
        id = pdb_name + '_' + chain_id + '_' + mut_info
        wt_aa_short=mut_info[0]
        mut_aa_short=mut_info[-1]
        foo=str(mut_info).replace(wt_aa_short,'').replace(mut_aa_short,'')
        loc=int(foo)
        if wt_aa_short not in amino_acid_map.values() or mut_aa_short not in amino_acid_map.values():
            error_obj.Something_Wrong(Prepare_for_Pred.__name__)
            return False
    except:
        error_obj.Something_Wrong(Prepare_for_Pred.__name__)
        return False

    true_loc = Get_True_Loc(loc, wt_aa_short, pdb_path,chain_id)
    if true_loc is False:
        error_obj.Something_Wrong(Prepare_for_Pred.__name__,'Something wrong in PDB file')
        return False
    wt_pdb_name = pdb_name
    wt_pdb_path = w_pdb_path + wt_pdb_name + '.pdb'
    Clean_PDB_by_Rosetta(pdb_path,w_pdb_path,clean_path,wt_pdb_name)
    if not os.path.exists(wt_pdb_path):
        error_obj.Something_Wrong(Prepare_for_Pred.__name__, 'PDB can not be cleaned, may only have CA')
        return False
    if not Fetch_Fasta_from_PDB(wt_pdb_path,wt_pdb_name,raw_fasta_path):
        error_obj.Something_Wrong(Prepare_for_Pred.__name__)
        return False
    wt_fasta_path=raw_fasta_path+wt_pdb_name+'.fasta'
    raw_seq_dict=Read_Seq_from_Fasta(wt_fasta_path)
    if raw_seq_dict==False:
        error_obj.Something_Wrong(Prepare_for_Pred.__name__)
        return False
    # if len(list(raw_seq_dict.keys()))>4:
    #     error_obj.Something_Wrong(Prepare.__name__, 'PDB has so many chains')
    #     return False
    # for key in raw_seq_dict.keys():
    #     if len(raw_seq_dict[key])>800:
    #         error_obj.Something_Wrong(Prepare.__name__, 'some chains in PDB have so many aa')
    #         return False
    mut_seq_dict = {}
    count = 0
    is_ok = True
    for key in raw_seq_dict.keys():
        seq_list = list(raw_seq_dict[key])
        for i in range(len(seq_list)):
            count += 1
            if count == true_loc:
                if not seq_list[i] == wt_aa_short:
                    is_ok = False
                seq_list[i] = mut_aa_short
        mut_seq_dict[key] = ''.join(seq_list)
    if not is_ok:
        error_obj.Something_Wrong(Prepare_for_Pred.__name__,'There may be HETATM during chain')
        return False
    mut_pdb_name=id
    if not Make_Fasta_from_Seq(mut_seq_dict, id, m_fasta_path):
        error_obj.Something_Wrong(Prepare_for_Pred.__name__)
        return False
    mut_fasta_path=m_fasta_path + id + '.fasta'
    mut_pdb_path=''
    wt_pssm_path=''
    mut_pssm_path=''
    wt_psi_blast_path=''
    mut_psi_blast_path=''
    wt_blastp_path=''
    mut_blastp_path=''
    if not os.path.exists(table_path):
        os.mkdir(table_path)
    is_file_existed=os.path.exists(table_path+res_table_name)
    with open(table_path+res_table_name,'a') as table:
        if not is_file_existed:
            table.write('id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,pH,temperature,ddg\n')
        if is_file_existed:
            table.write('\n')
        table.write(id+','+wt_aa_short+','+mut_aa_short+','+str(loc)+','+str(true_loc)+','+wt_pdb_name+','+wt_pdb_path+','+mut_pdb_name+','+mut_pdb_path+','+wt_fasta_path+','+mut_fasta_path+','+wt_pssm_path+','+mut_pssm_path+','+wt_psi_blast_path+','+mut_psi_blast_path+','+wt_blastp_path+','+mut_blastp_path+','+str(pH)+','+str(temperature)+','+str(0))
    return True

def Add_Reverse_Data(table_path,table_name):
    backup_lines = []
    with open(table_path+table_name,'r') as table:
        lines=table.readlines()
        if lines[0]!='id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,pH,temperature,ddg\n':
            error_obj.Something_Wrong(Add_Reverse_Data.__name__)
            exit(1)
        backup_lines.append(lines[0])
        for line in lines[1:]:
            item_list = str(line.replace('\n','')).split(',')
            id = item_list[0]
            wt_aa_short = item_list[1]
            mut_aa_short = item_list[2]
            loc = item_list[3]
            t_loc = item_list[4]
            wt_pdb_name = item_list[5]
            wt_pdb_path = item_list[6]
            mut_pdb_name = item_list[7]
            mut_pdb_path = item_list[8]
            wt_fasta_path = item_list[9]
            mut_fasta_path = item_list[10]
            wt_pssm_path = item_list[11]
            mut_pssm_path = item_list[12]
            wt_psi_blast_path = item_list[13]
            mut_psi_blast_path = item_list[14]
            wt_blastp_path = item_list[15]
            mut_blastp_path = item_list[16]
            pH = item_list[17]
            temperature = item_list[18]
            ddg = item_list[19]
            if line.find('\n')==-1:
                line=line+'\n'
            backup_lines.append(line)
            new_id=id.split('_')[0]+'_'+id.split('_')[1]+'_'+mut_aa_short+loc+wt_aa_short
            l=f'{new_id},{mut_aa_short},{wt_aa_short},{loc},{t_loc},{mut_pdb_name},{mut_pdb_path},{wt_pdb_name},{wt_pdb_path},{mut_fasta_path},{wt_fasta_path},{mut_pssm_path},{wt_pssm_path},{mut_psi_blast_path},{wt_psi_blast_path},{mut_blastp_path},{wt_blastp_path},{pH},{temperature},{str(-float(ddg))}\n'
            backup_lines.append(l)
    with open(table_path+table_name,'w') as w_table:
        for line in backup_lines:
            w_table.write(line)










def Fetch_Fasta_from_PDB(pdb_path, fasta_name, fasta_path):
    if os.path.exists(fasta_path)!=True:
        os.mkdir(fasta_path)
    files = os.listdir(fasta_path)
    fasta_names = []
    for file in files:
        fasta_names.append(file.split('.')[0])
    if fasta_name in fasta_names:
        return True
    pdb = PDBParser(QUIET=True)
    if not os.path.exists(pdb_path):
        error_obj.Is_Not_Existed(Fetch_Fasta_from_PDB.__name__, pdb_path)
        return False
    seq_dict={}
    structure = pdb.get_structure('foo', pdb_path)
    for chains in structure:
        for chain in chains:
            if chain.id not in seq_dict.keys():
                seq_dict[chain.id]=''
            for residue in chain:
                if residue.resname not in amino_acid_map.keys():
                    error_obj.Something_Wrong(Fetch_Fasta_from_PDB.__name__,'Amino acid is out of range')
                    return False
                seq_dict[chain.id]+=amino_acid_map[residue.resname]
    with open(fasta_path + fasta_name + '.fasta', 'w') as fasta:
        for key in seq_dict.keys():
            fasta.write('>' + fasta_name + '_' + key + '\n' + seq_dict[key])
            if key!=list(seq_dict.keys())[len(seq_dict.keys())-1]:
                fasta.write('\n')
    return True



def Read_Seq_from_Fasta(fasta_path):
    if os.path.exists(fasta_path)!=True:
        error_obj.Is_Not_Existed(Read_Seq_from_Fasta.__name__,fasta_path)
        return False
    query_seqres = SeqIO.parse(fasta_path, 'fasta')
    seq_dict={}
    for chain in query_seqres:
        seq_dict[str(chain.id).split('_')[len(str(chain.id).split('_'))-1]]=str(chain.seq)
    return seq_dict

def Read_Seq_from_AA_List(seq_dict:dict,aa_list:list):
    '''
    :purpose: By list of all AA obj, to generate a dict of sequence, this dict is divided by chain ID
    :param seq_dict: Input a output dict for receiving sequence info
    :param aa_list: Input a list of all AA obj
    :return: True/False
    :process: Read aa_list, divided by chain ID to record sequence
    '''
    try:
        for aa in aa_list:
            if aa.Chain_ID not in seq_dict.keys():
                seq_dict[aa.Chain_ID] = ''
            seq_dict[aa.Chain_ID] += aa.Type_short
        return True
    except:
        return False



def Make_Fasta_from_Seq(seq_dict:dict,fasta_name,fasta_path):
    if os.path.exists(fasta_path)!=True:
        os.mkdir(fasta_path)
    files = os.listdir(fasta_path)
    fasta_names = []
    for file in files:
        fasta_names.append(file.split('.')[0])
    if fasta_name in fasta_names:
        return True
    with open(fasta_path + fasta_name + '.fasta', 'w') as fasta:
        for key in seq_dict.keys():
            fasta.write('>' + fasta_name + '_' + key + '\n' + seq_dict[key])
            if key!=list(seq_dict.keys())[len(seq_dict.keys())-1]:
                fasta.write('\n')
    return True


def Fetch_PDB(pdb_name, pdb_path):
    if os.path.exists(pdb_path) == False:
        os.mkdir(pdb_path)
    files=os.listdir(pdb_path)
    pdbs_names=[]
    for file in files:
        pdbs_names.append(file.split('.')[0])
    if pdb_name in pdbs_names:
        return True
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
                      ' AppleWebKit/537.36 (KHTML, like Gecko)'
                      ' Chrome/99.0.4844.51 Safari/537.36Accept:'
                      ' text/html,application/xhtml+xml,application/xml;'
                      'q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,'
                      'application/signed-exchange;v=b3;q=0.9',
              }
    res=requests.get('https://files.rcsb.org/download/'+pdb_name+'.pdb',headers)
    if res.status_code!=200:
        error_obj.Request_Error(Fetch_PDB.__name__,pdb_name)
        return False
    else:
        with open(pdb_path+pdb_name+'.pdb', 'wb') as pdb:
            pdb.write(res.content)
        return True




def Get_Reasearched_Amino_Acid(amino_acid:Researched_Amino_Acid, pdb_name,pdb_path,loc_mutation, amino_acid_short_for_test):
    '''
    :purpose: Obtain info of WT/MUT and other AA
    :param amino_acid:Input a Researched_AA obj, which is waiting to assign value. This obj to record AA info
    :param pdb_name: Input a PDB file name
    :param pdb_path: Input a PDB file path
    :param loc_mutation: Input true location (after clean, may not equal to location in raw data) of AA position
    :param amino_acid_short_for_test: Input AA type for testing
    :return: True/False
    :process: By Biopthon PDBParser to fill AA info as well as central coordinate
    '''
    pdb = PDBParser(QUIET=True)
    if not os.path.exists(pdb_path):
        error_obj.Is_Not_Existed(Get_Reasearched_Amino_Acid.__name__,pdb_path)
        return False
    structure = pdb.get_structure(pdb_name,pdb_path)
    residue_l=[]
    chain_id=''
    count=0
    for chains in structure:
        for chain in chains:
            for residue in chain:
                count+=1
                residue_l.append(residue)
                if count==loc_mutation:
                    chain_id=chain.id
    try:
        residue_l[loc_mutation-1]
    except:
        error_obj.Something_Wrong(Get_Reasearched_Amino_Acid.__name__,'Out of array range')
        return False

    if residue_l[loc_mutation-1].resname not in amino_acid_map.keys():
        error_obj.Something_Wrong(Get_Reasearched_Amino_Acid.__name__,'Amino acid is out of range')
        return False
    if amino_acid_map[residue_l[loc_mutation-1].resname]!=amino_acid_short_for_test:
        error_obj.Something_Wrong(Get_Reasearched_Amino_Acid.__name__,'residue')
        return False

    amino_acid.Type=residue_l[loc_mutation-1].resname
    amino_acid.Type_short=amino_acid_map[residue_l[loc_mutation-1].resname]
    amino_acid.Num=loc_mutation
    amino_acid.Chain_ID=chain_id
    num_atom=0
    X_temp=0
    Y_temp=0
    Z_temp=0
    for atom in residue_l[loc_mutation-1]:
        a=Atom()
        a.Atom_Name=atom.name
        a.Atom_Full_Name=atom.fullname
        a.Element=atom.element
        a.X=atom.coord[0]
        X_temp+=a.X
        a.Y=atom.coord[1]
        Y_temp+=a.Y
        a.Z=atom.coord[2]
        Z_temp+=a.Z
        amino_acid.Atom_List.append(a)
        num_atom+=1
    amino_acid.Central_X=X_temp/num_atom
    amino_acid.Central_Y=Y_temp/num_atom
    amino_acid.Central_Z=Z_temp/num_atom
    return True


def Get_All_Amino_Acid(return_list:list,pdb_name,pdb_path):
    '''
    :purpose: Obtain info of all AA in protein, recorded into a list with Research_AA obj
    :param return_list: Input a list to receive all AA obj
    :param pdb_name: Input a PDB file name
    :param pdb_path: Input a PDB file path
    :return: True/False
    :process: By biopython to get all AA basic info and call Get_Researched_AA function to generate each AA obj
    '''
    pdb = PDBParser(QUIET=True)
    if not os.path.exists(pdb_path):
        error_obj.Is_Not_Existed(Get_All_Amino_Acid.__name__, pdb_path)
        return False
    chain_l = []
    structure = pdb.get_structure(pdb_name, pdb_path)
    for chains in structure:
        for chain in chains:
            chain_l.append(chain)

    index=0
    count=0
    aa_count=0
    not_count=0
    for chain in chain_l:
        for residue in chain:
            count+=1
            index+=1
            research_amino_acid=Researched_Amino_Acid()
            if residue.resname not in amino_acid_map.keys():
                error_obj.Something_Wrong(Get_All_Amino_Acid.__name__, str(count)+' amino acid is out of range')
                not_count+=1
                continue
            if not Get_Reasearched_Amino_Acid(research_amino_acid,pdb_name,pdb_path,index,amino_acid_map[residue.resname]):
                error_obj.Something_Wrong(Get_All_Amino_Acid.__name__,'foo')
                return False
            return_list.append(research_amino_acid)
            aa_count+=1
    if aa_count+not_count!=count:
        error_obj.Something_Wrong(Get_All_Amino_Acid.__name__,'Fail to get')
        return False
    return True


# def Run_HBPlus(return_list:list,hbplus_path,pdb_path,):
#     os.system(hbplus_path + 'hbplus ' + pdb_path)
#     files = os.listdir('./')
#     is_line = False
#     for file in files:
#         if file.split('.')[len(file.split('.')) - 1] == 'hb2':
#             with open(file, 'r') as hb2:
#                 lines = hb2.readlines()
#                 for line in lines:
#                     if is_line == True:
#                         temp_bond = Bond()
#                         temp_list = line.split()
#                         pre_temp = temp_list[0].split('-')[0]
#                         pre_temp = pre_temp.replace(pre_temp[0], '')
#                         for s in pre_temp:
#                             if s == '0':
#                                 pre_temp = pre_temp[1:]
#                             else:
#                                 break
#                         temp_num = int(pre_temp)
#                         temp_bond.AA_Num_1 = temp_num
#                         post_temp = temp_list[0].split('-')[1]
#                         if post_temp not in amino_acid_map.keys():
#                             error_obj.Something_Wrong(Run_HBPlus.__name__,post_temp)
#                             continue
#                         temp_bond.AA_1=amino_acid_map[post_temp]
#                         temp_bond.Atom_1=temp_list[1]
#
#
#                         pre_temp = temp_list[2].split('-')[0]
#                         pre_temp = pre_temp.replace(pre_temp[0], '')
#                         for s in pre_temp:
#                             if s == '0':
#                                 pre_temp = pre_temp[1:]
#                             else:
#                                 break
#                         temp_num = int(pre_temp)
#                         temp_bond.AA_Num_2 = temp_num
#                         post_temp = temp_list[2].split('-')[1]
#                         if post_temp not in amino_acid_map.keys():
#                             error_obj.Something_Wrong(Run_HBPlus.__name__, post_temp)
#                             continue
#                         temp_bond.AA_2 = amino_acid_map[post_temp]
#                         temp_bond.Atom_2=temp_list[3]
#
#                         temp_bond.Bond_Type='HB'
#
#                         return_list.append(temp_bond)
#
#                     if line.find('n') != -1 and line.find('s') != -1 and line.find('type') != -1 and line.find(
#                             'dist') != -1 and line.find('angle') != -1:
#                         is_line = True
#
#     for file in files:
#         if file.split('.')[len(file.split('.'))-1]=='hb2' or file=='hbdebug.dat':
#             os.remove(file)
#     return True



def Run_Prolego(pdb_path,hd_list:list,main_loc):
    '''
    :purpose: By Run_HD_Cluster function to compute HD Cluster info
    :param pdb_path: Input a PDB file path
    :param hd_list: Input an output list to save HD Cluster obj
    :param main_loc: Input the main location to meet requirement of definite path
    :return: number of HD Cluster in protein
    '''
    l=Run_HD_Cluser(pdb_path,main_loc)
    for hd_cluster in l:
        hd_list.append(hd_cluster)
    return len(l)

def Compute_AA_Categories(aa_list:list,pct_dict:dict,num_dict:dict):
    '''
    :purpose: Compute AA categories situation
    :param aa_list: Input a list of all AA obj
    :param pct_dict: Input an output dict of AA categories percentage from aa_list
    :param num_dict: Input an output dict of AA categories number from aa_list
    :return: None
    '''
    count=0
    uncharged_polar=0
    positively_charged_polar=0
    negatively_charged_polar=0
    nonpolar=0
    aromatic=0
    aliphatic=0
    heterocyclic=0
    sulfur_containing=0

    for aa in aa_list:
        count+=1
        if aa.Type_short in uncharged_polar_aa:
            uncharged_polar+=1
        if aa.Type_short in positively_charged_polar_aa:
            positively_charged_polar+=1
        if aa.Type_short in negatively_charged_polar_aa:
            negatively_charged_polar+=1
        if aa.Type_short in nonpolar_aa:
            nonpolar+=1
        if aa.Type_short in aromatic_aa:
            aromatic+=1
        if aa.Type_short in aliphatic_aa:
            aliphatic+=1
        if aa.Type_short in heterocyclic_aa:
            heterocyclic+=1
        if aa.Type_short in sulfur_containing_aa:
            sulfur_containing+=1
    if num_dict!={}:
        num_dict['uncharged_polar']=uncharged_polar
        num_dict['positively_charged_polar']=positively_charged_polar
        num_dict['negatively_charged_polar']=negatively_charged_polar
        num_dict['nonpolar']=nonpolar
        num_dict['aromatic']=aromatic
        num_dict['aliphatic']=aliphatic
        num_dict['heterocyclic']=heterocyclic
        num_dict['sulfur_containing']=sulfur_containing

    if pct_dict!={}:
        pct_dict['uncharged_polar'] = uncharged_polar/count
        pct_dict['positively_charged_polar'] =positively_charged_polar/count
        pct_dict['negatively_charged_polar'] =negatively_charged_polar/count
        pct_dict['nonpolar'] =nonpolar/count
        pct_dict['aromatic'] =aromatic/count
        pct_dict['aliphatic'] =aliphatic/count
        pct_dict['heterocyclic'] =heterocyclic/count
        pct_dict['sulfur_containing'] =sulfur_containing/count

def Judge_AA_Categories(aa:Researched_Amino_Acid):
    amino_acid_categories_map = {'uncharged_polar': 0, 'positively_charged_polar': 0, 'negatively_charged_polar': 0,
                                 'nonpolar': 0,'aliphatic':0,
                                 'aromatic': 0, 'heterocyclic': 0, 'sulfur_containing': 0}
    if aa.Type_short in uncharged_polar_aa:
        amino_acid_categories_map['uncharged_polar']=1
    if aa.Type_short in positively_charged_polar_aa:
        amino_acid_categories_map['positively_charged_polar']=1
    if aa.Type_short in negatively_charged_polar_aa:
        amino_acid_categories_map['negatively_charged_polar']=1
    if aa.Type_short in nonpolar_aa:
        amino_acid_categories_map['nonpolar']=1
    if aa.Type_short in aromatic_aa:
        amino_acid_categories_map['aromatic']=1
    if aa.Type_short in heterocyclic_aa:
        amino_acid_categories_map['heterocyclic']=1
    if aa.Type_short in aliphatic_aa:
        amino_acid_categories_map['aliphatic']=1
    if aa.Type_short in sulfur_containing_aa:
        amino_acid_categories_map['sulfur_containing']=1
    return amino_acid_categories_map

def Run_Dssp(pdb_name,pdb_path,seq_dict_for_test:dict):
    '''
    :purpose: By DSSP to get buried/exposed aa info and secondary structure percentage info
    :param pdb_name: Input a name
    :param pdb_path: Input a PDB file path
    :param seq_dict_for_test: Input sequence for testing
    :return: list of percentage of buried and exposed and a dict of each ss percentage from protein
    :process: By biopython to call DSSP to compute, iterative count and compute info of percentage
    '''
    seq=''
    for key in seq_dict_for_test.keys():
        seq+=seq_dict_for_test[key]
    p = PDBParser(QUIET=True)
    structure = p.get_structure(pdb_name, pdb_path)
    model = structure[0]
    dssp=DSSP(model,pdb_path)
    count=0
    errot_c=0
    for key in dssp.keys():
        if seq[count]!=dssp[key][1]:
            errot_c+=1
        count+=1
    if errot_c!=0:
        error_obj.Something_Wrong(Run_Dssp.__name__)
        return False
    count=0
    buried=0
    exposed=0
    ss_num_dict = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
    pct_dict = {'H': -99.99, 'B': -99.99, 'E': -99.99, 'G': -99.99, 'I': -99.99, 'T': -99.99, 'S': -99.99, '-': -99.99}
    for key in dssp.keys():
        count+=1
        if dssp[key][3]>0.25:
            exposed+=1
        else:
            buried+=1
        ss_num_dict[dssp[key][2]] += 1
    buried_pct=buried/count
    exposed_pct=exposed/count
    pct_dict['H'] = ss_num_dict['H'] / count
    pct_dict['B'] = ss_num_dict['B'] / count
    pct_dict['E'] = ss_num_dict['E'] / count
    pct_dict['G'] = ss_num_dict['G'] / count
    pct_dict['I'] = ss_num_dict['I'] / count
    pct_dict['T'] = ss_num_dict['T'] / count
    pct_dict['S'] = ss_num_dict['S'] / count
    pct_dict['-'] = ss_num_dict['-'] / count
    return [buried_pct,exposed_pct,pct_dict]

def Devide_Res_of_DSSP_by_Layers(dssp_list:list,aa_list:list[Researched_Amino_Acid],pct_dict:dict):
    aa_l=[]
    for aa in aa_list:
        aa_l.append(aa.Num)
    dssp_l=[]
    for dssp in dssp_list:
        if dssp[0] in aa_l:
            dssp_l.append(dssp)
    count=0
    ss_num_dict = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
    buried = 0
    exposed = 0
    for d in dssp_l:
        count+=1
        ss_num_dict[d[2]]+=1
        if d[3]>0.25:
            exposed+=1
        else:
            buried+=1
    buried_pct = buried / count
    exposed_pct = exposed / count
    pct_dict['H'] = ss_num_dict['H'] / count
    pct_dict['B'] = ss_num_dict['B'] / count
    pct_dict['E'] = ss_num_dict['E'] / count
    pct_dict['G'] = ss_num_dict['G'] / count
    pct_dict['I'] = ss_num_dict['I'] / count
    pct_dict['T'] = ss_num_dict['T'] / count
    pct_dict['S'] = ss_num_dict['S'] / count
    pct_dict['-'] = ss_num_dict['-'] / count
    return [buried_pct, exposed_pct]

def Get_Res_of_DSSP(pdb_name,pdb_path,seq_dict_for_test:dict,aa:Researched_Amino_Acid):
    '''
    :purpose: Targeting on one site, compute RSA, if_buried_or_exposed, ss, psi and phi info
    :param pdb_name: Input a name
    :param pdb_path: Input a PDB file path
    :param seq_dict_for_test: Input sequence for testing
    :param aa: Input an AA obj of this site
    :return: Return a list including RSA, is_buried_or_exposed, ss info and psi/phi
    :process: By biopython to call DSSP to get a dssp list, match corresponding AA to get info
    '''
    rsa=0.0
    is_buried_or_exposed=0
    ss=0
    ss_char=''
    psi=0.0
    phi=0.0
    seq = ''
    for key in seq_dict_for_test.keys():
        seq += seq_dict_for_test[key]
    p = PDBParser(QUIET=True)
    structure = p.get_structure(pdb_name, pdb_path)
    model = structure[0]
    dssp = DSSP(model, pdb_path)
    count = 0
    errot_c = 0
    for key in dssp.keys():
        if seq[count] != dssp[key][1]:
            errot_c += 1
        count += 1
    if errot_c != 0:
        error_obj.Something_Wrong(Run_Dssp.__name__)
        return False
    for key in dssp.keys():
        if aa.Num == key[1][1] and aa.Type_short==dssp[key][1]:
            rsa=dssp[key][3]
            if rsa>0.25:
                is_buried_or_exposed=1
            ss_char=dssp[key][2]
            ss=secondary_structure_map[ss_char]
            phi=dssp[key][4]
            psi=dssp[key][5]
    return [rsa,is_buried_or_exposed,int(ss),str(ss_char),float(psi),float(phi)]



def Run_FoldX(foldx_path,foldx_name,pdb_path,wt_aa,mut_aa,loc,chain_id,raw_dict:dict,diff_dict:dict,temp_path,o_folder_name):
    '''
    :purpose: Compute FoldX energy terms
    :param foldx_path: Input FoldX path
    :param foldx_name: Input FoldX program name
    :param pdb_path: Input a PDB file path
    :param wt_aa: Input wt AA for making mutation file
    :param mut_aa: Input mut AA for making mutation file
    :param loc: Input true loc for making mutation file
    :param chain_id: Input chain id for making mutation file
    :param raw_dict: Input an output dict to save raw energy terms
    :param diff_dict: Input an output dict to save diff energy terms
    :param temp_path: Input TMP Path
    :param o_folder_name: Input outpath folder name
    :return: True/False
    :outpath: like TMP/foldx_res_ID
    :process: 1. Check if resource folder exist
              2. make mutation file
              3. make temp pdb file
              4. run FoldX and output in outpath
              5. Read results and fill output dict
    '''
    outpath=temp_path+o_folder_name+'/'
    if os.path.exists(outpath):
        shutil.rmtree(outpath)
    os.mkdir(outpath)
    if not os.path.exists('./molecules/'):
        return False
    with open(f'{outpath}/individual_list.txt','w') as txt:
        txt.write(wt_aa+chain_id+str(loc)+mut_aa+';')

    with open(pdb_path,'r') as pdb:
        with open(f'{outpath}/temp.pdb','w') as new_pdb:
            new_pdb.write(pdb.read())

    os.system(f'{foldx_path}{foldx_name} --command=BuildModel --pdb=temp.pdb --pdb-dir {outpath} --mutant-file={outpath}individual_list.txt --output-dir={outpath}')

    files=os.listdir(outpath)


    for file in files:
        if os.path.isdir(outpath+file):
            continue
        is_data_line=False
        if file.split('_')[0]=='Raw':
            with open(outpath+file,'r') as f:
                for line in f.readlines():
                    if is_data_line:
                        data=line.split()
                        for i in range(len(list(raw_dict.keys()))):
                            raw_dict[list(raw_dict.keys())[i]]=data[i]
                    if line.find('Pdb')!=-1 and line.find('total')!=-1 and line.find('energy')!=-1 and line.find('Backbone')!=-1 and line.find('Electrostatics')!=-1:
                        is_data_line=True
        is_data_line = False
        if file.split('_')[0]=='Dif':
            with open(outpath+file,'r') as f:
                for line in f.readlines():
                    if is_data_line:
                        data=line.split()
                        for i in range(len(list(diff_dict.keys()))):
                            diff_dict[list(diff_dict.keys())[i]]=data[i]
                    if line.find('Pdb')!=-1 and line.find('total')!=-1 and line.find('energy')!=-1 and line.find('Backbone')!=-1 and line.find('Electrostatics')!=-1:
                        is_data_line=True
    # Clean_Main_Directory()
    # shutil.rmtree('./molecules/')
    return True

def Remove_FoldX_Resource():
    if os.path.exists('./molecules/'):
        shutil.rmtree('./molecules/')


def Fetch_Chain_ID_from_Seq(loc:int,seq_dict:dict,wt_aa_for_test):
    '''
    :purpose: By true location (after PDB clean), to locate mutation in which chain
    :param loc: Input true AA location
    :param seq_dict: Input a dict of sequence
    :param wt_aa_for_test: Input a AA type for testing
    :return: Chain ID/False
    :process: Maintain a count variable to iterative match true location then record chain ID
    '''
    count=0
    for key in seq_dict.keys():
        for aa in seq_dict[key]:
            count+=1
            if count==loc:
                if aa!=wt_aa_for_test:
                    error_obj.Something_Wrong(Fetch_Chain_ID_from_Seq.__name__)
                    return False
                else:
                    return key

def Clean_Main_Directory():
    files=os.listdir('./')
    for file in files:
        if os.path.isdir(file):
            continue
        if file.split('.')[len(file.split('.'))-1]=='py':
            continue
        os.remove(file)


def Run_Rdikit(pdb_path,rdkit_path,rdkit_fdef_name,res_dict:dict,aa:Researched_Amino_Acid,cutoff:float):
    '''
    :purpose: By function of rdkit to get all pharmacophore count surrounding AA site by a cutoff distance
    :param pdb_path: Input a PDB file path
    :param rdkit_path: Input rdkit resource path
    :param rdkit_fdef_name: Input rdkit resource name
    :param res_dict: Input an output dict of all pharmacophore count
    :param aa: Input an AA obj to extract its central x,y,z
    :param cutoff: Input a cutoff distance
    :return: True/False
    '''
    res=Compute_Pharmacophore_with_Rdkit(pdb_path,rdkit_path,rdkit_fdef_name,aa.Central_X,aa.Central_Y,aa.Central_Z,cutoff)
    if res is False:
        return False
    try:
        for key in res_dict.keys():
            res_dict[key]=res[key]
        return True
    except:
        error_obj.Something_Wrong(Run_Prolego.__name__)
        return False

def Subtract_Dict(dict_1:dict,dict_2:dict,dict_res:dict):
    try:
        for key in dict_1.keys():
            dict_res[key]=dict_2[key]-dict_1[key]
        return True
    except:
        error_obj.Something_Wrong(Subtract_Dict.__name__)
        return False

def Sulfur_Count(aa_list:list[Researched_Amino_Acid]):
    count=0
    try:
        for aa in aa_list:
            for atom in aa.Atom_List:
                if 'S'== atom.Element:
                    count+=1
        return count
    except:
        error_obj.Something_Wrong(Sulfur_Count.__name__)
        return False




def Run_NMA(wt_pdb_path,mut_pdb_path,loc:int,NMA_path,NMA_app_name,temp_path,o_folder_name):
    '''
    :purpose: By Bio3D R script to compute NMA info
    :param wt_pdb_path: Input a WT PDB file path to pass into Rscript
    :param mut_pdb_path: Input a MUT PDB file path to pass into Rscript
    :param loc: Input true location to pass into Rscript
    :param NMA_path: Input R script location
    :param NMA_app_name: Input R script name
    :outpath: like TMP_Path/nma_res_ID/
    :return: Dict/False
    '''
    outpath=temp_path+o_folder_name+'/'
    if os.path.exists(outpath):
        shutil.rmtree(outpath)
    os.mkdir(outpath)
    try:
        os.system(f'Rscript {NMA_path}{NMA_app_name} {wt_pdb_path} {mut_pdb_path} {loc} {NMA_path} {outpath}')
        with open(f'{outpath}/r_output.txt') as output:
                div=output.readlines()[0].split()
                #Clean_Main_Directory()
                return {'wt_fluctuation_loc':float(div[0]),'mut_fluctuation_loc':float(div[1]),'rmsip':float(div[2])}
    except:
        #Clean_Main_Directory()
        return False

def Run_Psipred(seq_dict:dict,chain_id,name,psipred_path):
    with open('./temp.fasta','w') as seq:
        seq.write(f'>{name}_{chain_id}\n')
        seq.write(f'{seq_dict[chain_id]}')
    os.system(f'{psipred_path}runpsipred ./temp.fasta')
    res_l=[]
    with open('./temp.ss2','r') as res:
        lines=res.readlines()
        for line in lines:
            if line.find('#')==-1 and line!='\n':
                line=line.replace('\n','')
                div=line.split()
                res_l.append([div[0],div[1],div[2]])
    Clean_Main_Directory()
    return res_l

def Get_Res_from_Psipred(seq_dict:dict,chain,aa:Researched_Amino_Acid,out_overall_dict:dict,in_list_psipred):
    count_single_chain=Fetch_Single_Chain_Loc(aa.Num,seq_dict,chain)
    H_n=0
    C_n=0
    E_n=0
    out_ss_char=''
    out_ss_num=0
    for item_list in in_list_psipred:
        if item_list[2]=='H':
            H_n+=1
        elif item_list[2]=='C':
            C_n+=1
        elif item_list[2]=='E':
            E_n+=1
        if int(item_list[0])==count_single_chain:
            if item_list[1]!=seq_dict[chain][count_single_chain-1]:
                error_obj.Something_Wrong(Get_Res_from_Psipred.__name__)
                return False
            out_ss_char=item_list[2]
            out_ss_num=secondary_structure_map[item_list[2]]
    out_overall_dict['H'] = H_n / len(seq_dict[chain])
    out_overall_dict['C'] = C_n / len(seq_dict[chain])
    out_overall_dict['E'] = E_n / len(seq_dict[chain])
    return [out_ss_char,out_ss_num]




def Run_ANGLOR(seq_dict:dict,chain_id,anglor_path,loc:int):
    with open(anglor_path+'example/e01/seq.txt','w') as seq:
        seq.write(seq_dict[chain_id])
    os.system(f'{anglor_path}ANGLOR/ANGLOR.pl e01')
    loc_=Fetch_Single_Chain_Loc(loc,seq_dict,chain_id)
    psi=''
    phi=''
    with open(anglor_path+'example/e01/psi.txt','r') as psi:
        lines=psi.readlines()
        for line in lines:
            line=line.replace('\n','')
            div=line.split()
            if div[0]==str(loc_):
                psi=div[1]
    with open(anglor_path+'example/e01/phi.txt','r') as phi:
        lines=phi.readlines()
        for line in lines:
            line=line.replace('\n','')
            div=line.split()
            if div[0]==str(loc_):
                phi=div[1]
    temp_files=os.listdir(anglor_path+'example/e01/')
    for file in temp_files:
        os.remove(anglor_path+'example/e01/'+file)
    return [float(psi),float(phi)]


def Get_Surrounding_AA(central_aa:Researched_Amino_Acid,all_aa:list[Researched_Amino_Acid],cutoff:float):
    aa_list=[]
    central_x=central_aa.Central_X
    central_y = central_aa.Central_Y
    central_z = central_aa.Central_Z
    for aa in all_aa:
        x=aa.Central_X
        y=aa.Central_Y
        z=aa.Central_Z
        dis=Get_Distance(central_x,central_y,central_z,x,y,z)
        if dis<=cutoff:
            aa_list.append(aa)
    return aa_list


def Get_Distance(x1,y1,z1,x2,y2,z2):
    dis = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2))
    return dis


def Get_True_Loc(loc:int,aa_short,pdb_path,chain_id):
    try:
        with open(pdb_path) as pdb:
            if_successful=False
            lines=pdb.readlines()
            last_aa_index=''
            for line in lines:
                if line[0:4] != "ATOM" and line[0:6] != 'HETATM':
                    continue
                else:
                    last_aa_index=line[22:27]
                    break
            aa_buffer=[]
            count=0
            for line in lines:
                if line[0:4] != "ATOM" and line[0:6] != 'HETATM':
                    continue
                else:
                    aa_index=line[22:27]
                    if aa_index!=last_aa_index:
                        hasCA=False
                        hasN=False
                        hasC=False
                        aa = aa_buffer[0][17:20]
                        if aa in amino_acid_map.keys():
                            aa_s=amino_acid_map[aa]
                        else:
                            aa_s='HETATM'
                        chain = aa_buffer[0][21:22]
                        for line_ in aa_buffer:
                            atomname = line_[12:16]
                            occupancy = float(line_[55:60])
                            if atomname == " CA " and occupancy > 0.0:
                                hasCA = True
                            if atomname == " N  " and occupancy > 0.0:
                                hasN = True
                            if atomname == " C  " and occupancy > 0.0:
                                hasC = True
                            if line_[17:20]!=aa:
                                return False
                            if line_[21:22]!=chain:
                                return False
                        if hasCA and hasN and hasC:
                            count+=1
                        if aa_s ==aa_short and chain==chain_id:
                            try:
                                loc_temp=int(last_aa_index.replace(' ', ''))
                                if loc_temp==loc:
                                    if_successful = True
                                    aa_buffer = []
                                    break
                            except:
                                pass
                        last_aa_index = aa_index
                        aa_buffer=[]

                    aa_buffer.append(line)
            if aa_buffer!=[]:
                hasCA = False
                hasN = False
                hasC = False
                last_aa_index = aa_buffer[0][22:27]
                aa = aa_buffer[0][17:20]
                if aa in amino_acid_map.keys():
                    aa_s = amino_acid_map[aa]
                else:
                    aa_s = 'HETATM'
                chain = aa_buffer[0][21:22]
                for line_ in aa_buffer:
                    atomname = line_[12:16]
                    occupancy = float(line_[55:60])
                    if atomname == " CA " and occupancy > 0.0:
                        hasCA = True
                    if atomname == " N  " and occupancy > 0.0:
                        hasN = True
                    if atomname == " C  " and occupancy > 0.0:
                        hasC = True
                    if line_[17:20] != aa:
                        return False
                    if line_[21:22] != chain:
                        return False
                if hasCA and hasN and hasC:
                    count += 1
                if aa_s == aa_short and chain == chain_id:
                    try:
                        loc_temp = int(last_aa_index.replace(' ', ''))
                        if loc_temp == loc:
                            if_successful = True
                            aa_buffer = []
                    except:
                        pass

    except:
            print(line)
            return False
    if not if_successful:
        return False
    return count




def Fetch_Single_Chain_Loc(true_loc:int,seq_dict:dict,chain_id):
    count = 0
    count_single_chain = 0
    is_end = False
    for chain_ in seq_dict.keys():
        for char in seq_dict[chain_]:
            if chain_ == chain_id:
                count_single_chain += 1
            count += 1
            if true_loc == count:
                is_end = True
                break
        if is_end: break
    return count_single_chain




def Read_XLS(Raw_Dataset_File):
    rb = xlrd.open_workbook(Raw_Dataset_File)
    rs = rb.sheet_by_index(0)
    rows = rs.nrows
    column = rs.ncols
    Raw_Data_List = []
    for i in range(1, rows):
        list_ = []
        for j in range(column):
            list_.append(rs.cell_value(i, j))
        Raw_Data_List.append(list_)
    temp_list=[]
    for data_list in Raw_Data_List:
        unique=data_list[0]+'_'+data_list[1]+'_'+str(data_list[2])
        temp_list.append(unique)
    temp_set=set(temp_list)
    if len(temp_list)!=len(temp_set):
        error_obj.Something_Wrong(Read_XLS.__name__, 'has repeated data')
        exit(1)
    return Raw_Data_List








def Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,Control_List:list):
    if Control_List==[]:
        Control_List=[1,1,1,1,1,1,1,1,1,1,1,1]
        #[1,1,0,0,1,1,0,0,0,0]
    if Control_List[0]==1:
        files=os.listdir(Table_Path)
        for file in files:
            os.remove(Table_Path+file)
    if Control_List[1] == 1:
        files = os.listdir(Features_Table_Path)
        for file in files:
            os.remove(Features_Table_Path+file)
    if Control_List[2] == 1:
        files = os.listdir(Raw_PDB_Path)
        for file in files:
            os.remove(Raw_PDB_Path+file)
    if Control_List[3] == 1:
        files = os.listdir(WT_PDB_Path)
        for file in files:
            os.remove(WT_PDB_Path+file)
    if Control_List[4] == 1:
        files = os.listdir(MUT_PDB_Path)
        for file in files:
            os.remove(MUT_PDB_Path+file)
    if Control_List[5] == 1:
        files = os.listdir(WT_Fasta_Path)
        for file in files:
            os.remove(WT_Fasta_Path+file)
    if Control_List[6] == 1:
        files = os.listdir(MUT_Fasta_Path)
        for file in files:
            os.remove(MUT_Fasta_Path+file)
    if Control_List[7] == 1:
        files = os.listdir(WT_PSSM_Data_Path)
        for file in files:
            os.remove(WT_PSSM_Data_Path + file)
    if Control_List[8] == 1:
        files = os.listdir(MUT_PSSM_Data_Path)
        for file in files:
            os.remove(MUT_PSSM_Data_Path + file)
    if Control_List[9] == 1:
        files = os.listdir(WT_PSI_BLAST_Data_Path)
        for file in files:
            os.remove(WT_PSI_BLAST_Data_Path + file)
    if Control_List[10] == 1:
        files = os.listdir(MUT_PSI_BLAST_Data_Path)
        for file in files:
            os.remove(MUT_PSI_BLAST_Data_Path + file)
    if Control_List[11] == 1:
        files = os.listdir(WT_BLASTP_Data_Path)
        for file in files:
            os.remove(WT_BLASTP_Data_Path + file)
    if Control_List[12] == 1:
        files = os.listdir(MUT_BLASTP_Data_Path)
        for file in files:
            os.remove(MUT_BLASTP_Data_Path + file)




def Clean_with_Error(docker_container_name):
    Clean_Main_Directory()
    Remove_FoldX_Resource()
    if scripts.Global_Value.D_or_S=='D':
        Docker_Remove_Container(docker_container_name)


def Generate_Raw_Dataset_for_Pred(pdb_name,vari_info,chain,pH,T,pdb_path,table_path,table_name):
    header=['PDB','Variation','Chain','pH','T']
    if os.path.exists(table_path+table_name):
        os.remove(table_path+table_name)
    if vari_info!='all' and str(vari_info).find('*')==-1:
        book=xlwt.Workbook()
        sheet=book.add_sheet('sheet1')
        for i in range(len(header)):
            sheet.write(0,i,header[i])
        sheet.write(1,0,pdb_name)
        sheet.write(1,1,vari_info)
        sheet.write(1,2,chain)
        sheet.write(1,3,pH)
        sheet.write(1,4,T)
        book.save(table_path+table_name)
        return True
    if vari_info!='all' and str(vari_info).find('*')!=-1:
        if vari_info[-1]!='*':
            return False
        if vari_info[0] not in amino_acid_map.values():
            return False
        foo = list(vari_info)
        foo.pop(0)
        foo.pop(-1)
        loc = ''.join(foo)
        try:
            int(loc)
        except:
            return False
        book = xlwt.Workbook()
        sheet = book.add_sheet('sheet1')
        for i in range(len(header)):
            sheet.write(0, i, header[i])
        aa_list=[]
        for aa in amino_acid_map.values():
            aa_list.append(aa)
        prefix=str(vari_info).replace('*','')
        for i in range(len(aa_list)):
            v=prefix+aa_list[i]
            sheet.write(i + 1, 0, pdb_name)
            sheet.write(i + 1, 1, v)
            sheet.write(i + 1, 2, chain)
            sheet.write(i + 1, 3, pH)
            sheet.write(i + 1, 4, T)
        book.save(table_path + table_name)
        return True
    if vari_info == 'all':
        record_dict={}
        pdb = PDBParser(QUIET=True)
        structure = pdb.get_structure('foo', pdb_path)
        for chains in structure:
            for chain in chains:
                for residue in chain:
                    if residue.resname not in amino_acid_map.keys():
                        continue
                    id=residue.get_id()
                    num=id[1]
                    chain_id=chain.id
                    aa=amino_acid_map[residue.resname]
                    if chain_id not in record_dict.keys():
                        record_dict[chain_id]={}
                    record_dict[chain_id][num]=aa
        record_list=[]
        for chain_ in record_dict.keys():
            for num_ in dict(record_dict[chain_]).keys():
                for aa in amino_acid_map.keys():
                    foo=[]
                    foo.append(pdb_name)
                    v=record_dict[chain_][num_]+str(num_)+amino_acid_map[aa]
                    foo.append(v)
                    foo.append(chain_)
                    foo.append(pH)
                    foo.append(T)
                    record_list.append(foo)
        book = xlwt.Workbook()
        sheet = book.add_sheet('sheet1')
        for i in range(len(header)):
            sheet.write(0, i, header[i])
        for i in range(len(record_list)):
            for j in range(len(record_list[i])):
                sheet.write(i+1,j,record_list[i][j])
        book.save(table_path + table_name)
        return True






