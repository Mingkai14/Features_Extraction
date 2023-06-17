import shutil
import xlrd
import requests
import os
from Bio import SeqIO
from Bio.PDB import *

import Global_Value
from Error import error_obj
from Classes import *
from bin.rdkit_2023_3_1.rdkit_compute import Compute_Pharmacophore_with_Rdkit
from bin.Protlego.Hydrophobic_cluster import *
from math import sqrt,pow
from Rosetta import Clean_PDB_by_Rosetta
from Docker import Docker_Remove_Container



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

secondary_structure_map={'H':0,'E':1,'C':2}

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

def Prepare_Table(Raw_Data_List,Table_Path,Res_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path):
    for Raw_Data in Raw_Data_List:
        Raw_PDB_Num=Raw_Data[0]
        Mut_Info=Raw_Data[1]
        Chain_ID=Raw_Data[2]
        pH=Raw_Data[4]
        Temperature=Raw_Data[5]
        DDG=Raw_Data[3]
        if not Prepare(Table_Path,Clean_Path,Res_Table_Name,Raw_PDB_Num,Mut_Info,Chain_ID,pH,Temperature,DDG,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path):
            error_obj.Something_Wrong(__name__,Raw_Data[0]+'_'+Raw_Data[1])
            exit(1)
    return True

def Prepare(table_path,clean_path,res_table_name,raw_pdb_num,mut_info,chain_id,pH,temperature,ddg,raw_pdb_path,w_pdb_path,m_pdb_path,raw_fasta_path,m_fasta_path,wt_pssm_data_path,mut_pssm_data_path,wt_psi_blast_data_path,mut_psi_blast_data_path):
    try:
        id=raw_pdb_num+'_'+mut_info
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
    wt_pdb_name = raw_pdb_num
    wt_pdb_path = w_pdb_path + wt_pdb_name + '.pdb'
    Clean_PDB_by_Rosetta(raw_pdb_path+raw_pdb_num+'.pdb',wt_pdb_path,clean_path)
    if not Fetch_Fasta_from_PDB(wt_pdb_path,wt_pdb_name,raw_fasta_path):
        error_obj.Something_Wrong(Prepare.__name__)
        return False
    wt_fasta_path=raw_fasta_path+wt_pdb_name+'.fasta'
    raw_seq_dict=Read_Seq_from_Fasta(wt_fasta_path)
    if raw_seq_dict==False:
        error_obj.Something_Wrong(Prepare.__name__)
        return False
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
        error_obj.Something_Wrong(Prepare.__name__)
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
    if not os.path.exists(table_path):
        os.mkdir(table_path)
    is_file_existed=os.path.exists(table_path+res_table_name)
    with open(table_path+res_table_name,'a') as table:
        if not is_file_existed:
            table.write('id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,pH,temperature,ddg\n')
        if is_file_existed:
            table.write('\n')
        table.write(id+','+wt_aa_short+','+mut_aa_short+','+str(loc)+','+str(true_loc)+','+wt_pdb_name+','+wt_pdb_path+','+mut_pdb_name+','+mut_pdb_path+','+wt_fasta_path+','+mut_fasta_path+','+wt_pssm_path+','+mut_pssm_path+','+wt_psi_blast_path+','+mut_psi_blast_path+','+str(pH)+','+str(temperature)+','+str(ddg))
    return True

def Add_Reverse_Data(table_path,table_name):
    backup_lines = []
    with open(table_path+table_name,'r') as table:
        lines=table.readlines()
        if lines[0]!='id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,pH,temperature,ddg\n':
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
            pH = item_list[15]
            temperature = item_list[16]
            ddg = item_list[17]
            if line.find('\n')==-1:
                line=line+'\n'
            backup_lines.append(line)
            new_id=id.split('_')[0]+'_'+mut_aa_short+loc+wt_aa_short
            l=f'{new_id},{mut_aa_short},{wt_aa_short},{loc},{t_loc},{mut_pdb_name},{mut_pdb_path},{wt_pdb_name},{wt_pdb_path},{mut_fasta_path},{wt_fasta_path},{mut_pssm_path},{wt_pssm_path},{mut_psi_blast_path},{wt_psi_blast_path},{pH},{temperature},{str(-float(ddg))}\n'
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
                    continue
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
    l=Run_HD_Cluser(pdb_path,main_loc)
    for hd_cluster in l:
        hd_list.append(hd_cluster)
    return len(l)

def Compute_AA_Categories(aa_list:list,pct_dict:dict,num_dict:dict):
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
    for key in dssp.keys():
        count+=1
        if dssp[key][3]>0.25:
            exposed+=1
        else:
            buried+=1
    buried_pct=buried/count
    exposed_pct=exposed/count
    return [buried_pct,exposed_pct]

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
    rsa=0.0
    is_buried_or_exposed=0
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
    return [rsa,is_buried_or_exposed]



def Run_FoldX(foldx_path,foldx_name,pdb_path,wt_aa,mut_aa,loc,chain_id,raw_dict:dict,diff_dict:dict):
    shutil.copytree(f'{foldx_path}molecules/','./molecules/')
    with open('./individual_list.txt','w') as txt:
        txt.write(wt_aa+chain_id+str(loc)+mut_aa+';')

    with open(pdb_path,'r') as pdb:
        with open('./temp.pdb','w') as new_pdb:
            new_pdb.write(pdb.read())

    os.system(f'{foldx_path}{foldx_name} --command=BuildModel --pdb=temp.pdb --mutant-file=individual_list.txt')

    files=os.listdir('./')


    for file in files:
        if os.path.isdir(file):
            continue
        is_data_line=False
        if file.split('_')[0]=='Raw':
            with open(file,'r') as f:
                for line in f.readlines():
                    if is_data_line:
                        data=line.split()
                        for i in range(len(list(raw_dict.keys()))):
                            raw_dict[list(raw_dict.keys())[i]]=data[i]
                    if line.find('Pdb')!=-1 and line.find('total')!=-1 and line.find('energy')!=-1 and line.find('Backbone')!=-1 and line.find('Electrostatics')!=-1:
                        is_data_line=True
        is_data_line = False
        if file.split('_')[0]=='Dif':
            with open(file,'r') as f:
                for line in f.readlines():
                    if is_data_line:
                        data=line.split()
                        for i in range(len(list(diff_dict.keys()))):
                            diff_dict[list(diff_dict.keys())[i]]=data[i]
                    if line.find('Pdb')!=-1 and line.find('total')!=-1 and line.find('energy')!=-1 and line.find('Backbone')!=-1 and line.find('Electrostatics')!=-1:
                        is_data_line=True
    Clean_Main_Directory()
    shutil.rmtree('./molecules/')

def Fetch_Chain_ID_from_Seq(loc:int,seq_dict:dict,wt_aa_for_test):
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
    res=Compute_Pharmacophore_with_Rdkit(pdb_path,rdkit_path,rdkit_fdef_name,aa.Central_X,aa.Central_Y,aa.Central_Z,cutoff)
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




def Run_NMA(wt_pdb_path,mut_pdb_path,loc:int,NMA_path,NMA_app_name):
    try:
        os.system(f'Rscript {NMA_path}{NMA_app_name} {wt_pdb_path} {mut_pdb_path} {loc} {NMA_path}')
        with open('./r_output.txt') as output:
                div=output.readlines()[0].split()
                Clean_Main_Directory()
                return {'wt_fluctuation_loc':float(div[0]),'mut_fluctuation_loc':float(div[1]),'rmsip':float(div[2])}
    except:
        Clean_Main_Directory()
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
    count = 0
    with open(pdb_path) as pdb:
        lines=pdb.readlines()
        loc_temp=0
        for line in lines:
            if line.split()[0]=='ATOM':
                aa=line[17:20]
                if aa in amino_acid_map.keys():
                    aa_=amino_acid_map[aa]
                else:
                    continue
                loc_=int(line[22:27].replace(' ',''))
                chain=line[21:22]
                if loc_!=loc_temp:
                    count+=1
                    loc_temp=loc_
                if aa_==aa_short and loc_==loc and chain==chain_id:
                    break
    return count




def Fetch_Single_Chain_Loc(true_loc:int,seq_dict:dict,chain_id):
    count = 0
    count_single_chain = 0
    is_end = False
    for chain_ in seq_dict.keys():
        for char in seq_dict[chain_id]:
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
    return Raw_Data_List








def Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,Control_List:list):
    if Control_List==[]:
        Control_List=[1,1,1,1,1,1,1,1,1,1]
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




def Clean_with_Error(docker_container_name):
    Clean_Main_Directory()
    if Global_Value.D_or_S=='D':
        Docker_Remove_Container(docker_container_name)
