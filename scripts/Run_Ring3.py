from scripts.Error import error_obj
import os
from scripts.Classes import Ring_Bond
from scripts.Utils import amino_acid_map,Researched_Amino_Acid
import shutil

def Run_Ring(pdb_path,ring_bin_path,bond_list:list,temp_path,o_folder_name):
    '''
    :purpose: Compute different bonds number by Ring3
    :param pdb_path: Input a PDB file path
    :param ring_bin_path: Input location of executable file of Ring3
    :param bond_list: Input an output list to save ring bond obj
    :param temp_path: Input TMP path as saving path
    :param o_folder_name: should include unique id and WT/MUT info, like ring3_res_ID_WT
    :return: Bond Number counting dict/False
    :outpath:temp_path/o_folder_name/
    :process: Make saving path in TMP, call ring3 to run and read results to return
    '''
    if not os.path.exists(pdb_path):
        error_obj.Is_Not_Existed(Run_Ring.__name__,pdb_path)
        return False

    if not os.path.exists(ring_bin_path+'ring'):
        error_obj.Is_Not_Existed(Run_Ring.__name__,ring_bin_path+'ring')
        return False
    if os.path.exists(temp_path+o_folder_name+'/'):
        shutil.rmtree(temp_path+o_folder_name+'/')
    os.mkdir(temp_path+o_folder_name+'/')
    os.system(ring_bin_path+'ring -i '+pdb_path+' --out_dir '+temp_path+o_folder_name+'/')
    count_dict = {'HBOND': 0, 'SSBOND': 0, 'IONIC': 0, 'VDW': 0, 'PICATION': 0, 'PIPISTACK': 0}
    try:
        with open(temp_path+o_folder_name+'/'+str(pdb_path).split('/')[len(str(pdb_path).split('/'))-1]+'_ringEdges','r') as edges:
            lines=edges.readlines()
            for line in lines[1:]:
                l=line.replace('\n','').split('\t')
                bond=Ring_Bond()
                bond.Type=l[1].split(':')[0]
                count_dict[bond.Type]+=1
                bond.Chain_ID=l[0].split(':')[0]
                try:
                    bond.AA_1_Num=int(l[0].split(':')[1])
                except:
                    error_obj.Something_Wrong(Run_Ring.__name__,'int')
                    continue
                try:
                    bond.AA_1=amino_acid_map[l[0].split(':')[len(l[0].split(':'))-1]]
                except:
                    error_obj.Something_Wrong(Run_Ring.__name__,'out of aa range')
                    continue
                try:
                    bond.AA_2_Num = int(l[2].split(':')[1])
                except:
                    error_obj.Something_Wrong(Run_Ring.__name__, 'int')
                    continue
                try:
                    bond.AA_2 = amino_acid_map[l[2].split(':')[len(l[2].split(':')) - 1]]
                except:
                    error_obj.Something_Wrong(Run_Ring.__name__, 'out of aa range')
                    continue
                bond.Atom_1=l[6]
                bond.Atom_2=l[7]
                bond.Distance=l[3]
                bond.Angle=l[4]
                bond.Energy=l[5]
                bond_list.append(bond)
    except:
        error_obj.Something_Wrong(Run_Ring.__name__, 'open_edges')
        return False

    # shutil.rmtree(ring_bin_path+'res/')
    # shutil.rmtree(temp_path+o_folder_name+'/')

    return count_dict

def Devide_Res_of_Ring_by_Layers(ring_bond_list:list[Ring_Bond],layer_aa_list:list[Researched_Amino_Acid]):
    bond_list=[]
    aa_list=[]
    for aa in layer_aa_list:
        aa_list.append(aa.Num)
    for bond in ring_bond_list:
        if bond.AA_1_Num in aa_list or bond.AA_2_Num in aa_list:
            bond_list.append(bond)
    count_dict = {'HBOND': 0, 'SSBOND': 0, 'IONIC': 0, 'VDW': 0, 'PICATION': 0, 'PIPISTACK': 0}
    for bond in bond_list:
        count_dict[bond.Type]+=1
    return count_dict

def Judge_Bond_of_Ring(ring_bond_list:list[Ring_Bond],aa:Researched_Amino_Acid):
    '''
    :purpose: Judge input AA obj belongs to which bonds
    :param ring_bond_list: Input a ring3 bond list
    :param aa: Input a Researched AA obj
    :return: A dict record bool value to represent bond situation of input AA
    :process: Read ring bond list to match input AA and return bond info
    '''
    judge_dict={'HBOND': 0, 'SSBOND': 0, 'IONIC': 0, 'VDW': 0, 'PICATION': 0, 'PIPISTACK': 0, 'IAC': 0}
    for bond in ring_bond_list:
        if bond.AA_1_Num==aa.Num or bond.AA_2_Num==aa.Num:
            judge_dict[bond.Type]=1
    return judge_dict