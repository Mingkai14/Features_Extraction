from Error import error_obj
import os
from Classes import Ring_Bond
from Utils import amino_acid_map,Researched_Amino_Acid
import shutil

def Run_Ring(pdb_path,ring_bin_path,bond_list:list):
    if not os.path.exists(pdb_path):
        error_obj.Is_Not_Existed(Run_Ring.__name__,pdb_path)
        return False

    if not os.path.exists(ring_bin_path+'ring'):
        error_obj.Is_Not_Existed(Run_Ring.__name__,ring_bin_path+'ring')
        return False
    if os.path.exists(ring_bin_path+'res/'):
        error_obj.Something_Wrong(Run_Ring.__name__,'res_existed')
        return False
    os.mkdir(ring_bin_path+'res/')
    os.system(ring_bin_path+'ring -i '+pdb_path+' --out_dir '+ring_bin_path+'res/')
    count_dict = {'HBOND': 0, 'SSBOND': 0, 'IONIC': 0, 'VDW': 0, 'PICATION': 0, 'PIPISTACK': 0}
    try:
        with open(ring_bin_path+'res/'+str(pdb_path).split('/')[len(str(pdb_path).split('/'))-1]+'_ringEdges','r') as edges:
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

    shutil.rmtree(ring_bin_path+'res/')

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
    judge_dict={'HBOND': 0, 'SSBOND': 0, 'IONIC': 0, 'VDW': 0, 'PICATION': 0, 'PIPISTACK': 0, 'IAC': 0}
    for bond in ring_bond_list:
        if bond.AA_1_Num==aa.Num or bond.AA_2_Num==aa.Num:
            judge_dict[bond.Type]=1
    return judge_dict