# from protlego.builder.chimera import Chimera
from Classes import *
# from moleculekit.molecule import Molecule

# mol=Molecule('./Mut_PDB/1AYF_C95S.pdb')
# # mol.view()
# print(mol)

# chimera=Chimera('../Raw_PDB/1AYF.pdb')
# # chimera.view()
# # print(chimera)
# # clusters= chimera.compute_hydrophobic_clusters(chain='A B')
# clusters= chimera.compute_hydrophobic_clusters(chain='A')
# print(clusters[0])
# print(len(clusters))
# print(type(clusters[0]))
# print(clusters[0].residues)
# print(clusters[0].area)
# print(type(clusters[0].contacts))
# print(float(clusters[0].ratio_area_residue))
# print(clusters[0].ratio_contacts_residue)

# def Run_HD_Cluser_old(pdb):
#     clusters_list=[]
#     chimera = Chimera(pdb,validateElements=False)
#     clusters = chimera.compute_hydrophobic_clusters(chain='A')
#     for cluster in clusters:
#         hd_cluster=Protlego_Hydrophobic_Cluster()
#         hd_cluster.Residues_List=cluster.residues
#         hd_cluster.Area=cluster.area
#         hd_cluster.Contacts=cluster.contacts
#         hd_cluster.Ratio_area_residue=cluster.ratio_area_residue
#         hd_cluster.Ratio_contacts_residue=cluster.ratio_contacts_residue
#         clusters_list.append(cluster)
#     return clusters_list

from .views import get_structures

def Run_HD_Cluser(pdb,main_loc):
    clusters_list = []
    res_list=get_structures(pdb.replace('./',main_loc))
    for cluster in res_list:
        hd_cluster=Protlego_Hydrophobic_Cluster()
        hd_cluster.Residues_List=cluster['residues']
        hd_cluster.Area=cluster['area']
        hd_cluster.Contacts=cluster['contacts']
        hd_cluster.Ratio_area_residue=cluster['ratio_area_residue']
        hd_cluster.Ratio_contacts_residue=cluster['ratio_contacts_residue']
        hd_cluster.Each_Chain_ID=cluster['chains']
        clusters_list.append(hd_cluster)
    return clusters_list

def Devide_Res_of_HD_Cluster_by_Layers(hd_cluster_list:list[Protlego_Hydrophobic_Cluster],layer_aa_list:list[Researched_Amino_Acid]):
    aa_list=[]
    for aa in layer_aa_list:
        aa_list.append(aa.Num)
    count=0
    for cluster in hd_cluster_list:
        for aa_num in cluster.Residues_List:
            if aa_num in aa_list:
                count+=1
                break
    return count

def Get_Max_Area(hd_cluster_list:list[Protlego_Hydrophobic_Cluster],layer_aa_list:list[Researched_Amino_Acid]):
    aa_list=[]
    for aa in layer_aa_list:
        aa_list.append(aa.Num)
    cluster_list=[]
    for cluster in hd_cluster_list:
        for aa_num in cluster.Residues_List:
            if aa_num in aa_list:
                cluster_list.append(cluster)
                break
    max=cluster_list[0].Area
    for cluster in cluster_list:
        if cluster.Area>max:
            max=cluster.Area
    return max


def Judge_If_in_Cluster(aa:Researched_Amino_Acid,cluster_list:list[Protlego_Hydrophobic_Cluster]):
    is_in_Cluster=0
    max=0
    for cluster in cluster_list:
        if aa.Num in cluster.Residues_List:
            is_in_Cluster=1
            if cluster.Area>max:
                max=cluster.Area
    return [is_in_Cluster,max]



