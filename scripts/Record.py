from scripts.Error import error_obj
from scripts.Classes import *
import xlwt,xlrd,csv
import os




pass_list=['WT_Structure','MUT_Structure','WT_Amino_Acid_short','MUT_Amino_Acid_short','Loc_of_Mutation',
           'True_Loc_of_Mutation','WT_Sequence_path','MUT_Sequence_path','WT_Amino_Acid','MUT_Amino_Acid',
           'WT_Amino_Acid_List','MUT_Amino_Acid_List','WT_Amino_Acid_List_Layer1','WT_Amino_Acid_List_Layer2',
           'WT_Amino_Acid_List_Layer3','MUT_Amino_Acid_List_Layer1','MUT_Amino_Acid_List_Layer2',
           'MUT_Amino_Acid_List_Layer3','Pfam_Num','Cutoff1','Cutoff2','Cutoff3','WT_Seq','MUT_Seq',
           'Chain_ID_of_Mut','WT_PSSM_Path','WT_PSI_BLAST_Path','MUT_PSSM_Path','MUT_PSI_BLAST_Path','WT_Ring_Bond_List',
           'MUT_Ring_Bond_List','WT_HD_Cluster_List','MUT_HD_Cluster_List','Dssp_List','COILS_line','REM465_line','HOTLOOPS_line',
           'WT_Secondary_Structure_Char','MUT_Secondary_Structure_Char',
           'WT_Psi_Angle','MUT_Psi_Angle','Diff_Psi_Angle','WT_Phi_Angle','MUT_Phi_Angle','Diff_Phi_Angle','RMSD_WT_MUT',
           'WT_FoldX_Energy_Term_Dict','WT_Rosetta_Energy_Term_Dict','MUT_Rosetta_Energy_Term_Dict',
           'Overall_Num_Amino_Acid_Categories',
           'WT_Num_Pharmacophore_Categories','MUT_Num_Pharmacophore_Categories',
           'WT_Num_Pharmacophore_Categories_Layer1','MUT_Num_Pharmacophore_Categories_Layer1',
           'WT_Num_Pharmacophore_Categories_Layer2','MUT_Num_Pharmacophore_Categories_Layer2',
           'WT_Num_Pharmacophore_Categories_Layer3','MUT_Num_Pharmacophore_Categories_Layer3',
           'WT_Psipred_List','MUT_Psipred_List','WT_BLASTP_Path','MUT_BLASTP_Path',
           'Is_Mut_Co_Evo','Co_Evo_AA_Type','Is_Group_Co_Evo','Co_Evo_Group_Num'
           ]

expand_dict={'WT_FoldX_Energy_Term_Dict':'wt_foldx_',
             'Diff_FoldX_Energy_Term_Dict':'diff_foldx_',
             'WT_Rosetta_Energy_Term_Dict':'wt_rosetta_',
             'MUT_Rosetta_Energy_Term_Dict':'mut_rosetta_',
             'Diff_Rosetta_Energy_Term_Dict':'diff_rosetta_',
             'WT_Pct_Amino_Acid_Categories':'wt_pct_aa_c_',
             'WT_Num_Amino_Acid_Categories':'wt_num_aa_c_',
             'Layer1_Pct_Amino_Acid_Categories':'layer1_pct_aa_c_',
             'Layer1_Num_Amino_Acid_Categories':'layer1_num_aa_c_',
             'Layer2_Pct_Amino_Acid_Categories':'layer2_pct_aa_c_',
             'Layer2_Num_Amino_Acid_Categories':'layer2_num_aa_c_',
             'Layer3_Pct_Amino_Acid_Categories':'layer3_pct_aa_c_',
             'Layer3_Num_Amino_Acid_Categories':'layer3_num_aa_c_',
             'WT_Pct_Secondary_Structure':'overall_pct_ss_',
             'Layer1_Pct_Secondary_Structure':'layer1_pct_ss_',
             'Layer2_Pct_Secondary_Structure':'layer2_pct_ss_',
             'Layer3_Pct_Secondary_Structure':'layer3_pct_ss_',
             'WT_Num_Pharmacophore_Categories':'wt_num_pharm_c_',
             'MUT_Num_Pharmacophore_Categories':'mut_num_pharm_c_',
             'Diff_Num_Pharmacophore_Categories':'diff_num_pharm_c_',
             'WT_Num_Pharmacophore_Categories_Layer1':'layer1_wt_num_pharm_c_',
             'MUT_Num_Pharmacophore_Categories_Layer1':'layer1_mut_num_pharm_c_',
             'Diff_Num_Pharmacophore_Categories_Layer1':'layer1_diff_num_pharm_c_',
             'WT_Num_Pharmacophore_Categories_Layer2':'layer2_wt_num_pharm_c_',
             'MUT_Num_Pharmacophore_Categories_Layer2':'layer2_mut_num_pharm_c_',
             'Diff_Num_Pharmacophore_Categories_Layer2':'layer2_diff_num_pharm_c_',
             'WT_Num_Pharmacophore_Categories_Layer3':'layer3_wt_num_pharm_c_',
             'MUT_Num_Pharmacophore_Categories_Layer3':'layer3_mut_num_pharm_c_',
             'Diff_Num_Pharmacophore_Categories_Layer3':'layer3_diff_num_pharm_c_',
             'WT_AAindex1':'wt_aaindex1_',
             'Diff_AAindex1':'diff_aaindex1_',
             'Descri_AAindex2':'descri_aaindex2_',
             'Descri_AAindex3':'descri_aaindex3_'}


# the conditions at which the stability change is measured
class_1_name=['pH','Temperature']
# the environment surrounding the mutation site
class_2_name=['WT_Pct_Amino_Acid_Categories','WT_Pct_Buried_Residue','WT_Pct_Exposed_Residue',
              'WT_Pct_Secondary_Structure','WT_Pct_coils','WT_Pct_rem465',
              'WT_Pct_hotloop','WT_NMA_Fluctuation','WT_FoldX_Energy_Term_Dict','WT_Rosetta_Energy_Term_Dict','WT_Num_HBOND_Ring','WT_Num_SSBOND_Ring','WT_Num_IONIC_Ring','WT_Num_VDW_Ring',
              'WT_Num_PICATION_Ring','WT_Num_PIPISTACK_Ring','WT_Num_HD_Cluster_Protlego','WT_Max_HD_Cluster_Area','WT_Num_Pharmacophore_Categories',
              'WT_Num_Pharmacophore_Categories_Layer1','WT_Num_Pharmacophore_Categories_Layer2','WT_Num_Pharmacophore_Categories_Layer3','WT_HD_Cluster_Area',
              'WT_B_Factor','WT_Psi','WT_Phi','WT_RSA','WT_Secondary_Structure','WT_AAindex1']
# changes in the space surrounding the mutation site
class_3_name=['Diff_NMA_Fluctuation','Diff_FoldX_Energy_Term_Dict','Diff_Rosetta_Energy_Term_Dict','Diff_Num_HBOND_Ring',
              'Diff_Num_SSBOND_Ring','Diff_Num_IONIC_Ring','Diff_Num_VDW_Ring','Diff_Num_PICATION_Ring',
              'Diff_Num_PIPISTACK_Ring','Diff_Num_HD_Cluster_Protlego','Diff_Max_HD_Cluster_Area',
              'Diff_Num_Pharmacophore_Categories','Diff_Num_Pharmacophore_Categories_Layer1','Diff_Num_Pharmacophore_Categories_Layer2',
              'Diff_Num_Pharmacophore_Categories_Layer3','Diff_HD_Cluster_Area','Diff_B_Factor','Diff_Psi','Diff_Phi','Diff_RSA','Diff_AAindex1']
# the type of mutation
class_4_name=['Descri_HBOND','Descri_SSBOND','Descri_IONIC','Descri_VDW','Descri_PICATION','Descri_PIPISTACK',
              'Descri_HD_Cluster','SIFT_Score','Descri_Buried_or_Exposed','Descri_SS','Descri_AA','Descri_Uncharged_Polar',
              'Descri_Positively_Charged_Polar','Descri_Negatively_Charged_Polar','Descri_Nonpolar','Descri_Aliphatic',
              'Descri_Aromatic','Descri_Heterocyclic','Descri_Sulfur_Containing','Descri_AAindex2','Descri_AAindex3']
# evolutionary information
class_5_name=['WT_PSSM_Score','MUT_PSSM_Score','WT_PSSM_Score_F1','WT_PSSM_Score_F2','WT_PSSM_Score_F3','WT_PSSM_Score_F4',
              'WT_PSSM_Score_F5','WT_PSSM_Score_B1','WT_PSSM_Score_B2','WT_PSSM_Score_B3','WT_PSSM_Score_B4',
              'WT_PSSM_Score_B5','WT_PSSM_Score_Aver','MUT_PSSM_Score_F1','MUT_PSSM_Score_F2','MUT_PSSM_Score_F3',
              'MUT_PSSM_Score_F4','MUT_PSSM_Score_F5','MUT_PSSM_Score_B1','MUT_PSSM_Score_B2','MUT_PSSM_Score_B3',
              'MUT_PSSM_Score_B4','MUT_PSSM_Score_B5','MUT_PSSM_Score_Aver','Diff_PSSM_Score','Diff_PSSM_Score_Aver']



def Record_Feature_Table(Feature_Obj_List:list[Feature_Object],Folder_Path):
    temp_obj=Feature_Obj_List[0]
    attribute_list=[]
    attribute_list.append('ID')
    for name in class_1_name:
        if name in expand_dict.keys():
            prefix=expand_dict[name]
            temp_dict=dict(temp_obj.__dict__[name])
            for key in temp_dict:
                attribute_list.append(prefix+key)
        else:
            attribute_list.append(name)
    for name in class_2_name:
        if name in expand_dict.keys():
            prefix = expand_dict[name]
            temp_dict = dict(temp_obj.__dict__[name])
            for key in temp_dict:
                attribute_list.append(prefix + key)
        else:
            attribute_list.append(name)
    for name in class_3_name:
        if name in expand_dict.keys():
            prefix=expand_dict[name]
            temp_dict=dict(temp_obj.__dict__[name])
            for key in temp_dict:
                attribute_list.append(prefix+key)
        else:
            attribute_list.append(name)
    for name in class_4_name:
        if name in expand_dict.keys():
            prefix=expand_dict[name]
            temp_dict=dict(temp_obj.__dict__[name])
            for key in temp_dict:
                attribute_list.append(prefix+key)
        else:
            attribute_list.append(name)
    for name in class_5_name:
        if name in expand_dict.keys():
            prefix=expand_dict[name]
            temp_dict=dict(temp_obj.__dict__[name])
            for key in temp_dict:
                attribute_list.append(prefix+key)
        else:
            attribute_list.append(name)
    attribute_list.append('Experimental_DDG')
    attribute_list.append('Experimental_DDG_Classification')

    w_data_whole=[]
    for obj in Feature_Obj_List:
        w_data=[]
        data_dict=obj.__dict__
        w_data.append(data_dict['ID'])
        for name in class_1_name:
            if name in expand_dict.keys():
                temp_dict = dict(data_dict[name])
                for key in temp_dict:
                    w_data.append(temp_dict[key])
            else:
                w_data.append(data_dict[name])
        for name in class_2_name:
            if name in expand_dict.keys():
                temp_dict = dict(data_dict[name])
                for key in temp_dict:
                    w_data.append(temp_dict[key])
            else:
                w_data.append(data_dict[name])
        for name in class_3_name:
            if name in expand_dict.keys():
                temp_dict = dict(data_dict[name])
                for key in temp_dict:
                    w_data.append(temp_dict[key])
            else:
                w_data.append(data_dict[name])
        for name in class_4_name:
            if name in expand_dict.keys():
                temp_dict = dict(data_dict[name])
                for key in temp_dict:
                    w_data.append(temp_dict[key])
            else:
                w_data.append(data_dict[name])
        for name in class_5_name:
            if name in expand_dict.keys():
                temp_dict = dict(data_dict[name])
                for key in temp_dict:
                    w_data.append(temp_dict[key])
            else:
                w_data.append(data_dict[name])
        w_data.append(data_dict['Experimental_DDG'])
        w_data.append(data_dict['Experimental_DDG_Classification'])
        w_data_whole.append(w_data)


    with open(Folder_Path+'features_table.csv', 'a', newline='') as w:
        w_csv = csv.writer(w, dialect='excel')
        w_csv.writerow(attribute_list)

    for w_l in w_data_whole:
        with open(Folder_Path + 'features_table.csv', 'a', newline='') as w:
            w_csv = csv.writer(w, dialect='excel')
            w_csv.writerow(w_l)

    return True




