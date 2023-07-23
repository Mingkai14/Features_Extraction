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
             'Overall_Pct_Amino_Acid_Categories':'overall_pct_aa_c_',
             'Overall_Num_Amino_Acid_Categories':'overall_num_aa_c_',
             'Layer1_Pct_Amino_Acid_Categories':'layer1_pct_aa_c_',
             'Layer1_Num_Amino_Acid_Categories':'layer1_num_aa_c_',
             'Layer2_Pct_Amino_Acid_Categories':'layer2_pct_aa_c_',
             'Layer2_Num_Amino_Acid_Categories':'layer2_num_aa_c_',
             'Layer3_Pct_Amino_Acid_Categories':'layer3_pct_aa_c_',
             'Layer3_Num_Amino_Acid_Categories':'layer3_num_aa_c_',
             'Overall_Pct_Secondary_Structure':'overall_pct_ss_',
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
             'Diff_AAindex1':'aaindex1_',
             'Overall_AAindex2':'aaindex2_',
             'Overall_AAindex3':'aaindex3_'}





def Record_Feature_Table(Feature_Obj_List:list[Feature_Object],Folder_Path):
    try:
        files=os.listdir(Folder_Path)
        for file in files:
            os.remove(Folder_Path+file)

        name_dict=Feature_Obj_List[0].__dict__
        name_list=name_dict.keys()
        name_w=[]

        for name in name_list:
            if name in pass_list:
                continue
            if name in expand_dict.keys():
                temp_dict=name_dict[name]
                for temp_name in dict(temp_dict).keys():
                    if temp_name!='Pdb' and temp_name!='total_score':
                        name_w.append(expand_dict[name]+temp_name)
            else:
                name_w.append(name)

        with open(Folder_Path+'features_table.csv', 'a', newline='') as w:
            w_csv = csv.writer(w, dialect='excel')
            w_csv.writerow(name_w)



        count=1
        error_count=0
        for obj in Feature_Obj_List:
            assert isinstance(obj, Feature_Object)

            data_dict=obj.__dict__
            data_w=[]
            data_w.append(data_dict['ID'])
            for key in data_dict.keys():
                if key in pass_list:
                    continue
                if type(data_dict[key])==str:
                    error_count+=1
                    continue
                if type(data_dict[key]) != dict:
                    try:
                        int(data_dict[key])
                    except:
                        error_count += 1
                        continue
                if key in expand_dict.keys():
                    temp_dict=data_dict[key]
                    for temp_name in dict(temp_dict).keys():
                        if temp_name != 'Pdb' and temp_name!='total_score':
                            data_w.append(temp_dict[temp_name])
                else:
                    if data_dict[key]==-99:
                        data_dict[key]=0
                    data_w.append(data_dict[key])
            if len(data_w)!= len(name_w):
                error_obj.Something_Wrong(Record_Feature_Table.__name__)
                return False

            with open(Folder_Path + 'features_table.csv', 'a', newline='') as w:
                w_csv = csv.writer(w, dialect='excel')
                w_csv.writerow(data_w)

            count += 1
        return True
    except:
        error_obj.Something_Wrong(Record_Feature_Table.__name__)
        return False




