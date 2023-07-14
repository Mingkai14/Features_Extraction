import argparse
import joblib
import pandas as pd
from scripts.Error import error_obj
from scripts.Utils import amino_acid_map
import scripts.Global_Value
from scripts.Global_Value import *
from scripts.Docker import Docker_Init_Container,Docker_Remove_Container
from scripts.Init import Init
from scripts.Utils import *
from scripts.Feature_Extracting import *
from scripts.Run_Modeller import *
from scripts.MSA import *
from scripts.Record import Record_Feature_Table
from ml.XGBoostRegression import XGBoostRegression_Predict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--pdb_name', type=str, default='')
    parser.add_argument('--pdb_path', type=str, default='')
    parser.add_argument('--variation', type=str, default='')
    parser.add_argument('--chain', type=str, default='')
    parser.add_argument('--pH', type=float, default=0.0)
    parser.add_argument('--T', type=float, default=0.0)
    parser.add_argument('--db_folder_path', type=str, default='')
    parser.add_argument('--db_name', type=str, default='')
    parser.add_argument('--psiblast_threads_num', type=int, default=4)
    parser.add_argument('--container_type', type=str, default='D')
    parser.add_argument('--mode', type=str, default='whole')

    print('Processing input arguments')

    args = parser.parse_args()
    pdb_name=args.pdb_name
    pdb_path=args.pdb_path
    vari_info=args.variation
    chain=args.chain
    pH=args.pH
    T=args.T
    db_fold_path=args.db_folder_path
    db_name=args.db_name
    psi_threads_num=args.psiblast_threads_num
    container_type=args.container_type
    mode=args.mode

    if pdb_name=='' or pdb_path=='' or vari_info=='' or chain=='' or pH==0.0 or T==0.0 or db_fold_path=='' or db_name=='' or psi_threads_num<1 or psi_threads_num>30 or container_type not in ['D','S']:
        error_obj.Something_Wrong(__name__, 'Incomplete agruments')
        exit(1)


    if not os.path.exists(pdb_path):
        error_obj.Something_Wrong(__name__,'PDB path is not existed')
        exit(1)
    if str(pdb_path).split('.')[-1]!='pdb':
        error_obj.Something_Wrong(__name__,'PDB file is wrong')
        exit(1)
    WT=vari_info[0]
    MUT=vari_info[-1]
    if (WT not in amino_acid_map.values()) or (MUT not in amino_acid_map.values()):
        error_obj.Something_Wrong(__name__, 'Vari_info is wrong')
        exit(1)
    foo=list(vari_info)
    foo.pop(0)
    foo.pop(-1)
    loc=''.join(foo)
    try:
        int(loc)
    except:
        error_obj.Something_Wrong(__name__, 'Vari_info is wrong')
        exit(1)
    if not os.path.isdir(db_fold_path):
        error_obj.Something_Wrong(__name__,'DB_path is not existed')
        exit(1)

    if mode not in ['blast_only','model_only','whole']:
        error_obj.Something_Wrong(__name__,'Mode is wrong')
        exit(1)

    scripts.Global_Value.MSA_DB_Path = db_fold_path
    scripts.Global_Value.MSA_DB_Name = db_name
    scripts.Global_Value.Psi_Threads_Num = psi_threads_num
    scripts.Global_Value.D_or_S = container_type
    scripts.Global_Value.Mode = mode


    print(f'Your input arguments:\n--pdb_path:{pdb_path}\n--variation:{vari_info}\n--chain:{chain}\n--pH:{pH}\n--T:{T}\n--db_folder_path:{scripts.Global_Value.MSA_DB_Path}\n--db_name:{scripts.Global_Value.MSA_DB_Name}\n--psiblast_threads_num:{scripts.Global_Value.Psi_Threads_Num}\n--container_type:{scripts.Global_Value.D_or_S}\n--mode:{scripts.Global_Value.Mode}\n')


    if scripts.Global_Value.D_or_S=='D':
        print('Initing Docker')
        Docker_Init_Container(Docker_Container_Name,Docker_Image_ID)

    try:
        print('Initing configuration')
        Init()

        print('Clearing')
        files = os.listdir(Pred_Table_Path)
        for file in files:
            os.remove(Pred_Table_Path + file)
        Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,[1,1,1,1,1,1,1,1,1,1,1,1,1])

        print('Preparing task table')
        if not Prepare_for_Pred(Pred_Table_Path,Clean_Path,Pred_Table_Name,pdb_name,pdb_path,vari_info,chain,pH,T,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path):
            error_obj.Something_Wrong(__name__, 'Preparing task table failed')
            exit(1)

        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='model_only':
            print('Modelling MUT models')
            Prepare_MUT_Models(Pred_Table_Path,Pred_Table_Name,MUT_PDB_Path)

        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='blast_only':
            print('Preparing Blast files')
            Prepare_Blast_Files(Pred_Table_Path, Pred_Table_Name, WT_PSSM_Data_Path, MUT_PSSM_Data_Path,
                                WT_PSI_BLAST_Data_Path, MUT_PSI_BLAST_Data_Path, WT_BLASTP_Data_Path,
                                MUT_BLASTP_Data_Path, scripts.Global_Value.MSA_DB_Path,
                                scripts.Global_Value.MSA_DB_Name)

        if scripts.Global_Value.Mode == 'whole':
            Feature_Object_List = []
            print('Beginning features extraction')
            Feature_Extraction_for_Pred(Pred_Table_Path, Pred_Table_Name, Feature_Object_List)

            print('Recording features results')
            if not Record_Feature_Table(Feature_Object_List, Features_Table_Path):
                error_obj.Something_Wrong(__name__,'Recording failed')
                exit(1)

            print('Compute ddG with XGB model')
            XGBoostRegression_Predict(Features_Table_Path+Features_Table_Name,Model_Path)


    except:
        error_obj.Something_Wrong(__name__,'exception')
        Clean_with_Error(Docker_Container_Name)
        exit(1)

    if scripts.Global_Value.D_or_S=='D':
        print('Removing Docker container')
        Docker_Remove_Container(Docker_Container_Name)

    Remove_FoldX_Resource()
    Clean_Main_Directory()

    print('Have finished')
    exit(0)













