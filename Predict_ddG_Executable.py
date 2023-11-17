import argparse
import joblib
import pandas as pd
from scripts.Error import error_obj
from scripts.Utils import amino_acid_map
import scripts.Global_Value
from scripts.Global_Value import *
from scripts.Init import Init
from scripts.Utils import *
from scripts.Feature_Extracting_Pred import *
from scripts.Run_Modeller import *
from scripts.MSA import *
from scripts.Record import Record_Feature_Table
from ml.DDGWizard import XGBoostRegression_Predict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--pred_dataset_path', type=str, default='')
    parser.add_argument('--db_folder_path', type=str, default='')
    parser.add_argument('--db_name', type=str, default='')
    parser.add_argument('--if_reversed_data', type=int, default=0)
    parser.add_argument('--blast_process_num', type=int, default=1)
    parser.add_argument('--mode', type=str, default='whole')
    parser.add_argument('--process_num', type=int, default=1)

    current_directory = os.getcwd()
    script_directory = os.path.dirname(os.path.abspath(__file__))
    if current_directory!=script_directory:
        error_obj.Something_Wrong(__name__,"Please cd to top folder of program!!!")
        exit(1)

    print('Processing input arguments')

    args = parser.parse_args()
    pred_dataset_path=args.pred_dataset_path
    db_folder_path=args.db_folder_path
    db_name=args.db_name
    if_reversed_data=args.if_reversed_data
    blast_process_num=args.blast_process_num
    mode=args.mode
    process_num=args.process_num


    if pred_dataset_path=='' or db_folder_path=='' or db_name=='' or args.if_reversed_data not in [0,1] or blast_process_num<1 or blast_process_num>200 or process_num>200 or process_num<1:
        error_obj.Something_Wrong(__name__,'Check your arguments')
        exit(1)
    if str(pred_dataset_path).split('.')[-1]!='xls':
        error_obj.Something_Wrong(__name__,'Check format of your pred_dataset')
        exit(1)
    if not os.path.exists(pred_dataset_path):
        error_obj.Something_Wrong(__name__,'Your pred dataset is not existed')
        exit(1)
    if not os.path.isdir(db_folder_path):
        error_obj.Something_Wrong(__name__,'Your db folder is not existed')
        exit(1)
    if str(mode) not in ['blast_only','model_only','whole']:
        error_obj.Something_Wrong(__name__,'check your mode argument')
        exit(1)




    scripts.Global_Value.MSA_DB_Path = db_folder_path
    scripts.Global_Value.MSA_DB_Name = db_name
    scripts.Global_Value.Is_Use_Reverse_Data = if_reversed_data
    scripts.Global_Value.BLAST_Process_Num = blast_process_num
    scripts.Global_Value.Mode = mode
    scripts.Global_Value.Process_Num = process_num
    scripts.Global_Value.Is_Pred=1

    print(f'Your input arguments:\n--pred_dataset_path:{pred_dataset_path}\n--db_folder_path:{scripts.Global_Value.MSA_DB_Path}\n--db_name:{scripts.Global_Value.MSA_DB_Name}\n--if_reversed_data:{scripts.Global_Value.Is_Use_Reverse_Data}\n--blast_process_num:{scripts.Global_Value.BLAST_Process_Num}\n--mode:{scripts.Global_Value.Mode}\n--process_num:{scripts.Global_Value.Process_Num}\n')

    try:
        print('Initing configuration')
        Init()

        print('Reading pred dataset ')
        Pred_Data_List = Read_Pred_XLS(pred_dataset_path)

        print('Clearing')
        files = os.listdir(Pred_Table_Path)
        for file in files:
            os.remove(Pred_Table_Path + file)
        Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,[1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0])

        print('Preparing task table')
        Prepare_Table(Pred_Data_List,Pred_Table_Path,Pred_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path)



        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='model_only':
            print('Modelling MUT models')
            Prepare_MUT_Models(Pred_Table_Path,Pred_Table_Name,MUT_PDB_Path,process_num)

        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='blast_only':
            print('Preparing Blast files')
            Prepare_Blast_Files(Pred_Table_Path, Pred_Table_Name, WT_PSSM_Data_Path, MUT_PSSM_Data_Path,
                                WT_PSI_BLAST_Data_Path, MUT_PSI_BLAST_Data_Path, WT_BLASTP_Data_Path,
                                MUT_BLASTP_Data_Path, scripts.Global_Value.MSA_DB_Path,
                                scripts.Global_Value.MSA_DB_Name)

        if scripts.Global_Value.Mode=='whole' and scripts.Global_Value.Is_Use_Reverse_Data:
            print('Adding reverse task')
            Add_Reverse_Data(Pred_Table_Path,Pred_Table_Name)

        if scripts.Global_Value.Mode == 'whole':
            Feature_Object_List = []
            print('Beginning features extraction')
            Feature_Extraction(Pred_Table_Path, Pred_Table_Name, Feature_Object_List,scripts.Global_Value.Process_Num)

            print('Recording features results')
            if not Record_Feature_Table(Feature_Object_List, Features_Table_Path):
                error_obj.Something_Wrong(__name__,'Recording failed')
                exit(1)

            print('Computing ddG with XGB model')
            XGBoostRegression_Predict(Features_Table_Path+Features_Table_Name,Model_Path,Pred_Res_Path)


    except Exception as e:
        print(e)
        error_obj.Something_Wrong(__name__,'exception')
        Clean_with_Error(None)
        exit(1)


    Remove_FoldX_Resource()
    Clean_Main_Directory()


    print('Cleaning temporary folder in ./src/TMP/')
    import shutil
    shutil.rmtree(TMP_Path)
    os.mkdir(TMP_Path)

    print('Have finished')
    exit(0)













