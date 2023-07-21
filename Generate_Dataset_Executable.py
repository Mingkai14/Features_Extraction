from scripts.Feature_Extracting import *
from scripts.Run_Modeller import *
from scripts.MSA import *
from scripts.Record import Record_Feature_Table
from scripts.Global_Value import *
import scripts.Global_Value
from scripts.Init import Init
from scripts.Utils import *
import argparse
from scripts.Docker import Docker_Init_Container,Docker_Remove_Container


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--raw_dataset_path', type=str, default='')
    parser.add_argument('--db_folder_path', type=str, default='')
    parser.add_argument('--db_name', type=str, default='')
    parser.add_argument('--if_reversed_data', type=int, default=1)
    parser.add_argument('--psiblast_threads_num', type=int, default=4)
    parser.add_argument('--container_type', type=str, default='D')
    parser.add_argument('--mode',type=str,default='whole')
    parser.add_argument('--process_num', type=int, default=1)

    print('Processing input arguments')
    args = parser.parse_args()
    if args.raw_dataset_path=='' or args.db_folder_path=='' or args.db_name=='' or args.if_reversed_data not in [0,1] or args.psiblast_threads_num<1 or args.psiblast_threads_num>30:
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.raw_dataset_path).split('.')[len(str(args.raw_dataset_path).split('.'))-1]!='xls':
        error_obj.Something_Wrong(__name__)
        exit(1)
    if not os.path.exists(args.raw_dataset_path):
        error_obj.Something_Wrong(__name__)
        exit(1)
    if not os.path.isdir(args.db_folder_path):
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.container_type) not in ['D','S']:
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.mode) not in ['blast_only','model_only','whole']:
        error_obj.Something_Wrong(__name__)
        exit(1)

    scripts.Global_Value.Raw_Dataset_file=args.raw_dataset_path
    scripts.Global_Value.MSA_DB_Path=args.db_folder_path
    scripts.Global_Value.MSA_DB_Name=args.db_name
    scripts.Global_Value.Is_Use_Reverse_Data=args.if_reversed_data
    scripts.Global_Value.Psi_Threads_Num=args.psiblast_threads_num
    scripts.Global_Value.D_or_S = args.container_type
    scripts.Global_Value.Mode = args.mode
    scripts.Global_Value.Process_Num = args.process_num

    print(f'Your input arguments:\n--raw_dataset_path:{scripts.Global_Value.Raw_Dataset_file}\n--db_folder_path:{scripts.Global_Value.MSA_DB_Path}\n--db_name:{scripts.Global_Value.MSA_DB_Name}\n--if_reversed_data:{scripts.Global_Value.Is_Use_Reverse_Data}\n--psiblast_threads_num:{scripts.Global_Value.Psi_Threads_Num}\n--container_type:{scripts.Global_Value.D_or_S}\n--mode:{scripts.Global_Value.Mode}\n--process_num:{scripts.Global_Value.Process_Num}\n')

    if scripts.Global_Value.D_or_S=='D':
        print('Initing Docker')
        Docker_Init_Container(Docker_Container_Name,Docker_Image_ID)

    try:
        print('Initing configuration')
        Init()

        print('Reading raw dataset ')
        Raw_Data_List = Read_XLS(scripts.Global_Value.Raw_Dataset_file)

        print('Clearing folders')
        Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,[1,1,1,0,0,1,1,0,0,0,0,0,0])

        print('Preparing task table')
        Prepare_Table(Raw_Data_List,Table_Path,Res_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path)


        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='model_only':
            print('Modelling MUT models')
            Prepare_MUT_Models(Table_Path,Res_Table_Name,MUT_PDB_Path)

        if scripts.Global_Value.Mode=='whole' or scripts.Global_Value.Mode=='blast_only':
            print('Preparing Blast files')
            Prepare_Blast_Files(Table_Path,Res_Table_Name,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,scripts.Global_Value.MSA_DB_Path,scripts.Global_Value.MSA_DB_Name)

        if scripts.Global_Value.Mode=='whole' and scripts.Global_Value.Is_Use_Reverse_Data:
            print('Adding reverse task')
            Add_Reverse_Data(Table_Path,Res_Table_Name)

        if scripts.Global_Value.Mode=='whole':
            Feature_Object_List = []

            print('Beginning features extraction')
            Feature_Extraction(Table_Path,Res_Table_Name,Feature_Object_List,scripts.Global_Value.Process_Num)

            print('Recording features results')
            if not Record_Feature_Table(Feature_Object_List,Features_Table_Path):
                error_obj.Something_Wrong(__name__)
                exit(1)

    except:
        error_obj.Something_Wrong(__name__)
        Clean_with_Error(Docker_Container_Name)
        exit(1)

    if scripts.Global_Value.D_or_S=='D':
        print('Removing Docker container')
        Docker_Remove_Container(Docker_Container_Name)

    Remove_FoldX_Resource()
    Clean_Main_Directory()

    print('Have finished')
    exit(0)













