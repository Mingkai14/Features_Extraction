from Feature_Extracting import *
from Run_Modeller import *
from MSA import *
from Record import Record_Feature_Table
from Global_Value import *
import Global_Value
from Init import Init
from Utils import *
import argparse
from Docker import Docker_Init_Container,Docker_Remove_Container


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--raw_dataset_path', type=str, default='')
    parser.add_argument('--db_folder_path', type=str, default='')
    parser.add_argument('--db_name', type=str, default='')
    parser.add_argument('--if_reversed_data', type=int, default=1)
    parser.add_argument('--psiblast_threads_num', type=int, default=4)
    parser.add_argument('--container_type', type=str, default='D')

    args = parser.parse_args()
    if args.raw_dataset_path=='' or args.db_folder_path=='' or args.db_name=='' or args.if_reversed_data not in [0,1] or args.psiblast_threads_num<1 or args.psiblast_threads_num>30:
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.raw_dataset_path).split('.')[len(str(args.raw_dataset_path).split('.'))-1]!='xls':
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.container_type) not in ['D','S']:
        error_obj.Something_Wrong(__name__)
        exit(1)

    Global_Value.Raw_Dataset_file=args.raw_dataset_path
    Global_Value.MSA_DB_Path=args.db_folder_path
    Global_Value.MSA_DB_Name=args.db_name
    Global_Value.Is_Use_Reverse_Data=args.if_reversed_data
    Global_Value.Psi_Threads_Num=args.psiblast_threads_num
    Global_Value.D_or_S = args.container_type

    if Global_Value.D_or_S=='D':
        Docker_Init_Container(Docker_Container_Name,Docker_Image_ID)

    try:
        Init()

        Raw_Data_List = Read_XLS(Global_Value.Raw_Dataset_file)

        Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,[1,1,1,0,0,1,1,0,0,0,0])

        Prepare_Table(Raw_Data_List,Table_Path,Res_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path)

        Prepare_MUT_Models(Table_Path,Res_Table_Name,MUT_PDB_Path)

        Prepare_Blast_Files(Table_Path,Res_Table_Name,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,Global_Value.MSA_DB_Path,Global_Value.MSA_DB_Name)

        if Global_Value.Is_Use_Reverse_Data:
            Add_Reverse_Data(Table_Path,Res_Table_Name)

        Feature_Object_List = []

        Feature_Extraction(Table_Path,Res_Table_Name,Feature_Object_List)

        if not Record_Feature_Table(Feature_Object_List,Features_Table_Path):
            error_obj.Something_Wrong(__name__)
            Clean_with_Error(Docker_Container_Name)
            exit(1)
    except:
        Clean_with_Error(Docker_Container_Name)
        exit(1)
    if Global_Value.D_or_S=='D':
        Docker_Remove_Container(Docker_Container_Name)
    exit(0)













