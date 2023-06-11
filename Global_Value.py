import os


Is_Use_Reverse_Data=True

Psi_Threads_Num=4

Docker_Container_Name='myprof'
Docker_Image_ID='7d9fe6723898'

Raw_Dataset_file=''


Raw_PDB_Path='./src/Raw_PDB/'
WT_PDB_Path='./src/WT_PDB/'
MUT_PDB_Path='./src/Mut_PDB/'
WT_Fasta_Path='./src/WT_Fasta/'
MUT_Fasta_Path='./src/MUT_Fasta/'
WT_PSSM_Data_Path='./src/PSSM_Data/WT/'
MUT_PSSM_Data_Path='./src/PSSM_Data/MUT/'
WT_PSI_BLAST_Data_Path='./src/PSI_BLAST_Data/WT/'
MUT_PSI_BLAST_Data_Path='./src/PSI_BLAST_Data/MUT/'

Table_Path='./src/Data_Table/'
Res_Table_Name='data_table.txt'

HBPlus_Path='./bin/hbplus/'

Ring_Path='./bin/ring-3.0.0/ring/bin/'

FoldX_Path='./bin/FoldX_5.0/'

FoldX_Name='foldx_20231231'

Rdkit_Path='./bin/rdkit_2023_3_1/'

Rdkit_Fdef_Name='BaseFeatures.fdef'


Features_Table_Path='./src/Features_Table/'


MSA_DB_Path=''
MSA_DB_Name=''


Prof_Path='./bin/Prof/'
Prof_Temp_Path='./bin/Prof/Middle_Files/'

Main_Location=os.path.dirname(os.path.abspath(__file__))+'/'

R_NMA_Path='./bin/R_NMA/'
R_NMA_App_Name='NMA.R'

DisEMBL_Path='./bin/DisEMBL_1_4/'

BLAST_Path='./bin/ncbi_blast_2_13_0+/bin/'
Caps_Path='./bin/caps_2_0/'

WT_MSA_Path='./src/WT_MSA/'
SIFT_Path='./bin/sift6_2_1/'

Psipred_Path='./bin/psipred_v4/'


ANGLOR_Path='./bin/ANGLOR_source/'

NR_Filter_Path='./src/nr_filter/nr.filter'

Rosetta_Path='./bin/rosetta/'
Rosetta_Bin_Path='./bin/rosetta/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/'
Rosetta_DB_Path='./bin/rosetta/rosetta_bin_linux_2021.16.61629_bundle/main/database/'

Clean_Path='./bin/clean/'

AAIndex1_Path='./src/AAindex/aaindex1'
AAIndex2_Path='./src/AAindex/aaindex2'
AAIndex3_Path='./src/AAindex/aaindex3'
