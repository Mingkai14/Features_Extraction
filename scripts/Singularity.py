import os

def Singularity_Make_Dir_on_Home(sif_path,path_name):
    os.system(f'singularity exec {sif_path} mkdir ~/{path_name}/')

def Singularity_Delete_Dir_on_Home(sif_path,path_name):
    os.system(f'singularity exec {sif_path} rm -rf ~/{path_name}/')

def Singularity_Copy_File(src_path,out_path):
    os.system(f'cp {src_path} {out_path}')

def Singularity_Run_Cmd(sif_path,cmd):
    os.system(f'singularity exec {sif_path} {cmd}')