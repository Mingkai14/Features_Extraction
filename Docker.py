import os

def Docker_Init_Container(container_name,image_id):
    os.system(f'docker rm -f {container_name} && exit 0')
    os.system(f'docker run -d -it --name {container_name} {image_id} /bin/bash && exit 0')

def Docker_Make_Dir(container_name,path):
    os.system(f'docker exec {container_name} mkdir {path}')

def Docker_Delete_Dir(container_name,path):
    os.system(f'docker exec {container_name} rm -rf {path}')

def Docker_Import_File(container_name,src_path,out_path):
    os.system(f'docker cp {src_path} {container_name}:{out_path}')

def Docker_Export_File(container_name,src_path,out_path):
    os.system(f'docker cp {container_name}:{src_path} {out_path}')

def Docker_Remove_Container(container_name):
    os.system(f'docker rm -f {container_name}')

def Docker_Run_Cmd(container_name,cmd):
    os.system(f'docker exec {container_name} {cmd}')