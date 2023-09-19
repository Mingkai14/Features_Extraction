import os
from scripts.Error import error_obj
from scripts.Global_Value import *
from scripts.Utils import Amino_Acid_Encode,SS_Encode
import scripts.AAindex

def Init():
    Init_for_SIFT()
    Amino_Acid_Encode()
    Init_for_FoldX()
    SS_Encode()
    scripts.AAindex.Init_AAindex(AAIndex1_Path,AAIndex2_Path,AAIndex3_Path)

def Init_for_SIFT():
    if not os.path.isdir(SIFT_Path):
        error_obj.Something_Wrong(Init_for_SIFT.__name__)
        return False
    backup_lines=[]
    with open(SIFT_Path+'bin/'+'SIFT_for_submitting_fasta_seq.csh','r') as config:
        lines=config.readlines()
        for line in lines:
            if line.find('setenv')!=-1 and line.find('NCBI')!=-1:
                abs_blast_path=BLAST_Path.replace('./',Main_Location)
                line=f'setenv NCBI {abs_blast_path}\n'
            if line.find('setenv')!=-1 and line.find('SIFT_DIR')!=-1 and line.find('tmpdir')==-1 and line.find('BLIMPS_DIR')==-1:
                abs_SIFT_path=SIFT_Path.replace('./',Main_Location)
                line=f'setenv SIFT_DIR {abs_SIFT_path}\n'
            backup_lines.append(line)
    with open(SIFT_Path+'bin/'+'SIFT_for_submitting_fasta_seq.csh','w') as w_config:
        for line in backup_lines:
            w_config.write(line)




def Init_for_FoldX():
    if not os.path.exists('./molecules/'):
        import shutil
        shutil.copytree(f'{FoldX_Path}molecules/','./molecules/')












