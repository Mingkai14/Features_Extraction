import os
from scripts.Error import error_obj
from scripts.Global_Value import *
from scripts.Utils import Amino_Acid_Encode,SS_Encode
import scripts.AAindex

def Init():
    Init_for_SIFT()
    # Init_for_Psipred()
    # Init_for_ANGLOR()
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

def Init_for_Psipred():
    if not os.path.isdir(Psipred_Path):
        error_obj.Something_Wrong(Init_for_SIFT.__name__)
        return False
    backup_lines = []
    with open(Psipred_Path+'runpsipred','r') as config:
        lines=config.readlines()
        for line in lines:
            if line.find('set')!=-1 and line.find('dbname')!=-1:
                abs_db_path=NR_Filter_Path.replace('./',Main_Location)
                line=f'set dbname = {abs_db_path}\n'
            if line.find('set')!=-1 and line.find('ncbidir')!=-1:
                abs_ncbi_path=Psipred_Path.replace('./',Main_Location)+'/blast-2.2.26/bin'
                line=f'set ncbidir = {abs_ncbi_path}\n'
            if line.find('set')!=-1 and line.find('execdir')!=-1:
                abs_execdir_path=Psipred_Path.replace('./',Main_Location)+'bin'
                line=f'set execdir = {abs_execdir_path}\n'
            if line.find('set')!=-1 and line.find('datadir')!=-1:
                abs_datadir_path=Psipred_Path.replace('./',Main_Location)+'data'
                line=f'set datadir = {abs_datadir_path}\n'
            backup_lines.append(line)
    with open(Psipred_Path+'runpsipred','w') as w_config:
        for line in backup_lines:
            w_config.write(line)

def Init_for_ANGLOR():
    if not os.path.isdir(ANGLOR_Path):
        error_obj.Something_Wrong(Init_for_SIFT.__name__)
        return False
    backup_lines = []
    with open(ANGLOR_Path+'library/bin/psipred24/runpsipred','r') as config:
        lines=config.readlines()
        for line in lines:
            if line.find('set')!=-1 and line.find('dbname')!=-1 and line.find('#')==-1:
                abs_db_path=NR_Filter_Path.replace('./',Main_Location)
                line=f'set dbname = {abs_db_path}\n'
            if line.find('set')!=-1 and line.find('ncbidir')!=-1 and line.find('#')==-1:
                abs_ncbi_path=ANGLOR_Path.replace('./',Main_Location)+'library/bin/blast/bin'
                line=f'set ncbidir = {abs_ncbi_path}\n'
            if line.find('set')!=-1 and line.find('execdir')!=-1 and line.find('#')==-1:
                abs_execdir_path=ANGLOR_Path.replace('./',Main_Location)+'library/bin/psipred24/bin'
                line=f'set execdir = {abs_execdir_path}\n'
            if line.find('set')!=-1 and line.find('datadir')!=-1 and line.find('#')==-1:
                abs_datadir_path=ANGLOR_Path.replace('./',Main_Location)+'library/bin/psipred24/data'
                line=f'set datadir = {abs_datadir_path}\n'
            backup_lines.append(line)
    with open(ANGLOR_Path+'library/bin/psipred24/runpsipred','w') as w_config:
        for line in backup_lines:
            w_config.write(line)


    backup_lines1=[]
    with open(ANGLOR_Path+'ANGLOR/ANGLOR.pl','r') as config1:
        lines = config1.readlines()
        for line in lines:
            if line.find('$datadir=')!=-1:
                abs_path=ANGLOR_Path.replace('./',Main_Location)+'example/$pdb'
                line=f'$datadir="{abs_path}";\n'
            if line.find('$libdir=')!=-1:
                abs_path=ANGLOR_Path.replace('./',Main_Location)+'library'
                line=f'$libdir="{abs_path}";\n'
            if line.find('$db=')!=-1:
                abs_path=NR_Filter_Path.replace('./',Main_Location)
                line=f'$db="{abs_path}";\n'
            if line.find('$work_dir=')!=-1:
                abs_path=Main_Location+'tmp/'
                line=f'$work_dir="{abs_path}";\n'
            backup_lines1.append(line)
    with open(ANGLOR_Path + 'ANGLOR/ANGLOR.pl', 'w') as w_config1:
        for line in backup_lines1:
            w_config1.write(line)

    backup_lines2 = []
    with open(ANGLOR_Path + 'library/bin/exp.pl', 'r') as config2:
        lines = config2.readlines()
        for line in lines:
            if line.find('system') != -1 and line.find('blastpgp') != -1:
                abs_path = NR_Filter_Path.replace('./', Main_Location)
                line = f'system("$librarydir/bin/blast/bin/blastpgp -i protein\.seq -Q protein\.mat3 -d {abs_path} -o protein\.out -e 1e-3 -h 1e-3 -M BLOSUM62 -j 3 -m 6");\n'
            backup_lines2.append(line)
    with open(ANGLOR_Path + 'library/bin/exp.pl', 'w') as w_config2:
        for line in backup_lines2:
            w_config2.write(line)



def Init_for_FoldX():
    if not os.path.exists('./molecules/'):
        import shutil
        shutil.copytree(f'{FoldX_Path}molecules/','./molecules/')












