import os
from Global_Value import WT_MSA_Path
from Utils import Clean_Main_Directory,Fetch_Single_Chain_Loc

def Share_Aligned_File(name,seq_dict:dict,chain_id,in_path,out_path=WT_MSA_Path):
    files=os.listdir(out_path)
    if name+'.aln.fasta' in files:
        return
    temp_lines=[]
    length=len(seq_dict[chain_id])
    temp_lines.append(f'>{name}_{chain_id}\n')
    temp_lines.append(f'{seq_dict[chain_id]}\n')
    with open(in_path,'r') as input:
        lines=input.readlines()
        for i in range(len(lines)):
            if lines[i].find('>')!=-1:
                if len(lines[i+1].replace('\n',''))==length:
                    temp_lines.append(lines[i])
                    temp_lines.append(lines[i+1])
    with open(out_path+name+'.aln.fasta','w') as output:
        for line in temp_lines:
            output.write(line)


def Make_Sub_file(wt_aa,mut_aa,loc:int):
    with open('temp.subst','w') as subset:
        subset.write(f'{wt_aa}{str(loc)}{mut_aa}')


def Clean_Files(files:list,path):
    for file in files:
        os.remove(path+file)

#bin/SIFT_for_submitting_fasta_seq.csh test/lacI.fasta <protein_database> test/lacI.subst
def Run_Sift_Useless(name,seq_dict:dict,chain_id,sift_path,wt_aa,mut_aa,loc:int,gap:int,db_path,db_name):
    Make_Sub_file(wt_aa,mut_aa,loc,gap)
    #os.system(f'export BLIMPS_DIR={bltmps_path} && {sift_bin_path}info_on_seqs {msa_path}{name}.aln.fasta ./temp.subst ./temp.SIFTprediction')
    sift_bin_path=sift_path+'bin/'
    with open('./temp.fasta','w') as fasta:
        fasta.write(f'>{name}_{chain_id}\n')
        fasta.write(f'{seq_dict[chain_id]}')
    os.system(f'{sift_bin_path}SIFT_for_submitting_fasta_seq.csh ./temp.fasta {db_path}{db_name} ./temp.subst')
    tmp_path=sift_path+'tmp/'
    files=os.listdir(tmp_path)
    predict_file=''
    for file in files:
        if file.split('.')[len(file.split('.'))-1]=='SIFTprediction':
            predict_file=file
    if predict_file=='':
        Clean_Files(files, tmp_path)
        Clean_Main_Directory()
        return False
    with open(tmp_path+predict_file,'r') as temp:
        line=temp.readlines()[0]
        div=line.split('\t')
        if div[1] not in ['TOLERATED','DELETERIOUS']:
            Clean_Files(files, tmp_path)
            Clean_Main_Directory()
            return False
        if div[1]=='TOLERATED':
            Clean_Files(files,tmp_path)
            Clean_Main_Directory()
            return 1
        elif div[1]=='DELETERIOUS':
            Clean_Files(files, tmp_path)
            Clean_Main_Directory()
            return -1
        else:
            Clean_Files(files, tmp_path)
            Clean_Main_Directory()
            return False

def Run_Sift(name,wt_aa,mut_aa,loc:int,sift_path,msa_path,seq_dict:dict,chain):
    loc_=Fetch_Single_Chain_Loc(loc,seq_dict,chain)
    Make_Sub_file(wt_aa, mut_aa, loc_)
    blimps_path=sift_path+'blimps/'
    sift_bin_path=sift_path+'bin/'
    os.system(f'export BLIMPS_DIR={blimps_path} && {sift_bin_path}info_on_seqs {msa_path}{name}.aln.fasta ./temp.subst ./temp.SIFTprediction')
    with open('./temp.SIFTprediction','r') as temp:
        lines=temp.readlines()
        line=''
        for l in lines:
            try:
                d=l.split('\t')
                if d[1] in ['TOLERATED','DELETERIOUS']:
                    line=l
                    break
            except:
                continue
        div=line.split('\t')
        if div[1] not in ['TOLERATED','DELETERIOUS']:
            Clean_Main_Directory()
            return False
        if div[1]=='TOLERATED':
            Clean_Main_Directory()
            return 1
        elif div[1]=='DELETERIOUS':
            Clean_Main_Directory()
            return -1
        else:
            Clean_Main_Directory()
            return False



