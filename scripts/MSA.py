import os
from scripts.Error import error_obj
from scripts.Utils import Fetch_Chain_ID_from_Seq,Read_Seq_from_Fasta,Fetch_Single_Chain_Loc
from scripts.Classes import Researched_Amino_Acid
from scripts.Global_Value import *
import scripts.Global_Value


data_list=[]
backup_lines=[]
def Prepare_Blast_Files(table_path,table_name,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,blast_db_path,blast_db_name):
    if os.path.exists(wt_pssm_path) == False:
        os.mkdir(wt_pssm_path)
    if os.path.exists(wt_psi_blast_path) == False:
        os.mkdir(wt_psi_blast_path)
    if os.path.exists(mut_pssm_path) == False:
        os.mkdir(mut_pssm_path)
    if os.path.exists(mut_psi_blast_path) == False:
        os.mkdir(mut_psi_blast_path)
    if os.path.exists(wt_blastp_path) == False:
        os.mkdir(wt_blastp_path)
    if os.path.exists(mut_blastp_path) == False:
        os.mkdir(mut_blastp_path)

    with open(table_path+table_name,'r') as table:
        lines=table.readlines()
        for line in lines:
            backup_lines.append(line)
        if lines[0]!='id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,pH,temperature,ddg\n':
            error_obj.Something_Wrong(Prepare_Blast_Files.__name__)
            exit(1)
        for line in lines[1:]:
            data_list.append(line.replace('\n',''))
    line_num = 1
    for data in data_list:
        item_list=str(data).split(',')
        if len(item_list)!=20:
            error_obj.Something_Wrong(Prepare_Blast_Files.__name__)
            exit(1)
        wt_fasta_path=item_list[9]
        mut_fasta_path=item_list[10]
        loc=int(item_list[3])
        t_loc=int(item_list[4])
        wt_aa=item_list[1]
        mut_aa=item_list[2]
        wt_pdb_name = item_list[5]
        mut_pdb_name = item_list[0]
        wt_seq_dict=Read_Seq_from_Fasta(wt_fasta_path)
        mut_seq_dict=Read_Seq_from_Fasta(mut_fasta_path)
        chain_id=Fetch_Chain_ID_from_Seq(t_loc,wt_seq_dict,wt_aa)

        wt_pssm_files = os.listdir(wt_pssm_path)
        mut_pssm_files = os.listdir(mut_pssm_path)
        wt_psi_blast_files = os.listdir(wt_psi_blast_path)
        mut_psi_blast_files = os.listdir(mut_psi_blast_path)
        wt_blastp_files = os.listdir(wt_blastp_path)
        mut_blastp_files = os.listdir(mut_blastp_path)
        wt_pssm_names = []
        mut_pssm_names = []
        wt_psi_blast_names = []
        mut_psi_blast_names = []
        wt_blastp_names = []
        mut_blastp_names = []

        for file in wt_pssm_files:
            wt_pssm_names.append(file.split('.')[0])
        for file in mut_pssm_files:
            mut_pssm_names.append(file.split('.')[0])
        for file in wt_psi_blast_files:
            wt_psi_blast_names.append(file.split('.')[0])
        for file in mut_psi_blast_files:
            mut_psi_blast_names.append(file.split('.')[0])
        for file in wt_blastp_files:
            wt_blastp_names.append(file.split('.')[0])
        for file in mut_blastp_files:
            mut_blastp_names.append(file.split('.')[0])

        with open('./temp.fasta','w') as temp_fasta:
            temp_fasta.write(f'>{wt_pdb_name}_{chain_id}\n')
            temp_fasta.write(wt_seq_dict[chain_id])

        if wt_pdb_name in wt_pssm_names and wt_pdb_name in wt_psi_blast_names and wt_pdb_name in wt_blastp_names:
            item_list[11]=wt_pssm_path+wt_pdb_name+'.pssm'
            item_list[13]=wt_psi_blast_path+wt_pdb_name+'.output'
            item_list[15]=wt_blastp_path+wt_pdb_name+'.output'
            if line_num!=len(backup_lines)-1:
                new_line=','.join(item_list)+'\n'
            else:
                new_line = ','.join(item_list)
            backup_lines[line_num]=new_line
        else:
            #run_blastpgp_2_6_0('./temp.fasta',blast_db_path,blast_db_name,wt_psi_blast_path,wt_pssm_path,wt_pdb_name)
            run_psiblast_2_13_0(BLAST_Path,'./temp.fasta',blast_db_path,blast_db_name,wt_psi_blast_path,wt_pssm_path,wt_pdb_name)
            item_list[11] = wt_pssm_path + wt_pdb_name + '.pssm'
            item_list[13] = wt_psi_blast_path + wt_pdb_name + '.output'
            run_blastp_2_13_0(BLAST_Path,'./temp.fasta',wt_blastp_path+wt_pdb_name+'.output',blast_db_path,blast_db_name,'6 sseqid sseq')
            item_list[15]=wt_blastp_path + wt_pdb_name + '.output'
            if line_num!=len(backup_lines)-1:
                new_line=','.join(item_list)+'\n'
            else:
                new_line = ','.join(item_list)
            backup_lines[line_num]=new_line
        os.remove('./temp.fasta')

        with open('./temp.fasta', 'w') as temp_fasta:
            temp_fasta.write(f'>{mut_pdb_name}_{chain_id}\n')
            temp_fasta.write(mut_seq_dict[chain_id])

        if mut_pdb_name in mut_pssm_names and mut_pdb_name in mut_psi_blast_names and mut_pdb_name in mut_blastp_names:
            item_list[12] = mut_pssm_path + mut_pdb_name + '.pssm'
            item_list[14] = mut_psi_blast_path + mut_pdb_name + '.output'
            item_list[16] = mut_blastp_path + mut_pdb_name + '.output'
            if line_num != len(backup_lines) - 1:
                new_line = ','.join(item_list) + '\n'
            else:
                new_line = ','.join(item_list)
            backup_lines[line_num] = new_line
        else:
            #run_blastpgp_2_6_0('./temp.fasta', blast_db_path, blast_db_name, mut_psi_blast_path, mut_pssm_path,mut_pdb_name)
            run_psiblast_2_13_0(BLAST_Path, './temp.fasta', blast_db_path, blast_db_name, mut_psi_blast_path, mut_pssm_path, mut_pdb_name)
            item_list[12] = mut_pssm_path + mut_pdb_name + '.pssm'
            item_list[14] = mut_psi_blast_path + mut_pdb_name + '.output'
            run_blastp_2_13_0(BLAST_Path,'./temp.fasta',mut_blastp_path+mut_pdb_name+'.output',blast_db_path,blast_db_name,'6 sseqid sseq')
            item_list[16] = mut_blastp_path + mut_pdb_name +'.output'
            if line_num != len(backup_lines) - 1:
                new_line = ','.join(item_list) + '\n'
            else:
                new_line = ','.join(item_list)
            backup_lines[line_num] = new_line
        os.remove('./temp.fasta')

        line_num+=1


    with open(table_path + table_name, 'w') as table:
        for line in backup_lines:
            table.write(line)
    return True


def run_blastpgp_2_6_0(fasta,db_path,db_name,psi_blast_path,pssm_path,pdb_name):
    os.system(f'blastpgp -i {fasta} -d {db_path+db_name} -j 3 -h 0.001 -o {psi_blast_path+pdb_name}.output -Q {pssm_path+pdb_name}.pssm')

def run_psiblast_2_13_0(blast_path,fasta,db_path,db_name,psi_blast_path,pssm_path,pdb_name):
    os.system(f'{blast_path}psiblast -query {fasta} -evalue .001 -inclusion_ethresh .002 -db {db_path+db_name} -num_iterations 3 -seg yes -out {psi_blast_path+pdb_name}.output -out_ascii_pssm {pssm_path+pdb_name}.pssm -outfmt 0 -num_threads {scripts.Global_Value.Psi_Threads_Num}')


pssm_map={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}


def find_pssm_score(pssm_path,aa_list:list[Researched_Amino_Acid],aa:Researched_Amino_Acid,seq_dict:dict,chain_id,windows=5):
    a_l=[]
    begin=aa.Num-windows
    end=aa.Num+windows
    #wrong

    res_list = []
    try:
        index_l=range(begin,end+1)

        for i in range(len(index_l)):
            if index_l[i]<1:
                res_list.append(0.0)

        for a in aa_list:
            if a.Num in index_l:
                a_l.append(a)


        for a in a_l:
            if a.Chain_ID!=chain_id:
                res_list.append(0.0)
            else:
                try:
                    data=read_pssm(pssm_path,a,seq_dict,chain_id)
                except:
                    res_list.append(0.0)
                    continue
                res_list.append(data)

        for i in range(len(index_l)):
            if index_l[i]>len(aa_list):
                res_list.append(0.0)


    except:
        error_obj.Something_Wrong(find_pssm_score.__name__)
        return False
    sum=0
    for data in res_list:
        sum+=data
    aver=sum/(2*windows+1)
    res_list.append(aver)
    return res_list

def read_pssm(pssm_path,aa:Researched_Amino_Acid,seq_dict:dict,chain):
    with open(pssm_path,'r') as pssm:
        loc=Fetch_Single_Chain_Loc(aa.Num,seq_dict,chain)
        lines=pssm.readlines()
        data=''
        for line in lines:
            line=line.replace('\n','')
            if line=='':
                continue
            div=line.split()
            try:
                int(div[0])
            except:
                continue
            if len(div)>=2 and int(div[0])==loc and div[1]==aa.Type_short:
                data=div[pssm_map[aa.Type_short]+2]
    return int(data)




def run_blastp_2_13_0(blast_path,fasta_path,out_path,db_path,db_name,blast_outfmt_info):
    os.system(f'{blast_path}blastp -query {fasta_path} -out {out_path} -db {db_path}{db_name} -outfmt \'{blast_outfmt_info}\' -evalue 1e-5 -num_threads {scripts.Global_Value.Psi_Threads_Num}')
#../ncbi-blast-2.13.0+/bin/blastp -query 1AAR_H68E.fasta -out res.fasta -db ~/blast_resource/uniref50 -outfmt '6 sseqid sseq'  -evalue 1e-5 -num_descriptions 100 -num_threads 4


