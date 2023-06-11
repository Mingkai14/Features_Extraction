from MSA import run_blastp_2_13_0
import os
from Utils import Clean_Main_Directory,amino_acid_num_map,Fetch_Single_Chain_Loc
from Run_Sift import Share_Aligned_File
def Trans_blast_2_fasta(blast_file,output_file,line_limit:int):
    count = 0
    with open(blast_file) as ori_blast:
        lines = ori_blast.readlines()
        with open(output_file, 'w') as new_blast:
            for line in lines:
                count += 1
                if count > line_limit:
                    break
                l = line.replace('\n', '')
                div = l.split('\t')
                new_blast.write(f'>{div[0]}\n')
                new_blast.write(f'{div[1]}\n')
    return output_file

def Run_Muscle(path,file,out_file):
    os.system(f'{path}muscle -in {file} -out {out_file}')
    return out_file
#./muscle -in ./blast.fasta -out blast_msa.fasta

def Run_Caps_2(caps_path,msa_fasta_path):
    with open(msa_fasta_path,'r') as fasta:
        temp=fasta.read()
        with open(caps_path+'sample/sample.fasta','w') as new_fasta:
            new_fasta.write(temp)
    os.system(f'{caps_path}caps -F {caps_path}sample/ --intra')
#./caps -F ./f/ --intra

pair_res=[]
group_res=[]
overlap_res=[]

def Identify_Blast_Out(name,seq_dict:dict,chain_id,in_path,out_path):
    length = len(seq_dict[chain_id])
    temp_lines = []
    with open(in_path,'r') as input:
        lines=input.readlines()
        for i in range(len(lines)):
            if lines[i].find('>')!=-1:
                if len(lines[i+1].replace('\n',''))==length:
                    temp_lines.append(lines[i])
                    temp_lines.append(lines[i+1])
    with open(out_path,'w') as output:
        for line in temp_lines:
            output.write(line)

def Read_Caps_2(path):
    files=os.listdir(path)
    outfile=''
    lines_cleaned=[]
    for file in files:
        if file.split('.')[len(file.split('.'))-1]=='out':
            outfile=file
    with open(path+outfile,'r') as out:
        lines=out.readlines()
        for line in lines:
            temp=line.replace('\n','').replace(' ','').replace('\t','')
            temp_count=0
            for c in temp:
                if c!='=' and c!='-':
                    temp_count+=1
            if temp_count!=0 or line=='\n':
                lines_cleaned.append(line)
    is_pair=False
    is_group=False
    is_overlap=False
    pair_list=[]
    group_list=[]
    overlap_list=[]
    for line in lines_cleaned:
        line=line.replace('\t\t','\t')
        if is_pair:
            if line!='\n' and line.find('Groups of coevolving Amino Acids')==-1:
                div=line.replace('\n','').split('\t')
                pair1=div[0].split('(')[0]
                pair2=div[1].split('(')[0]
                cor=div[4]
                pair_list.append([pair1,pair2,cor])
            if line=='\n' or line.find('Groups of coevolving Amino Acids')!=-1:
                is_pair=False
        if is_group:
            if line != '\n' and line.find('Overlapping groups of coevolving amino acids') == -1:
                line_ = line.split(':')[len(line.split(':')) - 1].replace(' ', '')
                div=line_.replace('\t\n','').split('\t')
                temp_l=[]
                for item in div:
                    if item!='':
                        temp_l.append(item.split('(')[0])
                group_list.append(temp_l)
            if line=='\n' or line.find('Overlapping groups of coevolving amino acids')!=-1:
                is_group=False
        if is_overlap:
            if line != '\n':
                line_=line.split(':')[len(line.split(':'))-1].replace(' ','')
                div = line_.replace('\t\n', '').split('\t')
                temp_l = []
                for item in div:
                    if item != '':
                        temp_l.append(item.split('(')[0])
                overlap_list.append(temp_l)
            if line=='\n':
                is_overlap=False
        if line.find('Position')!=-1 and line.find('AA1')!=-1 and line.find('Correlation')!=-1:
            is_pair=True
        if line.find('Groups of coevolving Amino Acids')!=-1:
            is_group=True
        if line.find('Overlapping groups of coevolving amino acids')!=-1:
            is_overlap=True
    global pair_res,group_res,overlap_res
    pair_res=pair_list
    group_res=group_list
    overlap_res=overlap_list

def Compute_Co_Evo(blast_path,seq_dict:dict,chain_id,name,db_path,db_name,caps_path,loc:int,seq_num:int):
    with open('./temp.fasta','w') as fasta:
        fasta.write(f'>{name}_{chain_id}\n')
        fasta.write(seq_dict[chain_id])
    run_blastp_2_13_0(blast_path,'./temp.fasta','./blast_out.blast',db_path,db_name,'6 sseqid sseq')
    Trans_blast_2_fasta('./blast_out.blast','./blast_out.fasta',seq_num)
    #Run_Muscle(caps_path,'./blast_out.fasta','./blast_out.aln.fasta')
    #Run_Caps_2(caps_path,'./blast_out.aln.fasta')
    Identify_Blast_Out(name,seq_dict,chain_id,'./blast_out.fasta','./blast_out.aln.fasta')
    Run_Caps_2(caps_path, './blast_out.aln.fasta')
    Read_Caps_2('./')
    loc_=loc
    loc_temp=Fetch_Single_Chain_Loc(loc_,seq_dict,chain_id)
    pair_l=[]
    for pair in pair_res:
        if str(loc_)==pair[0] or str(loc_)==pair[1]:
            pair_l.append(pair)
    max_v=-999999.9
    temp_l=[]
    aa_1=seq_dict[chain_id][loc_temp-1]
    loc__=0
    aa_2=''
    is_co=False
    for temp in pair_l:
        if float(temp[2])>max_v:
            max_v=float(temp[2])
            temp_l=temp
    if len(temp_l)!=0:
        is_co=True
        if temp_l[0]==str(loc_):
            loc__=int(temp_l[1])
            loc__temp=Fetch_Single_Chain_Loc(loc__,seq_dict,chain_id)
            aa_2 = seq_dict[chain_id][loc__temp - 1]
        else:
            loc__=int(temp_l[0])
            loc__temp = Fetch_Single_Chain_Loc(loc__, seq_dict, chain_id)
            aa_2 = seq_dict[chain_id][loc__temp - 1]
    is_in_group=False
    max_len=0
    temp_ll=[]
    for group in group_res:
        if str(loc_) in group:
            temp_ll.append(group)
    for group in overlap_res:
        if str(loc_) in group:
            temp_ll.append(group)
    if len(temp_ll)!=0:
        is_in_group=True
        for group in temp_ll:
            if len(group)>max_len:
                max_len=len(group)
    Share_Aligned_File(name,seq_dict,chain_id,'./blast_out.fasta')
    Clean_Main_Directory()
    res_l=[]
    if is_co:
        res_l.append(1)
    else:
        res_l.append(0)
    if aa_2=='':
        res_l.append(-1)
    else:
        res_l.append(amino_acid_num_map[aa_2])
    if is_in_group:
        res_l.append(1)
    else:
        res_l.append(0)
    res_l.append(max_len)
    return res_l





