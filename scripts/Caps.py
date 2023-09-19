from scripts.MSA import run_blastp_2_13_0
import os
from scripts.Utils import Clean_Main_Directory,amino_acid_num_map,Fetch_Single_Chain_Loc
# from Run_Sift import Share_Aligned_File
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

















