a_l=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
a_d={'0':['A','R','N','D','C','Q','E','G','H','I'],'1':['L','K','M','F','P','S','T','W','Y','V']}

aaindex1_list=[]
aaindex2_list=[]
aaindex3_list=[]




def Extract_Matrix_AAindex(aaindex_path):
    #res_list=[{'':{'AA':-99.99,}}]
    #res_list=['':{'A':-99.99,}]
    res_list=[]
    parts_list=[]
    with open(aaindex_path,'r') as aaindex:
        lines=aaindex.readlines()
        temp_lines=[]
        for line in lines:
            line=line.replace('\n','')
            temp_lines.append(line)
            if line=='//':
                parts_list.append(temp_lines)
                temp_lines=[]
    for part in parts_list:
        is_data=False
        data_l=[]
        key=''
        for line in part:
            if is_data==True and line.find('//')==-1:
                data_l.append(line)
            else:
                div=line.split()
                if div[0]=='D':
                    key=' '.join(div[1:])
            if line.find('M')!=-1 and line.find('rows')!=-1 and line.find('ARNDCQEGHILKMFPSTWYV')!=-1:
                is_data=True
        if data_l!=[]:
            d=Read_Matrix(data_l,a_l)
            res_list.append({key:d})
    return res_list

def Read_Matrix(lines:list[str],aa_list:list):
    matrix=[]
    for i in range(20):
        temp=[]
        for j in range(20):
            temp.append(-99.99)
        matrix.append(temp)
    lines_matrix=[]
    for line in lines:
        line=line.replace('\n','')
        div=line.split()
        lines_matrix.append(div)
    jj=1
    for i in range(20):
        for j in range(jj):
            if lines_matrix[i][j]=='-' or lines_matrix[i][j]=='NA':
                lines_matrix[i][j]=0.0
            matrix[i][j]=float(lines_matrix[i][j])
            if j!=i:
                matrix[j][i]=-float(lines_matrix[i][j])
        jj+=1
    res_dict={}
    for i in range(len(aa_list)):
        for j in range(len(aa_list)):
            res_dict[f'{aa_list[i]}{aa_list[j]}']=[i,j]
    for key in res_dict.keys():
        row=res_dict[key][0]
        col=res_dict[key][1]
        res_dict[key]=matrix[row][col]
    return res_dict


def Extract_Index_AAindex(aaindex_path):
    res_list=[]
    parts_list = []
    with open(aaindex_path, 'r') as aaindex:
        lines = aaindex.readlines()
        temp_lines = []
        for line in lines:
            line = line.replace('\n', '')
            temp_lines.append(line)
            if line == '//':
                parts_list.append(temp_lines)
                temp_lines = []
    for part in parts_list:
        is_data=False
        data_l=[]
        key=''
        for line in part:
            if is_data==True and line.find('//')==-1:
                data_l.append(line)
            else:
                div=line.split()
                if div[0]=='D':
                    key=' '.join(div[1:])
            if div[0]=='I' and line.find('A/L')!=-1 and line.find('R/K')!=-1:
                is_data=True
        if data_l!=[]:
            d=Read_Index(data_l,a_d)
            res_list.append({key:d})
    return res_list


def Read_Index(lines:list[str],aa_dict:dict):
    res_dict={}
    line_0=lines[0]
    div=line_0.split()
    for i in range(len(div)):
        if div[i]=='-' or div[i]=='NA':
            div[i]=0.0
        res_dict[aa_dict['0'][i]]=float(div[i])
    line_1=lines[1]
    div = line_1.split()
    for i in range(len(div)):
        if div[i]=='-' or div[i]=='NA':
            div[i]=0.0
        res_dict[aa_dict['1'][i]] = float(div[i])
    return res_dict



def Init_AAindex(aaindex1,aaindex2,aaindex3):
    global aaindex1_list,aaindex2_list,aaindex3_list
    aaindex1_list = Extract_Index_AAindex(aaindex1)
    aaindex2_list = Extract_Matrix_AAindex(aaindex2)
    aaindex3_list = Extract_Matrix_AAindex(aaindex3)

def Get_Mutation_Index_List_from_Matrix(mut_info:str,aaindex_list:list[dict],out_dict:dict):
    #res_list=[{'':{'AA':-99.99,}}]
    for d in aaindex_list:
        name=list(d.keys())[0]
        aa_dict=d[name]
        data=aa_dict[mut_info]
        out_dict[name]=data
def Get_Mutation_Index_List_from_Index(wt_aa:str,mut_aa:str,aaindex_list:list[dict],out_dict:dict):
    #res_list=[{'':{'A':-99.99,}}]
    wt_dict={}
    for d in aaindex_list:
        name=list(d.keys())[0]
        aa_dict=d[name]
        data=aa_dict[wt_aa]
        wt_dict[name]=data
    mut_dict={}
    for d in aaindex_list:
        name=list(d.keys())[0]
        aa_dict=d[name]
        data=aa_dict[mut_aa]
        mut_dict[name]=data
    for name in wt_dict.keys():
        out_dict[name]=mut_dict[name]-wt_dict[name]
    return wt_dict






