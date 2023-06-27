from pymol import cmd
import os

def Pymol_Remove_Water(pdb_path,out_path,is_replace=True):
    cmd.load(pdb_path,'A.pdb')
    cmd.remove('resn hoh')
    if is_replace:
        if os.path.exists(out_path):
            os.remove(out_path)
    if not os.path.exists(out_path):
        cmd.save(out_path,'A.pdb')
    cmd.delete('all')

def Pymol_Clean_Align_PDB_Pair(pdb_path1,pdb_path2,out_path1,out_path2,is_replace=True):
    cmd.load(pdb_path1,'A.pdb')
    cmd.remove('resn hoh')
    cmd.load(pdb_path2, 'B.pdb')
    res=cmd.align('A.pdb', 'B.pdb')
    if is_replace:
        if os.path.exists(out_path1):
            os.remove(out_path1)
        if os.path.exists(out_path2):
            os.remove(out_path2)
    if not os.path.exists(out_path1):
        cmd.save(out_path1,'A.pdb')
    if not os.path.exists(out_path2):
        cmd.save(out_path2,'B.pdb')
    cmd.delete('all')
    return res
