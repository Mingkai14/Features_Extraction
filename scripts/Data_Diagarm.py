'''

1. Raw Data(pdb_num, which, what to what) -> WT_Amino_Acid_short, MUT_Amino_Acid_short, Loc_of_Mutation
2. Raw Data -> WT_Structure(pdb)
3. WT_Structure -> WT_Sequence_path(fasta)
4. WT_Structure -> Difference_Gap
5. WT_Sequence_path,Difference_Gap -> MUT_Sequence_path
6. MUT_Sequence_path,WT_Structure,Raw Data -> MUT_Structure
7. WT_Structure,Difference_Gap,Loc_of_Mutation -> WT_Amino_Acid(property)
8. MUT_Structure,Difference_Gap,Loc_of_Mutation -> MUT_Amino_Acid
9. WT_Structure,Difference_Gap -> WT_Amino_Acid_List
10. WT_Amino_Acid_List -> WT_Seq
11. MUT_Amino_Acid_List -> MUT_Seq

'''
