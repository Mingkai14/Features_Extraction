
class Ring_Bond:
    def __init__(self):
        self.Type=''
        self.Chain_ID=''
        self.AA_1_Num=0
        self.AA_1=''
        self.AA_2_Num=0
        self.AA_2=''
        self.Atom_1=''
        self.Atom_2=''
        self.Distance=''
        self.Angle=''
        self.Energy=''

class Protlego_Hydrophobic_Cluster:
    def __init__(self):
        self.Residues_List=[]
        self.Area=-99.99
        self.Contacts=0
        self.Ratio_contacts_residue=-99.99
        self.Ratio_area_residue=-99.99
        self.Each_Chain_ID=[]


class PDB:
    def __init__(self):
        self.PDB_Name=''
        self.PDB_path=''

class Researched_Amino_Acid:
    def __init__(self):
        self.Chain_ID=''
        self.Num=0
        self.Type_short=''
        self.Type=''
        self.Atom_List=[]
        self.Central_X=-99
        self.Central_Y = -99
        self.Central_Z = -99
class Atom:
    def __init__(self):
        self.Atom_Name=''
        self.Atom_Full_Name=''
        self.Element=''
        self.X=-99
        self.Y=-99
        self.Z=-99


class Feature_Object:
    def __init__(self):
    #基本信息
        self.ID=''
        self.WT_Structure=PDB()
        self.MUT_Structure = PDB()
        self.WT_Amino_Acid_short=''
        self.MUT_Amino_Acid_short=''
        self.Loc_of_Mutation=0
        self.True_Loc_of_Mutation=0
        self.WT_Sequence_path=''
        self.MUT_Sequence_path=''
        self.WT_Amino_Acid=Researched_Amino_Acid()
        self.MUT_Amino_Acid=Researched_Amino_Acid()
        self.WT_Amino_Acid_List=[]
        self.MUT_Amino_Acid_List=[]
        self.WT_Amino_Acid_List_Layer1 = []
        self.WT_Amino_Acid_List_Layer2 = []
        self.WT_Amino_Acid_List_Layer3 = []
        self.MUT_Amino_Acid_List_Layer1 = []
        self.MUT_Amino_Acid_List_Layer2 = []
        self.MUT_Amino_Acid_List_Layer3 = []
        self.Pfam_Num=''
        self.Cutoff1=13.0
        self.Cutoff2 = 10.0
        self.Cutoff3 = 7.0
        self.Is_Beta = '0'



        self.WT_Seq={}
        self.MUT_Seq={}
        self.Chain_ID_of_Mut=''
        self.WT_PSSM_Path=''
        self.WT_PSI_BLAST_Path=''
        self.MUT_PSSM_Path = ''
        self.MUT_PSI_BLAST_Path = ''
        self.WT_BLASTP_Path=''
        self.MUT_BLASTP_Path=''

        # RMSD of alignment, by pymol
        self.RMSD_WT_MUT = -99.99

    #1. Conditions Features
        #experimental condition
        self.pH = -99.99
        self.Temperature = -99.99

    #2. Overall Environment Features
        # Amino Acid Categories Infomation, by my script
        self.WT_Pct_Amino_Acid_Categories = {'uncharged_polar': -99.99, 'positively_charged_polar': -99.99,
                                                  'negatively_charged_polar': -99.99, 'nonpolar': -99.99,
                                                  'aromatic': -99.99, 'aliphatic': -99.99, 'heterocyclic': -99.99,
                                                  'sulfur_containing': -99.99}
        self.WT_Num_Amino_Acid_Categories = {'uncharged_polar': 0, 'positively_charged_polar': 0,
                                                  'negatively_charged_polar': 0, 'nonpolar': 0, 'aromatic': 0,
                                                  'aliphatic': 0, 'heterocyclic': 0, 'sulfur_containing': 0}

        self.WT_Pct_Amino_Acid_Categories_Layer1 = {'uncharged_polar': -99.99, 'positively_charged_polar': -99.99,
                                                   'negatively_charged_polar': -99.99, 'nonpolar': -99.99,
                                                   'aromatic': -99.99, 'aliphatic': -99.99, 'heterocyclic': -99.99,
                                                   'sulfur_containing': -99.99}
        self.WT_Num_Amino_Acid_Categories_Layer1 = {'uncharged_polar': 0, 'positively_charged_polar': 0,
                                                   'negatively_charged_polar': 0, 'nonpolar': 0, 'aromatic': 0,
                                                   'aliphatic': 0, 'heterocyclic': 0, 'sulfur_containing': 0}

        self.WT_Pct_Amino_Acid_Categories_Layer2 = {'uncharged_polar': -99.99, 'positively_charged_polar': -99.99,
                                                   'negatively_charged_polar': -99.99, 'nonpolar': -99.99,
                                                   'aromatic': -99.99, 'aliphatic': -99.99, 'heterocyclic': -99.99,
                                                   'sulfur_containing': -99.99}
        self.WT_Num_Amino_Acid_Categories_Layer2 = {'uncharged_polar': 0, 'positively_charged_polar': 0,
                                                    'negatively_charged_polar': 0, 'nonpolar': 0, 'aromatic': 0,
                                                    'aliphatic': 0, 'heterocyclic': 0, 'sulfur_containing': 0}

        self.WT_Pct_Amino_Acid_Categories_Layer3 = {'uncharged_polar': -99.99, 'positively_charged_polar': -99.99,
                                                    'negatively_charged_polar': -99.99, 'nonpolar': -99.99,
                                                    'aromatic': -99.99, 'aliphatic': -99.99, 'heterocyclic': -99.99,
                                                    'sulfur_containing': -99.99}
        self.WT_Num_Amino_Acid_Categories_Layer3 = {'uncharged_polar': 0, 'positively_charged_polar': 0,
                                                    'negatively_charged_polar': 0, 'nonpolar': 0, 'aromatic': 0,
                                                    'aliphatic': 0, 'heterocyclic': 0, 'sulfur_containing': 0}


        # Buried_Residue, Exposed_Residue, Secondary_Structure Information, by DSSP
        self.WT_Pct_Buried_Residue = -99.99
        self.WT_Pct_Exposed_Residue = -99.99

        self.WT_Pct_Buried_Residue_Layer1 = -99.99
        self.WT_Pct_Exposed_Residue_Layer1 = -99.99

        self.WT_Pct_Buried_Residue_Layer2 = -99.99
        self.WT_Pct_Exposed_Residue_Layer2 = -99.99

        self.WT_Pct_Buried_Residue_Layer3 = -99.99
        self.WT_Pct_Exposed_Residue_Layer3 = -99.99


        self.Dssp_List = []
        self.WT_Psipred_List = []
        self.MUT_Psipred_List = []
        # self.Overall_Pct_Secondary_Structure = {'H': -99.99, 'E': -99.99, 'C': -99.99}

        self.WT_Pct_Secondary_Structure = {'H': -99.99, 'B': -99.99, 'E': -99.99, 'G': -99.99, 'I': -99.99,
                                                'T': -99.99, 'S': -99.99, '-': -99.99}

        self.WT_Pct_Secondary_Structure_Layer1 = {'H': -99.99, 'B': -99.99, 'E': -99.99, 'G': -99.99, 'I': -99.99,
                                                'T': -99.99, 'S': -99.99, '-': -99.99}

        self.WT_Pct_Secondary_Structure_Layer2 = {'H': -99.99, 'B': -99.99, 'E': -99.99, 'G': -99.99, 'I': -99.99,
                                                'T': -99.99, 'S': -99.99, '-': -99.99}

        self.WT_Pct_Secondary_Structure_Layer3 = {'H': -99.99, 'B': -99.99, 'E': -99.99, 'G': -99.99, 'I': -99.99,
                                                'T': -99.99, 'S': -99.99, '-': -99.99}

        # Coils, Rem465, and Hotloop Information, by disEMBL
        self.COILS_line = ''
        self.REM465_line = ''
        self.HOTLOOPS_line = ''

        self.WT_Pct_coils = -99.99
        self.WT_Whole_Length_coils = -99.99
        self.WT_Pct_rem465 = -99.99
        self.WT_Whole_Length_rem465 = -99.99
        self.WT_Pct_hotloop = -99.99
        self.WT_Whole_Length_hotloop = -99.99


    #3. Environmental Changes Features
        # NMA, by Bio3D
        self.WT_NMA_Fluctuation = -99.99
        self.MUT_NMA_Fluctuation = -99.99
        self.Diff_NMA_Fluctuation = -99.99
        self.Overall_Rmsip = -99.99

        #FoldX energy terms, by FoldX
        self.WT_FoldX_Energy_Term_Dict = {'Pdb': '', 'total energy': -99.99, 'Backbone Hbond': -99.99,
                                          'Sidechain Hbond': -99.99, 'Van der Waals': -99.99, 'Electrostatics': -99.99,
                                          'Solvation Polar': -99.99, 'Solvation Hydrophobic': -99.99,
                                          'Van der Waals clashes': -99.99, 'entropy sidechain': -99.99,
                                          'entropy mainchain': -99.99, 'sloop_entropy': -99.99, 'mloop_entropy': -99.99,
                                          'cis_bond': -99.99, 'torsional clash': -99.99, 'backbone clash': -99.99,
                                          'helix dipole': -99.99, 'water bridge': -99.99, 'disulfide': -99.99,
                                          'electrostatic kon': -99.99, 'partial covalent bonds': -99.99,
                                          'energy Ionisation': -99.99, 'Entropy Complex': -99.99}
        self.Diff_FoldX_Energy_Term_Dict = {'Pdb': '', 'total energy': -99.99, 'Backbone Hbond': -99.99,
                                            'Sidechain Hbond': -99.99, 'Van der Waals': -99.99, 'Electrostatics': -99.99,
                                            'Solvation Polar': -99.99, 'Solvation Hydrophobic': -99.99,
                                            'Van der Waals clashes': -99.99, 'entropy sidechain': -99.99,
                                            'entropy mainchain': -99.99, 'sloop_entropy': -99.99, 'mloop_entropy': -99.99,
                                            'cis_bond': -99.99, 'torsional clash': -99.99, 'backbone clash': -99.99,
                                            'helix dipole': -99.99, 'water bridge': -99.99, 'disulfide': -99.99,
                                            'electrostatic kon': -99.99, 'partial covalent bonds': -99.99,
                                            'energy Ionisation': -99.99, 'Entropy Complex': -99.99}




        #Bond Infomation, by Ring3
        #Whole protein
        self.WT_Ring_Bond_List=[]
        self.WT_Num_HBOND_Ring = -99
        self.WT_Num_SSBOND_Ring=-99
        self.WT_Num_IONIC_Ring=-99
        self.WT_Num_VDW_Ring=-99
        self.WT_Num_PICATION_Ring=-99
        self.WT_Num_PIPISTACK_Ring=-99

        self.MUT_Ring_Bond_List = []
        self.MUT_Num_HBOND_Ring = -99
        self.MUT_Num_SSBOND_Ring = -99
        self.MUT_Num_IONIC_Ring = -99
        self.MUT_Num_VDW_Ring = -99
        self.MUT_Num_PICATION_Ring = -99
        self.MUT_Num_PIPISTACK_Ring = -99

        self.Diff_Num_HBOND_Ring = -99
        self.Diff_Num_SSBOND_Ring = -99
        self.Diff_Num_IONIC_Ring = -99
        self.Diff_Num_VDW_Ring = -99
        self.Diff_Num_PICATION_Ring = -99
        self.Diff_Num_PIPISTACK_Ring = -99

        # Layer 1
        self.WT_Num_HBOND_Ring_Layer1 = -99
        self.WT_Num_SSBOND_Ring_Layer1 = -99
        self.WT_Num_IONIC_Ring_Layer1 = -99
        self.WT_Num_VDW_Ring_Layer1 = -99
        self.WT_Num_PICATION_Ring_Layer1 = -99
        self.WT_Num_PIPISTACK_Ring_Layer1 = -99

        self.MUT_Num_HBOND_Ring_Layer1 = -99
        self.MUT_Num_SSBOND_Ring_Layer1 = -99
        self.MUT_Num_IONIC_Ring_Layer1 = -99
        self.MUT_Num_VDW_Ring_Layer1 = -99
        self.MUT_Num_PICATION_Ring_Layer1 = -99
        self.MUT_Num_PIPISTACK_Ring_Layer1 = -99

        self.Diff_Num_HBOND_Ring_Layer1 = -99
        self.Diff_Num_SSBOND_Ring_Layer1 = -99
        self.Diff_Num_IONIC_Ring_Layer1 = -99
        self.Diff_Num_VDW_Ring_Layer1 = -99
        self.Diff_Num_PICATION_Ring_Layer1 = -99
        self.Diff_Num_PIPISTACK_Ring_Layer1 = -99

        # Layer 2
        self.WT_Num_HBOND_Ring_Layer2 = -99
        self.WT_Num_SSBOND_Ring_Layer2 = -99
        self.WT_Num_IONIC_Ring_Layer2 = -99
        self.WT_Num_VDW_Ring_Layer2 = -99
        self.WT_Num_PICATION_Ring_Layer2 = -99
        self.WT_Num_PIPISTACK_Ring_Layer2 = -99

        self.MUT_Num_HBOND_Ring_Layer2 = -99
        self.MUT_Num_SSBOND_Ring_Layer2 = -99
        self.MUT_Num_IONIC_Ring_Layer2 = -99
        self.MUT_Num_VDW_Ring_Layer2 = -99
        self.MUT_Num_PICATION_Ring_Layer2 = -99
        self.MUT_Num_PIPISTACK_Ring_Layer2 = -99

        self.Diff_Num_HBOND_Ring_Layer2 = -99
        self.Diff_Num_SSBOND_Ring_Layer2 = -99
        self.Diff_Num_IONIC_Ring_Layer2 = -99
        self.Diff_Num_VDW_Ring_Layer2 = -99
        self.Diff_Num_PICATION_Ring_Layer2 = -99
        self.Diff_Num_PIPISTACK_Ring_Layer2 = -99

        # Layer 3
        self.WT_Num_HBOND_Ring_Layer3 = -99
        self.WT_Num_SSBOND_Ring_Layer3 = -99
        self.WT_Num_IONIC_Ring_Layer3 = -99
        self.WT_Num_VDW_Ring_Layer3 = -99
        self.WT_Num_PICATION_Ring_Layer3 = -99
        self.WT_Num_PIPISTACK_Ring_Layer3 = -99

        self.MUT_Num_HBOND_Ring_Layer3 = -99
        self.MUT_Num_SSBOND_Ring_Layer3 = -99
        self.MUT_Num_IONIC_Ring_Layer3 = -99
        self.MUT_Num_VDW_Ring_Layer3 = -99
        self.MUT_Num_PICATION_Ring_Layer3 = -99
        self.MUT_Num_PIPISTACK_Ring_Layer3 = -99

        self.Diff_Num_HBOND_Ring_Layer3 = -99
        self.Diff_Num_SSBOND_Ring_Layer3 = -99
        self.Diff_Num_IONIC_Ring_Layer3 = -99
        self.Diff_Num_VDW_Ring_Layer3 = -99
        self.Diff_Num_PICATION_Ring_Layer3 = -99
        self.Diff_Num_PIPISTACK_Ring_Layer3 = -99



        # Hydrophobic cluster, by Protlego
        self.WT_HD_Cluster_List = []
        self.WT_Num_HD_Cluster_Protlego = -99
        self.WT_Num_HD_Cluster_Protlego_Layer1 = -99
        self.WT_Num_HD_Cluster_Protlego_Layer2 = -99
        self.WT_Num_HD_Cluster_Protlego_Layer3 = -99
        self.WT_Max_HD_Cluster_Area=-99.99

        self.MUT_HD_Cluster_List = []
        self.MUT_Num_HD_Cluster_Protlego = -99
        self.MUT_Num_HD_Cluster_Protlego_Layer1 = -99
        self.MUT_Num_HD_Cluster_Protlego_Layer2 = -99
        self.MUT_Num_HD_Cluster_Protlego_Layer3 = -99
        self.MUT_Max_HD_Cluster_Area = -99.99

        self.Diff_Num_HD_Cluster_Protlego = -99
        self.Diff_Num_HD_Cluster_Protlego_Layer1 = -99
        self.Diff_Num_HD_Cluster_Protlego_Layer2 = -99
        self.Diff_Num_HD_Cluster_Protlego_Layer3 = -99
        self.Diff_Max_HD_Cluster_Area = -99.99





        #Pharmacophore Information, by Rdkit
        #Whole protein
        self.WT_Num_Pharmacophore_Categories = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.MUT_Num_Pharmacophore_Categories = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                 'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.Diff_Num_Pharmacophore_Categories = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                  'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}

        #Layer 1
        self.WT_Num_Pharmacophore_Categories_Layer1 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.MUT_Num_Pharmacophore_Categories_Layer1 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                 'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.Diff_Num_Pharmacophore_Categories_Layer1 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                  'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        #Layer 2
        self.WT_Num_Pharmacophore_Categories_Layer2 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                       'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.MUT_Num_Pharmacophore_Categories_Layer2 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                        'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.Diff_Num_Pharmacophore_Categories_Layer2 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                         'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}

        # Layer 3
        self.WT_Num_Pharmacophore_Categories_Layer3 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                       'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.MUT_Num_Pharmacophore_Categories_Layer3 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                        'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}
        self.Diff_Num_Pharmacophore_Categories_Layer3 = {'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0,
                                                         'aromatics': 0, 'aliphatics': 0, 'hydrophobics': 0, 'ring': 0}




    #4. Desciption Features of Amino Acid Mutation
        #If join in bond, by Ring3
        self.Is_WT_HBOND=0
        self.Is_WT_SSBOND=0
        self.Is_WT_IONIC=0
        self.Is_WT_VDW=0
        self.Is_WT_PICATION=0
        self.Is_WT_PIPISTACK=0

        self.Is_MUT_HBOND = 0
        self.Is_MUT_SSBOND = 0
        self.Is_MUT_IONIC = 0
        self.Is_MUT_VDW = 0
        self.Is_MUT_PICATION = 0
        self.Is_MUT_PIPISTACK = 0

        self.Descri_HBOND=-99
        self.Descri_SSBOND = -99
        self.Descri_IONIC = -99
        self.Descri_VDW = -99
        self.Descri_PICATION = -99
        self.Descri_PIPISTACK = -99




        #If join in HD Cluster, by Protlego
        self.Is_WT_HD_Cluster=0
        self.WT_HD_Cluster_Area=0.0
        self.Is_MUT_HD_Cluster=0
        self.MUT_HD_Cluster_Area=0.0
        self.Diff_HD_Cluster_Area=0.0

        self.Descri_HD_Cluster = -99




        #B-factor, by Prof
        self.WT_B_Factor=-99.99
        self.MUT_B_Factor=-99.99
        self.Diff_B_Factor=-99.99


        #If harmful mutation, by SIFT
        self.SIFT_Score=0



        # Backbone_Torsional_Angle, by DSSP
        self.WT_Psi = -99.99
        self.WT_Phi = -99.99
        self.MUT_Psi = -99.99
        self.MUT_Phi = -99.99
        self.Diff_Psi = -99.99
        self.Diff_Phi = -99.99



        #RSA,SASA,Buried,Exposed,Secondary Structure by DSSP/Psipred/PHD
        self.WT_RSA = -99
        self.MUT_RSA = -99
        self.Diff_RSA = -99

        self.WT_Is_Buried_or_Exposed = -99
        self.MUT_Is_Buried_or_Exposed = -99
        self.Descri_Buried_or_Exposed=-99

        self.WT_Secondary_Structure = -99
        self.MUT_Secondary_Structure = -99
        self.WT_Secondary_Structure_Char = ''
        self.MUT_Secondary_Structure_Char = ''
        self.Descri_SS=-99


        #Mutation Substitute Information, by my script #Not yet
        self.WT_AA_Type=-99
        self.MUT_AA_Type = -99
        self.Descri_AA = -99


        self.Is_WT_Uncharged_Polar=0
        self.Is_WT_Positively_Charged_Polar=0
        self.Is_WT_Negatively_Charged_Polar=0
        self.Is_WT_Nonpolar=0
        self.Is_WT_Aliphatic=0
        self.Is_WT_Aromatic=0
        self.Is_WT_Heterocyclic=0
        self.Is_WT_Sulfur_Containing=0

        self.Is_MUT_Uncharged_Polar = 0
        self.Is_MUT_Positively_Charged_Polar = 0
        self.Is_MUT_Negatively_Charged_Polar = 0
        self.Is_MUT_Nonpolar = 0
        self.Is_MUT_Aliphatic = 0
        self.Is_MUT_Aromatic = 0
        self.Is_MUT_Heterocyclic = 0
        self.Is_MUT_Sulfur_Containing = 0

        self.Descri_Uncharged_Polar = 0
        self.Descri_Positively_Charged_Polar = 0
        self.Descri_Negatively_Charged_Polar = 0
        self.Descri_Nonpolar = 0
        self.Descri_Aliphatic = 0
        self.Descri_Aromatic = 0
        self.Descri_Heterocyclic = 0
        self.Descri_Sulfur_Containing = 0

        #PSSM matrix score, by blast
        self.WT_PSSM_Score=-99
        self.MUT_PSSM_Score=-99


        self.WT_PSSM_Score_F1=-99
        self.WT_PSSM_Score_F2 = -99
        self.WT_PSSM_Score_F3 = -99
        self.WT_PSSM_Score_F4 = -99
        self.WT_PSSM_Score_F5 = -99
        self.WT_PSSM_Score_B1 = -99
        self.WT_PSSM_Score_B2 = -99
        self.WT_PSSM_Score_B3 = -99
        self.WT_PSSM_Score_B4 = -99
        self.WT_PSSM_Score_B5 = -99
        self.WT_PSSM_Score_Aver = -99.99

        self.MUT_PSSM_Score_F1 = -99
        self.MUT_PSSM_Score_F2 = -99
        self.MUT_PSSM_Score_F3 = -99
        self.MUT_PSSM_Score_F4 = -99
        self.MUT_PSSM_Score_F5 = -99
        self.MUT_PSSM_Score_B1 = -99
        self.MUT_PSSM_Score_B2 = -99
        self.MUT_PSSM_Score_B3 = -99
        self.MUT_PSSM_Score_B4 = -99
        self.MUT_PSSM_Score_B5 = -99
        self.MUT_PSSM_Score_Aver = -99.99

        self.Diff_PSSM_Score = -99
        self.Diff_PSSM_Score_Aver= -99.99

        #AA Physicochemical properties, by aaindex
        self.WT_AAindex1={}
        self.Diff_AAindex1={}
        self.Descri_AAindex2 = {}
        self.Descri_AAindex3 = {}

    #Labels
        #Regression label #Not yet
        self.Experimental_DDG=-99.99
        #Three classification label #Not yet
        self.Experimental_DDG_Classification=-99

