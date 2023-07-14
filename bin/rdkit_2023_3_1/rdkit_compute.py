from rdkit import Chem
from rdkit.Chem import ChemicalFeatures



def get_atom_coordinates(feature,protein):
    atom_indices = feature.GetAtomIds()
    atom_coordinates = [protein.GetConformer().GetAtomPosition(idx) for idx in atom_indices]
    return atom_coordinates

def Compute_Pharmacophore_with_Rdkit(pdb_path,rdkit_path,rdkit_fdef_name,central_x,central_y,central_z,cutoff):
    protein = Chem.MolFromPDBFile(pdb_path)
    if protein==None:
        return False
    try:
        fdefFile = rdkit_path+rdkit_fdef_name
        factory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
        features = factory.GetFeaturesForMol(protein)
        features_list=[]
    except:
        return False

    res_dict={'hb_acceptors': 0, 'hb_donors': 0, 'positives': 0, 'negatives': 0, 'aromatics': 0, 'aliphatics': 0,
     'hydrophobics': 0, 'ring': 0}

    if cutoff!=0.0:
        for feature in features:
            atom_c=get_atom_coordinates(feature,protein)
            x=atom_c[0][0]
            y=atom_c[0][1]
            z=atom_c[0][2]
            from scripts.Utils import Get_Distance
            dis=Get_Distance(central_x,central_y,central_z,x,y,z)
            if dis<cutoff:
                features_list.append(feature)
    else:
        features_list=features




    for feature in features_list:
        if feature.GetFamily() == 'Acceptor':
            res_dict['hb_acceptors'] += 1
        elif feature.GetFamily() == 'Donor':
            res_dict['hb_donors'] += 1
        elif feature.GetFamily() == 'Positive':
            res_dict['positives'] += 1
        elif feature.GetFamily() == 'Negative':
            res_dict['negatives'] += 1
        elif feature.GetFamily() == 'Aromatic':
            res_dict['aromatics'] += 1
        elif feature.GetFamily() == 'Aliphatic':
            res_dict['aliphatics'] += 1
        elif feature.GetFamily() == 'Hydrophobe':
            res_dict['hydrophobics'] += 1
        elif feature.GetFamily() == 'Ring':
            res_dict['ring'] += 1

    return res_dict

