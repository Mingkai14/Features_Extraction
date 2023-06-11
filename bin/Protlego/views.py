from .clusters import *
import urllib





def get_structures(pdb_path):
    """
    Computes hydrophobic clusters for specific query or pdb id
    """
    #Directly from URL

    mol = Molecule(pdb_path, validateElements=False)

    clusters = compute_hydrophobic_clusters(mol)
    dict_clusters = []
    for cluster in clusters:
        dict_clusters.append({
            "area":cluster.area,
            'residues': cluster.residues,
            'chains': cluster.chains,
            "contacts": cluster.contacts,
            "ratio_contacts_residue": cluster.ratio_contacts_residue,
            "ratio_area_residue": cluster.ratio_area_residue
        })
    return dict_clusters







def compute_hydrophobic_clusters(mol,
                                 sel: str = "protein and not backbone and noh and resname ILE VAL LEU",
                                 cutoff_area: float = 10):
    """
    :param chain: Chain in the PDB to compute the hydrophobic clusters. Examples: "A", "A B C". Default: "A"
    :param sel: VMD selection on which to compute the clusters. Default is every sidechain heavy atom ILE, VAL and LEU residues. "protein and not backbone and noh and resname ILE VAL LEU"
    :return: A representation for each cluster
    """
    clusters = None

    resids = mol.get("resid", sel="protein and name CA and resname ILE VAL LEU") #ILV cluster residues
    chains = mol.get("chain", sel="protein and name CA and resname ILE VAL LEU") # chains of ILV cluster residues
    dims = len(resids) # length of ILV clusters
    indices = mol.get("index", sel=f"{sel}") #the indices of the atoms in ILV clusters

    # get a dictionary of atom index and resid position in resids list
    atom_to_residposition = {}
    for index in indices:
        resid = mol.get("resid", sel=f"index {index}")[0]
        chain = mol.get("chain", sel=f"index {index}")[0]
        index_residue = [j for j, residue in enumerate(resids) if (residue == resid and chains[j] == chain) ][0]
        atom_to_residposition[index] = index_residue

    contacts = np.zeros((dims, dims))

    for index in indices:
        a = Atom(index, mol)
        if not a.neighbor_indices.any():
            continue
        contacts = fill_matrices(a, mol, contacts, indices, atom_to_residposition)

    graph = create_graph(contacts, resids, chains, cutoff_area=cutoff_area)
    comp, _ = label_components(graph)
    if comp.a.any():
        clusters = add_clusters(graph, comp)
    else:
        print("There are not residues in contact for this selection")

    return clusters



