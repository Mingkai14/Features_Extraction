### RING - Residue Interaction Network Generator (v3.0.0)
#### BioComputingUP - University of Padova



### Usage

You can execute the program and print the help with the following command 
or just `Ring -h` if the installation folder has been added to the 
$PATH environment variable:

```
~/.ring/bin/Ring -h
```

Print nodes and edges (interactions) for all chains 
(both intra- and inter-chain contacts) and the first model of the structure:
```
Ring -i 2m6z.cif
```

Consider all models:
```
Ring -i 2m6z.cif --all_models
```

Write to file:
```
mkdir results
Ring -i 2m6z.cif --out_dir results
```
It generates two output files with `<fileName>_ringEdges` and `<fileName>_ringNodes` suffix.

##### RING-MD
RING - Molecular Dynamics (MD) is a RING module that computes aggregated statistics over multi-state structure files like
NMR structural ensembles or molecular dynamics snapshots (provided as different models in a PDB/mmCIF file). 

The module can be executed by adding the `--md` flag:

```
Ring -i 2m6z.cif --out_dir results --md
```

It generates the standard output files inside the **results/** folder and 
create a `md` subdirectory with four different types of data: 
- `<fileName>_cm_<type>`, all the contact maps, one per model. The first column is the model number
- `<fileName>_gcm_<type>`, the global contact map, number of contacts across models 
- `<fileName>_gfreq_<type>`, how many times each node is in contact over all models  
- `<fileName>_tdcm_<type>`, how many contacts for each node (rows) and each model (columns)

Where `<type>` is the type of interaction: `HBOND, IAC, IONIC, PICATION, PIPISTACK, SSBOND, VDW`







