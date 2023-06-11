### RING - Residue Interaction Network Generator
#### BioComputingUP - University of Padova

## CHANGELOG

### 3.0.0 | 05/01/2022

#### Input
 - Add support for PDBx/mmCIF format
 - Deprecate PSIBLAST input. 
It was used for calculation of mutual information and conservation columns in the nodes output
 
#### Algorithm optimization 
- Update core libraries to C++17 standard
- Add support for MacOS
- Calculation scales linearly with the dimension of the structure, before it was exponential
- Implement of sparse matrices for multi-state structures. It save RAM for structures with a lot of states (models)
- No limit on the number of states (models) for multi-state structures
- Minor and major bugfix on PDB reading and loading, low resolution structures, 
missing atoms or residues, model numeration
- Refined interactions between residues, with a new  that can give more accurate results.
- New internal representation of residues with automatic connection between atoms based on 
definitions provided in an external file (PDB format)
- Improve hydrogen placement, added support for nucleic acid hydrogen atoms

#### Usage 
- Add "install" directive in the make file
- No need to set the VICTOR_ROOT environment variable
- “--multi_edge” option include bifurcated interactions
- “--md” flag generate RING-MD output 
- Progress bar when running the program

#### Output
- Probabilistic network output for multi-sate structures (frequency of interaction across states)
- Modified output format of the edges and nodes files. Both contain a column indicating the model
- Edges and nodes output files are now single for each input file, 
even if it has been executed on multiple models 
- Deprecate naming of output files. When "--out_dir" is provided, the output files will have the name of the 
input file followed by “_ringEdges” or “_ringNodes” suffix. Otherwise, RING will print on the standard output/error 
- For RING-MD setting "--out_dir" is mandatory, a new “md” subdirectory will be created


### 2.0.0 | 01/06/2016

- RING-MD implementation
- PSIBLAST parsing to calculate nodes conservation and mutual information
- DSSP implementation
- Add hydrogen atoms
- Minor and major bug fixes

## 1.0.0 | 01/04/2011 

- Initial release