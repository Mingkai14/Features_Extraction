# Features_Extraction
conda create python=3.10
chmod -R +x .
Modeller
cp mkdssp dssp
singularity/docker
rosetta
R Bio3D
conda 3.10, singularity, R
model_only - conda 
blast_only - conda
whole - conda, singularity, R



# Tip: 

  This software has two parts: 1. Extracting features from the raw data of ddg 2. predicting 
  ddg. Currently only completed the feature extraction part.

# Function:

  Pipeline of ddg features generating. Accept a raw data set, continuously call a series of software and scripts 
  to generate features data in a large scale. At present, there is no available specific software pipeline for 
  generating feature data of ddg prediction. This software collects and compares most of the software used in 
  other ddg predictors, meanwhile, add some new software, to complete feature extraction from raw data. 
  The generated features can not only be used to predict ddg, but also have described protein stability and 
  can be used to predict other indicators of protein stability (e.g. dTm).


# Environment preparation steps:

## 1. Download and unzip Rosetta
[License and Download | RosettaCommons](https://www.rosettacommons.org/software/license-and-download)
  Download Rosetta bin version from official website. I have used  
  rosetta_bin_linux_2021.16.61629_bundle version for test. 
 
  Then unzip Rosetta tgz file to Your_Path/Features_Extraction/bin/rosetta/.
  like: tar -zxvf rosetta_bin_linux_3.13_bundle.tgz -C Your_Path//Features_Extraction/bin/rosetta/
## 2. Configure container software
### a. use docker
[Install Docker Engine | Docker Documentation](https://docs.docker.com/engine/install/)
  Follow docker official instructions to install docker. 
  
  You must add your user to docker super privilege user group, which have set when docker  
  was installed.
  Do this: 
  sudo usermod -aG docker Your_User_Name
  
  And restart your linux system.
  
  Then load docker image in local folder.
  Do this: 
  docker load -i Your_Path/Features_Extraction/src/Prof_Source/myprof.tar
### b. use singularity 
[Quick Start — Singularity User Guide 3.7 documentation (sylabs.io)](https://docs.sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps)
 Follow singularity official instructions to install singularity 
### c. tip
 Regardless of docker or singularity, only need to configure it, program will automatically call it.

## 3. Create your Fasta database by Blast 
  like: makeblastdb -in uniref50 -dbtype prot -out Your_Path/Your_DB_Name -parse_seqids
  tip: recommend to use blast+ 2.13.0
## 4. Create conda virtual environment and install conda dependencies
  Conda yml file in Your_Path/Features_Extraction/src/requirements.yml.
  Do this: 
  (Do not forget change prefix in requirements.yml)
  conda env create -f requirements.yml python=3.10
  conda activate Features_Extraction
## 5. I use modeller to generate mutation structures, it was already installed by conda. But it also need a license.
[Registration (salilab.org)](https://salilab.org/modeller/registration.html)
  Enter in Modeller official website to register and get license. 
  
  Then modify Modeller conda config file to add license, which should be in: 
  Your_Conda_Path/envs/Feature_Extraction/lib/modeller-10.4/modlib/modeller/config.py
  Replace XXXX to your license, save and close.

  Then Software can be used 
## 6. Configure DSSP
  Do this: 
  whereis mkdssp
  cd Your_Path_of_mkdssp
  cp mkdssp dssp

## 7. Configure R
   You also need to configure R and install R package "Bio3D"
   Do this:
   R 
   > install.package("bio3d")

## 8. Make sure all file have run permission
  Do this:
  cd Your_Path/Features_Extraction
  chmod -R +x .
  

# Usage: 
  You must cd to local top folder to run and make sure you are in Features_Extraction virtual environment and finish environment configuration. Then run Generate_Dataset_Executable.py.
  Like:
  conda activate Features_Extraction
  cd Your_Path/Features_Extraction/
  
python Generate_Dataset_Executable.py 
--raw_dataset_path Your_Raw_Dataset/dataset.xls 
--db_folder_path Your_Path/blast_db_folder/ 
--db_name db_name 
--if_reversed_data 1 
--psiblast_threads_num 8 
--container_type D 
--mode whole 
--process_num 10

# Arguments:
  --raw_dataset_path represent your raw data path. 
  It must save as xls format. 
  The first row must have the following columns in order: PDB, Variation, Chain, ddG, pH, T. 
  A sample file is also in Your_Path/Features_Extraction/src/sample.xls.
  
  --db_folder_path represent folder path of your blast database.
  
  --db_name represent your blast database name.
  
  --if_reversed_data represent if reverse data. Since 2018, many ddg predictors start to use inverse data to train based 
  on theory that forward mutation ddg is the negative number of reverse mutation ddg. 
  Refer to this paper:
  [On the critical review of five machine learning-based algorithms for predicting protein stability changes upon mutation | Briefings in Bioinformatics | Oxford Academic (oup.com).](https://academic.oup.com/bib/article/22/1/601/5688895?searchresult=1)

  --psiblast_threads_num 4 represent threads number to run PSI-Blast, which is limited from 1-30.

  --container_type "D" represent use docker, "S" represent use singularity

  --mode model_only mean only generate mutation models, blast_only mean only generate blast file, whole mean complete process

  --process_num Multiple process number

  

# Output:
  After running, it will generate a csv file name Features_Table.csv in 
  Your_Path/Features_Extraction/src/Features_Table/


 

