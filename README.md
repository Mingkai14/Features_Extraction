

# Function:

  1. Pipeline of ddG features generating. Accept a raw data set, continuously call a series of software and scripts 
  to generate features data in a large scale. At present, there is no available specific software pipeline for 
  generating feature data of ddG prediction. This pipeline collects and compares the most of software used in 
  other ddg predictors, and meanwhile, adds some new software to extract features. The generated features can not only
  be used to predict ddG, but also have described protein stability and can be used to predict other indicators
  of protein stability (e.g. dTm).

  2. Prediction of ddG. According to features we generated, we have trained models (AdaBoost Regressor, Decision Tree Regression, SVM Regression, Linear 
  Regression, Random Forest Regression, XG Boost Regression) and implemented predicting function.


# Environment preparation steps:
## 1. Git clone this repository
Do this:

cd Your_Path/

git clone https://github.com/Mingkai14/Features_Extraction.git

## 2. Configure container evironment

  In order to fix problems that some software can only run on specific system, we use container to run these software. Our script will automatically call commands to run container. 
  Docker and Singularity are supported. You only need to configure one of both, program will call it.

  ### a. use docker
[Install Docker Engine | Docker Documentation](https://docs.docker.com/engine/install/)
  Follow docker official instructions to install docker. 
  
  You must add your user to docker super privilege user group, which should have set when docker  
  was installed. 
  
  Do this: 
  
  sudo usermod -aG docker Your_User_Name
  
  And restart your linux system.
  
  Then load docker image from software folder.
  
  Do this: 
  
  docker load -i Your_Path/Features_Extraction/src/Prof_Source/myprof.tar

  ### b. use singularity 
[Quick Start — Singularity User Guide 3.7 documentation (sylabs.io)](https://docs.sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps)
   Follow singularity official instructions to install singularity 

## 3. Create conda virtual environment and install conda dependencies
  Conda yml file in Your_Path/Features_Extraction/src/environment.yml.

  Your gcc version must be greater than 4.8.5 and conda version must be greater than 23.0.
  
  Do this: 
  
  Change prefix in environment.yml to your own location of conda envs folder.

  cd Your_Path/Features_Extraction/src/
  
  conda env create -f environment.yml
  
  conda activate Features_Extraction

## 4. Create your Blast database 

  You need to prepare your own Blast database. 

  You need to download fasta database file first. We used uniref50.fasta file.

  We recommend you use blast+ 2.13.0 version.And there is a existing blast+ program folder in our program. Recommend you direct use this. 

  Do this:

  cd Your_Path/Features_Extraction/bin/ncbi_blast_2_13_0+/bin/

  chmod -R +x .

  ./makeblastdb -in Your_Path/uniref50.fasta -dbtype prot -out Your_Path/Your_DB_Name -parse_seqids

## 5. Download and unzip FoldX Linux version
[Homepage | FoldX (crg.eu)](https://foldxsuite.crg.eu/)
 Register and download FoldX Linux version from official website.

 There is already a FoldX bin program in Your_Path/Features_Extraction/bin/FoldX_5.0/, but it may expire. You should 
 delete all files in there and unzip your FoldX 5.0 to this folder.


## 6. Configure Modeller
  Modeller was used to generate mutation structures, it was already installed by conda. But it also need a license.
  
[Registration (salilab.org)](https://salilab.org/modeller/registration.html)
  Enter in Modeller official website to register and get license. 
  
  Then modify Modeller conda config file to add license, which should be in: 
  
  Your_Conda_Path/envs/Feature_Extraction/lib/modeller-10.4/modlib/modeller/config.py
  
  Replace XXXX to your license, save and close.

  Then Modeller can be used. 

## 7. Configure DSSP
  DSSP was used to calculate RSA and secondary stuctures. Due to version issues, you must do operations below to make DSSP can be used.
  
  Do this: 

  conda activate Features_Extraction
  
  whereis mkdssp
  
  cd Your_Path_of_mkdssp
  
  cp mkdssp dssp

## 8. Configure R
   You also need to configure R and install R package "Bio3D"
   (I assume you have already installed R)
   
   Do this:
   
   R 
   
   > install.packages("bio3d")

## 9. Make sure all file have run permission
  Do this:
  
  cd Your_Path/Features_Extraction
  
  chmod -R +x .
  

# Usage: 
  You must cd to the top folder to run and make sure you are in Features_Extraction virtual environment and finish environment preparation. 

## 1. Generate_Dataset_Executable.py

### a. Description
  This python program aims to extract features from raw data. You need to provide a raw dataset, blast database for multiple sequence alignment and some 
  addional arguments.

  You can run program Like:
  
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
  --process_num 8

### b. Arguments:
  --raw_dataset_path represent your raw data path.
  
  It must save as xls format. 
  
  The first row must have these attributes: 
  
  PDB, Variation, Chain, ddG, pH, T 
  
  A sample file is in Your_Path/Features_Extraction/src/sample.xls.
  
  --db_folder_path represent folder path of your blast database.
  
  --db_name represent your blast database name.
  
  --if_reversed_data represent if reverse data. Since 2018, many ddg predictors start to use inverse data to train based 
  on a theory that forward mutation ddG equal with the negative number of reverse mutation ddG. 
  
  Refer to this paper:
  [On the critical review of five machine learning-based algorithms for predicting protein stability changes upon mutation | Briefings in Bioinformatics | Oxford Academic (oup.com).](https://academic.oup.com/bib/article/22/1/601/5688895?searchresult=1)

  --psiblast_threads_num represent threads number to run PSI-Blast, which is limited from 1-40.

  --container_type "D" represent use docker, "S" represent use singularity

  --mode "model_only" mean only generate mutation models, "blast_only" mean only generate blast output files, "whole" mean completely process.
  
  You can generate features by separate sections. "model_only" and "blast_only" will gernerate and save files required by "whole". When run "whole",
  if it find you already have required files, it won't run "model_only" and "blast_only" again. So you can continuously run "model_only", "blast_only",   
  "whole". This mechanism is to prevent the program from running too long. 

  --process_num Multiple process number, which is limited from 1-40

  

### c. Output:
  After running, it will generate a csv file name features_table.csv in: 
  
  Your_Path/Features_Extraction/src/Features_Table/

## 2. Train_and_Evaluate_Models_Executable.py

### a. Description
  This python program aims to train and evaluate models from dataset after feature extraction. Program will read features_table.csv from Features_Table/    folder, so you must run Generate_Dataset_Executable.py first and don't move its output file.

  You can run program Like:
  (I assume you already run Generate_Dataset_Executable.py first, and features_table.csv is generated)
  
  conda activate Features_Extraction

  cd Your_Path/Features_Extraction/

  python Train_and_Evaluate_Models_Executable.py 
  --model_type XGB
  --model_saving_path ./models/

### b. Arguments:

  --model_type represents which model you want to train, which must be in these options: "XGB", "SVR", "RF", "LR", "DTR", "ABR".

  --model_saving_path represents you want to save model file to which folder

### c. Output:
  Program will generate a model file and a StandardScaler file corresponding to your chosen model_type to your chosen folder. It also will print
  evaluation infomaton to console.

## 3. Predict_ddG_Executable.py
### a. Description
  This python program aims to predict ddG value. You need to provide a PDB file and your mutation infomation. It will extract features from your PDB file
  first, then according to XGB model in "./models/" to predict ddG value.

  You can run program Like:

  conda activate Features_Extraction

   cd Your_Path/Features_Extraction/

   python Predict_ddG_Executable.py
   --pdb_name Your_PDB_Name
   --pdb_path Your_PDB_File_Path/example.pdb
   --variation A100S
   --chain A
   --pH 7
   --T 25
   --db_folder_path Your_Path/blast_db_folder/ 
   --db_name db_name 
   --psiblast_threads_num 8
   --container_type D
   --mode whole

### b. Arguments:
  --pdb_name represents your PDB name, it will be as ID of this task
  
  --pdb_path represents your PDB file path
  
  --variation means which amino acid site you want to mutate, and what are WT amino acid and mutation amino acid. The WT amino acid and mutation
  site must be consistent. And mutation amino acid must be in 20 type of amino acid.
  
  --chain represents your mutation site in which chain
  
  --pH represents condition of pH
  
  --T represents condition of temperature
  
  --db_folder_path represent folder path of your blast database
  
  --db_name represent your blast database name
  
  --psiblast_threads_num represent threads number to run PSI-Blast, which is limited from 1-40

  --container_type "D" represent use docker, "S" represent use singularity

  --mode "model_only" mean only generate mutation models, "blast_only" mean only generate blast output files, "whole" mean completely process

### c. Output:
  It will print predicted ddG value.
  





  
  
 


