# DDGWizard

  DDGWizard is a software pipeline for prediction of the changes in protein thermostability (ΔΔG/ddG) upon point mutation, based on broader feature space and data science process. DDGWizard continuously calls a series of software to extract features, then carries out RFE (Recursive Feature Elimination) feature selection and XGBoost machine learning. DDWizard generates HRMs (Hypothetical Reverse Mutations) involved in training to increase prediction ability of forward and reverse mutations, and according to benchmarking, it has achieved superior or comparable predictive performance to state-of-the-art algorithms. DDGWizard supports multi-process handling to meet the needs of large-scale computations, and its generation of 1547 relevant features makes it an equally effective tool for protein thermodynamic characterization.  

DDGWizard has two parts: A. Prediction part for calculating ddG. B. Characterization part for generating feature set to describe protein thermodynamics.   

They have different useages:   


# A. Prediction Part

## 1. Environment preparation steps

### (1). Git clone this repository
Do this:  

cd Your_Path/  

git clone https://github.com/Mingkai14/Features_Extraction.git  


### (2). Create conda virtual environment and install conda dependencies
  Conda yml file in Your_Path/Features_Extraction/src/environment.yml.  

  Your gcc version must be greater than 4.8.5 and conda version must be greater than 23.0.  
  
  Do this:   
  
  Change prefix in environment.yml to your own location of conda envs folder.  

  cd Your_Path/Features_Extraction/src/  
  
  conda env create -f environment.yml  
  
  conda activate Features_Extraction  

### (3). Create your Blast database 

  You need to prepare your own Blast database.   

  You need to download fasta database file first. We used uniref50.fasta file (https://ftp.uniprot.org/pub/databases/uniprot/uniref/).  

  We recommend you use blast+ 2.13.0 version.And there is a existing blast+ program folder in our program. Recommend you directly use this.   

  Do this:  

  cd Your_Path/Features_Extraction/bin/ncbi_blast_2_13_0+/bin/ 

  chmod -R +x .  

  ./makeblastdb -in Your_Path/uniref50.fasta -dbtype prot -out Your_Path/Your_DB_Name -parse_seqids  


### (4). Configure Modeller
  Modeller was used to generate mutation structures, it was already installed by conda. But it also need a license.  
  
[Registration (salilab.org)](https://salilab.org/modeller/registration.html)
  Enter in Modeller official website to register and get license.   
  
  Then modify Modeller conda config file to add license, which should be in:   
  
  Your_Conda_Path/envs/Feature_Extraction/lib/modeller-10.4/modlib/modeller/config.py  
  
  Replace XXXX to your license, save and close.  

  Then Modeller can be used.   

### (5). Configure DSSP
  DSSP was used to calculate RSA and secondary stuctures. Due to version issues, you must do operations below to make DSSP can be used.
  
  Do this: 

  conda activate Features_Extraction
  
  whereis mkdssp
  
  cd Your_Path_of_mkdssp
  
  cp mkdssp dssp

### (6). Configure R
   You also need to configure R and install R package "Bio3D"
   (I assume you have already installed R)
   
   Do this:
   
   R 
   
   > install.packages("bio3d")

### (7). Make sure all files have run permission
  Do this:
  
  cd Your_Path/Features_Extraction
  
  chmod -R +x .
  

## 2.Usage:

### (1). Tips
  You must cd to the top folder to run and make sure you are in Features_Extraction virtual environment and finish environment preparation.

  DDGWizard itself supports multiprocessing. We recommend utilizing our built-in multiprocessing fuction. Avoid running DDGWizard multiple times in the same folder, as conflicts may arise when the program matches files. If you genuinely need to implement multiprocessing or multithreading for running DDGWizard by your self, please make copies of the DDGWizard folder. Ensure that each instance of the DDGWizard program running in different processes/threads resides in a separate folder.

### (2). Predict_ddG_Executable.py

#### a. Description
  This python program aims to predict ddG.
#### b. Example

  You can run program Like:  
  
  conda activate Features_Extraction  
  
  cd Your_Path/Features_Extraction/  
  
  python Predict_ddG_Executable.py  
  --pred_dataset_path your_dataset.xls   
  --db_folder_path Your_Path/blast_db_folder/   
  --db_name db_name   
  --if_reversed_data 1   
  --blast_process_num 4    
  --mode whole   
  --process_num 4  

#### c. Arguments:
  --pred_dataset_path  
  Provide a xls file path, the file include data you want to predict.  
  
  File must be xls format and it has several attributes:  
  
  | Name | PDB_File_Path | Variation | Chain | pH | T |
  | 1SHG | /.../.../1SHG.pdb | Y57H |  A | 7 | 25 |
  
    <span style="font-size:small">Name (Identify this protein with a name consisting of fewer than 8 characters, and duplication is allowed)</span>  
    
    <span style="font-size:small">PDB_File_Path (The file path of the PDB protein that you need to predict. This must be an absolute path.)</span>  
    
    <span style="font-size:small">Variation (Specify the mutated amino acid, the mutation site number consistent with the PDB file, and the desired mutated amino acid. like: Y57H)</span>  
    
    <span style="font-size:small">Chain (Specify the chain number mutated amino acid located, consistent with the PDB file)</span>  
    
    <span style="font-size:small">pH (Specify pH)</span>  
    
    <span style="font-size:small">T (Specify temperature)</span>  

  There is a sample file in Your_Path/Features_Extraction/src/sample_pred.xls  
  
  --db_folder_path   
  Provide folder path of your blast database.   
  
  --db_name   
  Provide your blast database name.  
  
  --if_reversed_data   
  Provide 0 or 1, indicate if you also want to predict reverse mutations.    

  --blast_process_num 4  
  Provide a number less than 200 and greater than 0. DDGWizard will run blast in multi-process.  

  --mode   
  Default is "whole". DDGWizard will run complete processes. "model_only" mean only generate mutation structures. "blast_only" mean only run blast. When there is a large amount of data to be processed, this mechanism allows the task to be completed in segments.  

  --process_num   
  Provide a number less than 200 and greater than 0. DDGWizard will calculate data in multi-process.  

  

### c. Output:
There will be a output xls file in Your_Path/Features_Extraction/Pred_Res/, recording prediction results. 




# B. Characterization Part
## 1. Environment preparation steps
### (1)(2)(3)(4)(5)(6)(7)
  Complete same steps with Prediction part first.
### (8) Configure container evironment
  Characterization part also run some software in specific linux system. To solve platform compatibility, we use container to run these software. Our script will automatically call commands to run container. Docker and Singularity are supported. You only need to configure one of both.

  #### a. use docker
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

  #### b. use singularity 
[Quick Start — Singularity User Guide 3.7 documentation (sylabs.io)](https://docs.sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps)
   Follow singularity official instructions to install singularity. You don't need to additional configure using singularity.

## 2.Usage:




  





  
  
 


