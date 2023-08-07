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

# Environment preparation:

  This software try to be easily to migrate to different systems. It don’t need to install additional specific software, 
  except Rosetta, because it’s too big. So it need to download and unzip Rosetta to specific folder.
  
  And there is a software Prof which can only install in ubuntu 18.04, so I used docker to containerize it and made a light docker image in local folder.  Just need to install docker and load this image,  built-in scripts will automatically call this software.
  
  It also need to provide a Fasta database path, which must be created by Blast. 
  
  All other executable software and scripts have saved in local folder. Please follow next steps to prepare environment.

# Environment preparation steps:

## 1. Download and unzip Rosetta
[License and Download | RosettaCommons](https://www.rosettacommons.org/software/license-and-download)
  Download Rosetta bin version from official website. I have used  
  rosetta_bin_linux_2021.16.61629_bundle version for test. 
 
  Then unzip Rosetta tgz file to Your_Path/Features_Extraction/bin/rosetta/.
  like: tar -zxvf rosetta_bin_linux_3.13_bundle.tgz -C Your_Path//Features_Extraction/bin/rosetta/
## 2. Install and config docker
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
## 3. Create your Fasta database by Blast 
  like: makeblastdb -in uniref50 -dbtype prot -out Your_Path/Your_DB_Name -parse_seqids
## 4. Create conda virtual environment and install conda dependencies
  Conda yml file in Your_Path/Features_Extraction/src/requirements.yml.
  Do this: 
  (Do not forget change prefix in requirements.yml)
  conda env create -f requirements.yml
  conda activate Features_Extraction
## 5. I use modeller to generate mutation structures, it was already installed by conda. But it also need a license.
[Registration (salilab.org)](https://salilab.org/modeller/registration.html)
  Enter in Modeller official website to register and get license. 
  
  Then modify Modeller conda config file to add license, which should be in: 
  Your_Conda_Path/envs/Feature_Extraction/lib/modeller-10.4/modlib/modeller/config.py
  Replace XXXX to your license, save and close.

  Then Software can be used 

# Usage: 
  You must cd to local top folder to run and make sure you are in Features_Extraction virtual environment. Then run Generate_Dataset_Executable.py.
  Like:
  conda activate Features_Extraction
  cd Your_Path/Features_Extraction/
  
  python3 Generate_Dataset_Executable.py
 --raw_dataset_path Your_Raw_Dataset.xls 
--db_folder_path Your_DB_Path/ 
--db_name Your_DB_Name 
--if_reversed_data 1 
--psiblast_threads_num 4

# Arguments:
  --raw_dataset_path represent your raw data path. 
  It must save as xls format. The first row must have the following columns in order: PDB, Variation, Chain, ddG, pH, T. 
  A sample file is also in Your_Path/Features_Extraction/src/sample.xls.
  
  --db_folder_path represent folder path to your Fasta database.
  
  --db_name represent your Fasta database name.
  
  --if_reversed_data represent if reverse data. Since 2018, many ddg predictors start to use inverse data to train based 
  on theory that forward mutation ddg is the negative number of reverse mutation ddg. 
  Refer to this paper:
  [On the critical review of five machine learning-based algorithms for predicting protein stability changes upon mutation | Briefings in Bioinformatics | Oxford Academic (oup.com).](https://academic.oup.com/bib/article/22/1/601/5688895?searchresult=1)

  --psiblast_threads_num 4 represent threads number to run PSI-Blast, which is limited from 1-30.

# Output:
  After running, it will generate a csv file name Features_Table.csv in 
  Your_Path/Features_Extraction/src/Features_Table/


 
