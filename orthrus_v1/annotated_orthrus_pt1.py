# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# + [markdown] colab_type="text" id="view-in-github"
# <a href="https://colab.research.google.com/github/yc386/orthrus_metaproteomics/blob/main/annotated_orthrus_pt1_v1.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# + [markdown] id="VpisBTGyhi5v"
# <img src='https://drive.google.com/uc?export=view&id=19rmmQI1H2nIqgU598WROTcUNhOUoXcBP' width='400px' align='right'>
#
# # **Readme**
#
# ---
# Orthrus 🐾 is a hybrid, two-software pipeline that integrates [Casanovo](https://github.com/Noble-Lab/casanovo) (an AI transformer) with [Sage](https://github.com/lazear/sage) (a fast database search engine with advanced features like retention time alignment and machine learning-based rescoring).
#
# Designed to handle large search space and difficulties of selecting databases in metaproteomics and palaeoproteomics, Orthrus leverages de novo sequencing to define sample-specific databases, and uses probability ranking and conventional database searching to control FDRs (false discovery rates).
#
# Orthrus can be run online using Google Colab 🥳, or locally via Anaconda 🐍.
#
# ---

# + [markdown] id="Fva_FXk0htYl"
# # **Please note**❗️
# *   Before walking the dog, please change the runtime type to GPU (A100, L4, or T4. A100 most efficient but T4 is free)
# *   If you would like to connect your Google Drive, click the folder image 🗂️ on the left and mount the drive.
# *  Click `File` (top left) and save a copy in Drive or Github

# + [markdown] id="-K0wYqUruZmx"
# # Run `Orthrus`

# + cellView="form" id="Cv26MhBMhbT5"
#@title Add inputs -> click `Runtime` -> `Run all`
#@markdown **_De novo_ peptide sequencing algorithm inputs**
algorithm = "instanovo" #@param ["instanovo", "casanovo"]
#@markdown - use the drop-down menu to choose the de novo sequencing algorithm

folder_path="path/to/data"#@param {type:"string"}
#@markdown - a folder contains single or multiple `.mzML` or `.mgf` files for the de novo sequencing algorithm (`Instanovo` or `Casanovo`). Please check only _ (underscore) and no other special characters or space in a file name.
file_type="mzML" #@param ["mzML", "mgf"]
#@markdown - use the drop-down menu to choose the instrument file type

use_default = True #@param {type:"boolean"}
#@markdown **Advanced Options (ignored if using default settings)**

checkpoint = "path/to/model.ckpt" #@param {type:"string"}
#@markdown - path to a checkpoint `.ckpt` for a de novo peptide sequencing model
config = "path/to/config.yaml" #@param {type:"string"}
#@markdown - a `.yaml` configuration file for Casanovo

#@markdown **Inputs for converting Casanovo results to a `.fasta`**
use_SwissProt = True #@param {type:"boolean"}
#@markdown - use the latest, reviewed SwissProt form the UniProt FTP
database_path=""#@param {type:"string"}
#@markdown - path to a database (`.fasta`) for shortlisting proteins based on de novo results

# + cellView="form" colab={"base_uri": "https://localhost:8080/"} id="F9pGGlDUimMo" outputId="a8e81df3-8963-4180-d01e-6e99b37cbf75"
#@title install dependencies & modules

# %%time

# !pip install biopython mokapot

if algorithm == "instanovo":
  # !pip install instanovo
elif algorithm == "casanovo":
  # !pip install casanovo
else:
  raise ValueError("Invalid algorithm name")

import os

if not os.path.isfile("Orthrus_READY"):
  print("installing conda...")
  os.system("wget -qnc https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh")
  os.system("bash Mambaforge-Linux-x86_64.sh -bfp /usr/local")
  os.system("touch Orthrus_READY")
  os.system(f"conda install -c bioconda -c conda-forge sage-proteomics -y -q")


import glob
import pandas as pd
import numpy as np
import re
import json
from pyteomics import mztab
from numpy import string_
from joblib import Parallel, delayed
from itertools import chain
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import requests
import gzip
import shutil

'''
parse a .mztab file using pyteomics
add naked sequences (without PTMs) & the sequence length
input=path to .mztab file
output=pandas dataframe

'''
pattern = re.compile(r'(.\d*\.?\d+)')

def prep_mztab(mztab_path):
  m = mztab.MzTab(mztab_path)
  df = m.spectrum_match_table
  if df is None or df.empty:
    raise ValueError(f"{mztab_path} is empty")
  if 'sequence' not in df.columns:
    raise KeyError(f"'sequence' column is missing in the file: {mztab_path}")
  df.reset_index(drop=True)
  df1 = df.assign(sequence_naked=df['sequence'].str.replace(pattern, '', regex=True))
  df2= df1.assign(nAA=df1['sequence_naked'].str.len())
  df3=df2.sort_values(by='sequence_naked').drop_duplicates(subset='sequence_naked', keep="first").reset_index(drop=True)
  return df3

'''
parse a .fasta file using biopython
add UniProt ID e.g. P02754
input=path to .fasta file
output=pandas dataframe

'''

def fasta_to_df(fasta_file):

  data = []

  for record in SeqIO.parse(fasta_file, "fasta"):

    protein_id = record.id
    description = record.description
    sequence = str(record.seq)
    if not sequence:
      raise ValueError(f"Record with ID '{protein_id}' has no sequence in the fasta file.")

    data.append((protein_id, description, sequence))

  df = pd.DataFrame(data, columns=["Protein_ID", "Description", "Sequence"])
  df1=df.assign(UniProt_ID=df['Protein_ID'].str.split('|').str[1])

  return df1


'''
filter a Casanovo output file based on the maximum value below 0
search_engine_score[1] is a score assigned to each prediction by Casanovo, max=1,
if negative then outside the mass tolerance
input=pandas dataframe
output=pandas dataframe

'''

def casa_filter (df):
  np_array = df['search_engine_score[1]'].to_numpy()
  max_below_zero = np_array[np_array < 0].max()
  df1=df[df['search_engine_score[1]']>=max_below_zero]
  return df1


#prepare overlapping sequence tags for string matching
def get_seq_tags (sequence, k):
  return set(sequence[i:i+k] for i in range(len(sequence) - k + 1))


'''
match de novo-based tags with database tags
I=L in a reference database
inputs=path to fasta, filtered casanovo output dataframe, tag size=k, chunk size=10000 for processing
output=pandas dataframe

'''

def matching_count_v5 (df, df1, k, chunk_size=10000):

  sequence_set = get_seq_tags(''.join(chain.from_iterable(df1['sequence_naked'].astype(str))), k)
  print(f"📝 {len(sequence_set)} tags regenerated. Starting matching...")
  result_df = pd.DataFrame()
  for start in range(0, len(df), chunk_size):
    chunk = df.iloc[start:start+chunk_size].copy()
    chunk['seq_tags'] = chunk['Sequence'].astype(str).str.replace('I', 'L').apply(lambda x: get_seq_tags(x, k))
    chunk['matched_count'] = chunk['seq_tags'].apply(lambda seq_tags: len(seq_tags & sequence_set))
    chunk = chunk.assign(matched=chunk['matched_count'].apply(lambda x: 1 if x >= 2 else 0))
    result_df = pd.concat([result_df, chunk], ignore_index=True)
  total_matches=result_df['matched_count'].sum()
  print(f"Completed! {total_matches} matched 👍 ")
  return result_df

#get tryptic peptides per database entry
def count_tryptic_peptides(sequence):
  pattern=r'(?<=[KR])'

  peptides = re.split(pattern, sequence)

  filtered_peptides = [peptide for peptide in peptides if len(peptide) >= 6]

  return len(filtered_peptides)

#prepare a dataframe for NB classification
def prep_Bayes (df):
  print('🧑‍💻 Start Bayes probabilistic ranking...')
  df1=df.assign(length=df['Sequence'].astype(str).str.len(),
                 tryptic_count=df['Sequence'].apply(count_tryptic_peptides),
                 tag_count=df['seq_tags'].apply(len))
  df2=df1.assign(SAF=df1['matched_count']/df1['length'],
                 try_ratio=df1['tryptic_count']/df1['tag_count']
                 )
  return df2


'''
ranking matched proteins based on SAF and try_ratio
values normalised before NB classification
most probable matches (>= 95%) are shortlised

'''
def get_bayes_ranking_test (df, threshold=0.95):
  m=prep_Bayes(df)
  required_columns = {'SAF', 'try_ratio', 'matched'}
  if not required_columns.issubset(m.columns):
    missing_cols = required_columns - set(m.columns)
    raise ValueError(f"Missing columns in DataFrame: {missing_cols}")
  m1=m[m['tag_count']>0]
  X = m1[['SAF', 'try_ratio']].to_numpy()
  y = m1['matched'].to_numpy()
  scaler = MinMaxScaler()
  X_scaled = scaler.fit_transform(X.reshape(-1, 1)).reshape(*X.shape)
  X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=7)
  gnb = GaussianNB()
  gnb.fit(X_train, y_train)
  y_pred = gnb.predict(X_test)
  accuracy = accuracy_score(y_test, y_pred)
  precision = precision_score(y_test, y_pred)
  f1=f1_score(y_test, y_pred)
  print(f"✅ Gaussian Naive Bayes model ▶️ accuracy:{accuracy:.4f}, precision:{precision:.4f}, f1:{f1:.4f}")
  whole_pred=gnb.predict(X_scaled)
  class_probabilities = gnb.predict_proba(X_scaled)
  m2=m1.assign(pred=class_probabilities[:, 1])
  m3=m2[m2['pred']>=threshold]
  return m3


#combine previous functions together to output a shortlisted .fasta

def matching_ranking_to_fasta_mztab_v5(mztab_path, fasta_df):
  p=prep_mztab(mztab_path)
  p1=casa_filter(p)
  k=int(p1['nAA'].median())
  m=matching_count_v5 (fasta_df, p1, k, chunk_size=10000)
  m1=get_bayes_ranking_test (m)
  seq_records = []
  for index, row in m1.iterrows():
    header_id = f"{row['Description']}"
    sequence = Seq(row['Sequence'])
    description = ""
    seq_record = SeqRecord(sequence, id=header_id, description=description)
    seq_records.append(seq_record)

  output_fasta_filepath = mztab_path.replace('.mztab', '_matched.fasta')

  with open(output_fasta_filepath, 'w') as output_file:
    SeqIO.write(seq_records, output_file, 'fasta')
  print(f"🎊 Number of protein entries in the output fasta: {m1.shape[0]}")



#generate a de novo-first, experiment-specific .fasta for each input
def process_all_mztab_files_v2 (folder_path, database_path):
    mztab_filepaths = glob.glob(f"{folder_path}/*.mztab")
    print(f"🗂️ {len(mztab_filepaths)} file(s) collecting from {folder_path}...")
    fas=fasta_to_df(database_path)
    fasta_df=pd.DataFrame.from_dict(fas)
    print(f"⬆️ {database_path} loaded")
    print(f"📤 No. of proteins in the reference fasta: {fasta_df.shape[0]}")

    for mztab_filepath in mztab_filepaths:
      print(f"🚀 Processing file: {mztab_filepath}")
      matching_ranking_to_fasta_mztab_v5(mztab_filepath, fasta_df)

def process_all_csv_files(folder_path, database_path):
  csv_filepaths = glob.glob(f"{folder_path}/*.csv")
  print(f"🗂️ {len(csv_filepaths)} file(s) collecting from {folder_path}...")
  fas=fasta_to_df(database_path)
  fasta_df=pd.DataFrame.from_dict(fas)
  print(f"⬆️ {database_path} loaded")
  print(f"📤 No. of proteins in the reference fasta: {fasta_df.shape[0]}")

  for csv_filepath in csv_filepaths:
    print(f"🚀 Processing file: {csv_filepath}")
    matching_ranking_to_fasta(csv_filepath, fasta_df)

def process_all_files(folder_path, database_path, algorithm):
  if algorithm == "instanovo":
    process_all_csv_files(folder_path, database_path)
  elif algorithm == "casanovo":
    process_all_mztab_files_v2(folder_path, database_path)
  else:
    raise ValueError("Invalid algorithm name")



# +
#@title Run _de novo_ peptide sequencing algorithm

folder = glob.glob(f"{folder_path}/*.{file_type}")

if algorithm == "instanovo":
  for instrument_file in folder:
    output_path=instrument_file.replace(f".{file_type}", ".csv")
    if use_default:
      # ! curl -LRO https://github.com/instadeepai/InstaNovo/releases/download/1.0.0/instanovo_extended.ckpt
      # ! python -m instanovo.transformer.predict data_path={instrument_file} model_path="instanovo_extended.ckpt" denovo=True output_path={output_path} 
    else:
      # TODO add config
      # ! python -m instanovo.transformer.predict data_path={instrument_file} model_path={checkpoint} denovo=True output_path={output_path}
elif algorithm == "casanovo":
  for instrument_file in folder:
    output_path=instrument_file.replace(f".{file_type}", ".mztab")
    if use_default:
      # ! casanovo sequence {instrument_file} -v info -o {output_path}
    else:
      # ! casanovo sequence {instrument_file} -m {checkpoint} -c {config} -v info -o {output_path}

# + cellView="form" colab={"base_uri": "https://localhost:8080/"} id="CvDZIzfwrxbc" outputId="5114c954-2520-4fa8-e909-555fb99cf97f"
#@title Convert `Casanovo` results to .fasta per experiment

if use_SwissProt:
  url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
  output_file = "uniprot_sprot.fasta.gz"
  decompressed_file = "uniprot_sprot.fasta"
  response = requests.get(url, stream=True)
  if response.status_code == 200:
    with open(output_file, 'wb') as f:
      shutil.copyfileobj(response.raw, f)
    print(f"{output_file} downloaded successfully.")
  else:
    print(f"Failed to download {output_file}, status code: {response.status_code}")
  with gzip.open(output_file, 'rb') as f_in:
    with open(decompressed_file, 'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)
  sprot_path="uniprot_sprot.fasta"
  process_all_files(folder_path, sprot_path, algorithm)

else:
  process_all_files(folder_path, database_path, algorithm)

