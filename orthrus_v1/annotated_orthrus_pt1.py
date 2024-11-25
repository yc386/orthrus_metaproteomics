# %% [markdown]
# <a href="https://colab.research.google.com/github/yc386/orthrus_metaproteomics/blob/main/annotated_orthrus_pt1_v1.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# %% [markdown]
# <img src='https://drive.google.com/uc?export=view&id=19rmmQI1H2nIqgU598WROTcUNhOUoXcBP' width='400px' align='right'>
#
# # **Readme**
#
# ---
# Orthrus 🐾 is a hybrid, two-software pipeline that integrates de novo peptide sequencing (you have the choice between [InstaNovo](https://github.com/instadeepai/instanovo) and [Casanovo](https://github.com/Noble-Lab/casanovo)) with [Sage](https://github.com/lazear/sage) (a fast database search engine with advanced features like retention time alignment and machine learning-based rescoring).
#
# Designed to handle large search space and difficulties of selecting databases in metaproteomics and palaeoproteomics, Orthrus leverages de novo sequencing to define sample-specific databases, and uses probability ranking and conventional database searching to control FDRs (false discovery rates).
#
# Orthrus can be run online using Google Colab 🥳, or locally via Anaconda 🐍.
#
# ---

# %% [markdown]
# # **Please note**❗️
# *   Before walking the dog, please change the runtime type to GPU (A100, L4, or T4. A100 most efficient but T4 is free)
# *   If you would like to connect your Google Drive, click the folder image 🗂️ on the left and mount the drive.
# *  Click `File` (top left) and save a copy in Drive or Github

# %% [markdown]
# # Run `Orthrus`

# %%
import os
import subprocess


def run_command(command: str):
    """Runs a shell command and streams its output in real-time."""
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,  # Enables the use of a string command
        text=True,  # Ensures output is in string format
    )

    # Stream the output
    while True:
        output = process.stdout.readline()  # type: ignore
        if output == "" and process.poll() is not None:
            break
        if output:
            print(output.strip())

    # Check for errors
    stderr = process.stderr.read()  # type: ignore
    if stderr:
        print(f"Error: {stderr.strip()}")

    # Return the exit code
    return process.returncode


if os.getenv("COLAB_RELEASE_TAG"):
    print("Running in Colab")
    run_command("pip install -q condacolab")
    import condacolab

    condacolab.install()
    print(
        "Your session has restarted. You can ignore the 'Your session crashed for an unknown reason' warning."
    )
else:
    print("NOT running in Colab")
    print("Adding 'orthrus' conda environment to notebook")
    run_command("python -m ipykernel install --user --name=orthrus")

# %% [markdown]
# If running on Colab, please restart your runtime and run the cell below

# %%
import os
import urllib.request

if os.getenv("COLAB_RELEASE_TAG"):
    # TODO: update this  to the correct repository and branch once merged
    url = "https://raw.githubusercontent.com/BioGeek/orthrus_metaproteomics/refs/heads/aichor/environment.yml"
    urllib.request.urlretrieve(url, "environment.yml")
    run_command("conda env update -n base -f environment.yml")

# %%
#@title Add inputs -> click `Runtime` -> `Run all`
#@markdown **_De novo_ peptide sequencing algorithm inputs**
algorithm = "instanovo"  #@param ["instanovo", "casanovo"]
#@markdown - use the drop-down menu to choose the de novo sequencing algorithm

folder_path = "./data/PXD027613/mzML"  #@param {type:"string"}
#@markdown - a folder contains single or multiple `.mzML` or `.mgf` files for the de novo sequencing algorithm (`Instanovo` or `Casanovo`). Please check only _ (underscore) and no other special characters or space in a file name.
file_type = "mzML"  #@param ["mzML", "mgf"]
#@markdown - use the drop-down menu to choose the instrument file type

use_default = True  #@param {type:"boolean"}
#@markdown **Advanced Options (ignored if using default settings)**

checkpoint = "path/to/model.ckpt"  #@param {type:"string"}
#@markdown - path to a checkpoint `.ckpt` for a de novo peptide sequencing model
config = "path/to/config.yaml"  #@param {type:"string"}
#@markdown - a `.yaml` configuration file for Casanovo

#@markdown **Inputs for converting Casanovo results to a `.fasta`**
use_SwissProt = True  #@param {type:"boolean"}
#@markdown - use the latest, reviewed SwissProt form the UniProt FTP
database_path = ""  #@param {type:"string"}
#@markdown - path to a database (`.fasta`) for shortlisting proteins based on de novo results

# %%
#@title install dependencies & modules

import glob
import pandas as pd
import re
from pyteomics import mztab
from itertools import chain
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, precision_score, f1_score
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import requests
import gzip
import shutil
import s3fs

pattern = re.compile(r"(.\d*\.?\d+)")


def prep_mztab(mztab_path):
    """Parse a .mztab file using pyteomics.

    add naked sequences (without PTMs) & the sequence length
    input=path to .mztab file
    output=pandas dataframe
    """
    m = mztab.MzTab(mztab_path)
    df = m.spectrum_match_table
    if df is None or df.empty:
        raise ValueError(f"{mztab_path} is empty")
    if "sequence" not in df.columns:
        raise KeyError(f"'sequence' column is missing in the file: {mztab_path}")
    df.reset_index(drop=True)
    df1 = df.assign(sequence_naked=df["sequence"].str.replace(pattern, "", regex=True))
    df2 = df1.assign(nAA=df1["sequence_naked"].str.len())
    df3 = (
        df2.sort_values(by="sequence_naked")
        .drop_duplicates(subset="sequence_naked", keep="first")
        .reset_index(drop=True)
    )
    return df3


def prep_csv(csv_path):
    """Parse a .csv file using pandas.

    add naked sequences (without PTMs) & the sequence length
    """
    df = pd.read_csv(csv_path)
    if df is None or df.empty:
        raise ValueError(f"{csv_path} is empty")
    df = df.rename(columns={"preds": "sequence"})
    if "sequence" not in df.columns:
        raise KeyError(f"'sequence' column is missing in the file: {csv_path}")
    df.reset_index(drop=True)
    df1 = df.assign(sequence_naked=df["sequence"].str.replace(pattern, "", regex=True))
    df2 = df1.assign(nAA=df1["sequence_naked"].str.len())
    df3 = (
        df2.sort_values(by="sequence_naked")
        .drop_duplicates(subset="sequence_naked", keep="first")
        .reset_index(drop=True)
    )
    return df3


def fasta_to_df(fasta_file):
    """Parse a .fasta file using biopython.

    add UniProt ID e.g. P02754
    input=path to .fasta file
    output=pandas dataframe
    """
    data = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id
        description = record.description
        sequence = str(record.seq)
        if not sequence:
            raise ValueError(
                f"Record with ID '{protein_id}' has no sequence in the fasta file."
            )

        data.append((protein_id, description, sequence))

    df = pd.DataFrame(data, columns=["Protein_ID", "Description", "Sequence"])
    df1 = df.assign(UniProt_ID=df["Protein_ID"].str.split("|").str[1])

    return df1


def filter_casanovo(df):
    """Filter a de novo output file based on the maximum value below 0.

    search_engine_score[1] is a score assigned to each prediction by Casanovo, max=1,
    if negative then outside the mass tolerance
    """
    np_array = df["search_engine_score[1]"].to_numpy()
    max_below_zero = np_array[np_array < 0].max()
    df1 = df[df["search_engine_score[1]"] >= max_below_zero]
    return df1


def filter_instanovo(df):
    """TODO!"""
    return df


# prepare overlapping sequence tags for string matching
def get_seq_tags(sequence, k):
    """Generate overlapping sequence tags of size k."""
    return set(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def matching_count(df, df1, k, chunk_size=10000):
    """Match de novo-based tags with database tags.

    I=L in a reference database
    inputs=path to fasta, filtered casanovo output dataframe, tag size=k, chunk size=10000 for processing
    output=pandas dataframe
    """
    sequence_set = get_seq_tags(
        "".join(chain.from_iterable(df1["sequence_naked"].astype(str))), k
    )
    print(f"📝 {len(sequence_set)} tags regenerated. Starting matching...")
    result_df = pd.DataFrame()
    for start in range(0, len(df), chunk_size):
        chunk = df.iloc[start : start + chunk_size].copy()
        chunk["seq_tags"] = (
            chunk["Sequence"]
            .astype(str)
            .str.replace("I", "L")
            .apply(lambda x: get_seq_tags(x, k))
        )
        chunk["matched_count"] = chunk["seq_tags"].apply(
            lambda seq_tags: len(seq_tags & sequence_set)
        )
        chunk = chunk.assign(
            matched=chunk["matched_count"].apply(lambda x: 1 if x >= 2 else 0)
        )
        result_df = pd.concat([result_df, chunk], ignore_index=True)
    total_matches = result_df["matched_count"].sum()
    print(f"Completed! {total_matches} matched 👍 ")
    return result_df


# get tryptic peptides per database entry
def count_tryptic_peptides(sequence):
    """Count tryptic peptides in a protein sequence."""
    pattern = r"(?<=[KR])"

    peptides = re.split(pattern, sequence)

    filtered_peptides = [peptide for peptide in peptides if len(peptide) >= 6]

    return len(filtered_peptides)


# prepare a dataframe for NB classification
def prep_Bayes(df):
    """Prepare a dataframe for Naive Bayes classification."""
    print("🧑‍💻 Start Bayes probabilistic ranking...")
    df1 = df.assign(
        length=df["Sequence"].astype(str).str.len(),
        tryptic_count=df["Sequence"].apply(count_tryptic_peptides),
        tag_count=df["seq_tags"].apply(len),
    )
    df2 = df1.assign(
        SAF=df1["matched_count"] / df1["length"],
        try_ratio=df1["tryptic_count"] / df1["tag_count"],
    )
    return df2


def get_bayes_ranking_test(df, threshold=0.95):
    """Ranking matched proteins based on SAF and try_ratio.

    values normalised before NB classification
    most probable matches (>= 95%) are shortlised
    """
    m = prep_Bayes(df)
    required_columns = {"SAF", "try_ratio", "matched"}
    if not required_columns.issubset(m.columns):
        missing_cols = required_columns - set(m.columns)
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")
    m1 = m[m["tag_count"] > 0]
    X = m1[["SAF", "try_ratio"]].to_numpy()
    y = m1["matched"].to_numpy()
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X.reshape(-1, 1)).reshape(*X.shape)
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.2, random_state=7
    )
    gnb = GaussianNB()
    gnb.fit(X_train, y_train)
    y_pred = gnb.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    print(
        f"✅ Gaussian Naive Bayes model ▶️ accuracy:{accuracy:.4f}, precision:{precision:.4f}, f1:{f1:.4f}"
    )
    # whole_pred = gnb.predict(X_scaled)
    class_probabilities = gnb.predict_proba(X_scaled)
    m2 = m1.assign(pred=class_probabilities[:, 1])
    m3 = m2[m2["pred"] >= threshold]
    return m3


# combine previous functions together to output a shortlisted .fasta


def matching_ranking_to_fasta_mztab(mztab_path, fasta_df):
    """Generate a fasta file based on the matched proteins."""
    denovo_df = prep_mztab(mztab_path)
    denovo_df = filter_casanovo(denovo_df)
    filestem = os.path.splitext(mztab_path)[0]
    return matching_ranking_to_fasta(denovo_df, fasta_df, filestem)


def matching_ranking_to_fasta_csv(csv_path, fasta_df):
    """Generate a fasta file based on the matched proteins."""
    denovo_df = prep_csv(csv_path)
    denovo_df = filter_instanovo(denovo_df)
    filestem = os.path.splitext(csv_path)[0]
    return matching_ranking_to_fasta(denovo_df, fasta_df, filestem)


def upload_to_bucket(output_path):
    """Upload results to a bucket."""
    # Only applicable when running on https://aichor.ai/
    if "AICHOR_OUTPUT_PATH" in os.environ:
        s3_endpoint = "https://storage.googleapis.com"  # os.environ["S3_ENDPOINT"]
        s3_key = os.environ["AWS_ACCESS_KEY_ID"]
        s3_secret_key = os.environ["AWS_SECRET_ACCESS_KEY"]
        s3 = s3fs.S3FileSystem(
            client_kwargs={"endpoint_url": s3_endpoint},
            key=s3_key,
            secret=s3_secret_key,
        )
        bucket_path = (
            f"{os.environ['AICHOR_OUTPUT_PATH']}{os.path.basename(output_path)}"
        )
        with open(output_path, "r") as local_file, s3.open(
            bucket_path, mode="w"
        ) as bucket_file:
            bucket_file.write(local_file.read())
        print(f" 🪣 Results uploaded to {bucket_path}")


def matching_ranking_to_fasta(denovo_df, fasta_df, filestem):
    """Generate a fasta file based on the matched proteins."""
    k = int(denovo_df["nAA"].median())
    m = matching_count(fasta_df, denovo_df, k, chunk_size=10000)
    m1 = get_bayes_ranking_test(m)
    seq_records = []
    for _index, row in m1.iterrows():
        header_id = f"{row['Description']}"
        sequence = Seq(row["Sequence"])
        description = ""
        seq_record = SeqRecord(sequence, id=header_id, description=description)
        seq_records.append(seq_record)

    output_fasta_filepath = f"{filestem}_matched.fasta"

    with open(output_fasta_filepath, "w") as output_file:
        SeqIO.write(seq_records, output_file, "fasta")
    print(f"🎊 Number of protein entries in the output fasta: {m1.shape[0]}")
    upload_to_bucket(output_fasta_filepath)


# generate a de novo-first, experiment-specific .fasta for each input
def process_all_mztab_files(folder_path, database_path):
    """Process all .mztab files in a folder."""
    mztab_filepaths = glob.glob(f"{folder_path}/*.mztab")
    print(f"🗂️ {len(mztab_filepaths)} file(s) collecting from {folder_path}...")
    fas = fasta_to_df(database_path)
    fasta_df = pd.DataFrame.from_dict(fas)
    print(f"⬆️ {database_path} loaded")
    print(f"📤 No. of proteins in the reference fasta: {fasta_df.shape[0]}")

    for mztab_filepath in mztab_filepaths:
        print(f"🚀 Processing file: {mztab_filepath}")
        matching_ranking_to_fasta_mztab(mztab_filepath, fasta_df)


def process_all_csv_files(folder_path, database_path):
    """Process all .csv files in a folder."""
    csv_filepaths = glob.glob(f"{folder_path}/*.csv")
    print(f"🗂️ {len(csv_filepaths)} file(s) collecting from {folder_path}...")
    fasta_df = fasta_to_df(database_path)
    print(f"⬆️ {database_path} loaded")
    print(f"📤 No. of proteins in the reference fasta: {fasta_df.shape[0]}")

    for csv_filepath in csv_filepaths:
        print(f"🚀 Processing file: {csv_filepath}")
        matching_ranking_to_fasta_csv(csv_filepath, fasta_df)


def process_all_files(folder_path, database_path, algorithm):
    """Process all files in a folder."""
    if algorithm == "instanovo":
        process_all_csv_files(folder_path, database_path)
    elif algorithm == "casanovo":
        process_all_mztab_files(folder_path, database_path)
    else:
        raise ValueError("Invalid algorithm name")


# %%
#@title Run _de novo_ peptide sequencing algorithm

folder = glob.glob(f"{folder_path}/*.{file_type}")


if algorithm == "instanovo":
    print("🔍 Running InstaNovo...")
    for instrument_file in folder:
        print(f"🚀 Processing file: {instrument_file}")
        base, ext = instrument_file.rsplit(".", 1)
        output_path = f"{base}_{algorithm}.csv"
        if use_default:
            if not os.path.isfile("instanovo_extended.ckpt"):
                print("⏬ Downloading InstaNovo checkpoint")
                run_command(
                    "curl -LRO https://github.com/instadeepai/InstaNovo/releases/download/1.0.0/instanovo_extended.ckpt"
                )
            if not os.path.isfile(output_path):
                run_command(
                    f"python -m instanovo.transformer.predict data_path={instrument_file} model_path='instanovo_extended.ckpt' denovo=True output_path={output_path}"
                )
        else:
            # TODO add config
            if not os.path.isfile(output_path):
                run_command(
                    f"python -m instanovo.transformer.predict data_path={instrument_file} model_path={checkpoint} denovo=True output_path={output_path}"
                )
elif algorithm == "casanovo":
    print("🔍 Running Casanovo...")
    for instrument_file in folder:
        print(f"🚀 Processing file: {instrument_file}")
        base, ext = instrument_file.rsplit(".", 1)
        output_path = f"{base}_{algorithm}.mztab"
        if use_default:
            run_command(f"casanovo sequence {instrument_file} -v info -o {output_path}")
        else:
            run_command(
                f"casanovo sequence {instrument_file} -m {checkpoint} -c {config} -v info -o {output_path}"
            )
else:
    raise ValueError("Invalid algorithm name")
upload_to_bucket(output_path)


# %%
#@title Convert de novo results to .fasta per experiment

if use_SwissProt:
    url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
    output_file = "uniprot_sprot.fasta.gz"
    sprot_path = "uniprot_sprot.fasta"
    if not os.path.isfile(sprot_path):
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(output_file, "wb") as f:
                shutil.copyfileobj(response.raw, f)
            print(f"{output_file} downloaded successfully.")
        else:
            print(
                f"Failed to download {output_file}, status code: {response.status_code}"
            )
        with gzip.open(output_file, "rb") as f_in:
            with open(sprot_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    process_all_files(folder_path, sprot_path, algorithm)

else:
    process_all_files(folder_path, database_path, algorithm)

# %%
