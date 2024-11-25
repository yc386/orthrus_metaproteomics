# %% [markdown]
# <a href="https://colab.research.google.com/github/yc386/orthrus_metaproteomics/blob/main/annotated_orthrus_pt2.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# %% [markdown]
# # **Orthrus pt2- [`Sage`](https://github.com/lazear/sage) + [`Mokapot`](https://github.com/wfondrie/mokapot)**
#
# Please note: change to **TPU** runtime if RAM usage is expected to be high (due to files/PTMs/databases)

# %%
# @title Add inputs for `SAGE` -> click `Runtime` -> `Run all`

# @markdown **Parameters for `SAGE`**
peak_folder = "test_data/PXD027613/mzML"  # @param {type:"string"}
file_type = "mzML"  # @param ["mzML", "mgf"]
# @markdown - use the drop-down menu to choose the instrument file type

# @markdown **_De novo_ peptide sequencing algorithm**
algorithm = "instanovo"  # @param ["instanovo", "casanovo"]

# @markdown **Option 1: `SAGE` PTM plus**
# @markdown - Default `Sage` contains
use_PTM_plus = True  # @param {type:"boolean"}
max_variable_mods = 3  # @param {type:"number"}
missed_cleavages = 2  # @param {type:"number"}
AA_1 = "M"  # @param ["None", "[","]","A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AA_1_mod = 15.9949  # @param {type:"number"}
AA_2 = "P"  # @param ["None", "[","]","A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AA_2_mod = 15.9949  # @param {type:"number"}
AA_3 = "N"  # @param ["None", "[","]","A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AA_3_mod = 0.984016  # @param {type:"number"}
AA_4 = "Q"  # @param ["None", "[","]","A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AA_4_mod = 0.984016  # @param {type:"number"}
AA_5 = "None"  # @param ["None", "[","]","A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
# @markdown - [ = n-terminal
AA_5_mod = 42.010565  # @param {type:"number"}

# @markdown **Option 2: `Mokapot` rescoring**
use_Mokapot = True  # @param {type:"boolean"}
# @markdown - machine learning-based rescoring per experiment or across experiments
joint_modelling = False  # @param {type:"boolean"}
# @markdown - a joint model for low abundance samples
default_Percolator = False  # @param {type:"boolean"}
# @markdown - Python implementation of the Percolator SVM model

# %%
# @title install dependencies

import os
import shutil
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


def organise_files(directory, algorithm):
    """Organise mztabs and instrument files in the same folder."""
    if not os.path.isdir(directory):
        print(f"The directory {directory} does not exist.")
        return

    mzml_files = glob.glob(os.path.join(directory, "*.mzML"))

    for mzml_file in mzml_files:
        base_name = os.path.splitext(os.path.basename(mzml_file))[0]
        core_name = base_name.split("_")[
            0
        ]  # Extract the part before the de novo algorithm

        new_folder_path = os.path.join(directory, base_name)

        if not os.path.exists(new_folder_path):
            os.makedirs(new_folder_path)
            print(f"Created folder: {new_folder_path}")
        else:
            print(f"Folder already exists: {new_folder_path}")

        new_mzml_path = os.path.join(new_folder_path, os.path.basename(mzml_file))
        if not os.path.exists(new_mzml_path):
            shutil.move(mzml_file, new_folder_path)
            print(f"Moved {mzml_file} to {new_folder_path}")
        else:
            print(f".mzML file already exists in the destination: {new_mzml_path}")

        fasta_filename = f"{core_name}_{algorithm}_matched.fasta"
        fasta_file = os.path.join(directory, fasta_filename)

        if os.path.exists(fasta_file):
            new_fasta_path = os.path.join(new_folder_path, fasta_filename)
            if not os.path.exists(new_fasta_path):
                shutil.move(fasta_file, new_folder_path)
                print(f"Moved {fasta_file} to {new_folder_path}")
            else:
                print(
                    f".fasta file already exists in the destination: {new_fasta_path}"
                )
        else:
            print(f"No matching .fasta file found for {base_name}")


def get_sage_config(
    json_file_path,
    peak_folder,
    static_mods,
    new_mods,
    missed_cleavages,
    min_len,
    max_len,
    max_variable_mods,
    output_config_path,
):
    """Update the Sage config file with the new parameters."""
    with open(json_file_path, "r") as file:
        json_data = json.load(file)

        peak_files = glob.glob(peak_folder)
        print(f"🗂️ {len(peak_files)} file(s) collected from {peak_folder}")
        json_data["mzml_paths"] = peak_files
        json_data["database"]["static_mods"] = static_mods
        json_data["database"]["variable_mods"] = new_mods
        json_data["database"]["enzyme"]["missed_cleavages"] = missed_cleavages
        json_data["database"]["enzyme"]["min_len"] = min_len
        json_data["database"]["enzyme"]["max_len"] = max_len
        json_data["database"]["max_variable_mods"] = max_variable_mods
        json_data["database"]["decoy_tag"] = "rev_"
        json_data["database"]["generate_decoys"] = True

    with open(output_config_path, "w") as f:
        json.dump(json_data, f, indent=4)


# %%
# @title Run Sage
import glob
import json


organise_files(peak_folder, algorithm)
folder_path = peak_folder

if use_PTM_plus:
    AAs = [AA_1, AA_2, AA_3, AA_4, AA_5]
    mods = [AA_1_mod, AA_2_mod, AA_3_mod, AA_4_mod, AA_5_mod]
    PTMs = {}
    for AA, mod in zip(AAs, mods):
        if AA != "None":
            PTMs[AA] = [mod]
    big_folder = glob.glob(f"{folder_path}/*")
    for folder in big_folder:
        if not os.path.isdir(folder):
            continue
        mzml_files = glob.glob(f"{folder}/*.{file_type}")
        peak_path = mzml_files[0]
        output_json = peak_path.replace(f".{file_type}", ".json")
        json_file_path = "config.json"
        missed_cleavages = missed_cleavages
        min_len = 6
        max_len = 30
        max_variable_mods = max_variable_mods
        static_mods = {"C": 57.021464}
        get_sage_config(
            json_file_path,
            peak_path,
            static_mods,
            PTMs,
            missed_cleavages,
            min_len,
            max_len,
            max_variable_mods,
            output_json,
        )
        fasta_files = glob.glob(f"{folder}/*.fasta")
        fasta_path = fasta_files[0]
        run_command(
            f"sage {output_json} --fasta {fasta_path} --write-pin --output_directory {folder}"
        )

else:
    big_folder = glob.glob(f"{folder_path}/*")
    for folder in big_folder:
        if not os.path.isdir(folder):
            continue
        mzml_files = glob.glob(f"{folder}/*.{file_type}")
        peak_path = mzml_files[0]
        output_json = peak_path.replace(f".{file_type}", ".json")
        json_file_path = "config.json"
        missed_cleavages = 2
        min_len = 6
        max_len = 30
        max_variable_mods = 5
        static_mods = {"C": 57.021464}
        new_mods = {"M": [15.994915], "N": [0.984016], "Q": [0.984016]}
        get_sage_config(
            json_file_path,
            peak_path,
            static_mods,
            new_mods,
            missed_cleavages,
            min_len,
            max_len,
            max_variable_mods,
            output_json,
        )
        fasta_files = glob.glob(f"{folder}/*.fasta")
        fasta_path = fasta_files[0]
        run_command(
            f"sage {output_json} --fasta {fasta_path} --write-pin --output_directory {folder}"
        )

# %%
# @title Brew Mokapot

import mokapot
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV
import numpy as np

# XGBoost schema from Fondrie & Noble 2021). A non-linear XGBoost seems to be better for rescoring open search results.

grid = {
    "scale_pos_weight": np.logspace(0, 2, 3),
    "max_depth": [1, 3, 6],
    "min_child_weight": [1, 10, 100],
    "gamma": [0, 1, 10],
}
xgb_mod = GridSearchCV(
    XGBClassifier(),
    param_grid=grid,
    n_jobs=1,
    cv=3,
    scoring="roc_auc",
)


def get_all_pin_files(folder_path):
    """Get all .pin files in the folder."""
    psm_files = []
    for root, _dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".pin"):
                full_path = os.path.join(root, file)
                psm_files.append(full_path)
    return psm_files


if use_Mokapot:
    folder_path = peak_folder

    if joint_modelling:
        psm_files = get_all_pin_files(peak_folder)
        if default_Percolator:
            svm = mokapot.PercolatorModel()
            psm_list = mokapot.read_pin(psm_files)
            results, models = mokapot.brew(psm_list, svm)
            result_files = results.to_txt(peak_folder)
        else:
            mod = mokapot.Model(xgb_mod)
            psm_list = mokapot.read_pin(psm_files)
            results, models = mokapot.brew(psm_list, mod)
            result_files = results.to_txt(peak_folder)

    else:
        big_folder = sorted(glob.glob(f"{folder_path}/*"))
        for folder in big_folder:
            if not os.path.isdir(folder):
                continue
            print(f"Processing folder: {folder}")
            pin_files = glob.glob(f"{folder}/*.pin")
            if not pin_files:
                print(f"No .pin files found in {folder}. Skipping...")
                continue
            pin = pin_files[0]
            if default_Percolator:
                svm = mokapot.PercolatorModel()
                psm_list = mokapot.read_pin(pin)
                results, models = mokapot.brew(psm_list, svm)
                result_files = results.to_txt(folder)
            else:
                mod = mokapot.Model(xgb_mod)
                psm_list = mokapot.read_pin(pin)
                results, models = mokapot.brew(psm_list, mod)
                result_files = results.to_txt(folder)

else:
    print("Mokapot not brewed")
