# Orthrus: an AI-powered, cloud-ready, and open-source hybrid approach for metaproteomics
---
![orthrus_workflow_v2](https://github.com/user-attachments/assets/9203545f-751b-4c49-b9b0-5c8a6e41de9b)

## Quick start
For cloud execution:<br/>
Click orthrus_cloud_stable_v100 folder, open .ipynb notebooks in Colab, and just follow the notebook instructions from there!
- please note: a Github account will be needed for authorisation. Go to [Github](https://github.com/) and get one if you haven't >_< 

For local execution:<br/>
```Python
git clone https://github.com/yc386/orthrus_metaproteomics.git
cd orthrus_metaproteomics/orthrus_local_runner
mamba env create -f environment.yaml
conda activate orthrus_metaproteomics
python walking_orthrus_locally_stable_v100.py --help
```

basic usage
```Python
python walking_orthrus_locally_stable_v100.py --folder_path path_to_folder --file_type mgf --use_SwissProt \
--sage_path path_to_sage_binary --json_file_path path_to_sage_json \
--missed_cleavages 1 --max_variable_mods 2 \
--static_CAM --aas P N Q --mods 15.994915 0.984016 0.984016 \
--default_Percolator --joint_modelling

```

## FAQs
1. Do I need a Google Colab subscription?<br>Without any Colab+, general CPU runtime, T4 GPU, and TPU should still be accessible.
2. Do I need a Google account? <br> May make life easier. Colab can also access files in your Google drive (with permission).
3. Help! Still unsure how to walk Orthrus. <br> Open an issue or get in touch with me ([preprint](https://www.biorxiv.org/content/10.1101/2024.11.15.623814v1))
