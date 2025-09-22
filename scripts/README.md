# README for scripts

This folder contains scripts that are supposed to run sequentially. They are organised in three steps:

# Step 000
Follow the instructions in `docs/install_ChEMBL.md` to create a local copy of the ChEMBL DB. Currently this is set to ChEMBL_36 (latest version), but if a new release appears simply download the new one and change the DB name in `default.py`.
A folder inside `config` (named `chembl_activities`) will be created containing unmodified data extracted directly from ChEMBL tables.

# Step 001
Include information in the ChEMBL tables automatically using an LLM for processing activity comments and others. A series of modified files will be saved in `config`, each having a column that represents the outcome of the assay/activity (1: Active, -1: Not Active, 0: Unresolved or Direction: higher value == higher activity --> 1, lower value == higher activity --> -1, inconclusive --> 0 ).
For this step, it is necessary to download the Gemma Llamafile from [Hugging Face](https://huggingface.co/Mozilla/gemma-2-9b-it-llamafile).

# Step 002
Uses the LLM processed data to assign an outcome to each activity in ChEMBL (which are summarised in the Activities table).