# Requirements

Conda environment file: environment.yml

Required packages:

* rdkit
* typer
* tqdm
* pandas
* numpy
* (opencadd-klifs)[https://opencadd.readthedocs.io/en/latest/databases_klifs.html]

# Key step 1: download structures from KLIFS

`python main.py download-all-klifs-structures klifs_data.csv structure_data/`

Optionally, one may start with just first 100 structures.

`python main.py download-all-klifs-structures klifs_data.csv structure_data/ --smoke-test`

# Blockization example

`python main.py blockization-example klifs_data.csv`

# Interaction dataframe example

`python main.py interaction-example klifs_data.csv example_interaction_df.csv`

