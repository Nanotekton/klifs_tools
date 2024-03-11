import typer
from klifs_tools import klifs_structure_iterator, get_subpocket_to_residue_map, download_structure, get_subpocket_crd_and_insert_marks, get_block_crds, assign_block_center, KinFragLibSubPockets, get_interaction_data
from opencadd.databases.klifs import setup_remote
from rdkit import Chem
import pandas as pd
import traceback

app = typer.Typer()

@app.command()
def download_all_klifs_structures(summary_file:str, directory:str, smoke_test:bool=False):
    '''donwloads all KinFragLib-related structures from KLIFS'''

    result = []
    session = setup_remote()

    for structure_dict in klifs_structure_iterator(session=session):
        subpocket_residues = get_subpocket_to_residue_map(structure_dict['StructureId'], session=session)
        if not subpocket_residues:
            continue
        try:
            pdb_path, mol_path = download_structure(structure_dict['StructureId'], directory, session=session)
        except:
            print(traceback.format_exc())
            continue
        if mol_path == None:
            continue
        structure_dict['pdb_path'] = pdb_path
        structure_dict['mol_path'] = mol_path
        
        subpocket_crd = get_subpocket_crd_and_insert_marks(subpocket_residues, pdb_path)
        for k, v in subpocket_crd.items():
            for name, val in zip('xyz', v):
                structure_dict[f'{k}.{name}'] = val


        result.append(structure_dict)

        if smoke_test and len(result)>=100:
            break

    print(f'downloaded: {len(result)} structures')
    pd.DataFrame(result).to_csv(summary_file, index=False, sep=';')


@app.command()
def blockization_example(summary_file:str, index:int=0):
    d = pd.read_csv(summary_file, sep=';')
    mol_path = d.iloc[index]['mol_path']
    print(d['pdb_path'].iloc[0])
    pockets = {}
    for k in KinFragLibSubPockets:
        v = []
        for x in 'xyz':
            v.append(d[f'{k}.{x}'].iloc[index])
        pockets[k] = v

    with open(mol_path, 'r') as f:
        mol = Chem.MolFromMolBlock(f.read(), removeHs=False)

    print(Chem.MolToSmiles(Chem.RemoveHs(mol)), '\n')

    smiles_and_crd = get_block_crds(mol)
    for smiles, crd in smiles_and_crd:
        print(smiles, assign_block_center(pockets, crd))


@app.command()
def interaction_example(summary_file:str, output_file:str, index:int=0):
    d = pd.read_csv(summary_file, sep=';')
    mol_path = d.iloc[index]['mol_path']
    print(d['pdb_path'].iloc[0])
    structure_id = d.iloc[0]['StructureId']

    interaction = get_interaction_data(structure_id)
    print(interaction)
    interaction.to_csv(output_file, sep=';', index=False)

if __name__ == "__main__":
    app()
