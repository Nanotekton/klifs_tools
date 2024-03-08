from rdkit import Chem
import json
from collections import Counter
import pandas as pd
import numpy as np
from opencadd.databases.klifs import setup_remote
import os
import traceback
import time
from datetime import timedelta
from tqdm.auto import tqdm
from smiles_tools import  modularize_ligand

#===== KinFragLib Definitions ================
KinFragLibSubPockets = {
    'AP':[15, 46, 51, 75],
    'SE':[51],
    'FP':[10, 51, 72, 81],
    'GA':[17, 45, 81],
    'B1':[28, 38, 43, 81],
    'B2':[18, 24, 70, 83]
     }
KinFragLibResidues = sum(KinFragLibSubPockets.values(), [])

KinFragLibSubPocketsElements = {
    'AP':'C',
    'SE':'N',
    'FP':'O',
    'GA':'S',
    'B1':'P',
    'B2':'P'
     }

EXO_ALLOWED = [('SE', 'SE'), ('B1', 'B1'), ('B2', 'B2')]

KinFragLibConn = [
    ('SE', 'FP'),
    ('SE', 'AP'),
    ('AP', 'FP'),
    ('AP', 'GA'),
    ('FP', 'GA'),
    ('GA', 'B1'),
    ('GA', 'B2'),
    ('B1', 'B2')
]

KinFragLibAdj = {}
for a,b in KinFragLibConn:
    for center, other in [(a,b), (b,a)]:
        if center not in KinFragLibAdj:
            KinFragLibAdj[center] = []
        if other not in KinFragLibAdj[center]:
            KinFragLibAdj[center].append(other)

def check_allowed_connectivity(keys):
    to_check = [x for x in keys]
    leafs, remaining = to_check[:1], to_check[1:]
    visited = [x for x in leafs]
    while leafs:
        new_leafs = []
        for leaf in leafs:
            for nei in KinFragLibAdj[leaf]:
                if nei in visited:
                    continue
                visited.append(nei)
                if nei in remaining:
                    new_leafs.append(nei)
                    remaining.remove(nei)
        leafs = new_leafs
        # print(leafs, remaining)

    return remaining==[]


def is_allowed_pair(pocket1, pocket2):
    return (pocket1, pocket2) in KinFragLibConn or (pocket2, pocket1) in KinFragLibConn
#==========================

def check_pocket_is_complete(pocket):
    '''Checks whether all residues in KinFragLib definition are present in the structure'''
    all_klifs_residues_present_in_df = all([x in pocket['residue.klifs_id'].values for x in KinFragLibResidues])
    all_klifs_residues_present_in_pdb = not pocket.loc[pocket['residue.klifs_id'].isin(KinFragLibResidues), 'residue.id'].isna().any()
    all_klifs_residues_are_ok = not (pocket.loc[pocket['residue.klifs_id'].isin(KinFragLibResidues), 'residue.id']=='_').any()
    return all([all_klifs_residues_present_in_df, all_klifs_residues_present_in_pdb, all_klifs_residues_are_ok])

def get_subpocket_to_residue_map(structure_klifs_id, session=None):
    if session is None:
        session = setup_remote()
    
    pocket = session.pockets.by_structure_klifs_id(structure_klifs_id)
    if not check_pocket_is_complete(pocket):
        return {}
    
    result = {}
    for subpocket, subpocket_klifs_residues in KinFragLibSubPockets.items():
        mask = pocket['residue.klifs_id'].isin(subpocket_klifs_residues)
        residues_in_pdb = [int(x) for x in pocket.loc[mask, 'residue.id'].values]
        result[subpocket] = residues_in_pdb
    return result

def download_structure(structure_klifs_id, directory, session=None, ligand_name='ligand'):
    '''Downloads PDB and SDF (ligand) files of a given structure'''
    if not os.path.isdir(directory):
        os.popen(f'mkdir {directory}').read()

    target_dir = f'{directory}/Klifs_{structure_klifs_id}'
    if not os.path.isdir(target_dir):
        os.popen(f'mkdir {target_dir}').read()

    if session is None:
        session = setup_remote()

    mol = session.coordinates.to_rdkit(structure_klifs_id, entity='ligand', compute2d=False)
    try:
        block = Chem.MolToMolBlock(mol, confId=0)
    except:
        return None, None

    pdb_path = str(session.coordinates.to_pdb(structure_klifs_id, target_dir))
    mol_path = f'{target_dir}/{ligand_name}.sdf'
    with open(mol_path, 'w') as f:
        f.write(block)

    return pdb_path, mol_path


def parse_pdb_line(line):
    '''converts PDB line (section ATOM/HETATM} into dictionary'''
    result = dict(
        section = line[:6].strip(),
        num = int(line[6:11].strip()),
        atom = line[11:17].strip(),
        resn = line[17:20].strip(),
        chain = line[21],
        resid = int(line[22:26]),
        x = float(line[30:38]),
        y = float(line[38:46]),
        z = float(line[46:54]),
        occ = float(line[54:60]),
        beta = float(line[60:66])
    )
    return result


def make_pdb_line(idx, atname, resname, resid, x, y, z, section='ATOM', chain='A', occ=1.0, beta=0.0, element=''):
    return f'{section:<6s}{idx:>5d}  {atname:<4s}{resname:3s} {chain:1s}{resid:>4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{beta:6.2f}{"":10s}{element:>2s}\n'


def read_ca_center(residues, pdblines):
    '''returns geometric center of residues' CAs '''
    #HETATM 2606  N19 MI1 A1001      -0.689  17.510  43.856  1.00  0.00           N
    #HETATM 2710  N34 Q3E A1401      -3.848  13.419  49.416  1.00  0.00
    #ATOM   2409  CA  GLY A 946      16.674  10.665  48.134  1.00  0.00
    #012345678901234567890123456789012345678901234567890123456789012345678901234567890
    #          1         2         3         4         5         6         7
    vectors = []
    for line in pdblines:
        if line[:4]!='ATOM':continue
        atom_type = line[11:17].strip()
        if atom_type!='CA':continue
        resid = int(line[22:26].strip())
        if resid not in residues: continue
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        vectors.append(np.array([x,y,z]))
    assert len(vectors)==len(residues), f'too many chains?, {residues}'
    return np.mean(vectors, axis=0)

def simplify_structure_record(structure_row):
    '''Puts selected fields into a dictionary'''
    structure_id = structure_row['structure.klifs_id']
    ligand_id = structure_row['ligand.expo_id']
    pdb_code = structure_row['structure.pdb_id']

    return dict(DFG=structure_row['structure.dfg'], ResName=ligand_id, PDB_code=pdb_code, StructureId=structure_row['structure.klifs_id'])

def klifs_structure_iterator(batchsize:int=100, use_tqdm=True, session=None):
    '''Collects all KLIFS ligands and then iterates over corresponding structures'''
    if session is None:
        session = setup_remote()
    
    all_ligands = session.ligands.all_ligands()

    def try_batch(batch):
        try:
            result = session.structures.by_ligand_expo_id(batch).drop_duplicates('structure.klifs_id')
            result = (x[1] for x in result.iterrows())
        except:
            result = []
        return result
        
    ligand_iterator = all_ligands['ligand.expo_id'].values
    if use_tqdm:
        tqdm.pandas()
        ligand_iterator =  tqdm(ligand_iterator)

    batch = []
    for ligand_id in ligand_iterator:
        batch.append(ligand_id)
        if len(batch)%batchsize==0:
            for structure_record in try_batch(batch):
                yield simplify_structure_record(structure_record)
            batch = []
    if batch:
        for structure_record in try_batch(batch):
            yield simplify_structure_record(structure_record)


def get_subpocket_crd_and_insert_marks(subpocket_to_residue_map, pdb_file_path):
    '''Creates special residue - KFL - with atoms depicting subpockets. Returns coordinates of subpockets and writes modified PDB file'''
    with open(pdb_file_path, 'r') as f:
        pdblines = f.readlines()
    idx = None
    for i,x in enumerate(pdblines):
        if 'ATOM' in x or 'HETATM' in x:
            idx = i
    assert idx
    try:
        last_line = parse_pdb_line(pdblines[idx])
    except:
        print(pdblines[idx])
        raise
    idx += 1
    num, resid = last_line['num'], last_line['resid']
    resid += 1
    num += 1
    added_lines = []

    subpocket_crd = {}

    for key, residues in subpocket_to_residue_map.items():
        try:
            vector = read_ca_center(residues, pdblines)
        except:
            print(residues, pdb_file_path)
            raise
        x, y, z = vector
        subpocket_crd[key] = vector
        added_lines.append(make_pdb_line(num, key, 'KFL', resid, x, y, z, section='HETATM', element=KinFragLibSubPocketsElements[key]))
    
    pdblines = pdblines[:idx] + added_lines + pdblines[idx:]
    with open(pdb_file_path.replace('.pdb', '_with_subpockets.pdb'), 'w') as f:
        f.write(''.join(pdblines))

    return subpocket_crd


#========== Assignment ===============

def get_block_crds(mol_with_crd):
    module_data = modularize_ligand(mol=mol_with_crd, add_isotope_mapping=True)
    crd = mol_with_crd.GetConformer(0).GetPositions()
    result = []

    for block_smiles, block_indices in zip(module_data['blocks'], module_data['block_maps']):
        block_indices = [i-1 for i in block_indices]
        block_crd = crd[block_indices]
        result.append((block_smiles, block_crd))

    return result

# #assign block
def assign_block_center(subpocket_crd, block_crd, th=8):
    if len(block_crd.shape)==2:
        center = block_crd.mean(axis=0)
    else:
        center = block_crd
    assert len(center)==3
    results = []
    for pocket, crd in subpocket_crd.items():
        dist = crd - center
        dist = dist.dot(dist)**0.5
        results.append((pocket, dist))
    results.sort(key=lambda x:x[1])
    pocket, dist = results[0]
    if dist>th:
        pocket = 'X'
    return pocket, dist


def get_interaction_data(structure_klifs_id, session=None):
    if session==None:
        session = setup_remote()
    interactions = session.interactions.by_structure_klifs_id(structure_klifs_id)
    structure = session.pockets.by_structure_klifs_id(structure_klifs_id)
    if len(interactions)==0 or len(structure)==0:
        return None
    
    fgp = interactions['interaction.fingerprint'].iloc[0]
    assert len(fgp)/7==len(structure)
    
    all_interactions = session.interactions.interaction_types['interaction.name'].values

    assert len(all_interactions) == 7

    result = []
    i = 0
    for _, row in structure.iterrows():
        klifs_id, pdb_id = row['residue.klifs_id'], row['residue.id']
        
        block = fgp[i*7:(i+1)*7]
        assert len(block)==7
        for name, flag in zip(all_interactions, block):
            if flag=='1':
                result.append({'residue_klifs_id':klifs_id, 'residue_pdb_id':pdb_id, 'interatcion_type':name})

        i += 1

    return pd.DataFrame(result)



