from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import BRICS
from itertools import chain, product as cartesian_product, permutations, combinations
import re
re_patt_for_submol = re.compile('\[99[A-Z,a-z][a-z]?[H]?[0-9]?\]')
re_patt=re.compile('\[[0-9]+\*\]')

def neutralize_nitrogen(mol):
    for a in mol.GetAtoms():
        totalHs = a.GetTotalNumHs()
        if all([a.GetSymbol() == 'N', a.GetFormalCharge()==1, totalHs>0]):
            a.SetFormalCharge(0)
            a.SetNumExplicitHs(totalHs-1)

def clean_smiles(smi):
    try:
        smi = smi.replace('[N](=O)O', '[N+](=O)[O-]').replace('[N+H2]', 'N')
        if '+' in smi:
            m = Chem.MolFromSmiles(smi)
            neutralize_nitrogen(m)
            smi = Chem.MolToSmiles(m)
        smi = Chem.CanonSmiles(smi)
    except:
        smi = 'ERROR'
    return smi

def demap(smi):
    mol = Chem.MolFromSmiles(smi)
    for a in mol.GetAtoms():
        a.SetIsotope(0)
        a.SetAtomMapNum(0)
    smi = Chem.CanonSmiles(Chem.MolToSmiles(mol))
    return smi

def get_clean_idx_to_map(cleaned_mol, mapped_smiles):
    mapped_mol = Chem.MolFromSmiles(mapped_smiles)
    match = mapped_mol.GetSubstructMatch(cleaned_mol)
    result = {}
    for cleaned_idx, mapped_idx in enumerate(match):
        mapped_val = mapped_mol.GetAtomWithIdx(mapped_idx).GetIsotope()
        result[cleaned_idx] = mapped_val
    return result

def prepare_mapped_smiles(mapped_smiles):
    mol = Chem.MolFromSmiles(mapped_smiles.replace('At', 'Br'))
    for a in mol.GetAtoms():
        if a.GetSymbol()=='H':
            a.SetIsotope(0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))

def prepare_inhibitor_data(mapped):
    mapped_smiles = prepare_mapped_smiles(mapped)
    mapped_smiles = Chem.CanonSmiles(mapped_smiles)
    clean_smiles = Chem.CanonSmiles(demap(mapped_smiles))
    inhibitor_data = {'mapped':mapped_smiles,
                      'clean':clean_smiles,
                      'clean_no_stereo': Chem.CanonSmiles(clean_smiles, useChiral=False),
                      'smarts': make_smarts_from_smiles(clean_smiles)}
    return inhibitor_data

def get_map_to_crd_idx(crd_mol, clean_mol, clean_idx_to_map):
    result = {}
    for clean_idx, crd_idx in enumerate(crd_mol.GetSubstructMatch(clean_mol)):
        if clean_idx not in clean_idx_to_map: 
            continue

        result[clean_idx_to_map[clean_idx]] = crd_idx
    return result

def get_coordinates_and_map(mol_path, smiles_data, debug=False):
    #inhibitor_data = {'mapped':mapped_inhibitor, 'clean':clean_inhibitor, 'smarts': make_smarts_from_smiles(clean_inhibitor)}
    prefix = ''

    clean_smiles = smiles_data['clean']
    clean_no_stereo = smiles_data['clean_no_stereo']
    smarts = smiles_data['smarts']
    
    with open(mol_path, 'r') as f:
        mol_crd = Chem.MolFromMolBlock(f.read())

    found = False
    for smiles, kind, status_val in [(clean_smiles, 'smi', 'exact'), (clean_no_stereo, 'smi', 'exact_no_stereo'), (smarts, 'sma', 'scaffold')]:
        mol = Chem.MolFromSmiles(smiles) if kind=='smi' else Chem.MolFromSmarts(smiles)
        if mol_crd.HasSubstructMatch(mol):
            status = prefix + status_val
            found = True
            break
    else:
        status = prefix + 'no_match'
    
    if found:
        clean_idx_to_map = get_clean_idx_to_map(mol, smiles_data['mapped'])
        map_to_crd_idx = get_map_to_crd_idx(mol_crd, mol, clean_idx_to_map)
    else:
        map_to_crd_idx = {}

    if not map_to_crd_idx or debug:
        print('MATCH???', smiles_data['mapped'], Chem.MolToSmiles(mol_crd))

    return mol_crd.GetConformer(0).GetPositions(), map_to_crd_idx, status

def get_isotopes(mol, ommit=['*']):
    result = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ommit:
            continue
        #n = a.GetAtomMapNum()
        n = a.GetIsotope()
        if n!=0:
            result.append(n)
    return result


def find_largest_intersecting_element(idx, list_of_sets, adjacency_dict):
    result, size = None, None
    this_element = list_of_sets[idx]
    for i, x in enumerate(list_of_sets):
        is_adjacent = any(x_atom in adjacency_dict[this_element_atom] for x_atom in x for this_element_atom in this_element) 
        if i==idx or not is_adjacent:
            continue
        x_size = len(x)
        if size is None or x_size>size:
            size = x_size
            result = i
    assert not (result is None), 'No intersection!'
    return result
            
def get_isotope_mapping_of_brics_components(mapped_mol):
    parts_smiles = [x for x in BRICS.BRICSDecompose(mapped_mol)]
    parts_smiles.sort()
    parts = [Chem.MolFromSmarts(re_patt.sub('[*]', x)) for x in parts_smiles]

    istotope_adj = {}
    for a in mapped_mol.GetAtoms():
        iso = a.GetIsotope()
        nei = [x.GetIsotope() for x in a.GetNeighbors()]
        istotope_adj[iso] = nei

    # print('ISOTOPE ADJ', istotope_adj)

    all_indices = set(get_isotopes(mapped_mol))
    block_indices = [set(get_isotopes(x)) for x in parts]
    indices_present = set.union(*block_indices)
    assert all_indices==indices_present, f'wrong fragments, {Chem.MolToSmiles(mapped_mol)}, {parts_smiles}'
    
    all_ok = False
    while not all_ok:
        if len(block_indices)<=1:
            break
        for test_idx, block in enumerate(block_indices):
            if len(block)>4:
                continue
            to_merge_idx = find_largest_intersecting_element(test_idx, block_indices, istotope_adj)
            merge_partner = block_indices[to_merge_idx]
            merger = block.union(merge_partner)
            block_indices = [merger] + [x for i,x in enumerate(block_indices) 
                                        if i not in [test_idx, to_merge_idx]]
            break
        else:
            all_ok = True
    return block_indices

def clean_mol(mol):
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
        atom.SetAtomMapNum(0)
        
def select_submol(smi, indices, mode='index'):
    functions = {'index':lambda atom: atom.GetIdx(),
                 'iso':lambda atom: atom.GetIsotope(),
                 'map':lambda atom: atom.GetAtomMapNum(),
                }
    func = functions[mode]
    mol = Chem.MolFromSmiles(smi)
    path = set()
    to_mark = []
    for bond in mol.GetBonds():
        atoms = bond.GetBeginAtom(), bond.GetEndAtom()
        if all(func(x) in indices for x in atoms):
            path.add(bond.GetIdx())
        nei_bonds = [[(atom, nei) for nei in atom.GetNeighbors()] for atom in atoms]
        for atom, nei in chain(*nei_bonds):
            if func(atom) not in indices or func(nei) in indices:
                continue
            to_mark.append(nei.GetIdx())
            bnd = mol.GetBondBetweenAtoms(atom.GetIdx(), nei.GetIdx())
            path.add(bnd.GetIdx())
    clean_mol(mol)
    for idx in to_mark:
        mol.GetAtomWithIdx(idx).SetIsotope(99)
    smi2 = Chem.MolToSmiles(Chem.PathToSubmol(mol, list(path)))
    return re_patt_for_submol.sub('C', smi2)

def modularize_ligand(smiles=None, mol=None, add_isotope_mapping=True):
    if mol ==None:
        assert smiles
        mol = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
    if add_isotope_mapping:
        for atom in mol.GetAtoms():
            atom.SetIsotope(atom.GetIdx()+1)
            atom.SetAtomMapNum(0)

    if not smiles:
        smiles = clean_smiles(Chem.MolToSmiles(Chem.RemoveHs(mol)))

    istotope_adj = {}
    for bond in mol.GetBonds():
        a,b = bond.GetBeginAtom().GetIsotope(), bond.GetEndAtom().GetIsotope()
        for x,y in [(a,b), (b,a)]:
            if x not in istotope_adj:
                istotope_adj[x] = []
            if y not in istotope_adj[x]:
                istotope_adj[x].append(y)

    mapped_block = Chem.MolToSmiles(mol)

    brics_indices = get_isotope_mapping_of_brics_components(Chem.RemoveHs(mol))
    brics_smiles = [select_submol(mapped_block, x, mode='iso') for x in brics_indices]
    
    brics_adj = {}
    for i, block1 in enumerate(brics_indices):
        if i not in brics_adj:
            brics_adj[i] = []
        for j, block2 in enumerate(brics_indices[:i]):
            if j not in brics_adj:
                brics_adj[j] = []
            connected = any(any(atom2_iso in istotope_adj[atom1_iso] for atom1_iso in block1)
                            for atom2_iso in block2)
            if connected:
                brics_adj[i].append(j)
                brics_adj[j].append(i)
    result = {'smiles':smiles, 'blocks':brics_smiles, 'block_adj':brics_adj,
              'mapped_smiles':mapped_block, 'block_maps':brics_indices}
    return result

