#!/usr/bin/env python
# Quick and dirty parallel use of standardize_mol.
import os
import sys
import concurrent.futures as cf

from ruse.rdkit.rdkit_utils import standardize_mol

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

infile = sys.argv[1]
outfile = sys.argv[2]

fix_azide = AllChem.ReactionFromSmarts('[N-:1]=[N:2]#[N+:3]>>[N+0:1]#[N+:2][N-:3]')
normer = rdMolStandardize.Normalizer()
uncharger = rdMolStandardize.Uncharger()
metal_disconnector = rdMolStandardize.MetalDisconnector()


def process_smiles(inline):
    insmi, mol_name = inline.split()
    mol = Chem.MolFromSmiles(insmi, sanitize=False)
    if mol is None or not mol or not mol.GetNumAtoms():
        return None, mol_name
    norm_mol = standardize_mol(mol, normer, uncharger, metal_disconnector, fix_azide)
    if norm_mol is None:
        return None, mol_name
    else:
        return Chem.MolToSmiles(norm_mol), mol_name


def main():
    with open(infile, 'r') as inf, open(outfile, 'w') as outf:
        with cf.ProcessPoolExecutor(max_workers=os.cpu_count() - 2) as pool:
            for out_smi, mol_name in pool.map(process_smiles, inf):
                if out_smi is not None:
                    outf.write(f'{out_smi} {mol_name}\n')
                else:
                    print(f'{mol_name} is None')


if __name__ == '__main__':
    main()
