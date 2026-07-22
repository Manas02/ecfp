"""
Generate the ECFP input CSV from a SMILES string using RDKit.

The CSV is a square integer matrix:
    diagonal      = each atom's initial identifier (hash of its 7 invariants)
    off-diagonal  = bond orders (1 2 3 4 = single double triple aromatic)

For a chiral molecule an extra final row is appended, giving an (n+1) x n
matrix, holding a per-atom parity code (0 none, 1 = CIP R, 2 = CIP S).

    python smiles_to_csv.py "CCCC(=O)N" butyramide.csv
    python smiles_to_csv.py "C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O" --no-stereo thalidomide.csv
"""

import argparse
import numpy as np
from rdkit import Chem

HASH_MODULUS = 2147483647  # 2**31 - 1
HASH_BASE = 31
HASH_SEED = 1

BOND_ORDER = {
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.AROMATIC: 4,
}


def hash_integer_sequence(values):
    accumulator = HASH_SEED
    for value in values:
        accumulator = (accumulator * HASH_BASE + int(value)) % HASH_MODULUS
    return accumulator


def atom_invariants(atom):
    hydrogens = atom.GetTotalNumHs()
    return (
        atom.GetDegree(),
        atom.GetTotalValence() - hydrogens,
        atom.GetAtomicNum(),
        int(round(atom.GetMass())),
        atom.GetFormalCharge(),
        hydrogens,
        int(atom.IsInRing()),
    )


def smiles_to_csv(smiles, csv_path, use_stereo=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"could not parse SMILES: {smiles}")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    atom_count = mol.GetNumAtoms()

    matrix = np.zeros((atom_count, atom_count), dtype=np.int64)
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        matrix[index, index] = hash_integer_sequence(atom_invariants(atom))
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        matrix[i, j] = matrix[j, i] = BOND_ORDER[bond.GetBondType()]

    if use_stereo:
        code = {"R": 1, "S": 2}
        parity = [
            code.get(atom.GetPropsAsDict().get("_CIPCode"), 0)
            for atom in mol.GetAtoms()
        ]
        if any(parity):
            matrix = np.vstack([matrix, np.array([parity], dtype=np.int64)])

    np.savetxt(csv_path, matrix, fmt="%d", delimiter=",")
    return matrix


def main():
    parser = argparse.ArgumentParser(description="SMILES -> ECFP adjacency-matrix CSV")
    parser.add_argument("smiles")
    parser.add_argument("csv_path", nargs="?", default="molecule.csv")
    parser.add_argument(
        "--no-stereo",
        action="store_true",
        help="omit the parity row (achiral fingerprint input)",
    )
    args = parser.parse_args()

    matrix = smiles_to_csv(args.smiles, args.csv_path, use_stereo=not args.no_stereo)
    rows, cols = matrix.shape
    chiral = rows > cols
    print(
        f"{args.smiles} -> {args.csv_path}  ({cols} atoms, "
        f"{'chiral, parity row appended' if chiral else 'no parity row'})"
    )


if __name__ == "__main__":
    main()
