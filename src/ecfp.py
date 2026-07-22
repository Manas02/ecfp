"""
Extended-Connectivity Fingerprints (Rogers & Hahn, JCIM 2010) NumPy impl.

Input is a CSV holding a bond-weighted adjacency matrix:
    diagonal      = each atom's initial identifier (hash of its 7 invariants)
    off-diagonal  = bond orders (1 2 3 4 = single double triple aromatic)

A chiral molecule adds one final row, giving an (n+1) x n matrix, holding a
per-atom parity code (0 none, 1/2 for R/S). A square matrix is achiral;
enantiomers differ only in that row. Radius R corresponds to ECFP_{2R}.
"""

import argparse

import numpy as np
from collections import Counter

HASH_MODULUS = 2147483647
HASH_BASE = 31
HASH_SEED = 1


def hash_integer_sequence(values):
    accumulator = HASH_SEED
    for value in values:
        accumulator = (accumulator * HASH_BASE + int(value)) % HASH_MODULUS
    return accumulator


def load_molecule(path):
    table = np.atleast_2d(np.loadtxt(path, delimiter=",", dtype=np.int64))
    rows, cols = table.shape
    if rows == cols:
        matrix = table
        parity = np.zeros(cols, dtype=np.int64)
    elif rows == cols + 1:
        matrix = table[:-1, :]
        parity = table[-1, :].copy()
    else:
        raise ValueError(f"bad shape {table.shape}")
    identifiers = [int(x) for x in np.diag(matrix)]
    bond_orders = matrix.copy()
    np.fill_diagonal(bond_orders, 0)
    return identifiers, bond_orders, [int(x) for x in parity]


def build_neighbor_lists(bond_orders):
    n = bond_orders.shape[0]
    return [
        [
            (other, int(bond_orders[atom, other]))
            for other in range(n)
            if bond_orders[atom, other] != 0
        ]
        for atom in range(n)
    ]


def number_the_bonds(bond_orders):
    n = bond_orders.shape[0]
    bond_id_of = {}
    for i in range(n):
        for j in range(i + 1, n):
            if bond_orders[i, j] != 0:
                bond_id_of[(i, j)] = len(bond_id_of)
    return bond_id_of


def update_one_atom(
    atom,
    iteration_number,
    current_identifiers,
    neighbor_lists,
    parity_codes,
    already_disambiguated,
    use_stereo,
):
    hash_input = [iteration_number, current_identifiers[atom]]

    neighbor_pairs = sorted(
        (bond_order, current_identifiers[neighbor])
        for neighbor, bond_order in neighbor_lists[atom]
    )
    for bond_order, neighbor_identifier in neighbor_pairs:
        hash_input.append(bond_order)
        hash_input.append(neighbor_identifier)

    stereo_flag_added = False
    if use_stereo and parity_codes[atom] != 0 and not already_disambiguated[atom]:
        neighbor_identifiers = [
            current_identifiers[neighbor] for neighbor, _ in neighbor_lists[atom]
        ]
        if len(set(neighbor_identifiers)) == len(neighbor_identifiers):
            hash_input.append(int(parity_codes[atom]))
            stereo_flag_added = True

    return hash_integer_sequence(hash_input), stereo_flag_added


def grow_covered_bonds(atom, previous_cover, neighbor_lists, bond_id_of):
    grown = set(previous_cover[atom])
    for neighbor, _ in neighbor_lists[atom]:
        grown |= previous_cover[neighbor]
        low, high = (atom, neighbor) if atom < neighbor else (neighbor, atom)
        grown.add(bond_id_of[(low, high)])
    return frozenset(grown)


def generate_all_features(path, radius, use_stereo):
    initial_identifiers, bond_orders, parity_codes = load_molecule(path)
    atom_count = len(initial_identifiers)
    neighbor_lists = build_neighbor_lists(bond_orders)
    bond_id_of = number_the_bonds(bond_orders)

    current_identifiers = list(initial_identifiers)
    covered_bonds = [frozenset() for _ in range(atom_count)]
    already_disambiguated = [False] * atom_count

    feature_records = [
        (current_identifiers[atom], 0, covered_bonds[atom], atom)
        for atom in range(atom_count)
    ]

    for iteration_number in range(1, radius + 1):
        new_identifiers = [0] * atom_count
        new_covered_bonds = [None] * atom_count
        for atom in range(atom_count):
            new_identifiers[atom], stereo_added = update_one_atom(
                atom,
                iteration_number,
                current_identifiers,
                neighbor_lists,
                parity_codes,
                already_disambiguated,
                use_stereo,
            )
            new_covered_bonds[atom] = grow_covered_bonds(  # type: ignore
                atom, covered_bonds, neighbor_lists, bond_id_of
            )
            if stereo_added:
                already_disambiguated[atom] = True
        current_identifiers = new_identifiers
        covered_bonds = new_covered_bonds
        for atom in range(atom_count):
            feature_records.append(
                (current_identifiers[atom], iteration_number, covered_bonds[atom], atom)  # type: ignore
            )

    return feature_records


def remove_structural_duplicates(feature_records):
    survivors = []
    best_for_cover = {}
    for identifier, iteration_number, cover, central_atom in feature_records:
        if len(cover) == 0:
            survivors.append((identifier, iteration_number, cover, central_atom))
            continue
        candidate_rank = (iteration_number, identifier)
        incumbent = best_for_cover.get(cover)
        if incumbent is None or candidate_rank < (incumbent[1], incumbent[0]):
            best_for_cover[cover] = (identifier, iteration_number, cover, central_atom)
    survivors.extend(best_for_cover.values())
    return survivors


def build_fingerprint(path, radius=2, use_stereo=True, use_counts=False):
    feature_records = generate_all_features(path, radius, use_stereo)
    survivors = remove_structural_duplicates(feature_records)
    identifiers = [record[0] for record in survivors]
    return Counter(identifiers) if use_counts else set(identifiers)


def main():
    parser = argparse.ArgumentParser(description="ECFP fingerprint report from a CSV")
    parser.add_argument("--csv_path")
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--no-stereo", action="store_true")
    args = parser.parse_args()

    use_stereo = not args.no_stereo
    fingerprint = sorted(build_fingerprint(args.csv_path, args.radius, use_stereo))
    stereo_flag = 1 if use_stereo else 0
    ids = " ".join(str(x) for x in fingerprint)
    print(f"py←{ids}")
    print(
        f"py ≡ ({args.radius} {stereo_flag}) ECFP '{args.csv_path}'   ⍝ 1 if identical"
    )


if __name__ == "__main__":
    main()
