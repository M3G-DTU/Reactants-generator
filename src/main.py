# Main function for the workflow of this module.
# 1 - Get the descriptors of two molecules.
# 2 - Filter the descriptors according to some criterias
# 3 - Create atom pairs based on the filtered descriptors
# 4 - Create new initial geometries for the reactants based on the atom pairs

# Import the necessary classes and functions
from Descriptors.Filtering.filter_descriptors import get_atom_pairs
from ArrangeMolecules.rearange import OrientMoleculeSphere
import argparse
from ase.io import read
import os

def from_AMS_out_to_xyz(file_path: str, output_path: str = "temp.xyz"):
    with open(file_path, "r", encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "SINGLE POINT CALCULATION" in line:
            start = i + 9
            break
    symbols = []
    xyz = []
    for m in range(start, len(lines)):
        # Stop at empty line
        if lines[m].strip() == "":
            break
        ap = lines[m].split()
        symbols.append(ap[1])
        xyz.append((float(ap[2]), float(ap[3]), float(ap[4])))
    # Make temporary xyz file
    with open(output_path, "w") as f:
        f.write(f"{len(symbols)}\n")
        f.write("Temporary xyz file\n")
        for symbol, (x, y, z) in zip(symbols, xyz):
            f.write(f"{symbol} {x} {y} {z}\n")

def argument_parser():
    parser = argparse.ArgumentParser(description='Generate reactant geometries based on descriptors.')
    parser.add_argument('molecule1', type=str, help='Path to the output file of molecule 1')
    parser.add_argument('molecule2', type=str, help='Path to the output file of molecule 2')
    parser.add_argument('--method', type=str, default='Hirshfeld', choices=['FMO', 'Hirshfeld'], help='Method to use for descriptor calculation and filtering')
    parser.add_argument('--filter_descriptor', type=str, default='dual', choices=['fukui', 'dual', 'ElNuc'], help='Descriptor to use for filtering atoms')
    parser.add_argument('--max_pairs_per_molecule', type=int, default=5, help='Maximum number of atom pairs to generate per molecule')
    parser.add_argument('--include_hydrogens', action='store_false', help='Whether to include hydrogens in the atom pairs', default=True)
    parser.add_argument('--output_dir', type=str, default='output', help='Directory to save the generated geometries')
    return parser.parse_args()

def main():
    args = argument_parser()
    atom_pairs = get_atom_pairs(args.molecule1, args.molecule2, method=args.method, filter_descriptor=args.filter_descriptor, max_pairs_per_molecule=args.max_pairs_per_molecule, include_hydrogens=args.include_hydrogens)
    print('Atom pairs (mol1_idx, mol2_idx):')
    print(atom_pairs)
    # Create temporary xyz files for the two molecules
    from_AMS_out_to_xyz(args.molecule1, "temp1.xyz")
    from_AMS_out_to_xyz(args.molecule2, "temp2.xyz")
    molecule1 = read("temp1.xyz")
    molecule2 = read("temp2.xyz")
    for atom1_idx, atom2_idx in atom_pairs:
        orientor = OrientMoleculeSphere(molecule1, molecule2, atom1_idx, atom2_idx)

        # Optimize orientation
        result = orientor.optimize()

        # Merge and save the new geometry
        merged_xyz = orientor.molecule1 + orientor.molecule2
        with open(f"{args.output_dir}/reactant_{orientor.molecule1[atom1_idx].symbol}{atom1_idx}_{orientor.molecule2[atom2_idx].symbol}{atom2_idx}.xyz", "w") as f:
            f.write(f"{len(merged_xyz)}\n")
            f.write(f"Reactant geometry for atom pair ({atom1_idx}, {atom2_idx})\n")
            for atom in merged_xyz:
                f.write(f"{atom.symbol[0]} {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}\n")
    # Remove temporary xyz files
    os.remove("temp1.xyz")
    os.remove("temp2.xyz")
if __name__ == "__main__":
    main()