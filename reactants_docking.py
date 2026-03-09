# Import relevant libraries
import os
import math
#import pkg_resources
#from pkg_resources import DistributionNotFound, VersionConflict
import argparse

from molrearr.face2face import *
from molrearr.drone_move import *
from molrearr.quality_check import *
from molrearr.optimization import *
from molrearr.read_xyz2list import *

# Argument Parser
def argument_parser():
    parser = argparse.ArgumentParser(description='A script to perform reactant docking based on the output of FDL or FMO calculations in AMS.')
    parser.add_argument('-o', '--out_files', type=str, required=True, nargs=2, 
                        help='the AMS output files')
    parser.add_argument('--method', type=str, required=True, choices=['FMO', 'Hirshfeld', 'Mulliken', 'Voronoi'],
                        help='the method used for partitioning (FMO, Hirshfeld, Mulliken, or Voronoi)')
    parser.add_argument('-xyz', '--xyz_files', type=str, nargs=2,
                        help='the XYZ files of the two molecules. If not provided, it is expected that the XYZ files with the same name as the output files but with .xyz extension are in the same directory as the output files')
    parser.add_argument('-hs', '--hotspot_distance', type=float, default=2.2, 
                        help='the distance threshold for identifying hotspots (in angstroms)')
    parser.add_argument('-md', '--minimal_distance', type=float, default=1.6,
                        help='the minimal distance between the two atoms for different fragments (in angstroms)')
    parser.add_argument('-nc', '--n_candidates', type=int, default=0,
                        help='the number of top candidates to keep (if 0, keep all candidates). Can also be a percentage (e.g., 0.5 to keep the top 50%% of candidates)')
    parser.add_argument('-ah', '--add_aditional_H', type=int, default=None,
                        help='number of additional H atoms to add. If None, no additional H atoms are added')
    parser.add_argument('-ecut', '--elec_cutoff', type=float, default=10.0, 
                        help='the cutoff value for elec')
    parser.add_argument('-ncut', '--nucl_cutoff', type=float, default=10.0,
                        help='the cutoff value for nucl interactions')
    parser.add_argument('--output_dir', type=str, default='Output', 
                        help='the directory to save the output files')
    parser.add_argument('-rd', '--restrict_distance', type=float, default=2.2,
                        help='the distance threshold for restricting movements (in angstroms)')
    parser.add_argument('-sfl', '--strict_filter_level', type=int, default=3,
                        help='the level of strict filter >=3, higher is stricter')
    parser.add_argument('-d', '--descriptor', type=str, default='dual', choices=['fukui', 'dual', 'ElNuc'],
                        help='the descriptor to use for filtering (fukui, dual or ElNuc)')
    parser.add_argument('--fukui_cutoff', type=float, default=0.02,
                        help='the cutoff value for the Fukui function when using the fukui descriptor (default: 0.02)')
    parser.add_argument('--info_file', type=str, default='candidates.info',
                        help='the name of the info file to save the candidates and settings information (default: candidates.info)')
    return parser.parse_args()

#functions
def molecule_info_fdl(out_file: str) -> dict:
    """Extracts condensed local electrophilicity and nucleophilicity information from the output file of FDL calculation.
    \nparams:
    \n    out_file: the AMS output file from FDL calculation
    \nreturns:
    \n    atom_properties: a dictionary containing relevant information from the AMS output"""
    # Locate relevant sections in the AMS output file
    locate_sections = ['Hirshfeld Partitioning', 
                       'Mulliken Partitioning', 
                       'Voronoï Partitioning',
                       ]
    # Create empty dictionary to store the extracted information
    atom_properties = {}
    # Add the section names as keys to the dictionary
    for section in locate_sections:
        atom_properties[section] = {}

    with open(out_file, "r", encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    
    # Store the the start index of each section
    start_indencies = {}

    for i, line in enumerate(lines):
        if any(section in line for section in locate_sections):
            for section in locate_sections:
                if section in line:
                    start_indencies[section] = i
    
    # Now we can extract all information for each section
    for section, start_index in start_indencies.items():
        # Go to the start of the relevant section
        for i in range(start_index, len(lines)):
            if "MAIN ATOMIC DESCRIPTORS: CANONICAL ENSEMBLE" in lines[i]:
                start = i + 3
                break
        # Save the headers
        headers = lines[start].split()
        # Run until multiple "-" are found, which indicates the end of the section
        for m in range(start + 1, len(lines)):
            if "----" in lines[m]:
                break
            # Else add information to the dictionary
            ap = lines[m].split()
            for header, value in zip(headers, ap):
                if header not in atom_properties[section]:
                    atom_properties[section][header] = []
                if header == 'Atom':
                    atom_properties[section][header].append(value)
                else:
                    atom_properties[section][header].append(float(value))
        # Now we locate the Local electrophilicity and nucleophilicity information
        for n in range(m, len(lines)):
            if "CONDENSED LOCAL ELECTROPHILICITY AND NUCLEOPHILICITY" in lines[n]:
                start = n + 5
                break
        headers = ['Electrophilicity', 'Nucleophilicity']
        for p in range(start, len(lines)):
            if "----" in lines[p]:
                break
            ap = lines[p].split()[1:]
            for header, value in zip(headers, ap):
                if header not in atom_properties[section]:
                    atom_properties[section][header] = []
                    
                atom_properties[section][header].append(float(value))
    # Remove the "Partitioning" part to make the keys more concise
    for section in locate_sections:
        section_key = section.replace(" Partitioning", "")
        atom_properties[section_key] = atom_properties.pop(section)
    return atom_properties

def molecule_info_fmo(out_file: str) -> dict:
    """Extracts frontier molecular orbital information from the output file of FMO calculation.
    \nparams:
    \n    out_file: the AMS output file from FMO calculation
    \nreturns:
    \n    fmo_properties: a dictionary containing relevant information from the AMS output"""
    atom_properties = {'FMO': {}}
    with open(out_file, "r", encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "ATOMIC DESCRIPTORS: CANONICAL ENSEMBLE" in line:
            start = i + 3
            break
    headers = lines[start].split()
    for m in range(start + 2, len(lines)):
        if "----" in lines[m]:
            break
        ap = lines[m].split()
        # We need to make it into the right format to be able to use it later
        ap[0] = f"{ap[1]}({ap[0]})"
        # Remove indencies 1 and 2
        ap.pop(1)
        ap.pop(1)
        for header, value in zip(headers, ap):
            if header not in atom_properties['FMO']:
                atom_properties['FMO'][header] = []
            if header == 'Atom':
                atom_properties['FMO'][header].append(value)
            else:
                atom_properties['FMO'][header].append(float(value))
    headers = ['Electrophilicity', 'Nucleophilicity']
    for n in range(m, len(lines)):
        if "CONDENSED LOCAL ELECTROPHILICITY AND NUCLEOPHILICITY" in lines[n]:
            start = n + 5
            break
    for p in range(start, len(lines)):
        if "----" in lines[p]:
            break
        ap = lines[p].split(':')[1].split()
        for header, value in zip(headers, ap):
            if header not in atom_properties['FMO']:
                atom_properties['FMO'][header] = []
            atom_properties['FMO'][header].append(float(value))
    return atom_properties

# filter
def get_sorted_indexes(directory: dict, key_name: str, reverse: bool = False) -> list:
    return sorted(
        range(len(directory[key_name])),
        key=lambda i: directory[key_name][i],
        reverse=reverse
    )

def elec_nuc_filter(atom_properties: dict, method: str = 'Hirshfeld', filter_type: str = 'elec') -> list:
    """Filters the atoms based on their electrophilicity and nucleophilicity values.
    \nparams:
    \n    atom_properties: a dictionary containing relevant information from the AMS output
    \n    method: the method used for partitioning (Hirshfeld, Mulliken, or Voronoï)
    \n    filter_type: the type of filter to apply (elec or nuc)
    \nreturns:
    \n    filtered_atoms: a list of dictionaries containing the filtered atoms and their properties"""
    
    filtered_atoms = []
    if filter_type == 'elec':
        # Check if atom is electrophilic
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'Electrophilicity', reverse=True)
        # We loop through each atom
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['Electrophilicity'] > 0 and atom['Electrophilicity'] + atom['Nucleophilicity'] > 0:
                filtered_atoms.append(atom)
    elif filter_type == 'nuc':
        # Check if atom is nucleophilic
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'Nucleophilicity', reverse=False)
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['Nucleophilicity'] < 0 and atom['Electrophilicity'] + atom['Nucleophilicity'] <= 0:
                filtered_atoms.append(atom)
    else:
        raise ValueError("Invalid filter type. Please choose 'elec' or 'nuc'.")
    return filtered_atoms

def dual_filter(atom_properties: dict, method: str = 'Hirshfeld', filter_type: str = 'elec') -> list:
    """Filters the atoms based on their dual descriptor values.
    \nparams:
    \n    atom_properties: a dictionary containing relevant information from the AMS output
    \n    method: the method used for partitioning (Hirshfeld, Mulliken, or Voronoï)
    \n    filter_type: the type of filter to apply (elec or nuc)
    \nreturns:
    \n    filtered_atoms: a list of dictionaries containing the filtered atoms and their properties"""
    
    filtered_atoms = []
    if filter_type == 'elec':
        # Check if atom is electrophilic
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'Electrophilicity', reverse=True)
        # We loop through each atom
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['f(2)'] > 0:
                filtered_atoms.append(atom)
    elif filter_type == 'nuc':
        # Check if atom is nucleophilic
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'Nucleophilicity', reverse=False)
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['f(2)'] < 0:
                filtered_atoms.append(atom)
    else:
        raise ValueError("Invalid filter type. Please choose 'elec' or 'nuc'.")
    return filtered_atoms

def fukui_filter(atom_properties: dict, method: str = 'Hirshfeld', filter_type: str = 'elec', fukui_cutoff: float = 0.02) -> list:
    """Filters the atoms based on their Fukui function values.
    \nparams:
    \n    atom_properties: a dictionary containing relevant information from the AMS output
    \n    method: the method used for partitioning (Hirshfeld, Mulliken, or Voronoï)
    \n    filter_type: the type of filter to apply (elec or nuc)
    \n    fukui_cutoff: the cutoff value for the Fukui function
    \nreturns:
    \n    filtered_atoms: a list of dictionaries containing the filtered atoms and their properties"""
    
    filtered_atoms = []
    if filter_type == 'elec':
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'f+', reverse=True)
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['f+'] > fukui_cutoff:
                filtered_atoms.append(atom)
    elif filter_type == 'nuc':
        sorted_indexes = get_sorted_indexes(atom_properties[method], 'f-', reverse=True) # Reverse as the fukui functions are not negative
        for i in sorted_indexes:
            atom = {key: atom_properties[method][key][i] for key in atom_properties[method]}
            if atom['f-'] > fukui_cutoff:
                filtered_atoms.append(atom)
    else:
        raise ValueError("Invalid filter type. Please choose 'elec' or 'nuc'.")
    return filtered_atoms

def keep_top_candidates(filtered_atoms: list, n_candidates=0, add_aditional_H=None) -> list:
    """Keeps only the top candidates based on filtered atoms.
    \nparams:
    \n    filtered_atoms: a list of dictionaries containing the filtered atoms and their properties
    \n    n_candidates: the number of top candidates to keep (if 0, keep all candidates). Can also be a percentage (e.g., 0.5 to keep the top 50% of candidates).
    \n    add_aditional_H: number of additional H atoms to add. If None, no additional H atoms are added.
    \nreturns:
    \n    chosen_candidates: a list of dictionaries containing the chosen candidates and their properties"""
    chosen_candidates = []
    # First we check if n_candidates is a percentage or an integer
    if isinstance(n_candidates, float) and 0 < n_candidates < 1:
        p_index = math.ceil(n_candidates * len(filtered_atoms))
        chosen_candidates = filtered_atoms[:p_index]
    elif isinstance(n_candidates, int) and n_candidates > 0:
        chosen_candidates = filtered_atoms[:n_candidates]
    elif n_candidates == 0:
        chosen_candidates = filtered_atoms
    else:
        raise ValueError("Invalid value for n_candidates. Please provide a non-negative integer or a float between 0 and 1.")
    if add_aditional_H:
        # Add aditional hydrogen atoms to the candidates
        # Create new list to store all H atoms
        H_atoms = []
        for atom in filtered_atoms:
            if atom['Atom'].split('(')[0].startswith('H'):
                H_atoms.append(atom)
        # Remove H atoms already in chosen candidates
        H_atoms = [atom for atom in H_atoms if atom not in chosen_candidates]
        # Add additional H atoms to the chosen candidates
        chosen_candidates.extend(H_atoms[:add_aditional_H])
    return chosen_candidates

def xyz_from_out(outfile: str) -> list:
    """Generates the XYZ list necessary for the rest of the script from a FDL or FMO calculation in AMS."""
    with open(outfile, "r", encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "* SINGLE POINT CALCULATION *" in line:
            start = i + 9
            break
    xyz = []
    atoms = []
    for m in range(start, len(lines)):
        # If line is empty, we have reached the end of the XYZ section
        if lines[m].strip() == "":
            break
        # Add the line to the XYZ lines
        Axyz = lines[m].split()[1:]
        # Make the x y and z coordinates into floats
        Axyz[1], Axyz[2], Axyz[3] = map(float, Axyz[1:])
        atoms.append(Axyz[0])
        xyz.append(Axyz[1:])
    # Now we create 
    xyz_list = [atoms, xyz, '', len(atoms)]
    return xyz_list

def find_potential_interactions(filtered_atoms_elec: list, filtered_atoms_nuc: list, elec_cutoff: float, nuc_cutoff: float) -> list:
    """Finds indices of potential interactions between electrophilic and nucleophilic atoms based on their values and cutoffs.
    \nparams:
    \n    filtered_atoms_elec: a list of dictionaries containing the filtered electrophilic atoms and their properties
    \n    filtered_atoms_nuc: a list of dictionaries containing the filtered nucleophilic atoms and their properties
    \n    elec_cutoff: the cutoff value for electrophilic interactions
    \n    nuc_cutoff: the cutoff value for nucleophilic interactions
    \nreturns:
    \n    interactions: a list of tuples containing the indices of potential interactions between electrophilic and nucleophilic atoms"""
    interactions = []
    for i, atom1 in enumerate(filtered_atoms_elec):
        elec_atom_index = int(atom1['Atom'].split('(')[1].split(')')[0])
        elec_atom_value = atom1['Electrophilicity']
        if i == 0:
            most_electrophilic_value = elec_atom_value
            most_nucleophilic_value = filtered_atoms_nuc[0]['Nucleophilicity']
        if elec_atom_value >= most_electrophilic_value - elec_cutoff:
            for atom2 in filtered_atoms_nuc:
                nuc_atom_index = int(atom2['Atom'].split('(')[1].split(')')[0])
                nuc_atom_value = atom2['Nucleophilicity']
                if nuc_atom_value <= most_nucleophilic_value + nuc_cutoff:
                    interactions.append((elec_atom_index, nuc_atom_index))
                else:
                    continue
        else:
            continue
    return interactions

def merge_xyz(filename: str, mol_1: list, mol_2: list) -> None:
    with open(filename,"w") as merge_file:
        merge_file.write(str(int(mol_1[-1])+int(mol_2[-1]))+'\n')
        merge_file.write('\n')
        #print(mol_1[-1])

        for index in range(int(mol_1[-1])):
            merge_file.write("{:<3}{:17.9f}{:17.9f}{:17.9f}\n".format(mol_1[0][index],float(mol_1[1][index][0]),float(mol_1[1][index][1]),float(mol_1[1][index][2])))
        for index in range(int(mol_2[-1])):
            merge_file.write("{:<3}{:17.9f}{:17.9f}{:17.9f}\n".format(mol_2[0][index],float(mol_2[1][index][0]),float(mol_2[1][index][1]),float(mol_2[1][index][2])))

def molecule_rearrangement(xyz_file1: str, xyz_file2: str, hotspot1_index: int, hotspot2_index: int, 
                           output_file: str, info_file: str, restrict_distance: float, strict_filter_level: int, minimal_distance: float) -> str:
    molecule_1=read_xyz2list(xyz_file1)
    # basic parameters of mol1 and mol2 in original position
    center1 = center_of_geometry(molecule_1)

    hot1_coordinates = hot_coordinate(molecule_1,hotspot1_index)

    # move to face to face first

    position1 = face2faceS(xyz_file1,xyz_file2,hotspot1_index,hotspot2_index,restrict_distance)

    # optimization
    molecule_2final = optimal_position(xyz_file1,position1,hotspot1_index,hotspot2_index, minimal_distance, strict_filter_level)
    molecule_2finallist = read_xyz2list(molecule_2final)

    hot_new_coordinates = hot_coordinate(molecule_2finallist,hotspot2_index)
    center_new = center_of_geometry(molecule_2finallist)
    hot1_to_hot_new = two_point_distance(hot1_coordinates,hot_new_coordinates)
    center1_to_center_new=two_point_distance(center1,center_new)

    #final quality check
    #block atoms
    find_block=point_between_hotspots(molecule_1,molecule_2finallist,hotspot1_index,hotspot2_index,strict_filter_level)
    block_number_info = f'block_number: {find_block[-2]}'
    vec_angle_info = f'vec_angle<{strict_filter_level*10}: {find_block[-1]}'
    #other atoms within 1.6 ang
    other_bonds=other_near_atoms(molecule_1,molecule_2finallist,hotspot1_index,hotspot2_index,minimal_distance)
    unused_bonds_info = f'unused_bonds: {other_bonds}'
    if hot1_to_hot_new <= center1_to_center_new:
        is_good_start = f'Good start-point!'
    else:
        is_good_start = f'Need to check!'
    # merge xyz
    merge_xyz(output_file,molecule_1,molecule_2finallist)
    # add onto the info file
    info_file += f'{output_file}\n'
    info_file += f'{block_number_info}\n'
    info_file += f'{vec_angle_info}\n'
    info_file += f'{unused_bonds_info}\n'
    info_file += f'{is_good_start}\n\n'
    return info_file

# Parse arguments
args = argument_parser()

out_file1=args.out_files[0]
out_file2=args.out_files[1]
if args.xyz_files is not None:
    xyz_file1=args.xyz_files[0]
    xyz_file2=args.xyz_files[1]
else:
    xyz_file1=out_file1.replace('.out','.xyz')
    xyz_file2=out_file2.replace('.out','.xyz')

# Set cutoff criteria
elec_cutoff = args.elec_cutoff
nucl_cutoff = args.nucl_cutoff
fukui_cutoff = args.fukui_cutoff
# output directory
out_directory = args.output_dir
isExist = os.path.exists(out_directory)
if not isExist:
   # Create a new directory
   os.makedirs(out_directory)

# Method selection
method = args.method
if method == 'Voronoi':
    method = 'Voronoï' # We need to change this for the molecule_info_fdl function, as it is called Voronoï in the AMS output file

if method == 'FMO':
    info_out1 = molecule_info_fmo(out_file1)
    info_out2 = molecule_info_fmo(out_file2)
elif method in ['Hirshfeld', 'Mulliken', 'Voronoï']:
    info_out1 = molecule_info_fdl(out_file1)
    info_out2 = molecule_info_fdl(out_file2)
else:
    raise ValueError("Invalid method. Please choose 'FMO', 'Hirshfeld', 'Mulliken', or 'Voronoï'.")

info_file_path = f'{out_directory}/{args.info_file}'

# Get rest of arguments
hotspot_distance = args.hotspot_distance
strict_filter_level = args.strict_filter_level
restrict_distance = args.restrict_distance
minimal_distance = args.minimal_distance
n_candidates = args.n_candidates
add_aditional_H = args.add_aditional_H

# Make filters and keep top candidates, check for additional H atoms if needed
descriptor = args.descriptor
if descriptor == 'fukui':
    filtered_elec1 = fukui_filter(info_out1, method=method, filter_type='elec', fukui_cutoff=fukui_cutoff)
    filtered_nuc1 = fukui_filter(info_out1, method=method, filter_type='nuc', fukui_cutoff=fukui_cutoff)
    filtered_elec2 = fukui_filter(info_out2, method=method, filter_type='elec', fukui_cutoff=fukui_cutoff)
    filtered_nuc2 = fukui_filter(info_out2, method=method, filter_type='nuc', fukui_cutoff=fukui_cutoff)
elif descriptor == 'ElNuc':
    filtered_elec1 = elec_nuc_filter(info_out1, method=method, filter_type='elec')
    filtered_nuc1 = elec_nuc_filter(info_out1, method=method, filter_type='nuc')
    filtered_elec2 = elec_nuc_filter(info_out2, method=method, filter_type='elec')
    filtered_nuc2 = elec_nuc_filter(info_out2, method=method, filter_type='nuc')
elif descriptor == 'dual':
    filtered_elec1 = elec_nuc_filter(info_out1, method=method, filter_type='elec')
    filtered_nuc1 = elec_nuc_filter(info_out1, method=method, filter_type='nuc')
    filtered_elec2 = elec_nuc_filter(info_out2, method=method, filter_type='elec')
    filtered_nuc2 = elec_nuc_filter(info_out2, method=method, filter_type='nuc')
else:
    raise ValueError("Invalid descriptor. Please choose 'fukui', 'ElNuc', or 'dual'.")

top_elec1 = keep_top_candidates(filtered_elec1, n_candidates=n_candidates, add_aditional_H=add_aditional_H)
top_nuc1 = keep_top_candidates(filtered_nuc1, n_candidates=n_candidates, add_aditional_H=add_aditional_H)
top_elec2 = keep_top_candidates(filtered_elec2, n_candidates=n_candidates, add_aditional_H=add_aditional_H)
top_nuc2 = keep_top_candidates(filtered_nuc2, n_candidates=n_candidates, add_aditional_H=add_aditional_H)

# Create settings to info file
info_file = f"This is the information file for the reactant docking between {out_file1} and {out_file2}"
info_file += "***Settings***\n\n"
info_file += f'Method : {method}\n'
info_file += f'Descriptor : {descriptor}\n'
if descriptor == 'fukui':
    info_file += f'Fukui cutoff : {fukui_cutoff}\n'
info_file += f'Hotspot distance : {hotspot_distance} angstrom\n'
info_file += f'Strict filter level : {strict_filter_level}\n'
info_file += f'Restrict distance : {restrict_distance} angstrom\n'
info_file += f'Minimal distance : {minimal_distance} angstrom\n'
if isinstance(n_candidates, float) and 0 < n_candidates < 1:
    info_file += f'Number of candidates kept : top {int(n_candidates*100)}%\n'
elif isinstance(n_candidates, int) and n_candidates > 0:
    info_file += f'Number of candidates kept : top {n_candidates}\n'
else:
    info_file += f'Number of candidates kept : all candidates\n'
if add_aditional_H:
    info_file += f'Additional H atoms added : top {add_aditional_H} H atoms in each filter\n'
else:
    info_file += f'Additional H atoms added : no additional H atoms added\n'



# Add to the info file the relevant candidates
info_file += '**Candidate Section**\n\n'
info_file += f'Candidates for {out_file1}:\n'
info_file += f'Electrophilic atoms (sorted high to low): {[atom["Atom"] for atom in top_elec1]}\n'
info_file += f'Nucleophilic atoms (sorted low to high): {[atom["Atom"] for atom in top_nuc1]}\n'
info_file += f'Candidates for {out_file2}:\n'
info_file += f'Electrophilic atoms (sorted high to low): {[atom["Atom"] for atom in top_elec2]}\n'
info_file += f'Nucleophilic atoms (sorted low to high): {[atom["Atom"] for atom in top_nuc2]}\n\n'

# Add cutoff information to the info file
info_file += f'Cutoff criteria:\n'
info_file += f'Electrophilic cutoff: {elec_cutoff}\n'
info_file += f'Nucleophilic cutoff: {nucl_cutoff}\n\n'

# Find potential interactions and produce intermediate XYZ files
potential_interactions_1_2 = find_potential_interactions(top_elec1, top_nuc2, elec_cutoff, nucl_cutoff)
potential_interactions_2_1 = find_potential_interactions(top_elec2, top_nuc1, elec_cutoff, nucl_cutoff) 
print("Potential interactions found")
# Add potential interactions to the info file
info_file += f'Potential interactions between {out_file1} and {out_file2}:\n'
info_file += f'\nFrom {out_file1} to {out_file2}:\n'
for interaction in potential_interactions_1_2:
    elec_atom = next(atom for atom in top_elec1 if int(atom['Atom'].split('(')[1].split(')')[0]) == interaction[0])
    nuc_atom = next(atom for atom in top_nuc2 if int(atom['Atom'].split('(')[1].split(')')[0]) == interaction[1])
    info_file += f'{elec_atom["Atom"]} (elec: {elec_atom["Electrophilicity"]}) - {nuc_atom["Atom"]} (nucl: {nuc_atom["Nucleophilicity"]})\n'
info_file += f'\nFrom {out_file2} to {out_file1}:\n'
for interaction in potential_interactions_2_1:
    elec_atom = next(atom for atom in top_elec2 if int(atom['Atom'].split('(')[1].split(')')[0]) == interaction[0])
    nuc_atom = next(atom for atom in top_nuc1 if int(atom['Atom'].split('(')[1].split(')')[0]) == interaction[1])
    info_file += f'{elec_atom["Atom"]} (elec: {elec_atom["Electrophilicity"]}) - {nuc_atom["Atom"]} (nucl: {nuc_atom["Nucleophilicity"]})\n'

# Now we run the molecule rearrangement for each potential interaction and save the intermediate XYZ files

info_file += f'\nIntermediate XYZ files produced for potential interactions between {out_file1} and {out_file2}:\n'

for a_index, b_index in potential_interactions_1_2:
    elec_atom = next(atom for atom in top_elec1 if int(atom['Atom'].split('(')[1].split(')')[0]) == a_index)
    nuc_atom = next(atom for atom in top_nuc2 if int(atom['Atom'].split('(')[1].split(')')[0]) == b_index)
    intermediate_xyz_name = f'{out_directory}/{out_file1.split("/")[-1].split(".")[0]}_{out_file2.split("/")[-1].split(".")[0]}_ElecNu_{a_index}{elec_atom["Atom"].split("(")[0]}_{b_index}{nuc_atom["Atom"].split("(")[0]}.xyz'
    info_file = molecule_rearrangement(xyz_file1, xyz_file2, a_index, b_index, intermediate_xyz_name, info_file, restrict_distance, strict_filter_level, minimal_distance)
print(f'Intermediate XYZ files for interactions from {out_file1} to {out_file2} have been produced and saved in {out_directory}.\n')
# Now for the other direction
for a_index, b_index in potential_interactions_2_1:
    elec_atom = next(atom for atom in top_elec2 if int(atom['Atom'].split('(')[1].split(')')[0]) == a_index)
    nuc_atom = next(atom for atom in top_nuc1 if int(atom['Atom'].split('(')[1].split(')')[0]) == b_index)
    intermediate_xyz_name = f'{out_directory}/{out_file2.split("/")[-1].split(".")[0]}_{out_file1.split("/")[-1].split(".")[0]}_ElecNu_{a_index}{elec_atom["Atom"].split("(")[0]}_{b_index}{nuc_atom["Atom"].split("(")[0]}.xyz'
    info_file = molecule_rearrangement(xyz_file2, xyz_file1, a_index, b_index, intermediate_xyz_name, info_file, restrict_distance, strict_filter_level, minimal_distance)
print(f'Intermediate XYZ files for interactions from {out_file2} to {out_file1} have been produced and saved in {out_directory}.\n')

info_file += f'End of the information file.'
# Save the info file
with open(info_file_path, "w") as f:
    f.write(info_file)
print(f'Information about candidates, cutoffs, potential interactions, and intermediate XYZ files has been saved in {info_file_path}.\n')