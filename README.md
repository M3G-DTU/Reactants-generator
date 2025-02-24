# Reactants generator
Automated reactant generator based on local reactivity of the molecule
Automated reactant generator based on local reactivity of the molecule need numpy, geatpy

usage:

put input xyzfiles and outfiles in a input folder

change settings in reactants_docking_setting.txt

then

python reactants_docking.py

you will find something in your target output folder

if you want to copy the scripts to new folder, please note following files are necessary

reactants_docking.py

**molecule_rearrangement.py ** reactants_docking_setting.txt

molearr folder

split_xyz.py is a test split script after your process

you can use it as

python split_xyz.py filename atom_number_of_first_molecule atom_number_of_second_molecule

example

python split_xyz.py filename.xyz 12 13

it will generate two files filename_split_1.xyz and filename_split_2.xyz
