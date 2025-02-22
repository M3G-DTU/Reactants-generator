# Code Documentation and Instructions for Running Locally

This script processes molecular properties from .out files to generate rearrangements based on electrophilicity and nucleophilicity. It reads settings from a configuration file, applies filters on atomic properties, and generates new molecular configurations using a Python script.

## Overview

The script does the following:

	1.	Reads molecular data from .out files and parses specific sections related to electrophilicity and nucleophilicity.
	2.	Filters and sorts atoms based on their properties.
	3.	Generates new molecular configurations using a secondary script molecule_rearrangement.py with the atoms that meet certain criteria.
	4.	Allows customization of how the hydrogen atoms are weighted during the filtering process.
	5.	Outputs scripts that create new XYZ files containing the rearranged molecules.

## Prerequisites

Before running the script, ensure the following:

	1.	Python 3.x is installed.
	2.	The necessary Python packages are installed (numpy, geatpy).
	3.	A secondary script molecule_rearrangement.py should be available in the working directory.
	4.	.out files with the required molecular data are available.
	5.	A configuration file named reactants_docking_setting.txt should be correctly formatted.

## Installation Steps

### Step 1: Install Python

Ensure you have Python 3.x installed on your machine. You can download and install it from the official website: https://www.python.org/downloads/.

### Step 2: Create and Activate a Virtual Environment (Optional but Recommended)

It’s recommended to create a Python virtual environment to isolate dependencies for this project. You can do this by running:
```python
python -m venv myenv
source myenv/bin/activate  # On Windows use myenv\Scripts\activate
```

### Step 3: Install Required Dependencies

The script uses two Python libraries: numpy and geatpy. Install them using the following commands:
```python
pip install numpy
pip install geatpy
```
Alternatively, you can add more dependencies to the dependencies list in the script if necessary.

### Step 4: Prepare the Input Files

The script requires:

	1.	../Input/*.xyz files with molecular coordinates.
    2.	../Input/*.out files with molecular data.
    3.  molearr folder
	4.	./reactants_docking_setting.txt, a setting file formatted as shown below:

```
Line 1: #xyz_file1
Line 2: <path to xyz_file1>
Line 3: #xyz_file2
Line 4: <path to xyz_file2>
Line 5: #out_file1
Line 6: <path to out_file1>
Line 7: #out_file2
Line 8: <path to out_file2>
Line 9: #set critera
Line 10: #elec_cutoff. Set >=10 will ignore this restriction.
Line 11: <electrophilicity cutoff>
Line 12: #nucl_cutoff. Set >10 will ignore this restriction.
Line 13: <nucleophilicity cutoff>
Line 14: #name for output directory (don't need to create by yourself)
Line 15: <output directory>
Line 16: #hotspots distances setting (angstorm)
Line 17: <distance>
Line 18: #keep top N (1,2,3...) or N%(please convert to float, like 0.5)? Type 0 to keep all candidates.
Line 19: <top candidate count>
Line 20: #quality check
Line 21: #strict filter level: N (int >=3), higher stricter.
Line 22: <filter level>
Line 23: #minimal distance between two other atoms (at different molecules), higher stricter. Recommandation: 1.6
Line 24: <minimal distance>
```

### Step 5: Running the Script

Run the script using Python:

```python
python reactants_docking.py
```

You will be prompted to decide whether to apply more weight to hydrogen atoms in the filtering process. Type Y for yes or N for no.

## Example Workflow

When the script runs, it will:

	1.	Parse the .out files to extract atomic properties.
	2.	Filter and sort atoms based on the electrophilicity and nucleophilicity criteria specified in the reactants_docking_setting.txt.
	3.	Generate new molecular rearrangements by invoking the molecule_rearrangement.py script and create output XYZ files in the output directory.

## Output Files

The script will generate .xyz files in the output directory specified in the settings file. Each file represents a new molecular arrangement, with naming conventions based on the rearranged molecules and atoms involved.

## Split molecules

split_xyz.py is a test split script after your process

you can use it as 

```python
python split_xyz.py filename atom_number_of_first_molecule atom_number_of_second_molecule 
```

Example:

```python
python split_xyz.py filename.xyz 12 13
```
It will generate two files filename_split_1.xyz and filename_split_2.xyz


## Error Handling

	•	If a required package is not installed, a DistributionNotFound or VersionConflict will be thrown. Ensure the dependencies listed at the start are installed correctly.
	•	If the output directory does not exist, the script will automatically create it.

By following these steps and ensuring the correct files are in place, the script should execute smoothly on your local machine.
