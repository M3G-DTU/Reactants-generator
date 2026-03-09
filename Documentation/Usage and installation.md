# Code Documentation and Instructions for Running Locally

This script processes molecular properties from .out files to generate rearrangements based on electrophilicity and nucleophilicity. It reads settings from a configuration file, applies filters on atomic properties, and generates new molecular configurations using a Python script.

## Overview

The script does the following:

	1.	Reads molecular data from .out files and parses specific sections related to electrophilicity and nucleophilicity.
	2.	Filters and sorts atoms based on their properties.
	3.	Generates new molecular configurations with the atoms that meet certain criteria.
	4.	Allows customization of how the hydrogen atoms are weighted during the filtering process.
	5.	Outputs scripts that create new XYZ files containing the rearranged molecules.
	6.  To see available options use "python reactants_docking.py --help"
 
## Prerequisites

Before running the script, ensure the following:

	1.	Python 3.5, 3.6, 3.7, 3.8, 3.9, or 3.10 is installed.
	2.	The necessary Python packages are installed (numpy, geatpy).
	3.	.out files with the required molecular data are available.

## Installation Steps

### Step 1: Install using Python

Ensure you have Python 3.5, 3.6, 3.7, 3.8, 3.9, or 3.10 installed on your machine. You can download and install it from the official website: https://www.python.org/downloads/.

#### Step 1.1: Create and Activate a Virtual Environment (Optional but Recommended)
It’s recommended to create a Python virtual environment to isolate dependencies for this project. You can do this by running:
```python
python -m venv myenv
source myenv/bin/activate  # On Windows use myenv\Scripts\activate
```
The script uses two Python libraries: numpy and geatpy. Install them using the following commands:
```python
pip install numpy
pip install geatpy
```
Alternatively, you can add more dependencies to the dependencies list in the script if necessary.

### Step 2: Install using Conda and Pip
If you prefer using Conda, you can create and activate a Conda environment with the following
commands:
```python
conda create --name geat_env python=3.10 pip
conda activate geat_env
pip install --upgrade pip
pip install -r requirements.txt
python -c "import numpy, geatpy; print(numpy.__version__, geatpy.__version__)"
```
This will create a Conda environment named `geat_env`, install the required packages, and verify the installation by printing the versions of numpy and geatpy.

### Step 3: Before running the script

	Make sure you have the inital .out files. 



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

The script will generate .xyz files in the output directory which can be specified using the '--output_dir' option. Each file represents a new molecular arrangement, with naming conventions based on the rearranged molecules and atoms involved.


## Error Handling

	•	If a required package is not installed, a DistributionNotFound or VersionConflict will be thrown. Ensure the dependencies listed at the start are installed correctly.
	•	If the output directory does not exist, the script will automatically create it.

By following these steps and ensuring the correct files are in place, the script should execute smoothly on your local machine.
