# Reactants Generator
Automated reactant generator for the reactivity analysis of organic flow battery electrolytes

## Prerequisites and Installation

As geatpy requires Python version 3.5, 3.6, 3.7, 3.8, 3.9, or 3.10, make sure you have one of these versions installed on your machine. You can download and install Python from the official website: https://www.python.org/downloads/.

```bash
python=3.10
geatpy==2.7.0
numpy==1.26.4
```

For more about installation see [Usage and installation](Documentation/Usage%20and%20installation.md)

## Usage


Aquire .out files for the reactants. It should be a Conceptual DFT calculation calculation with output files from the Amsterdam Modeling Suite (AMS) program. [https://www.scm.com/amsterdam-modeling-suite/]

## Available options
To see available options use:

```bash
python reactants_docking.py --help
```

The only required input is the .out files for the reactants using the '-o' or '--out_files' options. So the most basic command to run the script is:

```bash
python reactants_docking.py -o path/to/reactant1.out path/to/reactant2.out
```

The default output folder is set to "Output". If you want to change it, you can specify the output directory using the `--output_dir` option when running the script. For example:

```bash
python reactants_docking.py -o path/to/reactant1.out path/to/reactant2.out --output_dir path/to/my/outputfolder
```