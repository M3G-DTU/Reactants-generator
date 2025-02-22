# Documentation for Molecular Optimization and Docking Script Using Genetic Algorithm

This Python script performs the optimization of molecular docking using a genetic algorithm implemented through the geatpy library. The goal of the script is to align two molecules at specific “hotspots” by optimizing the orientation and ensuring no significant overlap or undesired interactions between atoms. The script reads molecular coordinates from XYZ files, computes angles for molecular rotation, and runs a genetic algorithm to find the optimal docking configuration.

## Overview

The key components of the script are:

	1.	Molecular Input: Reading molecular coordinates from XYZ files.
	2.	Hotspot Alignment: Positioning molecule 2 relative to molecule 1 based on hotspot atoms.
	3.	Genetic Algorithm Optimization: Using a genetic algorithm to find optimal rotation angles that minimize atomic distance and avoid blocking.
	4.	Quality Checks: Ensuring no significant overlap or interactions between atoms.
	5.	Final Docking Configuration: The script outputs the optimal orientation of molecule 2 relative to molecule 1.

## Prepare Input Files

	•	Molecular XYZ files: You need two XYZ files for the molecules you want to dock, e.g., benzo.xyz and new1.xyz.
	•	Settings file: The reactants_docking_setting.txt file defines various parameters, including filter levels and distance settings. An example format is as follows:
```
Line 1-20: <parameters>
Line 21: <strict filter level>
Line 23: <maximum allowed distance>
```

## Script Functions

1. angle2dist(mol1, mol2, hot1, hot2, anglelist)

	•	Purpose: Moves mol2 based on the specified angles and computes the minimum distance between the atoms of mol2 and the hotspot of mol1.
	•	Parameters:
	•	mol1: Filename of the first molecule (XYZ format).
	•	mol2: Filename of the second molecule (XYZ format).
	•	hot1: Hotspot atom index for molecule 1.
	•	hot2: Hotspot atom index for molecule 2.
	•	anglelist: List of rotation angles to apply to molecule 2.
	•	Returns: The minimum distance between any atom in mol2 and the hotspot in mol1.

2. block_due2angle(mol1, mol2, hot1, hot2, anglelist, factor)

	•	Purpose: Moves mol2 based on the specified angles and checks if any atoms block the hotspot-to-hotspot interaction.
	•	Parameters:
	•	mol1, mol2: Filenames of the two molecules (XYZ format).
	•	hot1, hot2: Hotspot atom indices for the two molecules.
	•	anglelist: List of rotation angles to apply to molecule 2.
	•	factor: Strictness level for the blockage check.
	•	Returns: The number of blocking atoms between the two hotspots.

3. geat_optimization(mol1, mol2, hot1, hot2)

	•	Purpose: Defines the optimization problem using a genetic algorithm and finds the best set of angles to align mol2 to mol1 based on minimizing distances and avoiding blockages.
	•	Parameters:
	•	mol1, mol2: Filenames of the two molecules (XYZ format).
	•	hot1, hot2: Hotspot atom indices for the two molecules.
	•	Returns: The optimal set of rotation angles to apply to mol2.

4. optimal_position(mol1, mol2, hot1, hot2)

	•	Purpose: Uses the optimized angles from geat_optimization() to compute the final optimal position of mol2 relative to mol1.
	•	Parameters:
	•	mol1, mol2: Filenames of the two molecules (XYZ format).
	•	hot1, hot2: Hotspot atom indices for the two molecules.
	•	Returns: The filename of the final XYZ file with the optimized molecular positions.

## Genetic Algorithm Details

The genetic algorithm (geat_optimization()) optimizes the following:

	•	Objective: Minimize the distance between atoms of mol2 and the hotspot of mol1.
	•	Constraints: Ensure no atoms block the interaction between the hotspots of the two molecules.

The decision variables are the five rotation angles applied to mol2, which are bounded by [-60°, 60°] for each angle. The optimization process evolves over 50 generations, with a population size of 20.


## Output

	•	Optimal Angles: The best set of rotation angles for aligning mol2 to mol1.
	•	Optimized XYZ File: The final configuration of the docked molecules, saved as a new XYZ file.

## Conclusion

This script is a powerful tool for molecular docking and optimization using a genetic algorithm. It can efficiently compute the best orientation of two molecules based on their geometric and chemical properties while ensuring no significant interactions or overlaps occur. By following the above steps, you should be able to run this script successfully on your molecular systems.