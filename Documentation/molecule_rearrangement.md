# Documentation for the Molecular Docking Script

This Python script facilitates the docking of two molecules by aligning them based on their “hotspots,” performing an optimization, and then merging the molecular coordinates into a new XYZ file. The script reads molecular data from XYZ files, computes geometric properties, and ensures quality checks on the final molecular configuration. Below are the detailed steps and instructions on how to set up and run this code.

## Overview

The script follows these main steps:

	1.	Reads molecular data from two XYZ files.
	2.	Identifies hotspots on both molecules and computes relevant geometric parameters.
	3.	Positions molecule 2 relative to molecule 1, so they are aligned “face-to-face” at their hotspots.
	4.	Optimizes the position of molecule 2.
	5.	Performs quality checks on the final configuration to ensure no significant overlaps or undesirable bonds.
	6.	Merges the two molecules into a single XYZ file.

## Script Workflow

1.	Reading Molecules:

    The script uses read_xyz2list() to load the molecular coordinates from the given XYZ files.

2.	Geometric Computations:
The following geometric operations are performed:

	•	Center of Geometry: center_of_geometry() computes the geometric center of each molecule.

	•	Hotspot Coordinates: hot_coordinate() extracts the coordinates of specific atoms (hotspots) from each molecule.

	•	Distance and Vector Calculations: The script uses two_point_distance() to compute distances between atoms or centers, and np.subtract() to compute vectors.
3.	Face-to-Face Positioning:

    The script moves the second molecule so that its hotspot is positioned face-to-face with the hotspot of the first molecule using face2faceS(). The distance between hotspots is constrained by a user-defined restriction.
4.	Optimization:

    After positioning, the script refines the alignment using optimal_position().
5.	Quality Check:

    The final position of the molecules is checked for quality using:

	•	point_between_hotspots(): Ensures no atoms obstruct the bonding between hotspots.

	•	other_near_atoms(): Checks for unwanted short bonds between atoms.

6.	Merging Molecules:

    Once the molecules are aligned and optimized, they are merged into a single XYZ file using the merge_xyz() function.

## Output

The script will generate the following output:

	•	Merged XYZ File: A new file containing both molecules in their final docked positions.
	•	Quality Check Information: Printed information on blocking atoms, bond lengths, and other checks related to the docking quality.

## Error Handling

	•	Ensure all necessary modules are available in the molrearr package.
	•	Verify the XYZ files and settings file are correctly formatted.
	•	Check for proper hotspot indices (should be integers and correspond to atom numbers in the XYZ files).

