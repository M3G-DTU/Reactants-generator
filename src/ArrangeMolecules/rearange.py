# Class which can apply different operations to move a molecule around another molecule to generate different conformations of the reactant complex. 
# This is used to generate different conformations of the reactant complex for the reaction prediction model.

import numpy as np
from scipy import optimize
# Import differential evolution
from scipy.optimize import differential_evolution
from ase import Atoms
from ase.io import read, write
from xtb.ase.calculator import XTB
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw, rdDetermineBonds
from ase.data import covalent_radii, atomic_numbers

class OrientMoleculeSphere:
    def __init__(self, molecule1: Atoms, molecule2: Atoms, atom_1: int, atom_2: int, seed: int = None, sphere_radius: float = 5.0):
        self.molecule1_init = molecule1  # The molecule around which the other molecule will be translated (the "anchor" molecule)
        self.molecule2_init = molecule2  # The molecule which will be translated around the anchor molecule
        self.atom_1 = atom_1        # Atom index of molecule 1 
        self.atom_2 = atom_2        # Atom index of molecule 2
        if seed is not None:
            np.random.seed(seed)
        
        #self.sphere_center = self.molecule1_init.positions.mean(axis=0) # Center of the sphere is the mean positions of the first molecule
        #self.sphere_radius = np.max(np.linalg.norm(self.molecule1_init.positions - self.sphere_center, axis=1)) - 1# + sphere_radius_buffer # Adding a buffer.
        
        self.sphere_center = self.molecule1_init[self.atom_1].position # Center of the sphere is the position of the atom of interest in molecule 1
        self.sphere_radius = sphere_radius

        # Place molecule 2 at a random position on the surface of the sphere around molecule 1
        random_theta = np.random.uniform(0, 2 * np.pi)
        random_phi = np.random.uniform(0, np.pi)
        random_x = self.sphere_radius * np.sin(random_phi) * np.cos(random_theta)
        random_y = self.sphere_radius * np.sin(random_phi) * np.sin(random_theta)
        random_z = self.sphere_radius * np.cos(random_phi)
        random_position = np.array([random_x, random_y, random_z]) + self.sphere_center
        translation_vector = random_position - self.molecule2_init[self.atom_2].position
        for atom in self.molecule2_init:
            atom.position += translation_vector
        # Now for the movable molecules
        self.molecule1 = self.molecule1_init.copy()
        self.molecule2 = self.molecule2_init.copy()
        # Add calculator to get charges for electrostatic repulsion
        #!TODO add other options. Not a lot of calculators support charge calculation, might want to find other objective functions.
        self.molecule1.calc = XTB(method='GFN2-xTB')
        self.molecule2.calc = XTB(method='GFN2-xTB')

    # Utility functions
    def create_rdkit_molecule(self, ase_molecule: Atoms):
        """Convert an ASE Atoms object to an RDKit molecule."""
        rdkit_molecule = Chem.RWMol()
        for atom in ase_molecule:
            rdkit_molecule.AddAtom(Chem.Atom(atom.number))
        rdDetermineBonds.DetermineConnectivity(rdkit_molecule)
        return rdkit_molecule
    
    def normalize_vector(self, vector):
        return vector / np.linalg.norm(vector)
    
    def angle_between_vectors(self, v1, v2):
        """Return the angle in radians between vectors 'v1' and 'v2'"""
        v1_u = self.normalize_vector(v1)
        v2_u = self.normalize_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    
    # Translation and rotation on sphere
    def step_on_sphere(self, point_on_sphere: list, center: list, step_theta: float, step_phi: float):
        """Translate a point on the surface of a sphere by given angles theta and phi. Units of angles is radians."""
        
        x, y, z = point_on_sphere
        a, b, c = center

        r = self.sphere_radius # Radius of the sphere
        theta = np.arccos((z - c) / r) - step_theta
        phi = np.arctan((y - b) / (x - a)) - step_phi

        x_new = a + r * np.sin(theta) * np.cos(phi)
        y_new = b + r * np.sin(theta) * np.sin(phi)
        z_new = c + r * np.cos(theta)

        return [x_new, y_new, z_new]
    
    def translate_on_sphere(self, molecule: Atoms, initial_point: list, center: list, step_theta: float, step_phi: float):
        """Translate the molecule by moving it along the surface of a sphere."""
        new_position = self.step_on_sphere(initial_point, center, step_theta, step_phi)
        translation_vector = np.array(new_position) - np.array(initial_point)
        
        # Translate the molecule by the calculated translation vector
        for atom in molecule:
            atom.position += translation_vector
        return molecule.positions

    # Rotation matrix
    def Rx(self, angle):
        """Rotation matrix for rotation around x-axis. Units of angles is radians."""
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[1, 0, 0],
                         [0, c, -s],
                         [0, s, c]])
    
    def Ry(self, angle):
        """Rotation matrix for rotation around y-axis. Units of angles is radians."""
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])

    def Rz(self, angle):
        """Rotation matrix for rotation around z-axis. Units of angles is radians."""
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])

    def rotation_matrix(self, alpha, beta, gamma):
        """Combined rotation matrix for rotations around x, y, and z axes. Units of angles is radians."""
        return self.Rz(gamma) @ self.Ry(beta) @ self.Rx(alpha)
    
    # Rotate molecule
    def rotate_molecule(self, molecule: Atoms, center: list, alpha: float, beta: float, gamma: float):
        R = self.rotation_matrix(alpha, beta, gamma)
        return (molecule - center) @ R.T + center

    # Utility functions for objective function
    def lennard_jones_repulsion(self, R1, R2, sigma=1):
        diff = R1[:,None,:] - R2[None,:,:]
        r = np.linalg.norm(diff, axis=-1)
        return np.sum((sigma / r)**12)

    def electrostatic_repulsion(self, R1, R2, q1, q2):
        diff = R1[:, None, :] - R2[None, :, :]
        r = np.linalg.norm(diff, axis=-1)
        qprod = np.abs(q1[:, None] * q2[None, :])  # force repulsion
        return np.sum(qprod / r)

    # Objective function to minimize
    def objective_combined(self, Vars):
        step_theta, step_phi, alpha, beta, gamma = Vars
        
        # --- Translation ---
        new_center = self.step_on_sphere(self.molecule2[self.atom_2].position, self.sphere_center, step_theta, step_phi)

        # Shift molecule 2 to the new center
        translation_vector = new_center - self.molecule2[self.atom_2].position
        mol2_position = self.molecule2.positions + translation_vector

        # --- Rotation ---
        mol2_position = self.rotate_molecule(mol2_position, new_center, alpha, beta, gamma)
        mol1_position = self.molecule1.positions

        # --- Repulsion LJ ---
        rep = self.lennard_jones_repulsion(mol1_position, mol2_position)

        # --- Electrostatic repulsion ---
        # Get charges
        #q1 = self.molecule1.get_charges()
        #q2 = self.molecule2.get_charges()
#
        #rep += self.electrostatic_repulsion(mol1_position, mol2_position, q1, q2)
        return rep

    # Optimization function
    def optimize(self):

        # Bounds for the decision variables
        upper_bound = [ np.pi/2,  np.pi,  np.pi,  np.pi/2,  np.pi]  
        lower_bound = [-np.pi/2, -np.pi, -np.pi, -np.pi/2, -np.pi]
        bounds = optimize.Bounds(lower_bound, upper_bound)

        # Initial guess
        initial_guess = [0.0, 0.0, 0.0, 0.0, 0.0]
        # if len(self.molecule2) > 6 and len(self.molecule1) > 6:
        # We want more initial diversity to avoid local minima
        if len(self.molecule2.positions) > 6 and len(self.molecule1.positions) > 6:
            result = optimize.differential_evolution(self.objective_combined, bounds=bounds, maxiter=100, tol=1e-6)
        else:
        # Optimize translation first
            result = optimize.minimize(self.objective_combined, 
                                             x0=initial_guess, bounds=bounds, 
                                             method='L-BFGS-B', 
                                             options={'maxiter': 100, 'ftol': 1e-6})
        
        # Apply the best translation to the molecule before optimizing rotation
        best_step_theta, best_step_phi, best_alpha, best_beta, best_gamma = result.x
        
        # Find the vector to apply the translation to molecule 2
        self.molecule2.positions = self.translate_on_sphere(self.molecule2, self.molecule2[self.atom_2].position, self.sphere_center, best_step_theta, best_step_phi)
        self.molecule2.positions = self.rotate_molecule(self.molecule2.positions, self.molecule2[self.atom_2].position, best_alpha, best_beta, best_gamma)
        
        return result

if __name__ == "__main__":
    # Example usage
    molecule1 = read('molecule1.xyz')  # Replace with actual file path
    molecule2 = read('molecule2.xyz')  # Replace with actual file path

    pair_of_atoms = [(8, 1), (23, 1), (1, 1), (13, 1), (4, 1), (9, 1), (24, 1), (16, 1), (20, 1), (26, 1), (18, 1), (22, 1)]
    traj_file = 'optimization_trajectory.xyz'
    # Remove existing trajectory file if it exists
    import os
    if os.path.exists(traj_file):
        os.remove(traj_file)
    for atom1_index, atom2_index in pair_of_atoms:
        atom1_index -= 1 # Convert to 0-based index
        atom2_index -= 1 # Convert to 0-based index
        orientor = OrientMoleculeSphere(molecule1, molecule2, atom1_index, atom2_index)
        result = orientor.optimize()
        best_step_theta, best_step_phi, best_alpha, best_beta, best_gamma = result.x
        merged_molecule = orientor.molecule1 + orientor.molecule2
        #write(f'optimized_merged_molecule_{atom1_index}_{atom2_index}.xyz', merged_molecule,)
        # Append to trajectory file
        with open(traj_file, 'a') as traj:
            traj.write(f"{len(merged_molecule)}\n")
            traj.write(f"Optimized merged molecule for atom pair ({atom1_index}, {atom2_index})\n")
            for atom in merged_molecule:
                traj.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")

    traj_file2 = 'conf.xyz'
    # Remove existing trajectory file if it exists
    if os.path.exists(traj_file2):
        os.remove(traj_file2)

    atom_pairs = [(3, 2), (3, 11), (3, 5), (3, 15), (3, 6), (3, 14), (3, 3), (3, 10), (3, 25), (3, 7), (3, 17), (3, 12), (3, 19)]
    molecule1 = read('molecule2.xyz')
    molecule2 = read('molecule1.xyz')
    for atom1_index, atom2_index in atom_pairs:
        # Convert to 0-based index
        atom1_index -= 1
        atom2_index -= 1

        orientor = OrientMoleculeSphere(molecule1, molecule2, atom1_index, atom2_index)

        result = orientor.optimize()
        best_step_theta, best_step_phi, best_alpha, best_beta, best_gamma = result.x

        merged_molecule = orientor.molecule1 + orientor.molecule2
        #write(f'optimized_merged_molecule_{atom1_index}_{atom2_index}.xyz', merged_molecule,)
        with open(traj_file2, 'a') as traj:
            traj.write(f"{len(merged_molecule)}\n")
            traj.write(f"Optimized merged molecule for atom pair ({atom1_index}, {atom2_index})\n")
            for atom in merged_molecule:
                traj.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")