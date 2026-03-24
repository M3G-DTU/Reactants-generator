# Class which can apply different operations to move a molecule around another molecule to generate different conformations of the reactant complex. 
# This is used to generate different conformations of the reactant complex for the reaction prediction model.

import numpy as np
import geatpy as ea
from scipy import optimize
from ase import Atoms
from ase.io import read, write

class TranslateMoleculeSphere:
    def __init__(self, molecule1: Atoms, molecule2: Atoms, atom_1: int, atom_2: int):
        self.molecule1_init = molecule1  # The molecule around which the other molecule will be translated (the "anchor" molecule)
        self.molecule2_init = molecule2  # The molecule which will be translated around the anchor molecule
        self.atom_1 = atom_1        # Atom index of molecule 1 
        self.atom_2 = atom_2        # Atom index of molecule 2
        self.count = 0
        
        #self.sphere_center = self.molecule1_init.positions.mean(axis=0) # Center of the sphere is the mean positions of the first molecule
        #self.sphere_radius = np.max(np.linalg.norm(self.molecule1_init.positions - self.sphere_center, axis=1)) - 1# + sphere_radius_buffer # Adding a buffer.
        
        self.sphere_center = self.molecule1_init[self.atom_1].position # Center of the sphere is the position of the atom of interest in molecule 1
        self.sphere_radius = 3.0

        # Place molecule 2 at a random position on the surface of the sphere around molecule 1
        # set seed
        np.random.seed(42)
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
    # Utility functions
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

    def center_rotate(self, molecule: Atoms, center: list, alpha: float, beta: float, gamma: float):
        """Rotate the molecule around a given center using Euler angles for x, y, z axes. Unit of angles is radians."""
        origin = -1 * np.array(center) # Translate to origin
        molecule.positions += origin
        rotation_matrix = self.Rz(gamma) @ self.Ry(beta) @ self.Rx(alpha) # Combined rotation matrix
        molecule.positions = molecule.positions @ rotation_matrix.T # Apply rotation
        molecule.positions -= origin # Translate back to original center
        

    def translate_on_sphere(self, molecule: Atoms, initial_point: list, center: list, step_theta: float, step_phi: float):
        """Translate the molecule by moving it along the surface of a sphere."""
        new_position = self.step_on_sphere(initial_point, center, step_theta, step_phi)
        translation_vector = np.array(new_position) - np.array(initial_point)
        
        # Translate the molecule by the calculated translation vector
        for atom in molecule:
            atom.position += translation_vector

    def drone_translate(self, step_theta: float, step_phi: float, alpha: float, beta: float, gamma: float):
        """Translate molecule2 around molecule1 by moving it along the surface of a sphere and applying rotations."""
        self.translate_on_sphere(self.molecule2, self.molecule2[self.atom_2].position, self.sphere_center, step_theta, step_phi)
        self.center_rotate(self.molecule2, self.molecule2[self.atom_2].position, alpha, beta, gamma)

    def locate_hindrance_atoms(self):
        distances = []
        atom1 = self.molecule1[self.atom_1]
        atom2 = self.molecule2[self.atom_2]
        for index, atom in enumerate(self.molecule1):
            if index == self.atom_1:  # Skip the atom of interest in molecule 1
                continue
            distance = np.linalg.norm(atom.position - atom2.position)
            distances.append(distance)
        for index, atom in enumerate(self.molecule2):
            if index == self.atom_2:  # Skip the atom of interest in molecule 2
                continue
            distance = np.linalg.norm(atom.position - atom1.position)
            distances.append(distance)
        return distances # List of distances for all atoms in both molecules

    def distances_to_atom_of_interest(self):
        distances = []
        atom2 = self.molecule2[self.atom_2]
        for index, atom in enumerate(self.molecule1):
            if index == self.atom_1:  # Skip the atom of interest in molecule 1
                continue
            distance = np.linalg.norm(atom.position - atom2.position)
            distances.append(distance)
        return distances # List of distances for all atoms in molecule 1 to the atom of interest in molecule 2

    def locate_hindrance_molecule2(self):
        distances = []
        atom1 = self.molecule1[self.atom_1]
        for index, atom in enumerate(self.molecule2):
            if index == self.atom_2:  # Skip the atom of interest in molecule 2
                continue
            distance = np.linalg.norm(atom.position - atom1.position)
            distances.append(distance)
        return distances # List of distances for all atoms in molecule 2 to the atom of interest in molecule 1

    def check_hindrance(self, filter_size: float = 1.0):
        """Check number of atom is between the two atoms of interest"""
        n, m = 0, 0

        vector_m1_to_m2 = self.molecule2[self.atom_2].position - self.molecule1[self.atom_1].position
        distance = np.linalg.norm(vector_m1_to_m2)
        for index, atom in enumerate(self.molecule1):
            # Check if distance between atom and molecule 2 is less than the distance between the two atoms of interest
            if index == self.atom_1:  # Skip the atom of interest in molecule 1
                continue
            
            if np.linalg.norm(atom.position - self.molecule2[self.atom_2].position) < distance:
                n += 1 
                # Angle check
                angle_of_interest = self.angle_between_vectors(atom.position - self.molecule2[self.atom_2].position, vector_m1_to_m2)
                if angle_of_interest < 10*filter_size or angle_of_interest > 180-10*filter_size:
                    m += 1
        return n, m # n is the number of atoms in molecule 1 that are between the two atoms of interest, m is the number of those atoms that are also within a certain angle threshold (i.e. more likely to be hindrance)
    
    def locate_closest_atoms(self, threshold: float = 1.6):
        """Locate number of atoms in both molecules that are within a certain distance threshold of each other. Exclude the atoms of interest."""
        count_close = 0

        for index1, atom1 in enumerate(self.molecule1):
            if index1 == self.atom_1:  # Skip the atom of interest in molecule 1
                continue
            for index2, atom2 in enumerate(self.molecule2):
                if index2 == self.atom_2:  # Skip the atom of interest in molecule 2
                    continue
                distance = np.linalg.norm(atom1.position - atom2.position)
                if distance < threshold:
                    count_close += 1

        return count_close

    def optimization_geat(self):
        @ea.Problem.single
        def objective_function(Vars, minimal_allowed_distance: float = 1.6, filter_size: float = 1.0):
            """Objective function to minimize the number of atoms between the two atoms of interest and maximize the number of close contacts."""
            # Parse variables
            step_theta, step_phi, alpha, beta, gamma = Vars
            # Reset molecule positions to initial state before applying transformations
            self.molecule1 = self.molecule1_init.copy() # Reset molecule 1 to initial position
            self.molecule2 = self.molecule2_init.copy() # Reset molecule 2 to initial position

            # Apply transformations based on the decision variables
            self.drone_translate(step_theta, step_phi, alpha, beta, gamma)

            # Calculate
            distances = self.locate_hindrance_atoms()
            weak_hindrance_count, strong_hindrance_count = self.check_hindrance(filter_size)
            # We want to check the squared distance to penalize closer contacts more heavily
            # We want to have the fewest number of atoms between the two atoms of interest, to create a cleaner reactant complex.
            min_distance = min(distances) if distances else float('inf')

            # Objective
            score = 0

            if distances:
                min_distance = min(distances)
                # We want to maximize the minimum distance between any atom in molecule 2 and the atom of interest in molecule 1, to reduce hindrance.
                score += min_distance**3  # Square the distance to penalize closer contacts more heavily
            else:
                min_distance = float('inf')
                score += 1000  # If there are no other atoms, give a high score
            
            # We want to penelize having the distance between the two atoms compared to distances of other atoms,
            # to encourage the two atoms of interest to be closer together than any other atoms, which would indicate a cleaner reactant complex.
            #score -= weak_hindrance_count  # Penalize having more atoms between
            #print(len(distances), min_distance, minimal_allowed_distance)
            #if min_distance < minimal_allowed_distance:
            #    score -= 1000  # Penalize heavily if the minimum distance is below the threshold, to avoid steric clashes.
            #else:
            #    self.count += 1
            # Constraints
            CV = np.array([
                min_distance - minimal_allowed_distance, # Ensure that the minimum distance between any atom in molecule 2 and the atom of interest in molecule 1 is above a certain threshold to avoid steric clashes
                1 - weak_hindrance_count, # Dont allow more than 1 atom in molecule 1 to be between the two atoms of interest
            #     (-1) * strong_hindrance_count, # Dont allow any atoms in molecule 2 to be within the angle threshold (i.e. strong hindrance)
                ]) 
            print(score)
            return score, CV
        
        # Define the optimization problem
        # Bounds for the decision variables: step_theta, step_phi, alpha, beta, gamma
        upper_bound = [ np.pi,  np.pi/2,  np.pi,  np.pi/2,  np.pi]  # Max rotation of 360 degrees and max translation of 5 units
        lower_bound = [-np.pi, -np.pi/2, -np.pi, -np.pi/2, -np.pi]  # Min rotation of 0 degrees and min translation of -5 units
        problem = ea.Problem(name='Molecule Arrangement Optimization',
                             M=1,  # Number of objectives
                             maxormins=[-1],  # Maximize the objective function
                             #maxormins=[1],  # Minimize the objective function
                             Dim=5,  # Number of decision variables (step_theta, step_phi, alpha, beta, gamma)
                             varTypes=[0, 0, 0, 0, 0],  # All variables are continuous
                             lb=lower_bound,  # Lower bounds for decision variables
                             ub=upper_bound,  # Upper bounds for decision variables
                             evalVars=objective_function)  # Objective function to evaluate

        # Run the optimization algorithm
        algorithm = ea.soea_SEGA_templet(problem,
                                                  ea.Population(Encoding='RI', NIND=50),  # Population size
                                                  MAXGEN=100,  # Maximum number of generations
                                                  logTras=1,    # Log every generation
                                                  trappedValue=1e-6,    # Threshold for convergence
                                                  maxTrappedCount=10)   # Maximum number of generations to wait for convergence
        res = ea.optimize(algorithm, seed=1, verbose=False, drawing=0, outputMsg=False, drawLog=False, saveFlag=False, dirName='result')
        return res


    def optimization_rotate(self):
        # Take the best solution from the geat optimization and rotate the molecule accordingly, to ensure an optimal arrangement of the reactant complex.
        @ea.Problem.single
        def objective_function(Vars):
            alpha, beta, gamma = Vars
            self.center_rotate(self.molecule2, self.molecule2[self.atom_2].position, alpha, beta, gamma)
            # Get the minimum distance between atom of interest in molecule 1 and any atom in molecule 2
            distances = self.locate_hindrance_molecule2()
            
            score = min(distances)**2 # We want to maximize the minimum distance between any atom in molecule 2 and the atom of interest in molecule 1, to reduce hindrance.
            return score, np.array([])  # We want to maximize the minimum distance between
        
        # Define the optimization problem
        upper_bound = [ np.pi,  np.pi/2,  np.pi]  # Max rotation of 360 degrees
        lower_bound = [-np.pi, -np.pi/2, -np.pi]  # Max rotation of 360 degrees
        problem = ea.Problem(name='Molecule Rotation Optimization',
                             M=1,  # Number of objectives
                             maxormins=[-1],  # Maximize the objective function
                             Dim=3,  # Number of decision variables (alpha, beta, gamma)
                             varTypes=[0, 0, 0],  # All variables are continuous
                             lb=lower_bound,  # Lower bounds for decision variables
                             ub=upper_bound,  # Upper bounds for decision variables
                             evalVars=objective_function)  # Objective function to evaluate
        # Run the optimization algorithm
        algorithm = ea.soea_SEGA_templet(problem,
                                                  ea.Population(Encoding='RI', NIND=50),  # Population size
                                                  MAXGEN=100,  # Maximum number of generations
                                                  logTras=1,    # Log every generation
                                                  trappedValue=1e-6,    # Threshold for convergence
                                                  maxTrappedCount=10)   # Maximum number of generations to wait for convergence
        res = ea.optimize(algorithm, seed=1, verbose=False, drawing=0, outputMsg=False, drawLog=False, saveFlag=False, dirName='result')
#        print(self.count)
        return res

    def optimization_scipy(self):
        # Use scipy optimization to optimize the arrangement of the reactant complex, as an alternative to geatpy.
        def objective_function(Vars):
            step_theta, step_phi, alpha, beta, gamma = Vars
            self.drone_translate(step_theta, step_phi, alpha, beta, gamma)
            score = 0
            # reduce steric hindrance by not allowing atoms too close to each other
            distances = self.locate_hindrance_atoms()
            sum_squared_distances_inv = sum(1/d**2 for d in distances if d > 0)  # Inverse of squared distances to penalize closer contacts more heavily
            score += sum_squared_distances_inv
            #if min(distances) < 1.6:
            #    score += 100
            return score
            

        # Bounds for the decision variables: step_theta, step_phi, alpha, beta, gamma
        upper_bound = [ np.pi,  np.pi/2,  np.pi,  np.pi/2,  np.pi]  
        lower_bound = [-np.pi, -np.pi/2, -np.pi, -np.pi/2, -np.pi] 
        bounds = optimize.Bounds(lower_bound, upper_bound)

        # Initial guess (can be random or based on some heuristic)
        initial_guess = [0.0, 0.0, 0.0, 0.0, 0.0]

        result = optimize.differential_evolution(objective_function, bounds=bounds, strategy='best1bin', maxiter=100, popsize=15, tol=1e-6, seed=42)
        return result

    def angles_to_unit_vector(self, theta: float, phi: float):
        """Convert spherical angles to a unit vector."""
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return np.array([x, y, z])
    
    def normalize_vector(self, vector):
        """Normalize a vector to have unit length."""
        norm = np.linalg.norm(vector)
        if norm == 0:
            return vector
        return vector / norm

    def angles_to_unit(self, theta, phi):
        return np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
            ])

    def softmin(self, x, alpha=20):
        return -np.log(np.sum(np.exp(-alpha * x))) / alpha

    def objective_direction(self, angles):
        theta, phi = angles
        u = self.angles_to_unit(theta, phi)

        projections = self.molecule2.positions @ u
        return -self.softmin(projections)  # maximize minimum projection

    def rotate_to_optimal_direction(self, angles):
        # rotation to molecule 2 around atom of interest in molecule 2
        theta, phi = angles
        u = self.angles_to_unit(theta, phi)
        projections = self.molecule2.positions @ u
        return projections
    
    def opt_scipy_2(self):
        # First we only do translation on spherere to find a good position, then we do a second optimization step to find the best rotation from that position.
        #def objective_function_translate(Vars):
        #    step_theta, step_phi = Vars
        #    self.translate_on_sphere(self.molecule2, self.molecule2[self.atom_2].position, self.sphere_center, step_theta, step_phi)
        #    distances = self.distances_to_atom_of_interest()
        #    score = min(distances)
        #    print((-1) * score)
        #    return (-1) * score
        # f(x) = Lennard-Jones Repulsion
        def objective_function_translate(Vars, p: float = 6.0):
            step_theta, step_phi = Vars
            step_theta = np.clip(step_theta, 1e-8, np.pi - 1e-8)  # Avoid exactly 0 or pi to prevent singularities

            # Define position of molecule 2 based on the current step on the sphere
            x = self.molecule2[self.atom_2].position
            self.translate_on_sphere(self.molecule2, self.molecule2[self.atom_2].position, self.sphere_center, step_theta, step_phi)
            distances = np.linalg.norm(self.molecule1.positions - x, axis=1)
            score = np.sum(1.0 / distances**p)  # Inverse of squared distances to penalize closer contacts more heavily
            return score

        #def objective_function_rotate(Vars):
        #    alpha, beta, gamma = Vars
        #    self.center_rotate(self.molecule2, self.molecule2[self.atom_2].position, alpha, beta, gamma)
        #    # We want to rotate the molecule such that the minimum distance between atom of interest in molecule 1, and other atoms in molecule 2 (not the atom of interest) is maximized, to reduce hindrance.
        #    distances = self.locate_hindrance_molecule2()
        #    # Lennard-Jones Repulsion
        #    score = np.sum(1/np.array(distances)**6)  # Inverse of squared distances to penalize closer contacts more heavily
        #    return score
        def objective_function_rotate(Vars):
            projections = self.rotate_to_optimal_direction(Vars)
            return -self.softmin(projections)  # maximize minimum projection
        # Bounds for the decision variables
        upper_bound_translate = [ np.pi/2,  np.pi]  
        lower_bound_translate = [-np.pi/2, -np.pi]
        bounds_translate = optimize.Bounds(lower_bound_translate, upper_bound_translate)
        upper_bound_rotate = [ 100, 100]#, 100]  
        lower_bound_rotate = [-100, -100]#, -100]
        bounds_rotate = optimize.Bounds(lower_bound_rotate, upper_bound_rotate)
        # Initial guess
        initial_guess_translate = [0.0, 0.0]
        initial_guess_rotate = [0.0, 0.0]#, 0.0]
        # Optimize translation first
        result_translate = optimize.minimize(objective_function_translate, x0=initial_guess_translate, bounds=bounds_translate, method='L-BFGS-B', options={'maxiter': 100, 'ftol': 1e-6})
        
        # Apply the best translation to the molecule before optimizing rotation
        best_step_theta, best_step_phi = result_translate.x
        
        # Find the vector to apply the translation to molecule 2
        self.translate_on_sphere(self.molecule2, self.molecule2[self.atom_2].position, self.sphere_center, best_step_theta, best_step_phi)

        # Then optimize rotation from the new position
        result_rotate = optimize.minimize(objective_function_rotate, x0=initial_guess_rotate, bounds=bounds_rotate, method='L-BFGS-B', options={'maxiter': 100, 'ftol': 1e-6})
        return result_translate, result_rotate

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
        translator = TranslateMoleculeSphere(molecule1, molecule2, atom1_index, atom2_index)
        result_translate, result_rotate = translator.opt_scipy_2()
        print(result_rotate)
        best_step_theta, best_step_phi = result_translate.x
        best_alpha, best_beta = result_rotate.x
        translator.step_on_sphere(translator.molecule2[translator.atom_2].position, translator.molecule1[translator.atom_1].position, best_step_theta, best_step_phi)
        # Apply the best rotation to molecule 2
        translator.rotate_to_optimal_direction(result_rotate.x)

        #translator.rotate_to_optimal_direction(result_rotate.x)
        #translator.drone_translate(best_step_theta, best_step_phi, best_alpha, best_beta, best_gamma)
        ## Optionally, perform a second optimization step to fine-tune the rotation
        #result_rotate = translator.optimization_rotate()
        #best_alpha, best_beta, best_gamma = result_rotate['Vars'][0]
        #translator.center_rotate(translator.molecule2, translator.molecule2[translator.atom_2].position, best_alpha, best_beta, best_gamma)
        merged_molecule = translator.molecule1 + translator.molecule2
        write(f'optimized_merged_molecule_{atom1_index}_{atom2_index}.xyz', merged_molecule,)
        # Append to trajectory file
        with open(traj_file, 'a') as traj:
            traj.write(f"{len(merged_molecule)}\n")
            traj.write(f"Optimized merged molecule for atom pair ({atom1_index}, {atom2_index})\n")
            for atom in merged_molecule:
                traj.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")