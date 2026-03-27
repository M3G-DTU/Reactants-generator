# Class for filtering descripters based on user-defined criteria

def get_atom_pairs(molecule1, molecule2, method='Hirshfeld', filter_descriptor='dual', max_pairs_per_molecule: int = 5, include_hydrogens: bool = False) -> list:
    # Get the descriptors for molecule 1
    if method == 'FMO':
        descriptor1 = DescriptorFMO(molecule1)
        descriptor2 = DescriptorFMO(molecule2)
    else:
        descriptor1 = DescriptorFDL(molecule1)
        descriptor2 = DescriptorFDL(molecule2)

    descriptors1 = descriptor1.extract_descriptors()
    descriptors2 = descriptor2.extract_descriptors()

    # Filter the descriptors to get the most electrophilic and nucleophilic atoms
    filter1 = DescriptorFilter(descriptors1, method=method, filter_descriptor=filter_descriptor)
    filter2 = DescriptorFilter(descriptors2, method=method, filter_descriptor=filter_descriptor)
    filtered_atoms_elec_1, filtered_atoms_nuc_1 = filter1.get_filtered_atoms()
    filtered_atoms_elec_2, filtered_atoms_nuc_2 = filter2.get_filtered_atoms()

    atom_pairs = []
    count = 0
    # Molecule 1 is the electrophile and molecule 2 is the nucleophile
    for atom_elec in filtered_atoms_elec_1:
        for atom_nuc in filtered_atoms_nuc_2:
            if include_hydrogens == False and (atom_elec['Atom'].startswith('H(') or atom_nuc['Atom'].startswith('H(')):
                continue

            if count >= max_pairs_per_molecule:
                break
            count += 1
            atom1_name, atom2_name = atom_elec['Atom'], atom_nuc['Atom']
            atom1_idx, atom2_idx = int(atom1_name.split('(')[1].rstrip(')')) - 1, int(atom2_name.split('(')[1].rstrip(')')) - 1
            atom_pairs.append((atom1_idx, atom2_idx))
    # Molecule 2 is the electrophile and molecule 1 is the nucleophile
    for atom_elec in filtered_atoms_elec_2:
        for atom_nuc in filtered_atoms_nuc_1:
            if include_hydrogens == False and (atom_elec['Atom'].startswith('H(') or atom_nuc['Atom'].startswith('H(')):
                continue
            
            if count >= max_pairs_per_molecule * 2:
                break
            count += 1
            atom1_name, atom2_name = atom_elec['Atom'], atom_nuc['Atom']
            atom1_idx, atom2_idx = int(atom1_name.split('(')[1].rstrip(')')) - 1, int(atom2_name.split('(')[1].rstrip(')')) - 1
            atom_pairs.append((atom1_idx, atom2_idx))
    return atom_pairs

def get_sorted_indexes(directory: dict, method: str, key_name: str, reverse: bool = False) -> list:
    return sorted(
        range(len(directory[method][key_name])),
        key=lambda i: directory[method][key_name][i],
        reverse=reverse
    )

class DescriptorFilter:
    def __init__(self, descriptors: dict, method: str = 'FMO', filter_descriptor: str = 'ElNuc', fukui_threshold: float = 0.05):
        self.descriptors = descriptors
        self.method = method
        self.filter_descriptor = filter_descriptor # Either 'fukui', 'dual', 'ElNuc'
        self.filtered_atoms_elec = []
        self.filtered_atoms_nuc = []
        self.fukui_threshold = fukui_threshold # Threshold for fukui descriptor, can be adjusted based on the specific use case
    # Different filtering function for different descriptor types
    # We check the descriptor type and apply the corresponding filtering function
    def filter_descriptors(self):
        if self.filter_descriptor == 'ElNuc':
            self.filter_electrophilicity()
            self.filter_nucleophilicity()
        elif self.filter_descriptor == 'fukui':
            self.filter_fukui_electrophilicity()
            self.filter_fukui_nucleophilicity()
        elif self.filter_descriptor == 'dual':
            self.filter_dual_electrophilicity()
            self.filter_dual_nucleophilicity()
        else:
            raise ValueError(f"Unknown filter descriptor: {self.filter_descriptor}")

    # For Electrophilicity and Nucleophilicity
    def filter_electrophilicity(self):
        # Get the sorted indexes based on the electrophilicity descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'Electrophilicity', reverse=True)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['Electrophilicity'] > 0 and atom['Nucleophilicity'] + atom['Electrophilicity'] > 0: # Criteria for electrophilic atoms
                self.filtered_atoms_elec.append(atom)

    def filter_nucleophilicity(self):
        # Get the sorted indexes based on the nucleophilicity descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'Nucleophilicity', reverse=True)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['Nucleophilicity'] < 0 and atom['Nucleophilicity'] + atom['Electrophilicity'] < 0: # Criteria for nucleophilic atoms
                self.filtered_atoms_nuc.append(atom)
    
    # Using Fukui function as descriptors
    def filter_fukui_electrophilicity(self):
        # Get the sorted indexes based on the fukui electrophilicity descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'f+', reverse=True)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['f+'] > self.fukui_threshold: # Criteria for electrophilic atoms based on fukui descriptor
                self.filtered_atoms_elec.append(atom)

    def filter_fukui_nucleophilicity(self):
        # Get the sorted indexes based on the fukui nucleophilicity descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'f-', reverse=True)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['f-'] > self.fukui_threshold: # Criteria for nucleophilic atoms based on fukui descriptor
                self.filtered_atoms_nuc.append(atom)

    # Using fukui dual descriptor
    def filter_dual_electrophilicity(self):
        # Get the sorted indexes based on the dual descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'f(2)', reverse=True)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['f(2)'] > 0: # Criteria for electrophilic atoms based on dual descriptor
                self.filtered_atoms_elec.append(atom)

    def filter_dual_nucleophilicity(self):
        # Get the sorted indexes based on the dual descriptor
        sorted_indexes = get_sorted_indexes(self.descriptors, self.method, 'f(2)', reverse=False)
        for i in sorted_indexes:
            atom = {key: self.descriptors[self.method][key][i] for key in self.descriptors[self.method]}
            if atom['f(2)'] < 0: # Criteria for nucleophilic atoms based on dual descriptor
                self.filtered_atoms_nuc.append(atom)

    # Return two sets of filtered atoms, one for electrophilic and one for nucleophilic atoms
    # Based on what criteria the user wants to apply, they can choose to filter based on electrophilicity, nucleophilicity, fukui descriptors or dual descriptors
    def get_filtered_atoms(self) -> list:
        self.filter_descriptors()
        return self.filtered_atoms_elec, self.filtered_atoms_nuc