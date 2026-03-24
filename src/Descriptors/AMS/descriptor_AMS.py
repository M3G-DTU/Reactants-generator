from src.Descriptors.base_descriptor import BaseDescriptor

class DescriptorFMO(BaseDescriptor):
    def __init__(self, file_path: str):
        super().__init__(file_path)
        
    # Add functions
    def extract_descriptors(self) -> dict:
        # Add FMO to the descriptors dictionary
        self.descriptors['FMO'] = {}
        with open(self.file_path, "r", encoding='utf-8', errors='replace') as f:
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
                if header not in self.descriptors['FMO']:
                    self.descriptors['FMO'][header] = []
                if header == 'Atom':
                    self.descriptors['FMO'][header].append(value)
                else:
                    self.descriptors['FMO'][header].append(float(value))
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
                if header not in self.descriptors['FMO']:
                    self.descriptors['FMO'][header] = []
                self.descriptors['FMO'][header].append(float(value))
        return self.descriptors
    
class DescriptorFDL(BaseDescriptor):
    def __init__(self, file_path: str):
        super().__init__(file_path)

    def extract_descriptors(self) -> dict:
        # Locate relevant sections in the output file
        locate_sections = ['Hirshfeld Partitioning',
                           'Mulliken Partitioning', 
                           'Voronoï Partitioning',
                           ]
        for section in locate_sections:
            self.descriptors[section] = {}

            with open(self.file_path, "r", encoding='utf-8', errors='replace') as f:
                lines = f.readlines()

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
                    if header not in self.descriptors[section]:
                        self.descriptors[section][header] = []
                    if header == 'Atom':
                        self.descriptors[section][header].append(value)
                    else:
                        self.descriptors[section][header].append(float(value))
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
                    if header not in self.descriptors[section]:
                        self.descriptors[section][header] = []

                    self.descriptors[section][header].append(float(value))
        # Remove the "Partitioning" part to make the keys more concise
        for section in locate_sections:
            section_key = section.replace(" Partitioning", "")
            self.descriptors[section_key] = self.descriptors.pop(section)
        return self.descriptors
    