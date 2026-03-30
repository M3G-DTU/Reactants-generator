from ..base_descriptor import BaseDescriptor
from morfeus import XTB, read_xyz

class DescriptorXTB(BaseDescriptor):
    def __init__(self, file_path: str):
        super().__init__(file_path)
        # we only take .xyz files as input for xTB
        if not self.file_path.endswith(".xyz"):
            raise ValueError("Only .xyz files are supported for xTB descriptor extraction.")
    # Add functions
    def extract_descriptors(self) -> dict:
        elements, coordinates = read_xyz(self.file_path)
        xtb_calc = XTB(elements, coordinates, method='2', charge=0, n_unpaired=0)
        a_list = []
        for i, atom in enumerate(elements):
            a_list.append(f"{atom}({i+1})")
        electrophilicity = xtb_calc.get_fukui('local_electrophilicity').values()
        nucleophilicity = xtb_calc.get_fukui('local_nucleophilicity').values()
        f_plus = xtb_calc.get_fukui('electrophilicity').values()
        f_minus = xtb_calc.get_fukui('nucleophilicity').values()
        f_dual = xtb_calc.get_fukui('dual').values()
        f_radical = xtb_calc.get_fukui('radical').values()
        self.descriptors['xTB'] = {'Atom': a_list,
                                   'Electrophilicity': list(electrophilicity),
                                   'Nucleophilicity': list(nucleophilicity),
                                   'f+': list(f_plus),
                                   'f-': list(f_minus),
                                   'f(2)': list(f_dual),
                                   'f0': list(f_radical)}
        return self.descriptors