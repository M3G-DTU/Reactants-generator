# This code is the base descriptor class, which is used to store the information of the descriptors and the methods to calculate them.
# The specific descriptor classes will inherit from this base class and implement the methods to calculate or obtain descriptors from different sources (e.g. AMS, xTB, etc.)

from abc import ABC, abstractmethod


class BaseDescriptor(ABC):
    """Base class for all descriptor classes.
    Each descriptor class should be able to generate the descriptor values for each atom in a molecule
    and store the descriptor values in a dictionary format."""
    
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.descriptors = {}
        
    @abstractmethod
    def extract_descriptors(self):
        """Extract the descriptor values from the loaded file and store them in the descriptors dictionary."""
        pass
    
    def get_descriptors(self, name: str) -> dict:
        """Get the descriptor values for a specific descriptor."""
        return self.descriptors.get(name, {})
    
    def get_all_descriptors(self) -> dict:
        """Get all descriptor values."""
        return self.descriptors
