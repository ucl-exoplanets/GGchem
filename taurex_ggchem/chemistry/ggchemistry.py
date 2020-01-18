from taurex.chemistry import Chemistry
import numpy as np

class GGChem(Chemistry):


    allowed_elements = ['el','H','He','Li','Be','B','C','N','O','F','Ne',
                        'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',
                        'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
                        'As','Se','Br','Kr','Rb','Sr','Y','Zr','W']


    def __init__(self, dispol_files = None, 
                 elements = ['H', 'C', 'N', 'O'],
                 abundances = [0.25,0.25,0.25,0.25],
                 equilibrium_condensation=False)
        
        self._elements = elements
        self._element_number = [self.allowed_elements.index(e) for e in elements if e is not 'el']
        self._charge = 'el' in self._elements

