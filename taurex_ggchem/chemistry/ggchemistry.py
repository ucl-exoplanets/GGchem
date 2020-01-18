from taurex.chemistry import Chemistry
import numpy as np
import taurex_ggchem.external.fort_ggchem as fchem 
class GGChem(Chemistry):


    allowed_elements = ['el','H','He','Li','Be','B','C','N','O','F','Ne',
                        'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',
                        'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
                        'As','Se','Br','Kr','Rb','Sr','Y','Zr','W']
    element_index_dict = dict(zip(allowed_elements[1:],range(1,len(allowed_elements))))

    def __init__(self, dispol_files = None, 
                 elements = ['H', 'C', 'N', 'O'],
                 abundances = [0.25,0.25,0.25,0.25],
                 equilibrium_condensation=False, dustchem_file):
        
        self._elements = elements
        self._element_number = [self.allowed_elements.index(e) for e in elements if e is not 'el']
        self._charge = 'el' in self._elements
        

        fchem.fort_ggchem.init_lean()
        fchem.parameters.elements = ' '.join([s.ljust(2) for s in self._elements if s is not 'el']).ljust(200)
        dispol = [s.ljust(200) for s in dispol_files]  

        elements = fchem.dust_data.elnam.tostring().decode('utf-8')  
        elements = [elements[i:i+2].strip() for i in range(0, len(elements ), 2)] 

        self._elements = elements
        self._abundances = np.array(abundances)

        fchem.fort_ggchem.init_taurex_chemistry(dispol_files,self._charge)

        _mols = fchem.fort_ggchem.copy_molecule_names(fchem.chemistry.nmole) 
        mols = np.lib.stride_tricks.as_strided(_mols, strides=(_mols.shape[1],1))

        mols = [s[0].decode('utf-8') for s in np.char.strip(mols.view('S20'))]

        self._molecules = mols
        self.restore_molecule_names()

        self.info('Loading dustchem file %s',dustchem_file)
        fchem.dust_data.dustchem_file = dustchem_file.ljust(200)   
        fchem.init_dustchem_taurex()


    def restore_molecule_names(self):
        import re
        elm = [s for s in self._elements if len(s)>1]
        elem_upp = [s.upper() for s in self._elements if len(s)>1]

        rep = dict(zip(elem_upp, elm))

        rep = dict((re.escape(k), v) for k, v in rep.items()) 
        pattern = re.compile("|".join(rep.keys()))

        self._molecules = [pattern.sub(lambda m: rep[re.escape(m.group(0))], s) for s in self._molecules]

    def initialize_chemistry(self, nlayers=100, temperature_profile=None,
                             pressure_profile=None, altitude_profile=None):
        
        fchem.structure.tgas[:nlayers] = temperature_profile
        fchem.structure.press[:nlayers] = pressure_profile*10 # to dyn/cm2
        fchem.structure.estruc[:nlayers,:] = fchem.dust_data.eps0

        self.info('Running GGChem equilibrium code...')

