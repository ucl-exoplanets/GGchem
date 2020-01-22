from taurex.chemistry import Chemistry
import numpy as np
import taurex_ggchem.external.fort_ggchem as fchem 
import os

selected = ['C', 'N', 'O', 'Na', 'Mg', 'Si', 'Fe', 'Al', 'Ca', 'Ti', 'S', 'Cl', 'K', 'Li', 'F', 'P', 'V', 'Cr', 'Mn', 'Ni', 'Zr', 'W']

class GGChem(Chemistry):


    allowed_elements = ['el','H','He','Li','Be','B','C','N','O','F','Ne',
                        'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',
                        'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
                        'As','Se','Br','Kr','Rb','Sr','Y','Zr','W']

    element_index_dict = dict(zip(allowed_elements[1:],range(0,len(allowed_elements))))

    def __init__(self, dispol_files = None, 
                 abundance_profile='solar',
                 fill_elements = ['H','He'],
                 fill_ratios = 0.5,
                 include_charge = True,
                 trace_elements = selected,
                 equilibrium_condensation=False, 
                 dustchem_file=None,
                 Tfast=1000):
        super().__init__(self.__class__.__name__)

        elements = fill_elements + trace_elements
        if include_charge:
            elements.append('el')
         

        self._base_data_path = os.path.join(os.path.abspath(os.path.dirname(fchem.__file__)),'data')

        self._charge = 'el' in elements

        self.info('Using GGChem chemistry')

        ignored_elements = [s.strip() for s in elements if s not in self.allowed_elements]

        if len(ignored_elements) > 0:
            self.info('Elements ignored: %s',ignored_elements)

        if self._charge:
            self.info('Ions activated')


        elements = [s.strip() for s in elements if s in self.allowed_elements]

        self._passed_elements = elements

        fchem.fort_ggchem.init_lean()

        self.setup_abundances(abundance_profile)


        fchem.parameters.elements = ' '.join([s.ljust(2) for s in elements]).ljust(200)

        if dispol_files is None:
            dispol_files = [os.path.join(self._base_data_path,s) for s in ['dispol_BarklemCollet.dat',
                                                                           'dispol_StockKitzmann_withoutTsuji.dat',
                                                                           'dispol_WoitkeRefit.dat']]
                                                            

        dispol = [s.ljust(200) for s in dispol_files]  


        self.info('Initializing GGchemistry')
        fchem.fort_ggchem.init_taurex_chemistry(dispol,self._charge)

        self.update_abundances()

        _atms = fchem.fort_ggchem.copy_atom_names(fchem.chemistry.nelm) 
        atms = np.lib.stride_tricks.as_strided(_atms, strides=(_atms.shape[1],1))

        atms = [s[0].decode('utf-8') for s in np.char.strip(atms.view('S2'))]

        self._elements = atms
        self._passed_elements = atms

        fchem.parameters.model_eqcond = equilibrium_condensation

        fchem.parameters.tfast = Tfast

        self.info('T-fast %s K',Tfast)
        self.info('Eq. Condensation: %s', equilibrium_condensation)

        _mols = fchem.fort_ggchem.copy_molecule_names(fchem.chemistry.nmole) 
        mols = np.lib.stride_tricks.as_strided(_mols, strides=(_mols.shape[1],1))

        mols = [s[0].decode('utf-8') for s in np.char.strip(mols.view('S20'))]

        self._molecules = mols
        self.restore_molecule_names()

        if dustchem_file is None:
            dustchem_file = os.path.join(self._base_data_path,'DustChem.dat')

        self.info('Loading dustchem file %s',dustchem_file)
        fchem.dust_data.dustchem_file = dustchem_file.ljust(200)   
        fchem.init_dustchem_taurex()
        self.update_abundances()

        self._setup_active_inactive()


    def setup_abundances(self, profile):
        

        self.info('Selected abundance profile %s',profile)
        abundance_file = os.path.join(self._base_data_path,'Abundances.dat')
        chosen = 4 # 4 - earth 5 - ocean 6 - solar 7- meteorite
        profile = profile.lower()

        if profile in ('earthcrust',):
            chosen = 4
        elif profile in ('ocean',):
            chosen = 5
        elif profile in ('sun', 'star','solar',):
            chosen = 6
        elif profile in ('meteor', 'meteorites',):
            chosen = 7

        elements = []
        abundances = []

        with open(abundance_file,'r') as f:
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()

            for line in f:
                split_line = line.split()
                el = split_line[2].strip()

                if el not in self.allowed_elements:
                    continue
                elements.append(el)
                ab = float(split_line[chosen])
                abundances.append(ab)

        self.update_abundances(elements=elements,abundances=abundances)
        
        self._abundances = np.array([a for e,a in zip(elements,abundances) if e in self._passed_elements])


    def get_element_index(self, element):
        return int(getattr(fchem.chemistry,element.strip().lower()))-1



    def _setup_active_inactive(self):

        active_mols = []
        inactive_mols =[]
        active_indices = []
        inactive_indices = []
        nmol = fchem.chemistry.nmole 

        for idx,m in enumerate(self._elements + self._molecules):
            if m in self._avail_active:
                active_mols.append(m)
                active_indices.append(idx)
            else:
                inactive_mols.append(m)
                inactive_indices.append(idx)              

        self.info('Active molecules %s',active_mols)
        self.info('Inactive molecules %s',inactive_mols)

        self._active_index = np.array(active_indices,dtype=np.integer)
        self._inactive_index = np.array(inactive_indices,dtype=np.integer)

        self._active_mols = active_mols
        self._inactive_mols = inactive_mols

    
    def restore_molecule_names(self):
        import re
        elm = [s for s in self._passed_elements if len(s)>1]
        elem_upp = [s.upper() for s in self._passed_elements if len(s)>1]

        rep = dict(zip(elem_upp, elm))

        rep = dict((re.escape(k), v) for k, v in rep.items()) 
        pattern = re.compile("|".join(rep.keys()))

        self._molecules = [pattern.sub(lambda m: rep[re.escape(m.group(0))], s) for s in self._molecules]


    def update_abundances(self,elements=None,abundances=None):
        
        fchem.dust_data.muh = 0.0

        if elements is None:
            elements = self._passed_elements
            abundances = self._abundances
        


        Habun = abundances[elements.index('H')]

        for mol,abundance in zip(elements,abundances):
            self.info('%s %s',mol,abundance)
            try:
                mol_idx = self.get_element_index(mol)
                #mol_idx = self.element_index_dict[mol]
                fchem.dust_data.eps0[mol_idx] = abundance/Habun
                fchem.dust_data.muh += fchem.dust_data.mass[mol_idx]*abundance/Habun
                self.info('%s %s',fchem.dust_data.mass[mol_idx]*abundance,fchem.dust_data.muh)
            except AttributeError:
                pass 

    def initialize_chemistry(self, nlayers=100, temperature_profile=None,
                             pressure_profile=None, altitude_profile=None):
        from taurex.constants import KBOLTZ

        fchem.structure.tgas[:nlayers] = temperature_profile
        fchem.structure.press[:nlayers] = pressure_profile*10 # to dyn/cm2

        fchem.structure.estruc[:nlayers,:] = fchem.dust_data.eps0
        self.info('Running GGChem equilibrium code...')
        nmol = fchem.chemistry.nmole 
        nelem = fchem.chemistry.nelm 
        mols = fchem.fort_ggchem.run_ggchem(nlayers,nelem+nmol)
        self._mols = mols
        self._vmr = mols/np.sum(mols,axis=0)
        #self._vmr = mols/1e-6 #m-3
        #self._vmr /= pressure_profile/(KBOLTZ*temperature_profile)
        self.info('Computing mu Profile')
        self.compute_mu_profile(nlayers)



    @property
    def activeGases(self):
        return self._active_mols

    @property
    def inactiveGases(self):
        return self._inactive_mols

    @property
    def activeGasMixProfile(self):
        return self._vmr[self._active_index]

    @property
    def inactiveGasMixProfile(self):
        return self._vmr[self._inactive_index]

    def ratio(self, element_a, element_b):
        element_a = element_a.strip()
        element_b = element_b.strip()

        element_a_idx = int(getattr(fchem.chemistry,element_a.lower()))-1
        element_b_idx = int(getattr(fchem.chemistry,element_b.lower()))-1

        return float(fchem.dust_data.eps0[element_a_idx]/fchem.dust_data.eps0[element_b_idx])

    @classmethod
    def input_keywords(cls):
        return ['ggchem', ]