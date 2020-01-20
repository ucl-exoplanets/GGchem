from taurex.chemistry import Chemistry
import numpy as np
import taurex_ggchem.external.fort_ggchem as fchem 


default_elements = [s.strip() for s in [' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na','Mg','Al','Si',' P',' S',
            'Cl','Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
            'Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
            ' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
            'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Th',' U']]

default_abundances = [9.271E-01,7.159E-02,1.077E-11,1.382E-11,2.305E-10,3.112E-04,8.895E-05,7.008E-04,3.279E-08,6.174E-05,
              2.168E-06,3.588E-05,2.771E-06,3.992E-05,2.816E-07,1.554E-05,2.811E-07,2.183E-06,1.275E-07,2.176E-06,
              1.109E-09,1.041E-07,9.783E-09,4.792E-07,2.268E-07,2.231E-05,8.456E-08,1.698E-06,1.372E-08,3.810E-08,
              7.148E-10,3.430E-09,0.000E+00,0.000E+00,0.000E+00,0.000E+00,4.373E-10,7.110E-10,1.401E-10,5.463E-10,
              5.364E-11,1.169E-10,6.163E-11,2.421E-11,3.512E-11,1.155E-11,6.650E-11,4.340E-11,9.446E-11,1.023E-11,
              0.000E+00,0.000E+00,0.000E+00,7.499E-11,9.072E-11,1.794E-11,3.557E-11,8.842E-12,2.591E-11,8.286E-12,
              4.099E-12,1.585E-11,7.839E-13,1.533E-11,0.000E+00,7.449E-12,1.475E-12,7.200E-12,7.121E-12,6.980E-12,
              0.000E+00,2.711E-11,6.691E-13,1.310E-11,1.296E-11,5.748E-11,6.325E-05,1.242E-10,6.096E-12,6.013E-11,
              5.962E-11,1.611E-12,5.234E-12]


selected = ['H', 'He', 'C', 'N', 'O', 'Na', 'Mg', 'Si', 'Fe', 'Al', 'Ca', 'Ti', 'S', 'Cl', 'K', 'Li', 'F', 'P', 'V', 'Cr', 'Mn', 'Ni', 'Zr', 'W']
sel_ab = [default_abundances[default_elements.index(s)] for s in selected if s is not 'el']

class GGChem(Chemistry):


    allowed_elements = ['el','H','He','Li','Be','B','C','N','O','F','Ne',
                        'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',
                        'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
                        'As','Se','Br','Kr','Rb','Sr','Y','Zr','W']

    element_index_dict = dict(zip(allowed_elements[1:],range(0,len(allowed_elements))))

    def __init__(self, dispol_files = None, 
                 elements = selected,
                 abundances = sel_ab,
                 equilibrium_condensation=False, dustchem_file=None,Tfast=1000,
                 derived_ratios=['C/O']):
        super().__init__(self.__class__.__name__)

        self._charge = 'el' in elements
        
        clean_elements = []
        clean_abundances = []


        t_elements = [s.strip() for s in elements if s is not 'el']
        self._passed_elements = t_elements
        

        for el,ab in zip(t_elements,abundances):
            if el in self.allowed_elements:
                clean_elements.append(el)
                clean_abundances.append(ab)


        if self._charge:
            clean_elements.append('el')

        fchem.fort_ggchem.init_lean()
        fchem.parameters.elements = ' '.join([s.ljust(2) for s in clean_elements]).ljust(200)
        dispol = [s.ljust(200) for s in dispol_files]  

        self._abundances = np.array(clean_abundances)
        self.update_abundances()
        fchem.fort_ggchem.init_taurex_chemistry(dispol,self._charge)

        _atms = fchem.fort_ggchem.copy_atom_names(fchem.chemistry.nelm) 
        atms = np.lib.stride_tricks.as_strided(_atms, strides=(_atms.shape[1],1))

        atms = [s[0].decode('utf-8') for s in np.char.strip(atms.view('S2'))]

        self._elements = atms
        # print(self._elements)
        # quit()
        fchem.parameters.model_eqcond = equilibrium_condensation

        fchem.parameters.tfast = Tfast

        _mols = fchem.fort_ggchem.copy_molecule_names(fchem.chemistry.nmole) 
        mols = np.lib.stride_tricks.as_strided(_mols, strides=(_mols.shape[1],1))

        mols = [s[0].decode('utf-8') for s in np.char.strip(mols.view('S20'))]

        self._molecules = mols
        self.restore_molecule_names()

        self.info('Loading dustchem file %s',dustchem_file)
        fchem.dust_data.dustchem_file = dustchem_file.ljust(200)   
        fchem.init_dustchem_taurex()
        self.update_abundances()

        self._setup_active_inactive()

    # def compute_mu_profile(self, nlayers):
    #     from taurex.util import get_molecular_weight
    #     from taurex.constants import AMU
    #     """
    #     Computes molecular weight of atmosphere
    #     for each layer

    #     Parameters
    #     ----------
    #     nlayers: int
    #         Number of layers
    #     """

    #     self.mu_profile = np.zeros(shape=(nlayers,))
    #     if self.activeGasMixProfile is not None:
    #         for idx, gasname in enumerate(self.activeGases):
    #             self.mu_profile += self.activeGasMixProfile[idx] * \
    #                 get_molecular_weight(gasname)
    #     if self.inactiveGasMixProfile is not None:
    #         for idx, gasname in enumerate(self.inactiveGases):
    #             print(gasname,get_molecular_weight(gasname)/AMU)
    #             self.mu_profile += self.inactiveGasMixProfile[idx] * \
    #                 get_molecular_weight(gasname)
    #             print(self.inactiveGasMixProfile[idx,-1])
    #             print(self.mu_profile[-1]/AMU)

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


    def update_abundances(self):
        
        fchem.dust_data.muh = 0.0

        Habun = self._abundances[self._passed_elements.index('H')]

        for mol,abundance in zip(self._passed_elements,self._abundances):
            self.info('%s %s',mol,abundance)
            try:
                mol_idx = int(getattr(fchem.chemistry,mol.lower()))-1
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