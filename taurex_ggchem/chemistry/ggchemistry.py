from taurex.chemistry import Chemistry
import numpy as np
import taurex_ggchem.fort_ggchem as fchem
import os
from taurex.core import fitparam
import typing as t
import pathlib
import numpy.typing as npt

# from taurex.util.fortran import FortranSafeCaller
import importlib.resources as ires
from taurex.util.fortran import SafeFortranCaller, FortranStopException
from taurex.exceptions import InvalidModelException

PathType = t.Union[str, pathlib.Path, os.PathLike]

selected: t.List[str] = [
    "H",
    "He",
    "C",
    "N",
    "O",
    "Na",
    "Mg",
    "Si",
    "Fe",
    "Al",
    "Ca",
    "Ti",
    "S",
    "Cl",
    "K",
    "Li",
    "F",
    "P",
    "V",
    "Cr",
    "Mn",
    "Ni",
    "Zr",
    "W",
]


class GGChemNonSafe(Chemistry):
    """GGChem chemistry class.

    This implements GGChem through spawning a Fortran process and calling the GGChem Fortran code.
    The reason for this is the ``STOP`` command that GGChem uses to halt the code when an error occurs
    which essentially kills the Python process. We use a FortranSafeCaller which spawns a process to
    catch this error and reinitialize the GGChem code. Its not perfect but stops the code from crashing.


    """

    allowed_elements: t.List[str] = [
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "W",
    ]

    element_index_dict: t.Dict[str, int] = dict(
        zip(allowed_elements[1:], range(0, len(allowed_elements)))
    )

    def __init__(
        self,
        dispol_files: t.Optional[t.List[PathType]] = None,
        abundance_profile: t.Optional[str] = None,
        selected_elements: t.Optional[t.List[str]] = selected,
        ratio_elements: t.Optional[t.List[str]] = None,
        ratios_to_O: t.Optional[npt.ArrayLike] = None,
        he_h_ratio: t.Optional[float] = 0.083,
        metallicity: t.Optional[float] = 1.0,
        include_charge: t.Optional[bool] = False,
        equilibrium_condensation: t.Optional[bool] = False,
        dustchem_file: t.Optional[PathType] = None,
        Tfast: t.Optional[float] = 1000,
        new_back_it: t.Optional[int] = 6,
        new_back_fac: t.Optional[float] = 1e5,
        new_pre_method: t.Optional[int] = 2,
        new_full_it: t.Optional[bool] = False,
        new_fast_level: t.Optional[int] = 1,
    ) -> None:
        super().__init__(self.__class__.__name__)

        self._init_dispol_files = dispol_files
        self._init_abundance_profile = abundance_profile
        self._init_selected_elements = selected_elements
        self._init_ratio_elements = ratio_elements
        self._init_ratios_to_O = ratios_to_O
        self._init_he_h_ratio = he_h_ratio
        self._init_metallicity = metallicity
        self._init_include_charge = include_charge
        self._init_equilibrium_condensation = equilibrium_condensation
        self._init_dustchem_file = dustchem_file
        self._init_Tfast = Tfast
        self._init_new_back_it = new_back_it
        self._init_new_back_fac = new_back_fac
        self._init_new_pre_method = new_pre_method
        self._init_new_full_it = new_full_it
        self._init_new_fast_level = new_fast_level

        self.reinitialize_ggchem()
        self.restore_molecule_names()
        self._setup_active_inactive()
        self.add_ratio_params()

    def reinitialize_ggchem(self) -> None:
        """Initialize GGchem for the first time or after a STOP error."""
        dispol_files = self._init_dispol_files
        abundance_profile = self._init_abundance_profile
        selected_elements = self._init_selected_elements
        ratio_elements = self._init_ratio_elements
        ratios_to_O = self._init_ratios_to_O
        he_h_ratio = self._init_he_h_ratio
        metallicity = self._init_metallicity
        include_charge = self._init_include_charge
        equilibrium_condensation = self._init_equilibrium_condensation
        dustchem_file = self._init_dustchem_file
        Tfast = self._init_Tfast
        new_back_it = self._init_new_back_it
        new_back_fac = self._init_new_back_fac
        new_pre_method = self._init_new_pre_method
        new_full_it = self._init_new_full_it
        new_fast_level = self._init_new_fast_level

        elements = [s.strip() for s in selected_elements]
        if "O" not in elements:
            elements.append("O")
        if "He" not in elements:
            elements.insert(0, "He")
        if "H" not in elements:
            elements.insert(0, "H")

        self._he_h_ratio = he_h_ratio
        self._metallicity = metallicity
        self._base_data_path = ires.files("taurex_ggchem") / "data"

        self.info("Using GGChem chemistry")

        ignored_elements = [
            s.strip() for s in elements if s not in self.allowed_elements
        ]

        if len(ignored_elements) > 0:
            self.info("Elements ignored: %s", ignored_elements)

        elements = [s.strip() for s in elements if s in self.allowed_elements]
        if include_charge:
            self._charge = include_charge
            elements.append("el")
            self.info("Ions activated")
        else:
            self._charge = False
        fchem.chemistry.newfullit = new_full_it
        fchem.chemistry.newbackit = new_back_it
        fchem.chemistry.newbackfac = new_back_fac
        fchem.chemistry.newfastlevel = new_fast_level
        fchem.chemistry.newpremethod = new_pre_method

        self._passed_elements = elements
        self._elements = elements
        fchem.parameters.verbose = 0
        self._initial_abundances = fchem.fort_ggchem.init_lean()

        fchem.parameters.initchem_info = False

        if dispol_files is None:
            dispol_files = [
                self._base_data_path / s
                for s in [
                    "dispol_BarklemCollet.dat",
                    "dispol_StockKitzmann_withoutTsuji.dat",
                    "dispol_WoitkeRefit.dat",
                ]
            ]
        else:
            dispol_files = [pathlib.Path(s) for s in dispol_files]

        for dispol in dispol_files:
            if not dispol.exists():
                raise FileNotFoundError(f"Dispol file {dispol} not found")

        fchem.parameters.elements = " ".join([s.ljust(2) for s in elements]).ljust(200)

        print("ELEMENTS", fchem.parameters.elements)
        self.info("Elements in system: %s", elements)
        dispol = [str(s).ljust(200) for s in dispol_files]

        self.info("Initializing GGchemistry")

        fchem.fort_ggchem.init_taurex_chemistry(dispol, self._charge)

        self.setup_abundances(abundance_profile)
        self.setup_ratios(ratio_elements, ratios_to_O)

        fchem.parameters.model_eqcond = equilibrium_condensation

        fchem.parameters.tfast = Tfast

        self.info("T-fast %s K", Tfast)
        self.info("Eq. Condensation: %s", equilibrium_condensation)

        nmole = fchem.chemistry.nmole

        _mols = fchem.fort_ggchem.copy_molecule_names(nmole)
        self.debug("Molecules %s", _mols)

        mols = [s.decode("utf-8").strip() for s in _mols]
        self.debug("Molecules decoded %s", mols)
        self._molecules = mols

        if dustchem_file is None:
            dustchem_file = self._base_data_path / "DustChem.dat"

        self.info("Loading dustchem file %s", dustchem_file)
        fchem.dust_data.dustchem_file = str(dustchem_file).ljust(200)
        print("Start dustchem")
        fchem.init_dustchem_taurex()
        print("Chem initialized")
        ndust = fchem.dust_data.ndust

        _dust = fchem.fort_ggchem.copy_dust_names(ndust)
        dust = [s.decode("utf-8").strip() for s in _dust]
        self.debug("Dust decoded %s", dust)
        self._condensates = dust

        # quit()
        # self.update_abundances()

    @property
    def condensates(self):
        """
        Returns a list of condensates in the atmosphere.

        Returns
        -------
        active : :obj:`list`
            List of condensates

        """

        return self._condensates

    @property
    def condensateMixProfile(self):
        """
        **Requires implementation**

        Should return profiles of shape ``(ncondensates,nlayers)``.
        """
        if len(self.condensates) == 0:
            return None
        else:
            return self._dust

    def setup_abundances(self, profile: npt.NDArray[np.float64]) -> None:

        if profile is not None:
            self.info("Selected abundance profile %s", profile)
            abundance_file = self._base_data_path / "Abundances.dat"
            chosen = 4  # 4 - earth 5 - ocean 6 - solar 7- meteorite
            profile = profile.lower()

            if profile in ("earthcrust",):
                chosen = 4
            elif profile in ("ocean",):
                chosen = 5
            elif profile in (
                "sun",
                "star",
                "solar",
            ):
                chosen = 6
            elif profile in (
                "meteor",
                "meteorites",
            ):
                chosen = 7

            elements = []
            abundances = np.zeros(shape=(len(self.allowed_elements),))

            with open(abundance_file, "r") as f:
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
                    self.info("%s %s", el, ab)
                    self._initial_abundances[self.get_element_index(el)] = ab

        self._initial_abundances /= self._initial_abundances[0]
        self.info("H relative abundances.....")
        for el, ab in zip(self.allowed_elements, self._initial_abundances):
            self.info("%s %s", el, ab)

        # raise

    def setup_ratios(
        self, ratio_elements: t.List[str], ratios_to_O: npt.ArrayLike
    ) -> None:
        ratios_to_O = np.array(ratios_to_O)
        self._metal_idx = [
            self.get_element_index(e)
            for e in self._elements
            if e
            not in (
                "H",
                "He",
                "el",
            )
        ]
        self._metal_elements = [
            e
            for e in self._elements
            if e
            not in (
                "H",
                "He",
                "el",
            )
        ]
        self.info(
            "Metals detected %s", list(zip(self._metal_idx, self._metal_elements))
        )
        O_idx = self.get_element_index("O")
        self._O_abund = self._initial_abundances[O_idx]
        self._ratios = self._initial_abundances[self._metal_idx]
        self._ratios /= self._O_abund

        if ratio_elements is not None:
            for relm, rat in zip(ratio_elements, ratios_to_O):
                index = self._metal_elements.index(relm)
                self._ratios[index] = rat

        self.info("Metallic ratios against O")
        for met, rat in zip(self._metal_elements, self._ratios):
            self.info("%s %s", met, rat)

    def get_element_index(self, element: str) -> int:
        elnums = fchem.chemistry.elnum
        index = getattr(fchem.chemistry, f"{element.strip().lower()}", None)
        return elnums[index - 1] - 1
        # return self.get_element_index(element)

    def _setup_active_inactive(self) -> None:

        active_mols = []
        inactive_mols = []
        active_indices = []
        inactive_indices = []
        nmol = fchem.chemistry.nmole

        clean_elements = [el if el != "el" else "e-" for el in self._elements]

        for idx, m in enumerate(clean_elements + self._molecules):
            if m in self.availableActive:
                active_mols.append(m)
                active_indices.append(idx)
            else:
                inactive_mols.append(m)
                inactive_indices.append(idx)

        self.info("Active molecules %s", active_mols)
        self.info("Inactive molecules %s", inactive_mols)

        self._active_index = np.array(active_indices, dtype=np.int64)
        self._inactive_index = np.array(inactive_indices, dtype=np.int64)

        self._active_mols = active_mols
        self._inactive_mols = inactive_mols

    def restore_molecule_names(self) -> None:
        """Restores molecule names from GGChem style."""
        import re

        elm = [s for s in self._passed_elements if len(s) > 1]
        elem_upp = [s.upper() for s in self._passed_elements if len(s) > 1]

        rep = dict(zip(elem_upp, elm))

        rep = dict((re.escape(k), v) for k, v in rep.items())
        pattern = re.compile("|".join(rep.keys()))

        self._molecules = [
            pattern.sub(lambda m: rep[re.escape(m.group(0))], s)
            for s in self._molecules
        ]

    def update_abundances(self) -> npt.NDArray[np.float64]:
        """Updates the abundances of the system."""
        abundances = self._initial_abundances
        ratios = self._ratios
        O_abund = self._O_abund * self._metallicity
        abundances[self._metal_idx] = ratios * O_abund
        abundances /= abundances[self.get_element_index("H")]
        abundances[self.get_element_index("He")] = self._he_h_ratio

        return abundances

    def initialize_chemistry(
        self,
        nlayers: t.Optional[int] = 100,
        temperature_profile: t.Optional[npt.NDArray[np.float64]] = None,
        pressure_profile: t.Optional[npt.NDArray[np.float64]] = None,
        altitude_profile: t.Optional[npt.NDArray[np.float64]] = None,
    ):
        from taurex.constants import KBOLTZ

        ab = self.update_abundances()
        self.info("Running GGChem equilibrium code...")
        nmol = fchem.chemistry.nmole
        nelem = fchem.chemistry.nelm
        ndust = fchem.dust_data.ndust
        mols = []
        dust = []
        for t, p in zip(temperature_profile, pressure_profile):
            m, d = fchem.fort_ggchem.run_ggchem(nelem + nmol, ndust, t, p * 10, ab)

            mols.append(m)
            dust.append(d)
        self._mols = np.stack(mols).T
        self._dust = np.stack(dust).T
        # self._vmr = mols/np.sum(mols,axis=0)
        # self._vmr = self._mols
        # self._vmr = mols/1e-6 #m-3
        self._vmr = self._mols  # /= pressure_profile/(KBOLTZ*temperature_profile)
        self.info("Computing mu Profile")
        self.compute_mu_profile(nlayers)

    @fitparam(
        param_name="metallicity",
        param_latex="Z",
        default_bounds=[0.2, 2.0],
        default_fit=False,
    )
    def metallicity(self):
        return self._metallicity

    @metallicity.setter
    def metallicity(self, value):
        self._metallicity = value

    def add_ratio_params(self):

        for idx, element in enumerate(self._metal_elements):
            if element == "O":
                continue
            param_name = f"{element}_O_ratio"
            param_tex = f"{element}/O"
            id_val = self._metal_elements.index(element)

            def read_mol(self, idx=id_val):
                return self._ratios[idx]

            def write_mol(self, value, idx=id_val):
                self._ratios[idx] = value

            fget = read_mol
            fset = write_mol

            bounds = [1.0e-12, 0.1]

            default_fit = False
            self.add_fittable_param(
                param_name, param_tex, fget, fset, "linear", default_fit, bounds
            )

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

    def ratio(self, element_a: float, element_b: float) -> float:
        element_a = element_a.strip()
        element_b = element_b.strip()

        element_a_idx = self.get_element_index(element_a)
        element_b_idx = self.get_element_index(element_b)
        return float(
            fchem.dust_data.eps0[element_a_idx] / fchem.dust_data.eps0[element_b_idx]
        )

    @classmethod
    def input_keywords(cls):
        return [
            "ggchem-unsafe",
        ]

    BIBTEX_ENTRIES = [
        r"""
        @article{ ggchem,
            author = {{Woitke, P.} and {Helling, Ch.} and {Hunter, G. H.} and {Millard, J. D.} and {Turner, G. E.} and {Worters, M.} and {Blecic, J.} and {Stock, J. W.}},
            title = {Equilibrium chemistry down to 100 K - Impact of silicates and phyllosilicates on the carbon to oxygen ratio},
            DOI= "10.1051/0004-6361/201732193",
            url= "https://doi.org/10.1051/0004-6361/201732193",
            journal = {A\&A},
            year = 2018,
            volume = 614,
            pages = "A1",
        }
        """
    ]
