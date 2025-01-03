incdir_numpy = run_command(py3, [
    '-c',
    '''import os, numpy as np
try:
    print(os.path.relpath(np.get_include()))
except:
    print(np.get_include())
    ''',
], check: true).stdout().strip()
inc_np = include_directories(incdir_numpy)
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

# Dependencies for Fortran
fortranobject_lib = static_library(
    '_fortranobject',
    fortranobject_c,
    c_args: numpy_nodepr_api,
    dependencies: py3_dep,
    include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
    link_with: fortranobject_lib,
    include_directories: [inc_np, inc_f2py],
)

command_args = [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@', '--lower', '--backend', 'meson']

# Add '-fdefault-real-8', '-fdefault-double-8' and -O3 to the compile args
# for the Fortran code

# Add the Fortran arguments to the command


ggchem_module = custom_target(
    'fort_ggchemmodule',
    output : ['fort_ggchemmodule.c', 'fort_ggchem-f2pywrappers2.f90', 'fort_ggchem-f2pywrappers.f'],
    input : '../glue/fort_ggchem.pyf',
    command: command_args
)


py3.extension_module('fort_ggchem',
  [

                        '../src16/database.f',
                        '../src16/is_nan.F',
                        '../src16/nasa_polynomial.f',
                        '../src16/ggchem.f',
                        '../src16/gauss16.f',
                        '../src16/gauss8.f',
                        '../src16/gauss_nm.f',
                        '../src16/massfrac.f',
                        '../src16/linpack_q.f90',
                        '../src16/slatec_routines.f',
                        '../src16/smchem16.f',
                        '../src16/smchem8.f',
                        '../src16/supersat.f',
                        '../src16/stindex.f',
                        '../src16/upper.f',
                        '../src16/nucleation.f',
                        '../src16/equil_cond.f',
                        '../glue/init_dustchem_taurex.f',
                        '../glue/datamod.f90',
                        '../glue/taurex_glue.f90',
                        
  ggchem_module
  ],
  c_args: numpy_nodepr_api,
  dependencies: fortranobject_dep,
  install : true,
  link_language: 'fortran',
  subdir: 'taurex_ggchem',
)

py3.install_sources(
    ['__init__.py',],
    pure: false,
    subdir: 'taurex_ggchem'
)

py3.install_sources(
    ['chemistry/__init__.py','chemistry/ggchemistry_safe.py', 'chemistry/ggchemistry.py'],
    pure: false,
    subdir: 'taurex_ggchem/chemistry'
)


py3.install_sources(
    ['../data/abundances.pro',
 '../data/EarthCrust.dat',
 '../data/DustChem_GGchem.dat',
 '../data/dispol_StockKitzmann.dat',
 '../data/DustChem_Lee.dat',
 '../data/dispol_StockKitzmann_withoutTsuji.dat',
 '../data/dispol_GGchem.dat',
 '../data/dispol_WoitkeRefit.dat',
 '../data/dispol_fast.dat',
 '../data/DustChem.dat',
 '../data/DustChem_SH90.dat',
 '../data/Abundances.dat',
 '../data/DustChem_SUPCRTBL.dat',
 '../data/Burcat_ref-elements.dat',
 '../data/dispol_Tsuji.dat',
 '../data/dispol_BarklemCollet.dat',
 '../data/dispol_SharpHuebner.dat'],
    pure: false,
    subdir: 'taurex_ggchem/data'
)

py3.install_sources(
    ['../data/JANAF/K2SiO3_cr.txt',
 '../data/JANAF/SiCl2.txt',
 '../data/JANAF/FeCl2.txt',
 '../data/JANAF/W2O6.txt',
 '../data/JANAF/Mg2Cl4.txt',
 '../data/JANAF/H2WO4.txt',
 '../data/JANAF/FeS_cr.txt',
 '../data/JANAF/FeS.txt',
 '../data/JANAF/MgCl.txt',
 '../data/JANAF/Caatom.txt',
 '../data/JANAF/SiH.txt',
 '../data/JANAF/S7.txt',
 '../data/JANAF/AlCl3_cr.txt',
 '../data/JANAF/He+.txt',
 '../data/JANAF/KF.txt',
 '../data/JANAF/H+.txt',
 '../data/JANAF/Na2O2H2.txt',
 '../data/JANAF/Na2SiO3_cr.txt',
 '../data/JANAF/FeCl3.txt',
 '../data/JANAF/SiCl3.txt',
 '../data/JANAF/CaO_l.txt',
 '../data/JANAF/LiO.txt',
 '../data/JANAF/MgO_cr.txt',
 '../data/JANAF/CaCl.txt',
 '../data/JANAF/Zr_l.txt',
 '../data/JANAF/Li2SiO3_l.txt',
 '../data/JANAF/KCl_cr.txt',
 '../data/JANAF/W_l.txt',
 '../data/JANAF/PN.txt',
 '../data/JANAF/Cl-.txt',
 '../data/JANAF/fit_set.dat',
 '../data/JANAF/SiCH3Cl3.txt',
 '../data/JANAF/Al2Cl6.txt',
 '../data/JANAF/Ti4O7_l.txt',
 '../data/JANAF/Alref.txt',
 '../data/JANAF/Ca+.txt',
 '../data/JANAF/Ni3S2_l.txt',
 '../data/JANAF/H2+.txt',
 '../data/JANAF/AlO.txt',
 '../data/JANAF/CaS_cr.txt',
 '../data/JANAF/CKN_l.txt',
 '../data/JANAF/TiO_l.txt',
 '../data/JANAF/Ni_l.txt',
 '../data/JANAF/Featom.txt',
 '../data/JANAF/LiN.txt',
 '../data/JANAF/KOH.txt',
 '../data/JANAF/H3PO4_l.txt',
 '../data/JANAF/Ti_cr.txt',
 '../data/JANAF/Li+.txt',
 '../data/JANAF/V2O5_cr.txt',
 '../data/JANAF/VO_cr.txt',
 '../data/JANAF/SiCl4.txt',
 '../data/JANAF/ZrSiO4_cr.txt',
 '../data/JANAF/P4S3_cr.txt',
 '../data/JANAF/SiH4.txt',
 '../data/JANAF/Mg2SiO4_l.txt',
 '../data/JANAF/Li2O2H2.txt',
 '../data/JANAF/Co_cr.txt',
 '../data/JANAF/HAlO.txt',
 '../data/JANAF/testdata.pro',
 '../data/JANAF/Ti_l.txt',
 '../data/JANAF/SiN.txt',
 '../data/JANAF/SiO.txt',
 '../data/JANAF/CO2.txt',
 '../data/JANAF/COS.txt',
 '../data/JANAF/Hatom.txt',
 '../data/JANAF/CaO.txt',
 '../data/JANAF/S2O.txt',
 '../data/JANAF/H-.txt',
 '../data/JANAF/Al+.txt',
 '../data/JANAF/H2O.txt',
 '../data/JANAF/Crref.txt',
 '../data/JANAF/WCl2.txt',
 '../data/JANAF/MgOH.txt',
 '../data/JANAF/MgF2_cr.txt',
 '../data/JANAF/AllCond.txt',
 '../data/JANAF/Li2TiO3_cr.txt',
 '../data/JANAF/S_l.txt',
 '../data/JANAF/C4N2.txt',
 '../data/JANAF/K2Cl2.txt',
 '../data/JANAF/Liref.txt',
 '../data/JANAF/KCl.txt',
 '../data/JANAF/Cuatom.txt',
 '../data/JANAF/Sref.txt',
 '../data/JANAF/el-.txt',
 '../data/JANAF/Mgref.txt',
 '../data/JANAF/Fe.txt',
 '../data/JANAF/AlH.txt',
 '../data/JANAF/VO_l.txt',
 '../data/JANAF/Zr_cr.txt',
 '../data/JANAF/Cl+.txt',
 '../data/JANAF/S2.txt',
 '../data/JANAF/PH.txt',
 '../data/JANAF/Mn_cr.txt',
 '../data/JANAF/H2SO4.txt',
 '../data/JANAF/P3N5_cr.txt',
 '../data/JANAF/H2-.txt',
 '../data/JANAF/F2.txt',
 '../data/JANAF/Kref.txt',
 '../data/JANAF/Na2Cl2.txt',
 '../data/JANAF/Si_cr,l.txt',
 '../data/JANAF/Fe(CO)5_l.txt',
 '../data/JANAF/Naref.txt',
 '../data/JANAF/Li3Cl3.txt',
 '../data/JANAF/LiH.txt',
 '../data/JANAF/TiOCl2.txt',
 '../data/JANAF/WCl.txt',
 '../data/JANAF/NaF.txt',
 '../data/JANAF/Li2SiO3_cr.txt',
 '../data/JANAF/Co.txt',
 '../data/JANAF/MgCl2.txt',
 '../data/JANAF/NH.txt',
 '../data/JANAF/CH2.txt',
 '../data/JANAF/MgS.txt',
 '../data/JANAF/NaCl_l.txt',
 '../data/JANAF/H2SO4_l.txt',
 '../data/JANAF/H3O+.txt',
 '../data/JANAF/CH3.txt',
 '../data/JANAF/VO.txt',
 '../data/JANAF/Mg2SiO4_cr.txt',
 '../data/JANAF/Ni.txt',
 '../data/JANAF/K.txt',
 '../data/JANAF/CrN_cr.txt',
 '../data/JANAF/CN.txt',
 '../data/JANAF/KCl_l.txt',
 '../data/JANAF/LiCl_cr.txt',
 '../data/JANAF/TiOCl.txt',
 '../data/JANAF/Fe(CO)5.txt',
 '../data/JANAF/Si(CH3)4.txt',
 '../data/JANAF/AlO2.txt',
 '../data/JANAF/LiOH.txt',
 '../data/JANAF/P4O10.txt',
 '../data/JANAF/Cl.txt',
 '../data/JANAF/fitplot.py',
 '../data/JANAF/Niref.txt',
 '../data/JANAF/Niatom.txt',
 '../data/JANAF/Mg.txt',
 '../data/JANAF/Al2O3_l.txt',
 '../data/JANAF/Janaffit.f',
 '../data/JANAF/Zn.txt',
 '../data/JANAF/Cref.txt',
 '../data/JANAF/Zn_cr.txt',
 '../data/JANAF/TiO.txt',
 '../data/JANAF/LiH_cr.txt',
 '../data/JANAF/Catom.txt',
 '../data/JANAF/MgClF.txt',
 '../data/JANAF/FeCl2_l.txt',
 '../data/JANAF/NaCN_l.txt',
 '../data/JANAF/N+.txt',
 '../data/JANAF/WO3_l.txt',
 '../data/JANAF/H.txt',
 '../data/JANAF/Na2.txt',
 '../data/JANAF/Li2Si2O5_cr.txt',
 '../data/JANAF/MgTi2O5_cr.txt',
 '../data/JANAF/NaCl_cr.txt',
 '../data/JANAF/AllJANAF.txt',
 '../data/JANAF/AlCl3.txt',
 '../data/JANAF/PH2corr.txt',
 '../data/JANAF/FeO_cr.txt',
 '../data/JANAF/Na2SO4.txt',
 '../data/JANAF/HCN.txt',
 '../data/JANAF/Zn_l.txt',
 '../data/JANAF/P_l.txt',
 '../data/JANAF/CH4.txt',
 '../data/JANAF/CKN.txt',
 '../data/JANAF/Alatom.txt',
 '../data/JANAF/P_cr.txt',
 '../data/JANAF/MgF2.txt',
 '../data/JANAF/HCO.txt',
 '../data/JANAF/Mgatom.txt',
 '../data/JANAF/SiO2.txt',
 '../data/JANAF/TiCl4_l.txt',
 '../data/JANAF/NO+.txt',
 '../data/JANAF/AlCl2.txt',
 '../data/JANAF/CaCl2_l.txt',
 '../data/JANAF/K2SO4.txt',
 '../data/JANAF/NH4Cl_cr.txt',
 '../data/JANAF/MgS_cr.txt',
 '../data/JANAF/TiCl.txt',
 '../data/JANAF/Naatom.txt',
 '../data/JANAF/OH-.txt',
 '../data/JANAF/O.txt',
 '../data/JANAF/Feref.txt',
 '../data/JANAF/Ca_cr,l.txt',
 '../data/JANAF/Siref.txt',
 '../data/JANAF/KOH_l.txt',
 '../data/JANAF/Li2TiO3_l.txt',
 '../data/JANAF/WO.txt',
 '../data/JANAF/Mnatom.txt',
 '../data/JANAF/C2H4.txt',
 '../data/JANAF/Al2O2.txt',
 '../data/JANAF/NaOH_l.txt',
 '../data/JANAF/CrN.txt',
 '../data/JANAF/HCl.txt',
 '../data/JANAF/MgTi2O5_l.txt',
 '../data/JANAF/V2O3_cr.txt',
 '../data/JANAF/N.txt',
 '../data/JANAF/CaO_cr.txt',
 '../data/JANAF/CaF2.txt',
 '../data/JANAF/SiC2.txt',
 '../data/JANAF/SiH3Cl.txt',
 '../data/JANAF/K2O2H2.txt',
 '../data/JANAF/Katom.txt',
 '../data/JANAF/LiCl.txt',
 '../data/JANAF/Fe_cr.txt',
 '../data/JANAF/MgSiO3_cr.txt',
 '../data/JANAF/CaO2H2.txt',
 '../data/JANAF/PH2.txt',
 '../data/JANAF/K2SiO3_l.txt',
 '../data/JANAF/Si_cr.txt',
 '../data/JANAF/Cr_l.txt',
 '../data/JANAF/C.txt',
 '../data/JANAF/Al2O3_cr.txt',
 '../data/JANAF/FeO_l.txt',
 '../data/JANAF/Li.txt',
 '../data/JANAF/K2SO4_l.txt',
 '../data/JANAF/Na.txt',
 '../data/JANAF/MgTiO3_l.txt',
 '../data/JANAF/Natom.txt',
 '../data/JANAF/Al_cr,l.txt',
 '../data/JANAF/Zr.txt',
 '../data/JANAF/N2O4.txt',
 '../data/JANAF/Fe_l.txt',
 '../data/JANAF/O2.txt',
 '../data/JANAF/PH3.txt',
 '../data/JANAF/O2ref.txt',
 '../data/JANAF/NaCl.txt',
 '../data/JANAF/AlCl.txt',
 '../data/JANAF/fit.gnu',
 '../data/JANAF/NaCN.txt',
 '../data/JANAF/Na2S_cr.txt',
 '../data/JANAF/ZrO2_cr.txt',
 '../data/JANAF/TiCl3.txt',
 '../data/JANAF/WO3.txt',
 '../data/JANAF/Cr.txt',
 '../data/JANAF/W.txt',
 '../data/JANAF/ZrO2_l.txt',
 '../data/JANAF/Mn.txt',
 '../data/JANAF/Al2O2+.txt',
 '../data/JANAF/CaF2_cr.txt',
 '../data/JANAF/CaCl2.txt',
 '../data/JANAF/H3PO4_cr.txt',
 '../data/JANAF/Al2O+.txt',
 '../data/JANAF/MgO.txt',
 '../data/JANAF/FeO2H2.txt',
 '../data/JANAF/V.txt',
 '../data/JANAF/CS.txt',
 '../data/JANAF/Al.txt',
 '../data/JANAF/Pref.txt',
 '../data/JANAF/WO2.txt',
 '../data/JANAF/TiCl2.txt',
 '../data/JANAF/AlO-.txt',
 '../data/JANAF/ZrO2.txt',
 '../data/JANAF/NaO-.txt',
 '../data/JANAF/Ni_cr.txt',
 '../data/JANAF/Fatom.txt',
 '../data/JANAF/NaOH.txt',
 '../data/JANAF/AlOH.txt',
 '../data/JANAF/FeF2_cr.txt',
 '../data/JANAF/C_cr.txt',
 '../data/JANAF/FeS2_cr.txt',
 '../data/JANAF/TiO2_cr.txt',
 '../data/JANAF/Satom.txt',
 '../data/JANAF/fit_comp.py',
 '../data/JANAF/Co_l.txt',
 '../data/JANAF/SiC_cr.txt',
 '../data/JANAF/SO2.txt',
 '../data/JANAF/Ti.txt',
 '../data/JANAF/SO3.txt',
 '../data/JANAF/Mg+.txt',
 '../data/JANAF/W_cr.txt',
 '../data/JANAF/S.txt',
 '../data/JANAF/Ca.txt',
 '../data/JANAF/Li2Cl2.txt',
 '../data/JANAF/TiO2.txt',
 '../data/JANAF/Cratom.txt',
 '../data/JANAF/N2H4_l.txt',
 '../data/JANAF/Ti4O7_cr.txt',
 '../data/JANAF/LiH_l.txt',
 '../data/JANAF/FeF2.txt',
 '../data/JANAF/SiCl.txt',
 '../data/JANAF/CO2-.txt',
 '../data/JANAF/Na+.txt',
 '../data/JANAF/S2Ocorr.txt',
 '../data/JANAF/C2H4O.txt',
 '../data/JANAF/F.txt',
 '../data/JANAF/C5.txt',
 '../data/JANAF/Oatom.txt',
 '../data/JANAF/W3O9.txt',
 '../data/JANAF/N2.txt',
 '../data/JANAF/Ni3S2_cr.txt',
 '../data/JANAF/WF.txt',
 '../data/JANAF/V2O4_cr.txt',
 '../data/JANAF/W3O8.txt',
 '../data/JANAF/Cu.txt',
 '../data/JANAF/P.txt',
 '../data/JANAF/AlOCl.txt',
 '../data/JANAF/TiCl4.txt',
 '../data/JANAF/AlO+.txt',
 '../data/JANAF/Tiatom.txt',
 '../data/JANAF/SiO2_l.txt',
 '../data/JANAF/Mn_l.txt',
 '../data/JANAF/TiO_cr.txt',
 '../data/JANAF/test.txt',
 '../data/JANAF/Zrref.txt',
 '../data/JANAF/AlOH-.txt',
 '../data/JANAF/Mnref.txt',
 '../data/JANAF/MgO_l.txt',
 '../data/JANAF/AlO2H.txt',
 '../data/JANAF/AlCl3_l.txt',
 '../data/JANAF/S_cr.txt',
 '../data/JANAF/TiC_cr.txt',
 '../data/JANAF/KHF2_cr.txt',
 '../data/JANAF/AlS.txt',
 '../data/JANAF/Na2SO4_l.txt',
 '../data/JANAF/Si.txt',
 '../data/JANAF/Na2S_l.txt',
 '../data/JANAF/MgOH+.txt',
 '../data/JANAF/SiO2_cr.txt',
 '../data/JANAF/HSO3F.txt',
 '../data/JANAF/CaCl2_cr.txt',
 '../data/JANAF/VO2.txt',
 '../data/JANAF/KO.txt',
 '../data/JANAF/Na_cr.txt',
 '../data/JANAF/LiOH_cr.txt',
 '../data/JANAF/LiCl_l.txt',
 '../data/JANAF/MgCl2_l.txt',
 '../data/JANAF/dGatom.pro',
 '../data/JANAF/ZnSO4_cr.txt',
 '../data/JANAF/HF2K_cr.txt',
 '../data/JANAF/Li2Si2O5_l.txt',
 '../data/JANAF/W4O12.txt',
 '../data/JANAF/CN-.txt',
 '../data/JANAF/AlF6Na3_cr.txt',
 '../data/JANAF/Cl2.txt',
 '../data/JANAF/Clatom.txt',
 '../data/JANAF/Na_l.txt',
 '../data/JANAF/SiS2_cr.txt',
 '../data/JANAF/FeO.txt',
 '../data/JANAF/Siatom.txt',
 '../data/JANAF/TiO2_l.txt',
 '../data/JANAF/FeCl.txt',
 '../data/JANAF/SiS2_l.txt',
 '../data/JANAF/Cr_cr.txt',
 '../data/JANAF/Al2O.txt',
 '../data/JANAF/NaOH_cr.txt',
 '../data/JANAF/Caref.txt',
 '../data/JANAF/MgAl2O4_cr.txt',
 '../data/JANAF/H2O_l.txt',
 '../data/JANAF/Cl2ref.txt',
 '../data/JANAF/AlOH+.txt',
 '../data/JANAF/NH3.txt',
 '../data/JANAF/Liatom.txt',
 '../data/JANAF/P4O6.txt',
 '../data/JANAF/N2H4.txt',
 '../data/JANAF/H2.txt',
 '../data/JANAF/HS.txt',
 '../data/JANAF/S8.txt',
 '../data/JANAF/TiC_l.txt',
 '../data/JANAF/SN.txt',
 '../data/JANAF/MgSiO3_l.txt',
 '../data/JANAF/P4S3_l.txt',
 '../data/JANAF/Mg3P2O8_cr.txt',
 '../data/JANAF/Curef.txt',
 '../data/JANAF/He.txt',
 '../data/JANAF/Tiref.txt',
 '../data/JANAF/MgO2H2.txt',
 '../data/JANAF/Al_cr.txt',
 '../data/JANAF/NH2.txt',
 '../data/JANAF/LiOH_l.txt',
 '../data/JANAF/MgTiO3_cr.txt',
 '../data/JANAF/WO2Cl2.txt',
 '../data/JANAF/Na2SiO3_l.txt',
 '../data/JANAF/HScorr.txt',
 '../data/JANAF/KF_cr.txt',
 '../data/JANAF/Patom.txt',
 '../data/JANAF/AlF3_cr.txt',
 '../data/JANAF/CoO_cr.txt',
 '../data/JANAF/MgAl2O4_l.txt',
 '../data/JANAF/H2S.txt',
 '../data/JANAF/CaS.txt',
 '../data/JANAF/Zratom.txt',
 '../data/JANAF/WO3_cr.txt',
 '../data/JANAF/FeS_l.txt',
 '../data/JANAF/SiS.txt',
 '../data/JANAF/HCO+.txt',
 '../data/JANAF/K+.txt',
 '../data/JANAF/AlO2-.txt',
 '../data/JANAF/NaF_cr.txt',
 '../data/JANAF/P4O10_cr.txt'],
    pure: false,
    subdir: 'taurex_ggchem/data/JANAF'
)

py3.install_sources(
    ['../data/SUPCRT92/slop16.dat',
 '../data/SUPCRT92/SharpHuebner.dat',
 '../data/SUPCRT92/extractSU.py',
 '../data/SUPCRT92/spronsbl.dat',
 '../data/SUPCRT92/StandardStates.pdf',
 '../data/SUPCRT92/NISTatom.dat',
 '../data/SUPCRT92/dG.pdf',
 '../data/SUPCRT92/DustChemSUPCRTBL.tex',
 '../data/SUPCRT92/DustChemSUPCRTBL.dat',
 '../data/SUPCRT92/make_SH90.py'],
    pure: false,
    subdir: 'taurex_ggchem/data/SUPCRT92'
)