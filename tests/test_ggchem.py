import taurex_ggchem


def test_ggchem():
    """Test GGChem class if it functions properly"""
    from taurex_ggchem import GGChem
    from taurex.constants import AMU

    import numpy as np
    import math

    gg = GGChem()
    tp = np.logspace(math.log10(2500), math.log10(2000), 100)

    pp = np.ones(100) * 1e5
    pp
    gg.initialize_chemistry(nlayers=100, temperature_profile=tp, pressure_profile=pp)
    assert np.all(gg.muProfile / AMU != 0.0)


def test_ggchem():
    """Test GGChem class if it functions properly"""
    from taurex_ggchem import GGChem
    from taurex.constants import AMU

    import numpy as np
    import math

    gg = GGChem(equilibrium_condensation=True)
    tp = np.ones(1) * 2500

    pp = np.ones(1) * 1e7
    pp
    gg.initialize_chemistry(nlayers=1, temperature_profile=tp, pressure_profile=pp)
    assert np.all(gg.muProfile / AMU != 0.0)


def test_ggchem_unsafe():
    """Test GGChem Non-Safe class if it functions properly"""
    from taurex_ggchem import GGChemNonSafe as GGChem
    from taurex.constants import AMU

    import numpy as np
    import math

    gg = GGChem(equilibrium_condensation=True)
    tp = np.ones(1) * 2500

    pp = np.ones(1) * 1e7
    pp
    gg.initialize_chemistry(nlayers=1, temperature_profile=tp, pressure_profile=pp)
    assert np.all(gg.muProfile / AMU != 0.0)
