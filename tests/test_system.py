import pytest

from pmd.core import EMC, PSP, SolventSystem, System
from pmd.core.Builder import Builder


@pytest.fixture
def test_data():
    test_builder = PSP('opls-lbcc')

    syst = System('*CC*',
                  builder=test_builder,
                  density=0.5,
                  natoms_total=500,
                  natoms_per_chain=100)

    solventsyst = SolventSystem(smiles='*CC*',
                                solvent_smiles='C1CCCCC1',
                                ru_nsolvent_ratio=0.1,
                                builder=test_builder,
                                density=0.8,
                                natoms_total=5000,
                                natoms_per_chain=150)
    return {
        'builder': test_builder,
        'system': syst,
        'solventsystem': solventsyst,
    }


def test_system_initialization(test_data):
    syst = test_data['system']
    assert syst.smiles == '*CC*'
    assert syst.end_cap_smiles == '*C'
    assert syst.builder == test_data['builder']
    assert syst.data_fname == 'data.lmps'


def test_system_update(test_data):
    syst = test_data['system']
    new_builder = EMC('pcff')
    syst.smiles = '*CC(*)CC'
    syst.end_cap_smiles = '*[H]'
    syst.builder = new_builder
    assert syst.smiles == '*CC(*)CC'
    assert syst.end_cap_smiles == '*[H]'
    assert syst.builder == new_builder


def test_solventsystem_initialization(test_data):
    solv_syst = test_data['solventsystem']
    assert solv_syst.smiles == '*CC*'
    assert solv_syst.solvent_group == 'molecule <= 62'
    assert solv_syst.builder == test_data['builder']


def test_solventsystem_update(test_data):
    solv_syst = test_data['solventsystem']
    solv_syst.smiles = '*CC(*)CC'
    assert solv_syst.smiles == '*CC(*)CC'


def test_system_exceptions(test_data):
    # PSP not installed
    with pytest.raises(ImportError):
        syst = test_data['system']
        syst.write_data()

    # Invalid builder provided
    with pytest.raises(ValueError):
        syst.builder = 'pcff'

    # Invalid force_field provided for the builder
    with pytest.raises(ValueError):
        syst.builder = EMC('xyz')

    # No system size option provided
    with pytest.raises(ValueError):
        system = System('*CC*',
                        builder=EMC('pcff'),
                        density=0.5,
                        natoms_per_chain=100)
        system.write_data()

    # 2 chain length options simultaneously provided
    with pytest.raises(ValueError):
        system = System('*CC*',
                        builder=EMC('pcff'),
                        density=0.5,
                        natoms_per_chain=100,
                        mw_per_chain=1000)
        system.write_data()

    # Using Builder class itself is not allowed
    bare_builder = Builder('pcff', ('pcff'))
    with pytest.raises(NotImplementedError):
        bare_builder.write_data('.', '*CC*', 0.5, 1000, 10, 10, 'data.lmps',
                                True)

    with pytest.raises(NotImplementedError):
        bare_builder.write_functional_form(False)

    with pytest.raises(NotImplementedError):
        bare_builder.write_solvent_data('.', '*CC*', 'CCO', 0.5, 1000, 10, 10,
                                        10, 'data.lmps', True)


@pytest.mark.parametrize(
    'system',
    [
        System('*CC*',
               builder=EMC('pcff'),
               density=0.5,
               natoms_total=500,
               natoms_per_chain=100),
        System('*CC*',
               builder=EMC('pcff'),
               density=0.5,
               nchains_total=5,
               ru_per_chain=5),
        System('*CC*',
               builder=EMC('pcff'),
               density=0.5,
               natoms_total=500,
               mw_per_chain=500),
        System('*CC*',
               builder=EMC('opls-aa'),
               density=0.5,
               natoms_total=500,
               natoms_per_chain=100),
        System('*CC*',
               builder=EMC('opls-aa'),
               density=0.5,
               nchains_total=5,
               ru_per_chain=5),
        System('*CC*',
               builder=EMC('opls-aa'),
               density=0.5,
               natoms_total=500,
               mw_per_chain=500),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='C1CCCCC1',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('pcff'),
                      density=0.5,
                      natoms_total=500,
                      natoms_per_chain=100),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='C1CCCCC1',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('pcff'),
                      density=0.5,
                      nchains_total=5,
                      ru_per_chain=5),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='C1CCCCC1',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('pcff'),
                      density=0.5,
                      natoms_total=500,
                      mw_per_chain=500),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='CCCCC',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('opls-aa'),
                      density=0.5,
                      natoms_total=500,
                      natoms_per_chain=100),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='CCCCC',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('opls-aa'),
                      density=0.5,
                      nchains_total=5,
                      ru_per_chain=5),
        SolventSystem(smiles='*CC*',
                      solvent_smiles='CCCCC',
                      ru_nsolvent_ratio=0.1,
                      builder=EMC('opls-aa'),
                      density=0.5,
                      natoms_total=500,
                      mw_per_chain=500),
    ],
)
def test_system_emc_write_data(tmp_path, system):
    d = tmp_path / 'result'
    system.write_data(d)
    assert len(list(tmp_path.iterdir())) == 1
