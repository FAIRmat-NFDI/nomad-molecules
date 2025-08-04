import json
import sqlite3
import logging
import pytest
import ase.build
from ase import Atoms
from unittest.mock import MagicMock
from structlog.testing import capture_logs

from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.metainfo import runschema
from nomad.client import normalize_all
from nomad.config import config
from nomad_molecules.normalizers import MoleculesNormalizerEntryPoint

# --------------------- ASE Atoms Fixtures ---------------------
@pytest.fixture
def H2O_CO2_molecule_group():
    """A combined system of two H₂O molecules plus one CO₂ molecule."""
    atom_co2 = Atoms(
        symbols=["C", "O", "O"],
        positions=[[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
        cell=[[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]],
        pbc=False
    )
    water1 = ase.build.molecule('H2O')
    water2 = ase.build.molecule('H2O')
    water2.translate([5, 0, 0])
    system = water1 + water2 + atom_co2
    system.set_cell([10, 10, 10])
    system.set_pbc(False)
    return system

@pytest.fixture
def one_d_chain_atoms():
    """1D carbon chain, for dimensionality tests."""
    return Atoms(
        symbols=["C","C","C","C"],
        positions=[[0,0,0],[1.5,0,0],[3.0,0,0],[4.5,0,0]],
        cell=[[6.0,0,0],[0,10.0,0],[0,0,10.0]],
        pbc=[True, False, False]
    )

@pytest.fixture
def simple_heavy_water_atoms():
    """Heavy water (D₂O) without PBC for skeleton‐match tests."""
    atoms = Atoms(
        symbols=["O","H","H"],
        positions=[[2.5,2.5,2.5],[3.257,3.086,2.5],[1.743,3.086,2.5]],
        cell=[[5.0,0,0],[0,5.0,0],[0,0,5.0]],
        pbc=False
    )
    masses = atoms.get_masses()
    masses[1] = masses[2] = 2.01410177811
    atoms.set_masses(masses)
    return atoms

# ----------------- PubChem Temporary DB Fixture -----------------
@pytest.fixture(autouse=True)
def temporary_db(tmp_path, monkeypatch):
    """
    Creates an offline‐basic PubChem DB with H₂O in it,
    and sets MOLID_MASTER_DB / MOLID_CACHE_DB / MOLID_MODE accordingly.
    """
    db = tmp_path / "pubchem_data_test.db"
    if not db.exists():
        conn = sqlite3.connect(db)
        c = conn.cursor()
        c.execute("""
            CREATE TABLE compound_data (
                id INTEGER PRIMARY KEY,
                InChIKey TEXT UNIQUE,
                InChI TEXT,
                SMILES TEXT,
                InChIKey14 TEXT,
                Formula TEXT
            )
        """)
        c.execute(
            "INSERT INTO compound_data (InChIKey, InChI, SMILES, InChIKey14, Formula) VALUES (?,?,?,?,?)",
            ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "InChI=1S/H2O/h1H2", "O", "XLYOFNOQVPJJNP", "H2O")
        )
        conn.commit()
        conn.close()

    DummyCfg = MoleculesNormalizerEntryPoint(
        molid_mode       = "offline-basic",
        molid_master_db  = str(db),
        molid_cache_db   = str(db),
        max_atoms        = 4,
        min_atoms        = 2)
    # class DummyCfg:
    #     molid_mode       = "offline-basic"
    #     molid_master_db  = str(db)
    #     molid_cache_db   = str(db)
    #     max_atoms        = 4
    #     min_atoms        = 2

    monkeypatch.setattr(type(config),
                        "get_plugin_entry_point",
                        lambda self, entry_point_id: DummyCfg)

@pytest.fixture
def logger():
    """Fresh MagicMock logger."""
    return MagicMock()


# ---------------- Helper to build a minimal NOMAD archive ----------------
def get_section_system(atoms: Atoms):
    """Wrap an ASE Atoms into a runsystem.System + runsystem.Atoms section."""
    system = runschema.system.System()
    system.atoms = runschema.system.Atoms(
        positions=atoms.get_positions() * 1e-10,
        labels=atoms.get_chemical_symbols(),
        lattice_vectors=atoms.get_cell() * 1e-10,
        periodic=atoms.get_pbc(),
    )
    return system

def create_archive(ase_atoms: Atoms):
    """Makes an EntryArchive containing exactly one system with our ASE atoms."""
    archive = EntryArchive()
    run = runschema.run.Run()
    archive.run.append(run)
    run.program = runschema.run.Program(name='VASP', version='4.6.35')
    system = get_section_system(ase_atoms)
    run.system.append(system)
    calc = runschema.calculation.Calculation()
    calc.system_ref = system
    run.calculation.append(calc)
    archive.metadata = EntryMetadata()
    archive.metadata.domain = 'dft'
    return archive

def add_H2O_CO2_groups(system):
    # H2O molecules group
    H2O_group = runschema.system.AtomsGroup(
        label='H2O_GROUP',
        type='molecule_group',
        index=0,
        composition_formula='H4O2',
        n_atoms=6,
        atom_indices=[0,1,2, 3,4,5],
    )
    system.atoms_group.append(H2O_group)

    #   -- individual H2O molecules under H2O_GROUP
    H2O1 = runschema.system.AtomsGroup(
        label='H2O_MOL',
        type='molecule',
        index=0,
        composition_formula='H2O1',
        n_atoms=3,
        atom_indices=[0,1,2],
    )
    H2O_group.atoms_group.append(H2O1)

    H2O2 = runschema.system.AtomsGroup(
        label='H2O_MOL',
        type='molecule',
        index=1,
        composition_formula='H2O1',
        n_atoms=3,
        atom_indices=[3,4,5],
    )
    H2O_group.atoms_group.append(H2O2)

    # CO2 molecules group
    co2 = runschema.system.AtomsGroup(
        label='CO2_MOL',
        type='molecule',
        index=0,
        composition_formula='C1O2',
        n_atoms=3,
        atom_indices=[6,7,8],
    )
    system.atoms_group.append(co2)


# ---------------------- Integration Tests ----------------------
def test_molecule_group(H2O_CO2_molecule_group, temporary_db, caplog):
    archive = create_archive(H2O_CO2_molecule_group)
    system = archive.run[0].system[0]
    add_H2O_CO2_groups(system)

    caplog.set_level(logging.WARNING)
    normalize_all(archive)

    topologies = archive.results.material.topology
    # Should produce 5 entries: original, H2O_GROUP, H2O_MOL×2, CO2_MOL
    assert len(topologies) == 4

    # --- original ---
    original = topologies[0]
    assert original.label == "original"
    assert original.method == "parser"
    assert getattr(original, "cheminformatics", None) is None

    # --- group entries have no cheminformatics ---
    assert topologies[1].label == "H2O_GROUP" and topologies[1].cheminformatics is None

    # --- H2O molecule full match ---
    h2o = topologies[2]
    assert h2o.label == "H2O_MOL"
    cf_h2o = h2o.cheminformatics
    assert cf_h2o.inchi_key == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    assert cf_h2o.inchi     == 'InChI=1S/H2O/h1H2'
    assert cf_h2o.smiles    == 'O'
    assert cf_h2o.match_type == 'full structure'

    # --- CO2 molecule skeleton (no DB entry) ---
    co2 = topologies[3]
    assert co2.label == "CO2_MOL"
    assert getattr(co2, "cheminformatics", None) is None


def test_one_d_chain_atoms(one_d_chain_atoms, caplog):
    archive = create_archive(one_d_chain_atoms)
    caplog.set_level(logging.WARNING)
    with capture_logs() as captured:
        logging.getLogger().setLevel(logging.WARNING)
        normalize_all(archive)

    assert any(e["event"] == "System is 1D, only 0D systems are processed. Skipping normalization and continue." for e in captured)

    topologies = archive.results.material.topology
    # origin + subsystem + conventional cell = 3 entries
    assert len(topologies) == 3
    assert topologies[1].dimensionality == "1D"
    assert topologies[2].dimensionality == "1D"
    for topo in topologies:
        assert getattr(topo, "cheminformatics", None) is None

def test_missing_DB(H2O_CO2_molecule_group, tmp_path, monkeypatch, caplog):
    # Point to a non‐existent DB
    missing_db = tmp_path / "nofile.db"

    DummyCfg = MoleculesNormalizerEntryPoint(
        molid_mode       = "offline-basic",
        molid_master_db  = str(missing_db),
        molid_cache_db   = str(missing_db),
        max_atoms        = 4,
        min_atoms        = 2)

    monkeypatch.setattr(type(config),
                        "get_plugin_entry_point",
                        lambda self, entry_point_id: DummyCfg)

    archive = create_archive(H2O_CO2_molecule_group)
    system = archive.run[0].system[0]
    add_H2O_CO2_groups(system)

    caplog.set_level(logging.ERROR)
    normalize_all(archive)

    errs = [r.getMessage() for r in caplog.records if r.levelno == logging.ERROR]
    assert any(
        json.loads(m)["event"] == f"Database file '{missing_db}' not found or inaccessible."
        for m in errs
    )

    # And no cheminformatics anywhere
    for topo in archive.results.material.topology:
        assert getattr(topo, "cheminformatics", None) is None


test_validate_atom_count_cases = [
    {
        "species": [8],
        "positions": [[2.5, 2.5, 2.5]],
        "periodic": [False, False, False],
        "min_atoms": 2,
        "max_atoms": 5,
        "expected_result": False,
        "log_message": "System has only 1 atoms; minimum required is 2.",
        "case_description": "Below minimum atom count"
    },
    {
        "species": [8, 1, 1, 6, 7],
        "positions": [[2.5, 2.5, 2.5], [0.5, 2.5, 2.5], [2.5, 3.5, 2.5], [1.5, 1.5, 1.5], [3.0, 3.0, 3.0]],
        "periodic": [False, False, False],
        "min_atoms": 2,
        "max_atoms": 4,
        "expected_result": False,
        "log_message": "System has 5 atoms; exceeds maximum allowed 4.",
        "case_description": "Exceeds maximum atom count"
    }
]

@pytest.mark.parametrize(
    "case",
    test_validate_atom_count_cases,
    ids=[c["case_description"] for c in test_validate_atom_count_cases]
)
def test_validate_atom_count(case, temporary_db, caplog):
    """Test validate_atom_count function covering edge and valid cases."""
    atoms_data = Atoms(
        symbols=case["species"],
        positions=case["positions"],
        pbc=case["periodic"]
    )

    archive = create_archive(atoms_data)
    caplog.set_level(logging.WARNING)
    normalize_all(archive)
    topologies = archive.results.material.topology
    assert len(topologies) == 1
    errs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any(case["log_message"] in m for m in errs)
    assert getattr(topologies[0], "cheminformatics", None) is None
    assert getattr(topologies[0], "building_block", None) is None

# This is not working because Nomad has no option to identify isotops yet
# def test_D2O_skeleton_match(simple_heavy_water_atoms):
#     archive = create_archive(simple_heavy_water_atoms)
#     normalize_all(archive)

#     topo = archive.results.material.topology[0]
#     cf = topo.cheminformatics
#     # Since database only has H₂O, this will do a 14‐char skeleton match
#     assert cf.inchi_key.startswith("XLYOFNOQVPJJNP-")
#     assert cf.match_type == 'skeleton (core)'
#     # No full InChI or SMILES on skeleton match
#     assert not hasattr(cf, "inchi")
#     assert not hasattr(cf, "smiles")