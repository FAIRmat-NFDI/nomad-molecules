import numpy as np
import pytest
import sqlite3
import os
from ase import Atoms
from runschema.system import Atoms as NomadAtoms
from unittest.mock import MagicMock, call
from matid.geometry import get_dimensionality
from nomad.datamodel.results import System
from nomad_plugin_molecules.normalizers.atoms_utils import (
    get_atoms_data,
    wrap_atoms,
    validate_atom_count,
    validate_dimensionality,
    query_molecule_database_util,
    generate_topology_util
)

ATOMS_REF_NOMAD = NomadAtoms(atomic_numbers=[8, 1, 1],
                             labels=['O', 'H', 'H'],
                             positions=[[0.0, 0.0, -6.140480000000001e-13],
                                        [7.644331800000001e-11, 0.0, 5.8917024e-11],
                                        [-7.644331800000001e-11, 0.0, 5.8917024e-11]],
                             periodic=[False, False, False])

# ====================== Pytest Fixtures ======================
@pytest.fixture
def logger():
    """Returns a fresh MagicMock logger for each test."""
    return MagicMock()

@pytest.fixture
def simple_water_atoms():
    """Returns a simple water molecule with PBC enabled for testing."""
    return Atoms(
        symbols=["O", "H", "H"],
        positions=[[2.5, 2.5, 2.5], [3.257, 3.086, 2.5], [1.743, 3.086, 2.5]],
        cell=[[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]],
        pbc=True
    )

@pytest.fixture
def co2_atoms():
    """Returns a simple CO2 molecule without PBC for testing."""
    return Atoms(
        symbols=["C", "O", "O"],
        positions=[[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
        pbc=False
    )

@pytest.fixture
def simple_heavy_water_atoms():
    """Returns a simple heavy water molecule (D2O) with PBC enabled for testing."""
    # start from ordinary H2O
    atoms = Atoms(
        symbols=["O", "H", "H"],
        positions=[[2.5, 2.5, 2.5], [3.257, 3.086, 2.5], [1.743, 3.086, 2.5]],
        cell=[[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]],
        pbc=True
    )
    # relabel H → D for display
    # atoms.set_chemical_symbols(["O", "D", "D"])
    # override the masses of the two D atoms to the deuterium mass (≈ 2.014 u)
    masses = atoms.get_masses()
    masses[1] = masses[2] = 2.01410177811
    atoms.set_masses(masses)
    return atoms

@pytest.fixture
def zero_d_atoms():
    """Returns a 0D water-like system for dimensionality tests."""
    return Atoms(
        symbols="H2O",
        positions=[[2.5, 2.5, 2.5], [0.5, 2.5, 2.5], [2.5, 3.5, 2.5]],
        cell=[10.0, 10.0, 10.0],
        pbc=[False, False, False]
    )

@pytest.fixture
def one_d_chain_atoms():
    """Returns a 1D carbon chain for dimensionality tests."""
    return Atoms(
        symbols="CCCC",
        positions=[[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0], [4.5, 0, 0]],
        cell=[[6.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]],
        pbc=[True, False, False]
    )

@pytest.fixture(scope="function")
def test_database(tmp_path):
    """Creates a temporary SQLite database and yields its path for tests."""
    db_file = tmp_path / "test_master.db"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE compound_data (
            id INTEGER PRIMARY KEY,
            InChIKey TEXT UNIQUE,
            InChIKey14 TEXT,
            Name TEXT,
            Formula TEXT
        )
    """)
    cursor.execute("INSERT INTO compound_data (InChIKey, InChIKey14, Name, Formula) VALUES (?, ?, ?, ?)",
                   ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "XLYOFNOQVPJJNP", "Water", "H2O"))
    # cursor.execute("INSERT INTO compound_data (InChIKey14, Name, Formula) VALUES (?, ?, ?)",("XLYOFNOQVPJJNP", "H₂O_small", "H2O"))
    conn.commit()
    conn.close()
    yield str(db_file)

# Helper to create a topology
def make_topology(indices=None, atoms=None, atoms_ref=None):
    topo = System()
    topo.indices = indices
    topo.atoms = atoms
    topo.atoms_ref = atoms_ref
    return topo

# ---------------------- Tests for get_atoms_data ----------------------
# Happy path tests
@pytest.mark.parametrize(
    "indices, atoms, atoms_ref, expected_indices, expected_source, log_info, log_warning", [
        # With explicit indices
        ([[0, 2]], None, None, [0, 2], "indices", "Topology contains indices.", None),
        # Atoms attribute present (preferred over atoms_ref)
        (None, ATOMS_REF_NOMAD, None, None, "atoms", None, None),
        # Only atoms_ref attribute present
        (None, None, ATOMS_REF_NOMAD, None, "atoms_ref", None, None)
    ]
)
def test_get_atoms_data_happy_cases(indices, atoms, atoms_ref, expected_indices, expected_source, log_info, log_warning, logger):
    topo = make_topology(indices=indices, atoms=atoms, atoms_ref=atoms_ref)
    result = get_atoms_data(topo, ATOMS_REF_NOMAD, logger)
    if expected_source == "indices":
        expected = ATOMS_REF_NOMAD.to_ase()[expected_indices]
        assert isinstance(result, Atoms)
        assert result == expected
        np.testing.assert_array_equal(result.get_atomic_numbers(), [8, 1])
        logger.info.assert_called_once_with("Topology contains indices.")
        logger.warning.assert_not_called()
    else:
        assert isinstance(result, Atoms)
        assert result == ATOMS_REF_NOMAD.to_ase()
        logger.info.assert_not_called()
        logger.warning.assert_not_called()

# Edge case: no atoms data
def test_get_atoms_data_no_atoms_logs_warning_and_returns_none(logger):
    topo = make_topology(indices=None, atoms=None, atoms_ref=None)
    result = get_atoms_data(topo, ATOMS_REF_NOMAD, logger)
    assert result is None,  "No atoms data found (neither atoms nor atoms_ref exist)."

def test_get_atoms_data_prefers_atoms_over_atoms_ref_when_both_exist(logger):
    atoms_nomad = NomadAtoms(atomic_numbers=[1], labels=["H"], positions=[[0,0,0]], periodic=[False,False,False])
    atoms_ref_nomad = NomadAtoms(atomic_numbers=[1,1], labels=["H","H"], positions=[[0,0,0],[1,1,1]], periodic=[False,False,False])
    topo = make_topology(indices=None, atoms=atoms_nomad, atoms_ref=atoms_ref_nomad)
    result = get_atoms_data(topo, atoms_ref_nomad, logger)
    assert result == atoms_nomad.to_ase()
    logger.info.assert_not_called()  # no "Topology contains indices" on this branch
    logger.warning.assert_not_called()

# ---------------------- wrap_atoms Tests ----------------------
# Happy path: PBC wrapping and non-PBC behavior
@pytest.mark.parametrize(
    "symbols, positions, cell, expected_positions, case_description", [
        pytest.param(
            ["H", "H"],
            [[0, 0, 0], [1, 1, 1]],
            [[1, 0, 0], [0, 0, 0], [0, 0, 1]],
            [[0, 0, 0], [1, 1, 1]],
            "Singular cell - Atoms should remain unchanged",
            id="singular_cell"
        ),
        pytest.param(
            ["O", "H", "H"],
            [[0.2, 2.5, 2.5], [4.2, 2.5, 2.5], [0.2, 3.5, 2.5]],
            [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]],
            [[3.0, 2.0, 2.5], [2.0, 2.0, 2.5], [3.0, 3.0, 2.5]],
            "Periodic boundary condition - Atom should wrap correctly",
            id="pbc_wrapping"
        )
    ]
)
def test_wrap_atoms(symbols, positions, cell, expected_positions, case_description):
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    wrapped_atoms = wrap_atoms(atoms)
    assert np.allclose(wrapped_atoms.get_positions(), expected_positions), f"Failed: {case_description}"

# ---------------------- validate_atom_count Tests ----------------------
# Edge cases and valid range
test_validate_atom_count_cases = [
    {
        "species": [8],
        "positions": [[2.5, 2.5, 2.5]],
        "periodic": [False, False, False],
        "min_atoms": 2,
        "max_atoms": 5,
        "expected_result": False,
        "log_message": "System has only 1 atoms, which is less than the minimum required 2. Skipping normalization and continue.",
        "case_description": "Below minimum atom count"
    },
    {
        "species": [8, 1, 1],
        "positions": [[2.5, 2.5, 2.5], [0.5, 2.5, 2.5], [2.5, 3.5, 2.5]],
        "periodic": [False, False, False],
        "min_atoms": 2,
        "max_atoms": 5,
        "expected_result": True,
        "log_message": None,
        "case_description": "Within valid atom count range"
    },
    {
        "species": [8, 1, 1, 6, 7],
        "positions": [[2.5, 2.5, 2.5], [0.5, 2.5, 2.5], [2.5, 3.5, 2.5], [1.5, 1.5, 1.5], [3.0, 3.0, 3.0]],
        "periodic": [False, False, False],
        "min_atoms": 2,
        "max_atoms": 4,
        "expected_result": False,
        "log_message": "System has 5 atoms, which exceeds the maximum allowed 4. Skipping normalization and continue.",
        "case_description": "Exceeds maximum atom count"
    }
]

@pytest.mark.parametrize(
    "case",
    test_validate_atom_count_cases,
    ids=[c["case_description"] for c in test_validate_atom_count_cases]
)
def test_validate_atom_count(case, logger):
    """Test validate_atom_count function covering edge and valid cases."""
    atoms_data = Atoms(
        symbols=case["species"],
        positions=case["positions"],
        pbc=case["periodic"]
    )
    result = validate_atom_count(
        atoms_data,
        case["min_atoms"],
        case["max_atoms"],
        logger
    )
    assert result == case["expected_result"], f"Failed: {case['case_description']}"
    if case["log_message"]:
        logger.warning.assert_called_once_with(case["log_message"])
    else:
        logger.warning.assert_not_called()

# ---------------------- validate_dimensionality Tests ----------------------
# Valid (0D) and invalid (1D) systems
@pytest.mark.parametrize(
    "atoms, expected_dimensionality, expected_result, log_message, case_description", [
        ("zero_d_atoms", 0, True, None, "Valid 0D system"),
        ("one_d_chain_atoms", 1, False, "System is 1D, only 0D systems are processed. Skipping normalization and continue.", "1D system should be skipped")
    ]
)
def test_validate_dimensionality(atoms, expected_dimensionality, expected_result, log_message, case_description, request, logger):
    """Test validate_dimensionality function with 0D and 1D examples."""
    atoms_obj = request.getfixturevalue(atoms)
    assert get_dimensionality(atoms_obj) == expected_dimensionality, "Dimensionality calculation mismatch!"
    result = validate_dimensionality(atoms_obj, logger)
    assert result == expected_result, f"Failed: {case_description}"
    if log_message:
        logger.warning.assert_called_once_with(log_message)
    else:
        logger.warning.assert_not_called()

# ---------------------- query_molecule_database_util Tests ----------------------
@pytest.mark.parametrize(
    "atoms_fixture, expected_inchikey, expected_data, full_match, log_messages, case_description", [
        ("simple_water_atoms", "XLYOFNOQVPJJNP-UHFFFAOYSA-N", [{"id": 1, "InChIKey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N", 'InChIKey14': "XLYOFNOQVPJJNP", "Name": "Water", "Formula": "H2O"}], True, ["Computed InChIKey: XLYOFNOQVPJJNP-UHFFFAOYSA-N"], "molecule_found"),
        ("co2_atoms", "CURLTUGMZLYLDI-UHFFFAOYSA-N", [], None, ["Computed InChIKey: CURLTUGMZLYLDI-UHFFFAOYSA-N", "No match found in local DB for CURLTUGMZLYLDI-UHFFFAOYSA-N"], "molecule_not_found"),
        ("simple_water_atoms", None, [], None, ["Database file 'non_existent.db' not found or inaccessible."], "invalid_db"),
        ("simple_heavy_water_atoms", "XLYOFNOQVPJJNP-ZSJDYOACSA-N", [{"id": 1, "InChIKey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N", 'InChIKey14': "XLYOFNOQVPJJNP", "Name": "Water", "Formula": "H2O"}], False, ["Computed InChIKey: XLYOFNOQVPJJNP-ZSJDYOACSA-N"], "molecule_found by InChIKey14")
    ],
    ids=["molecule_found", "molecule_not_found", "invalid_db", "molecule_found by InChIKey14"]
)
def test_query_molecule_database_util(atoms_fixture, expected_inchikey, expected_data, full_match, log_messages, case_description, request, test_database, logger):
    """Test query_molecule_database_util covering found, not found, and missing DB."""
    atoms_obj = request.getfixturevalue(atoms_fixture)
    db_file = test_database if case_description != "invalid_db" else "non_existent.db"
    result_inchikey, result_data, match = query_molecule_database_util(atoms_obj, db_file, logger)
    if case_description == "molecule_found by InChIKey14":
        assert result_inchikey[14] == expected_inchikey[14], f"Failed: {case_description} (InChIKey mismatch)"
    else:
        assert result_inchikey == expected_inchikey, f"Failed: {case_description} (InChIKey mismatch)"
    assert result_data == expected_data, f"Failed: {case_description} (Molecule data mismatch)"
    assert full_match == match, f"Failed: {case_description} (Match flag mismatch)"

    if case_description in ["molecule_found", "molecule_found by InChIKey14"]:
        logger.info.assert_called_once_with(f"Computed InChIKey: {expected_inchikey}")
    elif case_description == "molecule_not_found":
        expected_calls = [call(msg) for msg in log_messages]
        logger.info.assert_has_calls(expected_calls, any_order=False)
    else:
        logger.error.assert_called_once_with(log_messages[0])

# ---------------------- Tests for generate_topology_util ----------------------
@pytest.fixture
def empty_system():
    """Returns a new System instance with empty attributes."""
    return System()

# Edge case: no molecule data provided
def test_generate_topology_util_no_molecule_data(empty_system, logger):
    inchikey = "DUMMY-KEY"
    molecule_data = None  # No data provided
    result = generate_topology_util(empty_system, inchikey, molecule_data, logger)
    assert result == [], f"No molecule data found for InChIKey {inchikey}."

# Happy path: system attributes updated when none exist
def test_generate_topology_util_updates_system_and_archive(empty_system, logger):
    inchikey = "DUMMY-KEY"
    molecule_data = ["some_data"]  # Non-empty data
    result = generate_topology_util(empty_system, inchikey, molecule_data, logger)
    assert empty_system.label == "molecule"
    assert empty_system.method == "parser"
    assert empty_system.building_block == "molecule"
    assert result == empty_system

# Edge case: existing system attributes should be preserved
def test_generate_topology_util_preserves_existing_system_values(logger):
    system = System()
    system.label = "pre-set label"
    system.method = "user"
    system.building_block = "monomer"
    inchikey = "DUMMY-KEY"
    molecule_data = ["some_data"]
    result = generate_topology_util(system, inchikey, molecule_data, logger)
    assert system.label == "pre-set label"
    assert system.method == "user"
    assert system.building_block == "monomer"
    assert result == system
