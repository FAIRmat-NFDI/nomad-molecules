import sqlite3
import logging
import pytest
from nomad.datamodel import EntryArchive
from nomad_plugin_molecules.normalizers.molecule import MoleculeNormalizer

# -----------------------------------------------------------------------------
# Helper functions to create minimal NOMAD atoms data and archive
# -----------------------------------------------------------------------------
def create_archive():
    """
    Create a minimal NOMAD archive containing a water molecule.
    The atoms are set with periodic boundaries to trigger wrapping.
    """
    species=[8, 1, 1]
    labels=["O", "H", "H"]
    positions=[[2.5e-10, 2.5e-10, 2.5e-10],
                [3.257e-10, 3.086e-10, 2.5e-10],
                [1.743e-10, 3.086e-10, 2.5e-10]]
    periodic=[False, False, False]
    cell = [[5.0e-10, 0, 0], [0, 5.0e-10, 0], [0, 0, 5.0e-10]]
    data = {
        "run": [{
            "system": [{
                "atoms": {
                    "species": species,
                    "labels": labels,
                    "positions": positions,
                    "lattice_vectors": cell,
                    "periodic": periodic
                }
            }]
        }],
        "results":{
            "material":{
                "topology":
                [
                    {"atoms_ref": "#/run/0/system/0/atoms",
                    }
                ]
            }
        }
    }
    return EntryArchive().m_from_dict(data)


# -----------------------------------------------------------------------------
# Fixture to set up a temporary SQLite database with expected PubChem data.
# -----------------------------------------------------------------------------
@pytest.fixture
def temporary_db(tmp_path):
    """
    Create a temporary SQLite database file with one record for water.
    The table schema and record match what the normalizer expects.
    """
    db_file = tmp_path / "pubchem_data_FULL.db"
    conn = sqlite3.connect(str(db_file))
    cursor = conn.cursor()
    cursor.execute(
        "CREATE TABLE compound_data (id INTEGER PRIMARY KEY, " \
        "InChIKey TEXT UNIQUE, " \
        "InChI TEXT, " \
        "SMILES TEXT," \
        "Name TEXT, " \
        "Formula TEXT)"
    )
    # InChIKey for water, matching what is used in atoms_utils.
    water_inchikey = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    water_inchi = "InChI=1S/H2O/h1H2"
    water_smiles = "O"
    cursor.execute(
        "INSERT INTO compound_data (InChIKey, InChI, SMILES, Name, Formula) VALUES (?, ?, ?, ?, ?)",
        (water_inchikey, water_inchi, water_smiles, "Water", "H2O")
    )
    conn.commit()
    conn.close()
    return str(db_file)

# -----------------------------------------------------------------------------
# Integration Test for MoleculeNormalizer
# -----------------------------------------------------------------------------
def test_molecule_normalizer_integration(temporary_db):
    """
    Integration test for MoleculeNormalizer.

    This test verifies that:
      - A NOMAD archive containing a water molecule is processed end-to-end.
      - The real atoms_utils functions are used as-is.
      - The PubChem query successfully retrieves the water record from the temporary DB.
      - The generated topology (a System object) contains the expected attributes.
    """
    # Create a test NOMAD archive with a water molecule.
    archive = create_archive()

    # Set up a simple logger.
    logger = logging.getLogger("test_molecule_normalizer")
    logger.setLevel(logging.INFO)

    # Instantiate the MoleculeNormalizer.
    normalizer = MoleculeNormalizer(database_file=temporary_db)

    # Run the normalization process.
    topology = normalizer.normalize(archive, logger=logger)

    # Validate the topology result.
    # The generate_topology_util function returns a list with one topology entry (a System object)
    # that should have attributes: method, label, and building_block.

    assert isinstance(topology, list), "Topology output should be a list."
    assert len(topology) == 1, "Expected one topology entry."

    topo_entry = topology[0]
    # Check the expected attributes.
    assert hasattr(topo_entry, "method"), "Topology entry missing 'method' attribute."
    assert topo_entry.method == "parser", "Unexpected method in topology entry."
    assert hasattr(topo_entry, "label"), "Topology entry missing 'label' attribute."
    assert topo_entry.label == "molecule", "Unexpected label in topology entry."
    assert hasattr(topo_entry, "building_block"), "Topology entry missing 'building_block' attribute."
    assert topo_entry.building_block == "molecule", "Unexpected building_block in topology entry."
