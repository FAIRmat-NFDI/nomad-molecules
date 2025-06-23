import logging
import pytest
import ase.build
from ase import Atoms
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo import runschema
from nomad.normalizing import normalizers

# -----------------------------------------------------------------------------
# Helper functions to create minimal NOMAD atoms data and archive
# -----------------------------------------------------------------------------
def get_section_system(atoms: Atoms):
    if runschema:
        system = runschema.system.System()
        system.atoms = runschema.system.Atoms(
            positions=atoms.get_positions() * 1e-10,
            labels=atoms.get_chemical_symbols(),
            lattice_vectors=atoms.get_cell() * 1e-10,
            periodic=atoms.get_pbc(),
        )
        return system

def create_archive():
    archive = EntryArchive()
    run = runschema.run.Run()
    archive.run.append(run)
    run.program = runschema.run.Program(name='VASP', version='4.6.35')

    # System
    water1 = ase.build.molecule('H2O')
    sys = water1
    sys.set_cell([10, 10, 10])
    sys.set_pbc(False)
    system = get_section_system(sys)
    run.system.append(system)

    # Calculation
    calc = runschema.calculation.Calculation()
    calc.system_ref = system
    run.calculation.append(calc)
    return run_normalize(archive)

from nomad.utils import get_logger
def run_normalize(entry_archive: EntryArchive) -> EntryArchive:
    for normalizer_class in normalizers:
        normalizer = normalizer_class(entry_archive)
        print("normalizer.normalizer: ", normalizer.normalizer)
        normalizer.normalize(logger=get_logger(__name__))
    return entry_archive

# -----------------------------------------------------------------------------
# Integration Test for MoleculeNormalizer
# -----------------------------------------------------------------------------
def test_molecule_normalizer_integration():
    """
    Integration test for MoleculeNormalizer.

    This test verifies that:
      - A NOMAD archive containing a water molecule is processed end-to-end.
      - The real atoms_utils functions are used as-is.
      - The PubChem query successfully retrieves the water record from the real full offline pubchem DB.
      - The generated topology (a System object) contains the expected attributes.
    """

    logger = logging.getLogger("test_molecule_normalizer")
    logger.setLevel(logging.INFO)
    archive = create_archive()
    assert len(archive.results.material.topology) == 1
    topology = archive.results.material.topology[0]

    # Check topology
    assert hasattr(topology, "building_block"), "Topology entry missing 'building_block' attribute."
    assert topology.building_block == "molecule", "Unexpected building_block in topology entry."
    assert hasattr(topology, "label"), "Topology entry missing 'label' attribute."
    assert topology.label == "original", "Unexpected label in topology entry."
    assert hasattr(topology, "method"), "Topology entry missing 'method' attribute."
    assert topology.method == "parser", "Unexpected method in topology entry."

    # Check cheminformatics
    assert hasattr(topology, 'cheminformatics'), "Topology entry missing 'cheminformatics' attribute."
    cheminformatics = topology.cheminformatics
    # InChIkey
    assert hasattr(cheminformatics, 'inchi_key'), "Cheminfomatics entry missing inchikey' attribute"
    assert cheminformatics.inchi_key == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    # InChI
    assert hasattr(cheminformatics, 'inchi'), "Cheminfomatics entry missing 'inchi' attribute"
    assert cheminformatics.inchi == 'InChI=1S/H2O/h1H2'
    # SMILES
    assert hasattr(cheminformatics, 'smiles'), "Cheminfomatics entry missing 'smiles' attribute"
    assert cheminformatics.smiles == 'O'