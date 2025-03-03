import numpy as np
from ase import Atoms

from nomad.datamodel import EntryArchive
from nomad.datamodel.results import Material, System
from molid import query_pubchem_database


def create_atoms_object(atoms_data) -> Atoms:
    """Creates an ASE Atoms object from NOMAD atomic data."""
    # Convert positions (in meters) to angstroms
    atomic_positions = np.array(atoms_data.positions)
    atomic_positions_angstrom = atomic_positions * 1e10

    # Convert lattice vectors from meters to angstroms if available
    lattice_vectors = getattr(atoms_data, "lattice_vectors", None)
    if lattice_vectors is not None:
        lattice_vectors = np.array(lattice_vectors, dtype=float) * 1e10
    else:
        lattice_vectors = None

    atoms = Atoms(
        symbols=atoms_data.labels,
        positions=atomic_positions_angstrom.astype(float),
        cell=lattice_vectors,
        pbc=atoms_data.periodic
    )
    atoms = unwrap_atoms(atoms)
    return atoms


def unwrap_atoms(atoms: Atoms) -> Atoms:
    """Unwraps atoms positions using the minimum image convention."""
    positions = atoms.get_positions(copy=True)
    cell = atoms.get_cell()
    # Skip unwrapping if the cell is singular
    if abs(np.linalg.det(cell)) < 1e-8:
        return atoms

    ref = positions[0]
    for i in range(1, len(positions)):
        disp = positions[i] - ref
        try:
            frac_disp = np.linalg.solve(cell, disp)
        except np.linalg.LinAlgError:
            frac_disp = np.zeros_like(disp)
        disp -= cell.dot(np.round(frac_disp))
        positions[i] = ref + disp

    atoms.set_positions(positions)
    atoms.center()  # Optionally recenter the molecule.
    return atoms


def validate_atom_count(atoms: Atoms, min_atoms: int, max_atoms: int, logger) -> bool:
    """Check if the system has an acceptable number of atoms."""
    num_atoms = len(atoms)
    if num_atoms < min_atoms:
        logger.warning(
            f"System has only {num_atoms} atoms, which is less than the minimum required {min_atoms}. Skipping normalization."
        )
        print(f"System has only {num_atoms} atoms, which is less than the minimum required {min_atoms}. Skipping normalization.")
        return False
    if num_atoms > max_atoms:
        logger.warning(
            f"System has {num_atoms} atoms, which exceeds the maximum allowed {max_atoms}. Skipping normalization."
        )
        print(f"System has {num_atoms} atoms, which exceeds the maximum allowed {max_atoms}. Skipping normalization.")
        return False
    return True


def query_molecule_database_util(atoms: Atoms, logger):
    """Queries the local PubChem database using molid."""
    database_file = '../molecule_identification/OpenBabel/pubchem_data_FULL.db'
    try:
        inchikey, molecule_data = query_pubchem_database(atoms, database_file)
        logger.info(f"InChIKey: {inchikey}")
        return inchikey, molecule_data
    except Exception as e:
        logger.error(f"Error querying the local PubChem DB: {e}")
        return None, None


def generate_topology_util(archive, inchikey, molecule_data, logger) -> list:
    """Generates molecular topology data in a format compatible with NOMAD."""
    if not molecule_data:
        logger.warning(f"No molecule data found for InChIKey {inchikey}.")
        return []

    # Ensure the archive has the results/material section available.
    if not archive.results:
        archive.results = EntryArchive().results
    if not archive.results.material:
        archive.results.material = archive.results.m_create(Material)
    topology_container = archive.results.material.topology
    topology_entry = System(
        method='parser',
        label='molecule',
        building_block='molecule'
    )
    topology_container.append(topology_entry)
    print('topology_container after appending:', topology_container)
    return topology_container