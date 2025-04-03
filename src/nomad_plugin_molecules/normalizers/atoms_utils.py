import os
import numpy as np
from ase import Atoms
from matid.geometry import get_dimensionality
from nomad.datamodel.results import Material, Results
from molid import query_pubchem_database


def get_atoms_data(topology, topology_original, logger):
    """
    Retrieves atom data from a topology object.

    Args:
        topology: The topology object containing atom information.
        logger: Logger for logging messages.

    Returns:
        The atoms data if available; otherwise, None.
    """
    if topology.indices is not None:
        logger.info("Topology contains indices.")
        ase_atoms = topology_original.atoms_ref.to_ase()[topology.indices[0]]
        return ase_atoms
    elif topology.atoms is not None:
        return topology.atoms.to_ase()
    elif topology.atoms_ref is not None:
        return topology.atoms_ref.to_ase()
    else:
        logger.warning("No atoms data found (neither atoms nor atoms_ref exist).")
        return None


def wrap_atoms(atoms: Atoms) -> Atoms:
    """
    Unwraps atom positions using the minimum image convention.

    Args:
        atoms: ASE Atoms object.

    Returns:
        ASE Atoms object with unwrapped positions.
    """
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


def validate_atom_count(atoms_data: Atoms, min_atoms: int, max_atoms: int, logger) -> bool:
    """
    Checks if the number of atoms in the system is within acceptable limits.

    Args:
        atoms_data: ASE Atoms object.
        min_atoms: Minimum number of atoms required.
        max_atoms: Maximum number of atoms allowed.
        logger: Logger for logging messages.

    Returns:
        True if the atom count is acceptable; otherwise, False.
    """
    num_atoms = len(atoms_data)
    if num_atoms < min_atoms:
        logger.info(
            f"System has only {num_atoms} atoms, which is less than the minimum required {min_atoms}. Continue."
        )
        return False
    if num_atoms > max_atoms:
        logger.info(
            f"System has {num_atoms} atoms, which exceeds the maximum allowed {max_atoms}. Continue."
        )
        return False
    return True


def validate_dimensionality(atoms: Atoms, logger) -> bool:
    """
    Validates that the system is 0D (non-periodic) for processing.

    Args:
        atoms: ASE Atoms object.
        logger: Logger for logging messages.

    Returns:
        True if the system is 0D; otherwise, False.
    """

    dimensionality = get_dimensionality(atoms)
    if dimensionality != 0:
        logger.info(
            f"System is {dimensionality}D, only 0D systems are processed. Continue."
        )
        return False
    return True


def query_molecule_database_util(atoms: Atoms, database_file, logger):
    """
    Queries the local PubChem database using the provided atoms data.

    Args:
        atoms: ASE Atoms object.
        database_file: Path to the SQLite database file.
        logger: Logger for logging messages.

    Returns:
        Tuple of (InChIKey, molecule_data) if successful; otherwise, (None, None).
    """

    if not os.path.isfile(database_file):
        logger.error(f"Database file '{database_file}' not found or inaccessible.")
        return None, None

    try:
        inchikey, molecule_data = query_pubchem_database(atoms, database_file)
        logger.info(f"InChIKey: {inchikey}")
        return inchikey, molecule_data
    except Exception as e:
        logger.error(f"Error querying the local PubChem DB: {e}")
        return None, None


def generate_topology_util(system, inchikey, molecule_data, logger):
    """
    Generates molecular topology data in a format compatible with NOMAD.

    Args:
        system: The original topology/system object.
        inchikey: The InChIKey of the molecule.
        molecule_data: Data retrieved from the database.
        logger: Logger for logging messages.

    Returns:
        Updated system object with molecular topology data, or None if generation fails.
    """
    if not molecule_data:
        logger.info(f"No molecule data found for InChIKey {inchikey}.")
        return []

    # if system.building_block:
    #     logger.warning("building_block already set")

    if not system.label:
        system.label = 'molecule'
    if not system.method:
        system.method = 'parser'
    if not system.building_block:
        system.building_block = 'molecule'

    return system