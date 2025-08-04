import os

import numpy as np

from ase import Atoms
from matid.geometry import get_dimensionality
from molid.search.db_lookup import basic_offline_search
from molid.utils.conversion import atoms_to_inchikey

# threshold below which we treat a cell as singular
SINGULAR_CELL_DET_THRESHOLD = 1e-8


def get_atoms_data(topology, atoms_ref, logger):
    """
    Retrieves atom data from a topology object.

    Args:
        topology: The topology object containing atom information.
        atoms_ref: Reference to ASE atoms or indices.
        logger: Logger for logging messages.

    Returns:
        The ASE Atoms object if available; otherwise, None.
    """
    # If specific indices are provided, select the first block
    if getattr(topology, 'indices', None) is not None:
        logger.info("Topology contains indices.")
        try:
            # indices may be a nested list
            idx = topology.indices[0]
            return atoms_ref.to_ase()[idx]
        except Exception as e:
            logger.error(f"Failed to apply topology indices: {e}")
            return None
    # Otherwise, use the full Atoms
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
    """
    positions = atoms.get_positions(copy=True)
    cell = atoms.get_cell()
    # Skip unwrapping if the cell is singular
    if abs(np.linalg.det(cell)) < SINGULAR_CELL_DET_THRESHOLD:
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
    atoms.center()
    return atoms


def validate_atom_count(atoms_data: Atoms, min_atoms: int, max_atoms: int, logger) -> bool:
    """
    Checks if the number of atoms is within acceptable limits.
    """
    num_atoms = len(atoms_data)
    if num_atoms < min_atoms:
        logger.warning(
            f"System has only {num_atoms} atoms; minimum required is {min_atoms}."
        )
        return False
    if num_atoms > max_atoms:
        logger.warning(
            f"System has {num_atoms} atoms; exceeds maximum allowed {max_atoms}."
        )
        return False
    return True


def validate_dimensionality(atoms: Atoms, logger) -> bool:
    """
    Validates that the system is 0D (non-periodic).
    """
    dimensionality = get_dimensionality(atoms)
    if dimensionality != 0:
        logger.warning(
            f"System is {dimensionality}D, only 0D systems are processed. Skipping normalization and continue."
        )
        return False
    return True


def query_molecule_database_util(atoms: Atoms, database_file: str, logger):
    """
    Compute an InChIKey from the ASE Atoms, then do an offline-basic lookup.
    Returns (inchikey, results_list, full_match_flag).
    """
    if not os.path.isfile(database_file):
        logger.error(f"Database file '{database_file}' not found or inaccessible.")
        return None, [], None

    try:
        inchikey = atoms_to_inchikey(atoms)
        logger.info(f"Computed InChIKey: {inchikey}")
    except Exception as e:
        logger.error(f"Error computing InChIKey: {e}")
        return None, [], None

    try:
        results = basic_offline_search(database_file, inchikey)
        if not results or not results[0]:
            logger.info(f"No match found in local DB for {inchikey}")
            return inchikey, [], None
    except Exception as e:
        logger.error(f"Error querying local DB: {e}")
        return None, [], None

    try:
        full_match = (results[0]["InChIKey"] == inchikey)
        return inchikey, results, full_match
    except Exception as e:
        logger.error(f"Error processing DB results: {e}")
        return None, [], None


def generate_topology_util(system, inchikey, molecule_data, logger):
    """
    Generates molecular topology data in NOMAD format.
    """
    if molecule_data is None:
        logger.info(f"No molecule data found for InChIKey {inchikey}.")
        return system

    if not getattr(system, 'label', None):
        system.label = 'molecule'
    if not getattr(system, 'method', None):
        system.method = 'parser'
    if not getattr(system, 'building_block', None):
        system.building_block = 'molecule'

    return system