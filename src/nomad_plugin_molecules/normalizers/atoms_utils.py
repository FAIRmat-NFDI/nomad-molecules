import os
import numpy as np
from ase import Atoms
from matid.geometry import get_dimensionality
from nomad.datamodel.results import Material, Results
from molid.utils.conversion import atoms_to_inchikey
from molid.search.db_lookup  import basic_offline_search


def get_atoms_data(topology, atoms_ref, logger):
    """
    Retrieves atom data from a topology object.

    Args:
        topology: The topology object containing atom information.
        logger: Logger for logging messages.

    Returns:
        The atoms data if available; otherwise, None.
    """
    # import pdb; pdb.set_trace()
    if topology.indices is not None:
        logger.info("Topology contains indices.")
        atoms = atoms_ref.to_ase()[topology.indices[0]]
        return atoms
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
        logger.warning(
            f"System has only {num_atoms} atoms, which is less than the minimum required {min_atoms}. Skipping normalization and continue."
        )
        return False
    if num_atoms > max_atoms:
        logger.warning(
            f"System has {num_atoms} atoms, which exceeds the maximum allowed {max_atoms}. Skipping normalization and continue."
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
        logger.warning(
            f"System is {dimensionality}D, only 0D systems are processed. Skipping normalization and continue."
        )
        return False
    return True


def query_molecule_database_util(atoms: Atoms, database_file, logger):
    """
    Compute an InChIKey from the ASE Atoms, then do an offline-basic lookup.
    Returns (inchikey, results_list, full_match_flag).
    """
    if not os.path.isfile(database_file):
        logger.error(f"Database file '{database_file}' not found or inaccessible.")
        return None, [], None

    try:
        # 1) compute InChIKey
        inchikey = atoms_to_inchikey(atoms)
        logger.info(f"Computed InChIKey: {inchikey}")

        # 2) query the master DB for full‐ or 14‐char match
        results = basic_offline_search(database_file, inchikey)

        if not results or not results[0]:
            logger.info(f"No match found in local DB for {inchikey}")
            return inchikey, [], None

        # 3) detect whether it was a full‐Key match
        full_match = (results[0]["InChIKey"] == inchikey)
        return inchikey, results, full_match

    except Exception as e:
        logger.error(f"Error querying local PubChem DB: {e}")
        return None, [], None


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

    if not system.label:
        system.label = 'molecule'
    if not system.method:
        system.method = 'parser'
    if not system.building_block:
        system.building_block = 'molecule'

    return system