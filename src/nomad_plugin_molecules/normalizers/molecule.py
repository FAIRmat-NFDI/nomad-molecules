# normalizer.py (Updated to match example structure)

import numpy as np
from ase import Atoms
from molid import query_pubchem_database
from nomad.utils import get_logger
from nomad.normalizing.normalizer import Normalizer


class MolecularNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""

    normalizer_level = 1

    def __init__(self):
        super().__init__()

    def normalize(self, archive, logger) -> None:
        # super().normalize(archive, logger)
        if logger is None:
            logger = get_logger(__name__)
        self.logger = logger
        logger.info("Starting molecular normalization process.")

        if not archive.run:
            return

        try:
            atoms_data = archive.run[0].system[0].atoms
            atoms = self.create_atoms_object(atoms_data)
            inchikey, molecule_data = self.query_molecule_database(atoms)
            self.attach_molecule_data(archive, inchikey, molecule_data)
        except Exception as e:
            logger.error(f"Error in normalization: {e}", exc_info=True)

    def create_atoms_object(self, atoms_data) -> Atoms:
        """Creates an ASE Atoms object from NOMAD atomic data."""
        positions = np.array(atoms_data.positions)
        lattice_vectors = getattr(atoms_data, "lattice_vectors", None)
        atoms = Atoms(
            symbols=atoms_data.labels,
            positions=positions.astype(float),
            cell=np.array(lattice_vectors, dtype=float) if lattice_vectors else None,
            pbc=atoms_data.periodic
        )
        return atoms

    def query_molecule_database(self, atoms: Atoms):
        """Queries the local PubChem database using molid."""
        database_file = '../molecule_identification/OpenBabel/pubchem_data_FULL.db'
        try:
            inchikey, molecule_data = query_pubchem_database(atoms, database_file)
            self.logger.info(f"InChIKey: {inchikey}")
            return inchikey, molecule_data
        except Exception as e:
            self.logger.error(f"Error querying the local PubChem DB: {e}")
            return None, None

    def attach_molecule_data(self, archive, inchikey, molecule_data):
        """Attaches molecular data to the NOMAD archive."""
        from nomad.datamodel import EntryArchive
        from nomad.datamodel.results import Material, Topology

        if molecule_data:
            if not archive.results:
                archive.results = EntryArchive().results
            if not archive.results.material:
                archive.results.material = archive.results.m_create(Material)
            if not archive.results.material.topology:
                archive.results.material.topology = []
            topology_entry = archive.results.material.topology.m_create(Topology)
            topology_entry.building_block = 'molecule'
            self.logger.info(f"Molecular data attached: {inchikey}")
        else:
            self.logger.warning(f"No molecule data found for InChIKey {inchikey}.")
