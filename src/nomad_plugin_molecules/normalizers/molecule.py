# nomad_plugin_molecules/normalizers/molecule.py

import numpy as np
from ase import Atoms
from molid import query_pubchem_database
from nomad.utils import get_logger
# from nomad.normalizing.normalizer import Normalizer


class MolecularNormalizer:
    """Normalizer for molecular data extraction."""

    normalizer_level = 1

    def __init__(self, entry_archive, logger=None):
        self.entry_archive = entry_archive
        self.logger = logger or get_logger(__name__)

    def normalize(self) -> list:
        self.logger.info("Starting molecular normalization process.")

        if not self.entry_archive.run:
            return self.entry_archive

        try:
            # Get the atoms section from the archive
            atoms_data = self.entry_archive.run[0].system[0].atoms
            atoms = self.create_atoms_object(atoms_data)
            inchikey, molecule_data = self.query_molecule_database(atoms)
            # Optionally, attach detailed molecule data:
            # self.attach_molecule_data(inchikey, molecule_data)
            return self.generate_topology(inchikey, molecule_data)
        except Exception as e:
            self.logger.error(f"Error in normalization: {e}", exc_info=True)
            return self.entry_archive

    def create_atoms_object(self, atoms_data) -> Atoms:
        """Creates an ASE Atoms object from NOMAD atomic data."""
        # Get the positions as a numpy array (in meters)
        atomic_positions = np.array(atoms_data.positions)
        # Convert from meters to angstrom by multiplying by 1e10
        atomic_positions_angstrom = atomic_positions * 1e10

        lattice_vectors = getattr(atoms_data, "lattice_vectors", None)
        atoms = Atoms(
            symbols=atoms_data.labels,
            positions=atomic_positions_angstrom.astype(float),
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

    # def attach_molecule_data(self, inchikey, molecule_data):
    #     """Attaches molecular data to the NOMAD archive."""
    #     from nomad.datamodel.datamodel import EntryArchive
    #     from nomad.datamodel.results import Material
    #     from nomad.datamodel.results.materials import Topology


    #     if molecule_data:
    #         if not self.entry_archive.results:
    #             self.entry_archive.results = EntryArchive().results
    #         if not self.entry_archive.results.material:
    #             self.entry_archive.results.material = self.entry_archive.results.m_create(Material)
    #         if not self.entry_archive.results.material.topology:
    #             self.entry_archive.results.material.topology = []
    #         topology_entry = self.entry_archive.results.material.topology.m_create(Topology)
    #         topology_entry.building_block = 'molecule'
    #         self.logger.info(f"Molecular data attached: {inchikey}")
    #     else:
    #         self.logger.warning(f"No molecule data found for InChIKey {inchikey}.")

    def generate_topology(self, inchikey, molecule_data) -> list:
        """Returns molecular topology data in a format compatible with NOMAD."""
        if not molecule_data:
            self.logger.warning(f"No molecule data found for InChIKey {inchikey}.")
            return []

        from nomad.datamodel.datamodel import EntryArchive
        # For creating a new material section if needed.
        from nomad.datamodel.results import Material, System

        # Ensure the archive has the results/material section available.
        if not self.entry_archive.results:
            self.entry_archive.results = EntryArchive().results
        if not self.entry_archive.results.material:
            self.entry_archive.results.material = self.entry_archive.results.m_create(Material)
        # import pdb; pdb.set_trace()
        topology_container = self.entry_archive.m_xpath('results.material.topology')
        if not topology_container:
            topology_container = self.entry_archive.results.material.topology

        # Create a new System to represent the molecule topology.
        topology_entry = System(
            method='parser',
            label='molecule',
            building_block='molecule'
        )
        # Optionally, you can set additional properties (e.g. symmetry, material_id, etc.)
        # For example:
        # topology_entry.material_id = <your generated id>

        # Append the new system to the topology container.
        topology_container.append(topology_entry)

        return [topology_entry]


