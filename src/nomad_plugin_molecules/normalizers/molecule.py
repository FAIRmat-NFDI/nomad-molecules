import numpy as np
from ase import Atoms

from nomad.normalizing import Normalizer
from nomad.datamodel import EntryArchive
from nomad.datamodel.results import Material, System
from molid import query_pubchem_database


class MoleculeNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""
    def normalize(self, archive: EntryArchive, logger=None) -> list:
        self.logger.info("Starting molecular normalization process.")

        if not archive.run:
            return

        try:
            # Get the atoms section from the archive
            atoms_data = archive.run[0].system[0].atoms
            atoms = self.create_atoms_object(atoms_data)
            inchikey, molecule_data = self.query_molecule_database(atoms)
            return self.generate_topology(archive, inchikey, molecule_data)
        except Exception as e:
            self.logger.error(f"Error in normalization: {e}", exc_info=True)
            return

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

    def generate_topology(self, archive, inchikey, molecule_data) -> list:
        """Returns molecular topology data in a format compatible with NOMAD."""
        if not molecule_data:
            self.logger.warning(f"No molecule data found for InChIKey {inchikey}.")
            return []

        # Ensure the archive has the results/material section available.
        if not archive.results:
            archive.results = EntryArchive().results
        if not archive.results.material:
            archive.results.material = archive.results.m_create(Material)
        # import pdb; pdb.set_trace()
        topology_container = archive.m_xpath('results.material.topology')
        if not topology_container:
            topology_container = archive.results.material.topology

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


