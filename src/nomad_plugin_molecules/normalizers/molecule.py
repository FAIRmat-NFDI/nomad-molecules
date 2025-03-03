from nomad.normalizing import Normalizer
from nomad.datamodel import EntryArchive
from .atoms_utils import create_atoms_object, validate_atom_count, query_molecule_database_util, generate_topology_util


class MoleculeNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""

    def __init__(self, max_atoms=50, min_atoms=2, **kwargs):
        super().__init__(**kwargs)
        self.max_atoms = max_atoms
        self.min_atoms = min_atoms

    def normalize(self, archive: EntryArchive, logger=None) -> list:
        self.logger.info("Starting molecular normalization process.")

        if not archive.run:
            print('archive.run: ', archive.run)
            return
        # if not EntryArchive().results:
        #     print('EntryArchive().results: ', EntryArchive().results)
        #     return

        try:
            # Convert NOMAD atoms data into an ASE Atoms object.
            atoms_data = archive.run[0].system[0].atoms
            atoms = create_atoms_object(atoms_data)

            # Validate atom count (must be >1 and <= max_atoms)
            if not validate_atom_count(atoms, self.min_atoms, self.max_atoms, self.logger):
                return []

            # Query the local PubChem database.
            inchikey, molecule_data = query_molecule_database_util(atoms, self.logger)
            print('inchikey:', inchikey)

            # Generate and return topology information.
            topology = generate_topology_util(archive, inchikey, molecule_data, self.logger)
            print('Topology: ', topology)
            return topology
        except Exception as e:
            self.logger.error(f"Error in normalization: {e}", exc_info=True)
            print('Error (molecule):', e)
            return
