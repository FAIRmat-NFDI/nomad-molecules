from nomad.normalizing import Normalizer
from nomad.datamodel.results import Cheminformatics
from .atoms_utils import wrap_atoms, validate_atom_count, validate_dimensionality, query_molecule_database_util, generate_topology_util, get_atoms_data


class MoleculeNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""

    def __init__(self, max_atoms=100, min_atoms=2, database_file=None, **kwargs):
        super().__init__(**kwargs)
        self.max_atoms = max_atoms
        self.min_atoms = min_atoms
        self.database_file = database_file or '../molecule_identification/pubchem_data_FULL.db'

    def normalize(self, archive, logger=None) -> list:
        """
        Normalizes molecular data from the provided archive.

        Args:
            archive: The archive containing molecular data.
            logger: Optional logger for logging messages.

        Returns:
            List of processed topology/system objects.
        """
        logger.info("Starting molecular normalization process.")

        if not getattr(archive, "results", False):
            logger.info("No results found in archive. Exiting normalization.")
            return []
        elif not getattr(archive.results, "material", False):
            logger.info("No material found in archive. Exiting normalization.")
            return []

        try:
            topologies = archive.results.material.topology
        except AttributeError as e:
            logger.error("Topology data not found in archive.", exc_info=True)
            return []

        for topology in topologies:
            if topology.label == 'original':
                topology_original = topology

            if topology.label == 'conventional cell':
                logger.debug("Skipping 'conventional cell' topology.")
                continue

            ase_atoms = get_atoms_data(topology, topology_original, logger)

            if ase_atoms is None:
                continue

            # Validate atom count (must be >1 and <= max_atoms)
            if not validate_atom_count(ase_atoms, self.min_atoms, self.max_atoms, logger):
                continue

            if all(ase_atoms.pbc):
                ase_atoms = wrap_atoms(ase_atoms)

            if not validate_dimensionality(ase_atoms, logger):
                continue

            # Query the local PubChem database.
            inchikey, molecule_data, matched_full = query_molecule_database_util(ase_atoms, self.database_file, logger)
            if not matched_full:
                logger.info("Molecule identification based on the first block of the InChIKey (connectivity only, first 14 characters).")
            topology.cheminformatics = Cheminformatics(
                smiles=molecule_data[0][1],
                inchi_key = molecule_data[0][2],
                inchi=molecule_data[0][3],
            )

            if inchikey is None:
                continue

            topology = generate_topology_util(topology, inchikey, molecule_data, logger)

        return topologies
