from nomad.config import config
from nomad.normalizing import Normalizer
from nomad.datamodel.results import Cheminformatics
from molid.utils.settings import save_config

from .atoms_utils import (
    generate_topology_util,
    get_atoms_data,
    query_molecule_database_util,
    validate_atom_count,
    validate_dimensionality,
    wrap_atoms,
)

class MoleculesNormalizer(Normalizer):
    """Normalizer for molecular data extraction."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # This is not how its suppose to be but for now the only solution
        from simulationworkflowschema import load_modules
        load_modules()

    def normalize(self, archive, logger=None) -> list:
        """
        Normalizes molecular data from the provided archive.

        Args:
            archive: The archive containing molecular data.
            logger: Optional logger for logging messages.

        Returns:
            List of processed topology/system objects.
        """
        entry_point_id = "nomad_molecules.normalizers:moleculesnormalizer"
        # maybe easier to mock this function for testing with db_path and so on
        cfg = config.get_plugin_entry_point(entry_point_id)

        molid_mode = cfg.molid_mode
        if molid_mode == "offline-basic":
            database_file = cfg.molid_master_db
            save_config(master_db = cfg.molid_master_db, mode = molid_mode)
        else:
            database_file = cfg.molid_cache_db
            save_config(cache_db = cfg.molid_cache_db, mode = molid_mode)
        max_atoms = cfg.max_atoms
        min_atoms = cfg.min_atoms

        if logger:
            logger.info("Starting molecular normalization process.")
        # Quick existence checks
        results = getattr(archive, "results", None)
        if not results:
            if logger:
                logger.info("No results found in archive. Exiting normalization.")
            return

        material = getattr(results, "material", None)
        if not material:
            if logger:
                logger.info("No material found in archive. Exiting normalization.")
            return

        topologies = getattr(material, "topology", []) or []
        if not topologies:
            if logger:
                logger.info("No topology data found in archive. Exiting normalization.")
            return

        # Process each topology
        for topology in topologies:
            label = getattr(topology, "label", None)
            if label == 'conventional cell':
                if logger:
                    logger.debug("Skipping 'conventional cell' topology.")
                continue
            if label == 'original' and len(topologies) > 1:
                if logger:
                    logger.debug("Skipping 'original' topology in a composite system.")
                continue

            # Determine atoms reference
            atoms_ref = None
            if getattr(topology, "atoms_ref", None):
                atoms_ref = topology.atoms_ref
            elif getattr(topology, "atoms", None):
                atoms_ref = topology.atoms
            else:
                if logger:
                    logger.warning("No atoms or atoms_ref found in topology; skipping.")
                continue

            # Retrieve ASE Atoms
            ase_atoms = get_atoms_data(topology, atoms_ref, logger)
            if ase_atoms is None:
                continue

            # Validate atom count
            if not validate_atom_count(ase_atoms, min_atoms, max_atoms, logger):
                continue

            # Unwrap periodic coordinates if fully periodic
            if all(ase_atoms.pbc):
                ase_atoms = wrap_atoms(ase_atoms)

            # Ensure a 0D (non-periodic) molecule
            if not validate_dimensionality(ase_atoms, logger):
                continue

            # Query the local PubChem offline DB
            inchikey, molecule_data, matched_full = query_molecule_database_util(
                ase_atoms, database_file, logger
            )
            if inchikey is None or matched_full is None:
                continue

            # Attach cheminformatics metadata
            if matched_full:
                match_type = 'full structure'
                topology.cheminformatics = Cheminformatics(
                    smiles=molecule_data[0]["SMILES"],
                    inchi_key=molecule_data[0]["InChIKey"],
                    inchi=molecule_data[0]["InChI"],
                    match_type=match_type
                )
            else:
                if logger:
                    logger.info(
                        "Molecule identification based on the first block of the InChIKey "
                        "(connectivity only, first 14 characters)."
                    )
                match_type = 'skeleton (core)'
                topology.cheminformatics = Cheminformatics(
                    inchi_key=inchikey,
                    match_type=match_type
                )

            # Generate NOMAD-compatible topology data
            generate_topology_util(topology, inchikey, molecule_data, logger)
