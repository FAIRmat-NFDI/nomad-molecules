from nomad.normalizing.normalizer import Normalizer
from nomad.normalizing.common import ase_atoms_from_nomad_atoms
from nomad_plugin_molecules.common.molecule_information import get_molecule_information
from ase import Atoms


class MoleculeNormalizer(Normalizer):
    """This normalizers is used to add basic molecule information for systems
    stored in NOMAD.
    """
    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        topology = archive.results.materials.topology
        for system in topology:
            if system.atoms and system.atoms.n_atoms < 20:
                ase_atoms = ase_atoms_from_nomad_atoms(system.atoms)
                info = get_molecule_information(ase_atoms)
                system.inchi_key = info.get('inchi_key')
