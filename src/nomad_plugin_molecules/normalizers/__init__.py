from pydantic import Field
from nomad.config.models.plugins import NormalizerEntryPoint


class MoleculeNormalizerEntryPoint(NormalizerEntryPoint):
    def load(self):
        from nomad_plugin_molecules.normalizers.molecule import MoleculeNormalizer

        return MoleculeNormalizer(**self.dict())


moleculenormalizer = MoleculeNormalizerEntryPoint(
    name='MoleculeNormalizer',
    description='Normalizer that identifies molecules in results.material.topology and adds new molecule related information extracted from PubChem.',
)
