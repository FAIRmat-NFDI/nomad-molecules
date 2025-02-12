# from pydantic import Field
# from nomad.config.models.plugins import NormalizerEntryPoint

# class MolecularNormalizerEntryPoint(NormalizerEntryPoint):
#     def load(self):
#         # Lazy import inside the function to avoid circular import
#         from nomad_plugin_molecules.normalizers.molecule import MolecularNormalizer
#         return MolecularNormalizer()

# molecularnormalizer_entry_point = MolecularNormalizerEntryPoint(
#     name='MolecularNormalizer',
#     description='Normalizer that identifies molecules in results.material.topology and adds new molecule-related information extracted from PubChem.',
# )



