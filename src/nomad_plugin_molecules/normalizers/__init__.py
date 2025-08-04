from nomad.config.models.plugins import NormalizerEntryPoint
from pydantic import Field


class MoleculeNormalizerEntryPoint(NormalizerEntryPoint):
    # model_config = ConfigDict(extra='ignore')
    MOLID_MASTER_DB: str = Field(
        'pubchem_data_FULL.db',
        description='Path to the PubChem “master” offline database'
    )
    MOLID_CACHE_DB: str = Field(
        'user_cache.db',
        description='Path to the PubChem cache database'
    )
    MOLID_MODE: str = Field(
        'offline-basic',
        description='MolID search mode (offline-basic, offline-advanced, …)'
    )
    max_atoms: int = Field(
        100,
        description='Maximum number of atoms per molecule to attempt matching'
    )
    min_atoms: int = Field(
        2,
        description='Minimum number of atoms per molecule to attempt matching'
    )

    def load(self):
        # Lazy import inside the function to avoid circular import
        # noqa: PLC0415 — we *want* this lazy import
        from nomad_plugin_molecules.normalizers.molecule import MoleculeNormalizer

        return MoleculeNormalizer(**self.dict())


moleculenormalizer = MoleculeNormalizerEntryPoint(
    name='MoleculeNormalizer',
    description='Normalizer that identifies molecules in results.material.topology and ' \
    'adds new molecule-related information extracted from PubChem.',
    master_db='./pubchem_data_FULL.db',
    level=6,
    MOLID_CACHE_DB = "./user_cache.db",
    MOLID_MODE = "offline-basic",
    max_atoms = 100,
    min_atoms = 2
)
