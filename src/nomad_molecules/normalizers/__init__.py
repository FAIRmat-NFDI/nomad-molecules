from nomad.config.models.plugins import NormalizerEntryPoint
from pydantic import Field


class MoleculesNormalizerEntryPoint(NormalizerEntryPoint):
    # model_config = ConfigDict(extra='ignore')
    molid_master_db: str = Field(
        'pubchem_data_FULL.db',
        description='Path to the PubChem “master” offline database'
    )
    molid_cache_db: str = Field(
        'user_cache.db',
        description='Path to the PubChem cache database'
    )
    molid_mode: str = Field(
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
        from nomad_molecules.normalizers.molecules import MoleculesNormalizer

        return MoleculesNormalizer(**self.dict())


moleculesnormalizer = MoleculesNormalizerEntryPoint(
    name='MoleculesNormalizer',
    description='Normalizer that identifies molecules in results.material.topology and ' \
    'adds new molecule-related information extracted from PubChem.',
    level=6,
    molid_master_db='./pubchem_data_FULL.db',
    molid_cache_db = "./user_cache.db",
    molid_mode = "offline-basic",
    max_atoms = 100,
    min_atoms = 2
)
