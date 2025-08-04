import pytest

from nomad.datamodel.metainfo import runschema

@pytest.fixture(scope="session", autouse=True)
def preload_nomad_metainfo():
    # Ensures metainfo schema is loaded early
    _ = runschema.run.Run