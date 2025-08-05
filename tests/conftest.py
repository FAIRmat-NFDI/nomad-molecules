import pytest

from nomad.datamodel.metainfo import runschema
# This is not how its suppose to be but for now the only solution
from simulationworkflowschema import load_modules
load_modules()

@pytest.fixture(scope="session", autouse=True)
def preload_nomad_metainfo():
    # Ensures metainfo schema is loaded early
    _ = runschema.run.Run