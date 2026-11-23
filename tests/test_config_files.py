"""Every configuration file shipped in configFiles/ must load into the operadar Config dataclass."""
from pathlib import Path

import pytest

from operadar.load_config_file import load_configuration_file
from operadar.forward_operator import Config, load_config_instance

CONFIG_DIR = Path(__file__).resolve().parents[1] / "configFiles"
CONFIG_FILES = sorted(p.name for p in CONFIG_DIR.glob("*.py") if not p.name.startswith("_"))


@pytest.mark.parametrize("name", CONFIG_FILES)
def test_config_file_loads(name):
    conf = load_config_instance(load_configuration_file(name))
    assert isinstance(conf, Config)
    assert conf.radar_band
    assert isinstance(conf.sensitivity_threshold, bool)
