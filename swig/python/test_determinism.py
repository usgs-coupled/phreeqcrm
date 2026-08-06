import numpy as np
import pytest

from phreeqcrm import BMIPhreeqcRM

from constants import FilePaths


@pytest.fixture
def create_yaml_file():
    try:
        import WriteYamlFile as yrm
        yrm.WriteYamlFile()
    except Exception as e:
        pytest.skip(f"Cannot create yaml file: {e}")


def test_two_instance_determinism(create_yaml_file):
    # Two instances in one process have different heap layouts; bitwise
    # equality catches any address-ordered iteration feeding results.
    def run():
        bmi = BMIPhreeqcRM()
        bmi.initialize(FilePaths.YAML)
        bmi.SetTimeStep(0.2 * 86400.0)
        bmi.RunCells()
        c = np.array(bmi.GetConcentrations())
        bmi.finalize()
        return c

    a = run()
    b = run()
    assert np.array_equal(a, b)
