import numpy as np
from numpy.testing import assert_array_almost_equal

from phreeqcrm import BMIPhreeqcRM

from constants import FilePaths

def test_set_value():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z0 = model.get_value_ptr("Temperature")
    z1 = np.zeros_like(z0) - 1

    model.set_value("Temperature", z1)

    new_z = model.get_value_ptr("Temperature")

    assert new_z is z0
    assert new_z is not z1
    assert_array_almost_equal(new_z, z1)

def test_set_value_at_indices():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    model.set_value_at_indices("Temperature", np.array([0, 2, 4]), np.array([-1.0, -1.0, -1.0]))

    new_z = model.get_value_ptr("Temperature")
    assert_array_almost_equal(new_z.take((0, 2, 4)), -1.0)

def test_set_value_SelectedOutputOn():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z0 = model.get_value_ptr("SelectedOutputOn")
    assert str(z0.dtype) == "int32"
    assert z0[0] == 1

    z1 = np.zeros_like(z0)
    model.set_value("SelectedOutputOn", z1)

    new_z = model.get_value_ptr("SelectedOutputOn")

    assert new_z is z0
    assert new_z is not z1
    assert new_z == z1

def test_numpy_integers_compatibility():
    # if this fails it probably has something
    # to do with the numpy file pyfragments.swg
    # not being integrated into PhreeqcRM.i

    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z1 = np.full(1, 0, dtype=int)
    n = model.GetNthSelectedOutputUserNumber(z1[0])
