import numpy as np
from numpy.testing import assert_array_almost_equal

from phreeqcrm import BMIPhreeqcRM

def test_set_value():
    model = BMIPhreeqcRM()
    model.initialize("AdvectBMI_py.yaml")

    z0 = model.get_value_ptr("Temperature")
    z1 = np.zeros_like(z0) - 1

    model.set_value("Temperature", z1)

    new_z = model.get_value_ptr("Temperature")

    # assert new_z is z0    @todo
    assert new_z is not z1
    assert_array_almost_equal(new_z, z1)


def test_set_value_at_indices():
    model = BMIPhreeqcRM()
    model.initialize("AdvectBMI_py.yaml")

    # model.set_value_at_indices("Temperature", np.array([0, 2, 4]), np.array([-1, -1, -1]))
    model.set_value_at_indices("Temperature", np.array([0, 2, 4]), np.array([-1.0, -1.0, -1.0]))

    new_z = model.get_value_ptr("Temperature")
    assert_array_almost_equal(new_z.take((0, 2, 4)), -1.0)
