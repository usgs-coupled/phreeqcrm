import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less, assert_array_equal

from phreeqcrm import bmi_phreeqcrm

def test_get_initial_value():
    model = bmi_phreeqcrm()
    model.initialize()

    z0 = model.get_value_ptr("Porosity")

def test_get_value_copy():
    model = bmi_phreeqcrm()
    model.initialize()

    dest0 = np.empty(model.get_grid_size(0), dtype=float)
    dest1 = np.empty(model.get_grid_size(0), dtype=float)

    z0 = model.get_value("Porosity", dest0)
    z1 = model.get_value("Porosity", dest1)

    assert z0 is not z1
    assert_array_almost_equal(z0, z1)

def test_get_value_copy_int_scalar():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    dest0 = np.empty((1,), dtype=int)
    dest1 = np.empty((1,), dtype=int)

    z0 = model.get_value("ComponentCount", dest0)
    z1 = model.get_value("ComponentCount", dest1)

    assert z0 is not z1
    assert z0 is dest0
    assert z1 is dest1
    assert_array_equal(z0, z1)

def test_get_value_copy_float_scalar():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    dest0 = np.empty((1,), dtype=float)
    dest1 = np.empty((1,), dtype=float)

    z0 = model.get_value("Time", dest0)
    z1 = model.get_value("Time", dest1)

    assert z0 is not z1
    assert z0 is dest0
    assert z1 is dest1
    assert_array_almost_equal(z0, z1)

def test_get_value_copy_str():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    components = model.get_value_ptr("Components")
    dest0 = np.empty_like(components)
    dest1 = np.empty_like(components)

    z0 = model.get_value("Components", dest0)
    z1 = model.get_value("Components", dest1)

    assert z0 is not z1
    assert_array_equal(z0, z1)

def test_get_value_pointer():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    dest1 = np.empty(model.get_grid_size(0), dtype=float)

    z0 = model.get_value_ptr("Porosity")
    z1 = model.get_value("Porosity", dest1)

    assert z0 is not z1
    assert_array_almost_equal(z0.flatten(), z1)

    # this will throw if solutions aren't set for each cell (ie in initialize())
    for _ in range(10):
        model.update()
    assert z0 is model.get_value_ptr("Porosity")


def test_get_value_at_indices():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    temps = np.linspace(18.0, 30.0, model.get_grid_size(0))
    model.set_value("Temperature", temps)

    dest = np.empty(3, dtype=float)

    z0 = model.get_value_ptr("Temperature")
    z1 = model.get_value_at_indices("Temperature", dest, np.array([0, 2, 4]))
    # z1 = model.get_value_at_indices("Temperature", dest, [0, 2, 4]) @todo handle list

    assert_array_almost_equal(z0.take((0, 2, 4)), z1)

def test_value_size():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    z = model.get_value_ptr("Temperature")
    assert model.get_grid_size(0) == z.size


def test_value_nbytes():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    z = model.get_value_ptr("Temperature")
    assert model.get_var_nbytes("Temperature") == z.nbytes