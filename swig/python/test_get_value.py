import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less, assert_array_equal

from phreeqcrm import BMIPhreeqcRM

from constants import FilePaths

def test_get_initial_value():
    model = BMIPhreeqcRM()
    model.initialize()

    z0 = model.get_value_ptr("Porosity")

def test_get_value_copy():
    model = BMIPhreeqcRM()
    model.initialize()

    dest0 = np.empty(model.get_grid_size(0), dtype=float)
    dest1 = np.empty(model.get_grid_size(0), dtype=float)

    z0 = model.get_value("Porosity", dest0)
    z1 = model.get_value("Porosity", dest1)

    assert z0 is not z1
    assert_array_almost_equal(z0, z1)

def test_get_value_copy_int_scalar():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    dest0 = np.empty((1,), dtype=int)
    dest1 = np.empty((1,), dtype=int)

    z0 = model.get_value("ComponentCount", dest0)
    z1 = model.get_value("ComponentCount", dest1)

    assert z0 is not z1
    assert z0 is dest0
    assert z1 is dest1
    assert_array_equal(z0, z1)

def test_get_value_copy_float_scalar():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    dest0 = np.empty((1,), dtype=float)
    dest1 = np.empty((1,), dtype=float)

    z0 = model.get_value("Time", dest0)
    z1 = model.get_value("Time", dest1)

    assert z0 is not z1
    assert z0 is dest0
    assert z1 is dest1
    assert_array_almost_equal(z0, z1)

def test_get_value_copy_str():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    components = model.get_value_ptr("Components")
    dest0 = np.empty_like(components)
    dest1 = np.empty_like(components)

    z0 = model.get_value("Components", dest0)
    z1 = model.get_value("Components", dest1)

    assert z0 is not z1
    assert_array_equal(z0, z1)

def test_get_value_pointer():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    dest1 = np.empty(model.get_grid_size(0), dtype=float)

    z0 = model.get_value_ptr("Porosity")
    z1 = model.get_value("Porosity", dest1)

    assert z0 is not z1
    assert_array_almost_equal(z0.flatten(), z1)

    # this will throw if solutions aren't set for each cell (ie in initialize())
    for _ in range(10):
        model.update()
    assert z0 is model.get_value_ptr("Porosity")

def test_get_value_pointer2():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z0 = model.get_value_ptr("Porosity")
    assert(z0[0] != 0.999)
    z0[0] = 0.999

    z1 = model.get_value_ptr("Porosity")
    assert(z1[0] == 0.999)

    assert_array_almost_equal(z0, z1)

    assert(z1[0] != 0.888)
    z1[0] = 0.888
    assert(z0[0] == 0.888)

    assert_array_almost_equal(z0, z1)


def test_get_value_at_indices():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    temps = np.linspace(18.0, 30.0, model.get_grid_size(0))
    model.set_value("Temperature", temps)

    dest = np.empty(3, dtype=float)

    z0 = model.get_value_ptr("Temperature")
    z1 = model.get_value_at_indices("Temperature", dest, [0, 2, 4])

    assert_array_almost_equal(z0.take((0, 2, 4)), z1)

def test_value_size():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z = model.get_value_ptr("Temperature")
    assert model.get_grid_size(0) == z.size

def test_value_nbytes():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    z = model.get_value_ptr("Temperature")
    assert model.get_var_nbytes("Temperature") == z.nbytes

def test_all_pointable():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    vars = model.get_pointable_var_names()
    assert(isinstance(vars, tuple))

    for var in vars:
        print(var)
        z0 = model.get_value_ptr(var)
        assert z0 is not None

        assert(str(z0.dtype) == model.get_var_type(var))

        dest1 = np.empty_like(z0)
        z1 = model.get_value(var, dest1)

        assert z0 is not z1
        assert_array_almost_equal(z0.flatten(), z1)

        assert z0 is model.get_value_ptr(var)
