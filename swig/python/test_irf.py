from io import StringIO

import numpy as np
import yaml
from numpy.testing import assert_almost_equal, assert_array_equal, assert_array_less
##from six.moves import range

from phreeqcrm import BMIPhreeqcRM

from constants import FilePaths

def test_component_name():
    model = BMIPhreeqcRM()

    name = model.get_component_name()
    assert name[:13] == "BMI PhreeqcRM"
    #assert model.get_component_name() is name  @todo

def test_start_time():
    model = BMIPhreeqcRM()
    model.initialize()

    assert_almost_equal(model.get_start_time(), 0.0)


def test_end_time():
    model = BMIPhreeqcRM()
    model.initialize()

    # assert_almost_equal(model.get_end_time(), np.finfo("d").max) @todo

def test_initialize_defaults():
    model = BMIPhreeqcRM()
    model.initialize()

    # assert_almost_equal(model.get_current_time(), 0.0)
    # z0 = model.get_value_ptr("plate_surface__temperature")
    # assert_array_less(z0, 1.0)
    # assert_array_less(0.0, z0)


def test_initialize_from_file_like():
    config = StringIO(yaml.dump({"shape": [7, 5]}))
    model = BMIPhreeqcRM()
    # model.initialize(config)

    # ndim = model.get_grid_rank(0)
    # shape = np.empty(ndim, dtype=np.int32)

    # assert_array_equal(model.get_grid_shape(0, shape), (7, 5))


def test_initialize_from_file():
    import os
    import tempfile

    import yaml

    with tempfile.NamedTemporaryFile("w", delete=False) as fp:
        fp.write(yaml.dump({"shape": [7, 5]}))
        name = fp.name

    print(name)

    model = BMIPhreeqcRM()
    # model.initialize(name)

    # os.remove(name)

    # ndim = model.get_grid_rank(0)
    # shape = np.empty(ndim, dtype=np.int32)

    # assert_array_equal(model.get_grid_shape(0, shape), (7, 5))

def test_update():
    model = BMIPhreeqcRM()
    ## model.initialize() #  @todo
    model.initialize(FilePaths.YAML)

    for inc in range(10):
        model.update()
        assert_almost_equal(model.get_current_time(), (inc + 1) * model.get_time_step())

def test_update_until():
    model = BMIPhreeqcRM()
    model.initialize()

    model.update_until(10.1)
    assert_almost_equal(model.get_current_time(), 10.1)

def test_finalize():
    model = BMIPhreeqcRM()
    model.initialize()
    model.update()
    model.finalize()