import numpy as np
from numpy.testing import assert_array_equal

from phreeqcrm import bmi_phreeqcrm

import pytest

grid_id = 0
invalid_grid_id = 12345

def test_grid_var_names():
    model = bmi_phreeqcrm()
    model.initialize()

    names = model.get_input_var_names()
    assert names == ('Concentrations', 'DensityUser', 'FilePrefix', 'NthSelectedOutput', 'SaturationUser', 'Time', 'TimeStep', 'Porosity', 'Pressure', 'SelectedOutputOn', 'Temperature')

    names = model.get_output_var_names()
    assert names == ('ComponentCount', 'Components', 'Concentrations', 'DensityCalculated', 'ErrorString', 'FilePrefix', 'Gfw', 'GridCellCount', 'SaturationCalculated', 'SelectedOutput', 'SelectedOutputColumnCount', 'SelectedOutputCount', 'SelectedOutputRowCount', 'SolutionVolume', 'Time', 'TimeStep', 'CurrentSelectedOutputUserNumber', 'Porosity', 'Pressure', 'SelectedOutputOn', 'Temperature', 'Viscosity')

def test_grid_var_item_count():
    model = bmi_phreeqcrm()
    model.initialize()

    count = model.get_input_item_count()
    assert count == 11

    count = model.get_output_item_count()
    assert count == 22

def test_grid_var_units():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_var_units("Pressure") == "atm"
    ##assert model.get_var_units("Saturation") == "unitless"
    assert model.get_var_units("SolutionVolume") == "L"
    assert model.get_var_units("Viscosity") == "mPa s"

def test_grid_id():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_var_grid("Pressure") == grid_id

def test_grid_var_rank():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_grid_rank(grid_id) == 1

def test_grid_var_size():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_grid_size(grid_id) == 10

def test_grid_var_shape():
    model = bmi_phreeqcrm()
    model.initialize()
    ndim = model.get_grid_rank(0)
    shape = np.empty(ndim, dtype=np.int32)
    # assert_array_equal(model.get_grid_shape(grid_id, shape), (10))
    with pytest.raises(RuntimeError, match="Not Implemented"):
        model.get_grid_shape(grid_id, shape)

def test_grid_var_spacing():
    model = bmi_phreeqcrm()
    model.initialize()
    ndim = model.get_grid_rank(0)
    spacing = np.empty(ndim, dtype=np.int32)
    with pytest.raises(RuntimeError, match="Not Implemented"):
        model.get_grid_spacing(grid_id, spacing)

def test_grid_var_origin():
    model = bmi_phreeqcrm()
    model.initialize()
    ndim = model.get_grid_rank(0)
    origin = np.empty(ndim, dtype=np.int32)
    with pytest.raises(RuntimeError, match="Not Implemented"):
        model.get_grid_origin(grid_id, origin)

def test_grid_var_type():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_var_type("Pressure") == "float64"
    assert model.get_var_type("SaturationUser") == "float64"
    assert model.get_var_type("SolutionVolume") == "float64"
    assert model.get_var_type("Viscosity") == "float64"

def test_grid_type():
    model = bmi_phreeqcrm()
    model.initialize()
    assert model.get_grid_type(grid_id) == "points"