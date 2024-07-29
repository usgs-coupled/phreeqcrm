import os
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less
import pytest

from phreeqcrm import BMIPhreeqcRM, IRM_OK, State

from constants import FilePaths

ERROR_GET_VALUE_PTR_NOT_SUPPORTED = "get_value_ptr not supported for this variable."
ERROR_SET_VALUE_NOT_SUPPORTED     = "set_value not supported for this variable."

def test_main():
    # for debugging
    print(f"PYTHONPATH={os.getenv('PYTHONPATH')}")

def test_prereqs():
    # for debugging
    cwd = os.getcwd()
    print(f"Current working directory: {cwd}")

    yaml = FilePaths.YAML

    # Check if src file exists
    assert os.path.exists(yaml), f"{yaml} does not exist"

    # Check if src file size is greater than 0 bytes
    assert os.path.getsize(yaml) > 0, f"{yaml} is empty"

    database = FilePaths.DATABASE

    # Check if src file exists
    assert os.path.exists(database), f"{database} does not exist"

    # Check if src file size is greater than 0 bytes
    assert os.path.getsize(database) > 0, f"{database} is empty"

    # database must be in current working directory
    database_dest = os.path.join(cwd, os.path.basename(database))
    assert os.path.exists(database_dest), f"{database_dest} does not exist"

    pqi = FilePaths.PQI
    print(f"pqi: {pqi}")

    # Check if src file exists
    assert os.path.exists(pqi), f"{pqi} does not exist"

    # Check if file size is greater than 0 bytes
    assert os.path.getsize(pqi) > 0, f"{pqi} is empty"

    # pqi must be in current working directory
    pqi_dest = os.path.join(cwd, os.path.basename(pqi))
    assert os.path.exists(pqi_dest), f"{pqi_dest} does not exist"

def test_dtor():
    model = BMIPhreeqcRM()
    del model

def test_initialized():
    model = BMIPhreeqcRM()
    assert model._state == State.UNINITIALIZED

def test_finalize_not_initialized():
    model = BMIPhreeqcRM()

    assert model._state == State.UNINITIALIZED
    model.finalize()
    assert model._state == State.UNINITIALIZED

def test_get_components_is_tuple():
    model = BMIPhreeqcRM()

    # may want to make this throw to be consistent
    components = model.GetComponents()
    assert(isinstance(components, tuple))
    assert(len(components) == 0)

def test_get_grid_size_AdvectBMI():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    nxyz = model.get_grid_size(0)
    assert(nxyz == 40)

def test_get_value_ptr_is_ndarray():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)
    nxyz = model.get_grid_size(0)

    temperature = model.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(25.))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))

def test_get_temperature_ptr_is_writeable():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)
    nxyz = model.get_grid_size(0)

    # test WRITEABLE
    temperature = model.get_value_ptr("Temperature")
    assert(temperature.flags.writeable)

    # write
    assert(temperature[0] == pytest.approx(25.))
    temperature[0] = 0.0

    temperature = model.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(0.0))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 1] == pytest.approx(25.))

def test_get_value_ptr_ComponentCount_is_readonly():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # make sure get_value_ptr can be called before get_input_var_names
    component_count = model.get_value_ptr("ComponentCount")
    assert(isinstance(component_count, np.ndarray))
    assert(len(component_count) == 1)
    assert(component_count.size == 1)
    assert(component_count[0] == 8)
    assert(component_count.dtype == 'int32')
    assert(component_count.flags.writeable == False)

    # make sure ComponentCount is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        component_count[0] = 25
    assert(component_count[0] == 8)

def test_get_input_var_names_is_tuple():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    input_vars = model.get_input_var_names()
    assert(isinstance(input_vars, tuple))
    assert(isinstance(input_vars[0], str))

def test_get_output_var_names_is_tuple():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    output_vars = model.get_output_var_names()
    assert(isinstance(output_vars, tuple))
    assert(isinstance(output_vars[0], str))

def test_GetComponents_is_tuple():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    components = model.GetComponents()
    assert(isinstance(components, tuple))
    assert(len(components) == 8)

def test_get_value_ptr_Components_is_ndarray():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    components = model.get_value_ptr("Components")
    assert(isinstance(components, np.ndarray))
    assert(components.dtype == '<U6')
    assert(len(components) == 8)
    assert(components.size == 8)
    assert(isinstance(components[0], np.str_))
    assert(components[0] == 'H')
    assert(components[1] == 'O')

def test_get_value_ptr_Components_is_readonly():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    components = model.get_value_ptr("Components")

    # make sure Components is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        components[0] = 'XX'
    assert(components[0] == 'H')
    assert(components[1] == 'O')

def test_get_value_ptr_Porosity():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)
    nxyz = model.get_grid_size(0)

    porosity = model.get_value_ptr("Porosity")
    assert(isinstance(porosity, np.ndarray))
    assert(len(porosity) == nxyz)
    assert(porosity.size == nxyz)
    assert(porosity.dtype == 'float64')

def test_get_value_ptr_SolutionVolume():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)
    nxyz = model.get_grid_size(0)

    solution_volume = model.get_value_ptr("SolutionVolume")
    assert(isinstance(solution_volume, np.ndarray))
    assert(len(solution_volume) == nxyz)
    assert(solution_volume.size == nxyz)
    assert(solution_volume.dtype == 'float64')
    assert(solution_volume[0] == pytest.approx(0.2))
    assert(solution_volume[1] == pytest.approx(0.2))

def test_get_value_ptr_Concentrations():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)
    nxyz = model.get_grid_size(0)

    components = model.get_value_ptr("Components")

    concentrations = model.get_value_ptr("Concentrations")
    assert(isinstance(concentrations, np.ndarray))
    assert(len(concentrations) == nxyz * len(components))
    assert(concentrations.size == nxyz * len(components))
    assert(concentrations.dtype == 'float64')

#     # # density_list = [1.0] * nxyz
#     # # bmi.set_value("Density", density_list)

#     # # density_tuple = tuple([1.0] * nxyz)
#     # # bmi.set_value("Density", density_tuple)

#     # # density_ndarray = np.full((nxyz,), 1.0)
#     # # bmi.set_value("Density", density_ndarray)

#     # '''
#     # # Set properties
#     # status = phreeqc_rm.SetComponentH2O(False)
#     # phreeqc_rm.UseSolutionDensityVolume(False)

#     # # Open files
#     # status = phreeqc_rm.SetFilePrefix("SimpleAdvect_py")
#     # phreeqc_rm.OpenFiles()

#     # # Set concentration units
#     # status = phreeqc_rm.SetUnitsSolution(2)           # 1, mg/L; 2, mol/L; 3, kg/kgs
#     # status = phreeqc_rm.SetUnitsExchange(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

#     # # Set conversion from seconds to user units (days)
#     # time_conversion = 1.0 / 86400
#     # status = phreeqc_rm.SetTimeConversion(time_conversion)

#     # # Set initial porosity
#     # por = [0.2] * nxyz
#     # status = phreeqc_rm.SetPorosity(por)

#     # ################################################################################################
#     # # BMI
#     # ################################################################################################

#     # # Set properties
#     # bmi.set_value("ComponentH2O", False)                    # YAML "SetComponentH2O"
#     # bmi.set_value("UseSolutionDensityVolume", False)        # YAML "UseSolutionDensityVolume"

#     # # Open files
#     # bmi.set_value("FilePrefix", "SimpleAdvect_py")          # YAML "SetFilePrefix"
#     # ????

#     # # Set concentration units
#     # bmi.set_value("UnitsSolution", 2)                       # YAML "SetUnitsSolution"              # 1, mg/L; 2, mol/L; 3, kg/kgs
#     # bmi.set_value("UnitsExchange", 1)                       # YAML "SetUnitsExchange"              # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

#     # # Set conversion from seconds to user units (days)
#     # time_conversion = 1.0 / 86400
#     # ???status = phreeqc_rm.SetTimeConversion(time_conversion)

#     # # Set initial porosity
#     # por = [0.2] * nxyz
#     # status = phreeqc_rm.SetPorosity(por)
#     # por = np.full((nxyz,), 0.2)
#     # bmi.set_value("Porosity", por)

#     # '''

#     # # Set properties
#     # # c_h2o = np.full((1,), False)
#     # # bmi.set_value("ComponentH2O", c_h2o)                    # YAML "SetComponentH2O"

#     # # use_sol_dens_vol = np.full((1,), False)
#     # # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"

#     # file_prefix = np.full((1,), "prefix")
#     # bmi.set_value("FilePrefix", file_prefix)

#     # # file_prefix = np.full((nxyz,), "nxyz_prefix")
#     # # with pytest.raises(RuntimeError, match="dimension error in set_value: FilePrefix"):
#     # #     bmi.set_value("FilePrefix", file_prefix)

#     # # use_sol_dens_vol = np.full((1,), False)
#     # # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"

def test_get_value_ptr_ndarray_str():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    components = model.get_value_ptr("Components")
    comps = np.empty_like(components)

    comps = model.get_value("Components", comps)
    assert(comps[0] == 'H')
    assert(comps[1] == 'O')
    assert(comps[2] == 'Charge')
    assert(comps[3] == 'Ca')
    assert(comps[4] == 'Cl')
    assert(comps[5] == 'K')
    assert(comps[6] == 'N')
    assert(comps[7] == 'Na')

def test_SetComponentH2O():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    status = model.SetComponentH2O(True)
    assert(status == IRM_OK)

def test_get_value_FilePrefix():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # FilePrefix doesn't support get_value_ptr
    expected = "AdvectBMI_py"
    spaces = " " * len(expected)
    fileprefix = np.array([spaces])
    fileprefix = model.get_value("FilePrefix", fileprefix)

    assert(len(fileprefix[0]) == 12)
    assert(fileprefix[0] == expected)

def test_get_value_FilePrefix_fail():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # FilePrefix doesn't support get_value_ptr
    sz = 4
    expected = "AdvectBMI_py"
    spaces = " " * sz
    fileprefix = np.array([spaces])

    fileprefix = model.get_value("FilePrefix", fileprefix)

    # Current version truncates the string rather than
    # throwing "buffer too small"

    # with pytest.raises(RuntimeError, match="buffer too small"):
    #     fileprefix = model.get_value("FilePrefix", fileprefix)

    assert(len(fileprefix[0]) == 4)
    assert(fileprefix[0] == expected[:sz])

def test_get_value_FilePrefix_big_buffer():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # FilePrefix doesn't support get_value_ptr
    expected = "AdvectBMI_py"
    spaces = " " * 80
    fileprefix = np.array([spaces])
    fileprefix = model.get_value("FilePrefix", fileprefix)

    assert(fileprefix[0] == expected)


def test_get_value_ComponentCount():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # ComponentCount supports get_value_ptr @todo
    ##component_count = np.full((1,), 0)
    component_count = np.empty((1,), dtype=int)
    component_count = model.get_value("ComponentCount", component_count)

    assert(component_count[0] == 8)

def test_get_value_ptr_ComponentCount():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # ComponentCount supports get_value_ptr @todo
    component_count = np.empty((1,), dtype=int)
    component_count = model.get_value_ptr("ComponentCount")

    assert(component_count[0] == 8)

def test_get_value_ptr_failure():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # "Components" is currently excluded
    should_raise = ("CurrentSelectedOutputUserNumber", "DensityUser", "ErrorString", "FilePrefix", "NthSelectedOutput", "SaturationUser", "SelectedOutput", "SelectedOutputColumnCount", "SelectedOutputCount", "SelectedOutputHeadings", "SelectedOutputRowCount")

    for item in should_raise:
        print(f"testing {item}")
        with pytest.raises(RuntimeError, match=ERROR_GET_VALUE_PTR_NOT_SUPPORTED):
            arr = model.get_value_ptr(item)

def test_read_only_vars():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # "Components" is currently excluded
    readonly_list = model.get_readonly_var_names()

    for item in readonly_list:
        print(f"testing {item}")

        itemsize = model.get_var_itemsize(item)
        nbytes = model.get_var_nbytes(item)
        if item == 'ErrorString':
            assert itemsize == 0
            assert nbytes == 0
            continue
        dim = nbytes // itemsize
        vtype = model.get_var_type(item)
        dest = np.empty(dim, dtype=vtype)

        with pytest.raises(RuntimeError, match=ERROR_SET_VALUE_NOT_SUPPORTED):
            model.set_value(item, dest)

def test_get_value_Time():
    model = BMIPhreeqcRM()
    model.initialize(FilePaths.YAML)

    # Time supports get_value_ptr @todo
    time = np.empty((1,), dtype=float)
    time = model.get_value("Time", time)

    assert(time[0] == 0.0)