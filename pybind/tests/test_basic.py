import os
##import phreeqcrm as rm
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less
import pytest

from phreeqcrm import bmi_phreeqcrm, IRM_RESULT

def test_main():
    # for debugging
    print(f"PYTHONPATH={os.getenv('PYTHONPATH')}")

def test_dtor():
    model = bmi_phreeqcrm()
    del model

def test_initialized():
    model = bmi_phreeqcrm()
    assert(not model._initialized)

def test_finalize_not_initialized():
    model = bmi_phreeqcrm()

    model.finalize()
    assert(not model._initialized)

def test_get_components_is_tuple():
    model = bmi_phreeqcrm()

    # may want to make this throw to be consistent
    components = model.get_components()
    assert(isinstance(components, tuple))
    assert(len(components) == 0)

def test_get_grid_size_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.get_grid_size(0)

def test_get_input_var_names_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.get_input_var_names()

def test_get_output_var_names_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.get_output_var_names()

def test_get_value_ptr_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.get_value_ptr("Temperature")

def test_update_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.update()

def test_update_until_raises_uninitialized():
    model = bmi_phreeqcrm()
    with pytest.raises(RuntimeError, match="must call initialize first"):
        model.update_until(1.0)

def test_initialize_AdvectBMI():
    model = bmi_phreeqcrm()

    assert(not model._initialized)
    model.initialize("AdvectBMI_py.yaml")
    assert(model._initialized)

def test_get_grid_size_AdvectBMI():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    nxyz = model.get_grid_size(0)
    assert(nxyz == 40)

def test_get_value_ptr_is_ndarray():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")
    nxyz = model.get_grid_size(0)

    temperature = model.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(25.))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))

def test_get_temperature_ptr_is_writeable():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")
    nxyz = model.get_grid_size(0)

    # test WRITEABLE
    temperature = model.get_value_ptr("Temperature")
    assert(temperature.flags.writeable)

    # write
    temperature[0] = 0.0

    temperature = model.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(0.0))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 1] == pytest.approx(25.))

def test_get_ComponentCount_ptr_is_readonly():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # make sure get_value_ptr can be called before get_input_var_names
    component_count = model.get_value_ptr("ComponentCount")
    assert(isinstance(component_count, np.ndarray))
    assert(len(component_count) == 1)
    assert(component_count.size == 1)
    assert(component_count[0] == 8)
    assert(component_count.dtype == 'int32')

    # make sure ComponentCount is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        component_count[0] = 25
    assert(component_count[0] == 8)

def test_get_input_var_names_is_tuple():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    input_vars = model.get_input_var_names()
    assert(isinstance(input_vars, tuple))
    assert(isinstance(input_vars[0], str))

def test_get_output_var_names_is_tuple():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    output_vars = model.get_output_var_names()
    assert(isinstance(output_vars, tuple))
    assert(isinstance(output_vars[0], str))

def test_get_components_is_tuple():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    components = model.get_components()
    assert(isinstance(components, tuple))
    assert(len(components) == 8)

def test_get_Components_ptr_is_ndarray():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    components = model.get_value_ptr("Components")
    assert(isinstance(components, np.ndarray))
    assert(components.dtype == '<U6')
    assert(len(components) == 8)
    assert(components.size == 8)
    assert(isinstance(components[0], np.str_))
    assert(components[0] == 'H')
    assert(components[1] == 'O')

def test_get_Components_ptr_is_readonly():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    components = model.get_value_ptr("Components")

    # make sure Components is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        components[0] = 'XX'
    assert(components[0] == 'H')
    assert(components[1] == 'O')

def test_get_Porosity_ptr():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")
    nxyz = model.get_grid_size(0)

    porosity = model.get_value_ptr("Porosity")
    assert(isinstance(porosity, np.ndarray))
    assert(len(porosity) == nxyz)
    assert(porosity.size == nxyz)
    assert(porosity.dtype == 'float64')

def test_get_SolutionVolume_ptr():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")
    nxyz = model.get_grid_size(0)

    solution_volume = model.get_value_ptr("SolutionVolume")
    assert(isinstance(solution_volume, np.ndarray))
    assert(len(solution_volume) == nxyz)
    assert(solution_volume.size == nxyz)
    assert(solution_volume.dtype == 'float64')
    assert(solution_volume[0] == pytest.approx(0.2))
    assert(solution_volume[1] == pytest.approx(0.2))

def test_get_Concentrations_ptr():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")
    nxyz = model.get_grid_size(0)

    components = model.get_value_ptr("Components")

    concentrations = model.get_value_ptr("Concentrations")
    assert(isinstance(concentrations, np.ndarray))
    assert(len(concentrations) == nxyz * len(components))
    assert(concentrations.size == nxyz * len(components))
    assert(concentrations.dtype == 'float64')

    # # density_list = [1.0] * nxyz
    # # bmi.set_value("Density", density_list)

    # # density_tuple = tuple([1.0] * nxyz)
    # # bmi.set_value("Density", density_tuple)

    # # density_ndarray = np.full((nxyz,), 1.0)
    # # bmi.set_value("Density", density_ndarray)

    # '''
    # # Set properties
    # status = phreeqc_rm.SetComponentH2O(False)
    # phreeqc_rm.UseSolutionDensityVolume(False)

    # # Open files
    # status = phreeqc_rm.SetFilePrefix("SimpleAdvect_py")
    # phreeqc_rm.OpenFiles()

    # # Set concentration units
    # status = phreeqc_rm.SetUnitsSolution(2)           # 1, mg/L; 2, mol/L; 3, kg/kgs
    # status = phreeqc_rm.SetUnitsExchange(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

    # # Set conversion from seconds to user units (days)
    # time_conversion = 1.0 / 86400
    # status = phreeqc_rm.SetTimeConversion(time_conversion)

    # # Set initial porosity
    # por = [0.2] * nxyz
    # status = phreeqc_rm.SetPorosity(por)

    # ################################################################################################
    # # BMI
    # ################################################################################################

    # # Set properties
    # bmi.set_value("ComponentH2O", False)                    # YAML "SetComponentH2O"
    # bmi.set_value("UseSolutionDensityVolume", False)        # YAML "UseSolutionDensityVolume"

    # # Open files
    # bmi.set_value("FilePrefix", "SimpleAdvect_py")          # YAML "SetFilePrefix"
    # ????

    # # Set concentration units
    # bmi.set_value("UnitsSolution", 2)                       # YAML "SetUnitsSolution"              # 1, mg/L; 2, mol/L; 3, kg/kgs
    # bmi.set_value("UnitsExchange", 1)                       # YAML "SetUnitsExchange"              # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

    # # Set conversion from seconds to user units (days)
    # time_conversion = 1.0 / 86400
    # ???status = phreeqc_rm.SetTimeConversion(time_conversion)

    # # Set initial porosity
    # por = [0.2] * nxyz
    # status = phreeqc_rm.SetPorosity(por)
    # por = np.full((nxyz,), 0.2)
    # bmi.set_value("Porosity", por)

    # '''

    # # Set properties
    # # c_h2o = np.full((1,), False)
    # # bmi.set_value("ComponentH2O", c_h2o)                    # YAML "SetComponentH2O"

    # # use_sol_dens_vol = np.full((1,), False)
    # # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"

    # file_prefix = np.full((1,), "prefix")
    # bmi.set_value("FilePrefix", file_prefix)

    # # file_prefix = np.full((nxyz,), "nxyz_prefix")
    # # with pytest.raises(RuntimeError, match="dimension error in set_value: FilePrefix"):
    # #     bmi.set_value("FilePrefix", file_prefix)

    # # use_sol_dens_vol = np.full((1,), False)
    # # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"

def test_get_value_ndarray_str():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

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
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    status = model.SetComponentH2O(True)
    assert(status == IRM_RESULT.IRM_OK)

def test_get_value_FilePrefix():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # FilePrefix doesn't support get_value_ptr
    expected = "AdvectBMI_py"
    spaces = " " * len(expected)
    fileprefix = np.array([spaces])
    fileprefix = model.get_value("FilePrefix", fileprefix)

    assert(len(fileprefix[0]) == 12)
    assert(fileprefix[0] == expected)

def test_get_value_FilePrefix_fail():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # FilePrefix doesn't support get_value_ptr
    spaces = " " * 4
    fileprefix = np.array([spaces])

    with pytest.raises(RuntimeError, match="buffer too small"):
        fileprefix = model.get_value("FilePrefix", fileprefix)

def test_get_value_FilePrefix_big_buffer():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # FilePrefix doesn't support get_value_ptr
    expected = "AdvectBMI_py"
    spaces = " " * 80
    ##f = np.empty((1,), 80, dtype=char)
    fileprefix = np.array([spaces])
    fileprefix = model.get_value("FilePrefix", fileprefix)

    assert(fileprefix[0] == expected)


def test_get_value_ComponentCount():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # ComponentCount supports get_value_ptr @todo
    ##component_count = np.full((1,), 0)
    component_count = np.empty((1,), dtype=int)
    component_count = model.get_value("ComponentCount", component_count)

    assert(component_count[0] == 8)

def test_get_value_ptr_ComponentCount():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # ComponentCount supports get_value_ptr @todo
    ##component_count = np.full((1,), 0)
    component_count = np.empty((1,), dtype=int)
    component_count = model.get_value_ptr("ComponentCount")

    assert(component_count[0] == 8)

def test_get_value_ptr_failure():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # "Components" is currently excluded
    should_raise = ("CurrentSelectedOutputUserNumber", "DensityUser", "ErrorString", "FilePrefix", "NthSelectedOutput", "SaturationUser", "SelectedOutput", "SelectedOutputColumnCount", "SelectedOutputCount", "SelectedOutputHeadings", "SelectedOutputRowCount")

    for item in should_raise:
        print(f"testing {item}")
        with pytest.raises(RuntimeError, match="This variable does not support get_value_ptr"):
            arr = model.get_value_ptr(item)

def test_get_value_Time():
    model = bmi_phreeqcrm()
    model.initialize("AdvectBMI_py.yaml")

    # Time supports get_value_ptr @todo
    time = np.empty((1,), dtype=float)
    time = model.get_value("Time", time)

    assert(time[0] == 0.0)