import os
import phreeqcrm as rm
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less
import pytest

def test_main():
    # for debugging
    print(f"PYTHONPATH={os.getenv('PYTHONPATH')}")

    bmi = rm.bmi_phreeqcrm()
    assert(not bmi._initialized)

    #
    # test methods that don't require initialize call (no dependency on VarManager)
    #

    bmi.finalize()
    assert(not bmi._initialized)

    component_name = bmi.get_component_name()
    assert(component_name == "BMI PhreeqcRM")

    # may want to make this throw to be consistent
    components = bmi.get_components()
    assert(isinstance(components, tuple))
    assert(len(components) == 0)

    #
    # test methods that require initialize call (depends on VarManager)
    #

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.get_grid_size(0)

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.get_input_var_names()

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.get_output_var_names()

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.get_value_ptr("Temperature")

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.update()

    with pytest.raises(RuntimeError, match="must call initialize first"):
        bmi.update_until(1.0)

    #
    # initialize
    #

    bmi.initialize("AdvectBMI_py.yaml")
    assert(bmi._initialized)

    nxyz = bmi.get_grid_size(0)
    assert(nxyz == 40)

    temperature = bmi.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(25.))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))

    # test WRITEABLE
    temperature[0] = 0.0
    temperature = bmi.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(0.0))
    assert(temperature[1] == pytest.approx(25.))
    assert(temperature[nxyz - 2] == pytest.approx(25.))
    assert(temperature[nxyz - 1] == pytest.approx(25.))

    # test set value as numpy.array
    temperature = np.full((nxyz,), 20.0)
    bmi.set_value("Temperature", temperature)


    temperature = bmi.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(20.))
    assert(temperature[1] == pytest.approx(20.))
    assert(temperature[nxyz - 2] == pytest.approx(20.))
    assert(temperature[nxyz - 1] == pytest.approx(20.))

    # test WRITEABLE
    temperature[0] = 1.1
    temperature = bmi.get_value_ptr("Temperature")
    assert(isinstance(temperature, np.ndarray))
    assert(len(temperature) == nxyz)
    assert(temperature[0] == pytest.approx(1.1))
    assert(temperature[1] == pytest.approx(20.))
    assert(temperature[nxyz - 2] == pytest.approx(20.))
    assert(temperature[nxyz - 1] == pytest.approx(20.))


    # make sure get_value_ptr can be called before get_input_var_names
    component_count = bmi.get_value_ptr("ComponentCount")
    assert(isinstance(component_count, np.ndarray))
    assert(len(component_count) == 1)
    assert(component_count.size == 1)
    assert(component_count[0] == 8)
    assert(component_count.dtype == 'int32')

    # make sure ComponentCount is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        component_count[0] = 25
    assert(component_count[0] == 8)

    input_vars = bmi.get_input_var_names()
    assert(isinstance(input_vars, tuple))
    assert(isinstance(input_vars[0], str))

    output_vars = bmi.get_output_var_names()
    assert(isinstance(output_vars, tuple))
    assert(isinstance(output_vars[0], str))

    # components = bmi.get_components()
    # assert(len(components) == 8)

    # saturation_user = bmi.get_value_ptr("SaturationUser")
    # assert(isinstance(saturation_user, np.ndarray))
    # assert(len(saturation_user) == nxyz)

    components = bmi.get_value_ptr("Components")
    assert(isinstance(components, np.ndarray))
    assert(components.dtype == '<U6')
    assert(len(components) == 8)
    assert(components.size == 8)
    assert(isinstance(components[0], np.str_))
    assert(components[0] == 'H')
    assert(components[1] == 'O')

    # make sure Components is readonly
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        components[0] = 'XX'
    assert(components[0] == 'H')
    assert(components[1] == 'O')

    porosity = bmi.get_value_ptr("Porosity")
    assert(isinstance(porosity, np.ndarray))
    assert(len(porosity) == nxyz)
    assert(porosity.size == nxyz)
    assert(porosity.dtype == 'float64')

    solution_volume = bmi.get_value_ptr("SolutionVolume")
    assert(isinstance(solution_volume, np.ndarray))
    assert(len(solution_volume) == nxyz)
    assert(solution_volume.size == nxyz)
    assert(solution_volume.dtype == 'float64')
    assert(solution_volume[0] == pytest.approx(0.2))
    assert(solution_volume[1] == pytest.approx(0.2))

    concentrations = bmi.get_value_ptr("Concentrations")
    assert(isinstance(concentrations, np.ndarray))
    assert(len(concentrations) == nxyz * len(components))
    assert(concentrations.size == nxyz * len(components))
    assert(concentrations.dtype == 'float64')

    # density_list = [1.0] * nxyz
    # bmi.set_value("Density", density_list)

    # density_tuple = tuple([1.0] * nxyz)
    # bmi.set_value("Density", density_tuple)

    # density_ndarray = np.full((nxyz,), 1.0)
    # bmi.set_value("Density", density_ndarray)



    '''
    # Set properties
    status = phreeqc_rm.SetComponentH2O(False)
    phreeqc_rm.UseSolutionDensityVolume(False)

    # Open files
    status = phreeqc_rm.SetFilePrefix("SimpleAdvect_py")
    phreeqc_rm.OpenFiles()

    # Set concentration units
    status = phreeqc_rm.SetUnitsSolution(2)           # 1, mg/L; 2, mol/L; 3, kg/kgs
    status = phreeqc_rm.SetUnitsExchange(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

    # Set conversion from seconds to user units (days)
    time_conversion = 1.0 / 86400
    status = phreeqc_rm.SetTimeConversion(time_conversion)

    # Set initial porosity
    por = [0.2] * nxyz
    status = phreeqc_rm.SetPorosity(por)

    ################################################################################################
    # BMI
    ################################################################################################

    # Set properties
    bmi.set_value("ComponentH2O", False)                    # YAML "SetComponentH2O"
    bmi.set_value("UseSolutionDensityVolume", False)        # YAML "UseSolutionDensityVolume"

    # Open files
    bmi.set_value("FilePrefix", "SimpleAdvect_py")          # YAML "SetFilePrefix"
    ????

    # Set concentration units
    bmi.set_value("UnitsSolution", 2)                       # YAML "SetUnitsSolution"              # 1, mg/L; 2, mol/L; 3, kg/kgs
    bmi.set_value("UnitsExchange", 1)                       # YAML "SetUnitsExchange"              # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

    # Set conversion from seconds to user units (days)
    time_conversion = 1.0 / 86400
    ???status = phreeqc_rm.SetTimeConversion(time_conversion)

    # Set initial porosity
    por = [0.2] * nxyz
    status = phreeqc_rm.SetPorosity(por)
    por = np.full((nxyz,), 0.2)
    bmi.set_value("Porosity", por)

    '''

    # Set properties
    # c_h2o = np.full((1,), False)
    # bmi.set_value("ComponentH2O", c_h2o)                    # YAML "SetComponentH2O"

    # use_sol_dens_vol = np.full((1,), False)
    # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"

    file_prefix = np.full((1,), "prefix")
    bmi.set_value("FilePrefix", file_prefix)

    # file_prefix = np.full((nxyz,), "nxyz_prefix")
    # with pytest.raises(RuntimeError, match="dimension error in set_value: FilePrefix"):
    #     bmi.set_value("FilePrefix", file_prefix)

    # use_sol_dens_vol = np.full((1,), False)
    # bmi.set_value("UseSolutionDensityVolume", use_sol_dens_vol)        # YAML "UseSolutionDensityVolume"


    bmi.finalize()
    assert(not bmi._initialized)
