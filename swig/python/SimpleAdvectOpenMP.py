import phreeqcrm
import sys
import numpy as np

"""
    Functions that accept a constant vector reference
    (const std::vector<double>&) should be callable
    with any of the three Python data types, namely
    np.ndarray, list, and tuple.

    For example, the porosity can be initialized:

    # numpy.ndarray
    porosity = np.full((nxyz,), 0.2)

    # list (using list repetition)
    porosity = [0.2] * nxyz

    # list (using list comprehension)
    porosity = [0.2 for i in range(nxyz)]

    # tuple (using a generator expression)
    porosity = tuple(0.2 for i in range(nxyz))

    # And used:
    nxyz = 20
    nthreads = 3
    phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)
    status = phreeqc_rm.SetPorosity(porosity)

    # likewise for integers (const std::vector<int>&)

    Currently functions that take non-const vectors (ie OUT or INOUT)
    must use the SWIG wrapped vector class like this:

    c_dbl_vect = phreeqcrm.DoubleVector(nxyz * len(components))
    status = phreeqc_rm.GetConcentrations(c_dbl_vect)
"""

def SimpleAdvect():
    nxyz = 20
    nthreads = 3

    phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)

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
    #por = [0.2] * nxyz
    por = np.full((nxyz), 0.2)
    status = phreeqc_rm.SetPorosity(por)

    # Set cells to print chemistry when print chemistry is turned on
    #print_chemistry_mask = [1] * nxyz
    print_chemistry_mask = np.full((nxyz), 1)
    status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask)
    nchem = phreeqc_rm.GetChemistryCellCount()

    # --------------------------------------------------------------------------
    # Set initial conditions
    # --------------------------------------------------------------------------

    # Set printing of chemistry file
    status = phreeqc_rm.SetPrintChemistryOn(False, True, False)  # workers, initial_phreeqc, utility

    # Load database
    status = phreeqc_rm.LoadDatabase("phreeqc.dat")

    # Run file to define solutions and reactants for initial conditions, selected output
    status = phreeqc_rm.RunFile(True, True, True, "advect.pqi")

    # Clear contents of workers and utility
    input = "DELETE; -all"
    status = phreeqc_rm.RunString(True, False, True, input)

    # Determine number of components to transport
    ncomps = phreeqc_rm.FindComponents()

    # Get component information (as a tuple)
    components = phreeqc_rm.GetComponents()

    for comp in components:
        phreeqc_rm.OutputMessage(comp)
    phreeqc_rm.OutputMessage("\n")

    # Set array of initial conditions
    if 'numpy' in sys.modules:
        # this may require numpy to be linked in
        ic1 = np.full((nxyz * 7,), -1)
    else:
        ic1 = [-1] * nxyz * 7
    ic1 = [-1] * nxyz * 7     # Need to fix with numpy
    for i in range(nxyz):
        ic1[i]            =  1  # Solution 1
        ic1[nxyz + i]     = -1  # Equilibrium phases none
        ic1[2 * nxyz + i] =  1  # Exchange 1
        ic1[3 * nxyz + i] = -1  # Surface none
        ic1[4 * nxyz + i] = -1  # Gas phase none
        ic1[5 * nxyz + i] = -1  # Solid solutions none
        ic1[6 * nxyz + i] = -1  # Kinetics none
    status = phreeqc_rm.InitialPhreeqc2Module(ic1)

    # Initial equilibration of cells
    time = 0.0
    time_step = 0.0
    status = phreeqc_rm.SetTime(time)
    status = phreeqc_rm.SetTimeStep(time_step)
    status = phreeqc_rm.RunCells()

    # for now use std::vector<double> wrapper for [inout] arrays
    #c_dbl_vect = phreeqcrm.DoubleVector(nxyz * len(components), 0.0)
    #status = phreeqc_rm.GetConcentrations(c_dbl_vect)
    c_dbl_vect = phreeqc_rm.GetConcentrations()

    # --------------------------------------------------------------------------
    # Set boundary condition
    # --------------------------------------------------------------------------

    # for now use std::vector<double> wrapper for [inout] arrays
    #bc_conc_dbl_vect = phreeqcrm.DoubleVector()
    nbound = 1
    #bc1 = [0] * nbound                            # solution 0 from Initial IPhreeqc instance
    bc1 = np.full((1), 0)
    #status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc_dbl_vect, bc1)
    bc_conc_dbl_vect = phreeqc_rm.InitialPhreeqc2Concentrations(bc1)

    # --------------------------------------------------------------------------
    # Transient loop
    # --------------------------------------------------------------------------
    status = phreeqc_rm.SetTemperature([20.0] * nxyz)
    status = phreeqc_rm.SetPressure([2.0] * nxyz)

    time_step = 86400.0
    status = phreeqc_rm.SetTimeStep(time_step)

    nsteps = 10
    phreeqc_rm.SetScreenOn(True)
    for steps in range(nsteps):
        # Transport calculation here

        message = 'Beginning transport calculation              {} days\n'.format(phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion())
        phreeqc_rm.LogMessage(message)
        phreeqc_rm.ScreenMessage(message)
        
        message = '          Time step                          {} days\n'.format(phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion())
        phreeqc_rm.LogMessage(message)
        phreeqc_rm.ScreenMessage(message)

        simpleadvection(c_dbl_vect, bc_conc_dbl_vect, ncomps, nxyz, nbound)

        # Transfer data to PhreeqcRM for reactions
        print_selected_output_on = (steps == nsteps - 1)
        print_chemistry_on = (steps == nsteps - 1)
        status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on)
        status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, False, False)  # workers, initial_phreeqc, utility
        status = phreeqc_rm.SetConcentrations(c_dbl_vect)         # Transported concentrations
        time += time_step
        status = phreeqc_rm.SetTime(time)

        # Run cells with transported conditions
        message = 'Beginning reaction calculation               {} days\n'.format(time * phreeqc_rm.GetTimeConversion())
        phreeqc_rm.LogMessage(message)
        phreeqc_rm.ScreenMessage(message)
        status = phreeqc_rm.RunCells()

        # Transfer data from PhreeqcRM for transport
        #status = phreeqc_rm.GetConcentrations(c_dbl_vect)
        c_dbl_vect = phreeqc_rm.GetConcentrations()
        
    # Clean up
    status = phreeqc_rm.CloseFiles()
    status = phreeqc_rm.MpiWorkerBreak()
def simpleadvection(c, bc_conc, ncomps, nxyz, dim):
    """
    TODO
    """
    for i in range(nxyz - 1, 0, -1):
        for j in range(ncomps):
            c[j * nxyz + i] = c[j * nxyz + i - 1]              # component j
    
    # Cell zero gets boundary condition
    for j in range(ncomps):
        c[j * nxyz] = bc_conc[j * dim];                        # component j


if __name__ == '__main__':
    SimpleAdvect()