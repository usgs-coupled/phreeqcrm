import phreeqcrm
import numpy as np
import yamlphreeqcrm

# Note: If phreeqcrm is built with MPI support, mpi4py must be available;
# otherwise, the phreeqcrm module will fail to load.
if phreeqcrm.has_mpi():
    from mpi4py import MPI

class AdvectBMI(phreeqcrm.BMIPhreeqcRM):

    def __init__(self, nxyz):
        if phreeqcrm.has_mpi():
            super().__init__(nxyz, MPI.COMM_WORLD)
        else:
            super().__init__(nxyz)

        # Bogus conductivity field for Basic callback demonstration
        self.hydraulic_K = [0.0] * nxyz
        for i in range(nxyz):
            self.hydraulic_K[i] = i*2.0 

    def display_all_selected_output(self):

        nxyz = self.get_value_ptr("GridCellCount")[0]
        components = self.get_value_ptr("Components")
        density = self.get_value_ptr("DensityCalculated")
        volume = self.get_value_ptr("SolutionVolume")
        c = self.get_value_ptr("Concentrations")

        ncomps = len(components)

        n = self.get_scalar("SelectedOutputCount")
        for isel in range(n):
            self.set_scalar("NthSelectedOutput", isel)
            n_user = self.get_scalar("CurrentSelectedOutputUserNumber")

            print(f"Selected output sequence number: {isel}")
            print(f"Selected output user number:     {n_user}")

            # Get 2D array of selected output values
            cols = self.get_scalar("SelectedOutputColumnCount")
            rows = self.get_scalar("SelectedOutputRowCount")

            # Get headings
            headings = self.get_selected_output_headings()

            # Get selected output
            selected_out = self.get_selected_output()
            
            # Print results
            for i in range(rows//2):
                print("Cell number ", i)
                print(f"     Density:    {density[i]}")
                print(f"     Volume:     {volume[i]}")
                print(f"     Components: ")
                for j in range(ncomps):
                    print(f"{j}, {components[j]}: {c[j * nxyz + i]}")
                print(f"     Selected output: ")
                for j in range(cols):
                    print(f"{j}, {headings[j]}, {selected_out[j * nxyz + i]}")

    def display_results(self):
        # Use GetValue to extract exchange composition and pH
        # YAMLAddOutputVars can be used to select groups of
        # variables that are available through GetValue
        nxyz = self.get_scalar("GridCellCount")
        pH = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_ph", pH)
        Ca = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_total_molality_Ca", Ca)
        K = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_total_molality_K", K)
        Na = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_total_molality_Na", Na)
        Cl = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_total_molality_Cl", Cl)
        N = np.empty(nxyz, dtype=float)
        x = self.get_value("solution_total_molality_N", N)
        print("Cell         pH         Ca          K         Na         Cl        NO3")
        for i in range(int(nxyz / 2)):
            print(f"{i:4.0f}", f"{pH[i]:10.5f}", f"{Ca[i]:10.5f}", f"{K[i]:10.5f}", f"{Na[i]:10.5f}",f"{Cl[i]:10.5f}",f"{N[i]:10.5f}") 
        return
    
    def get_scalar(self, var_name):
        itemsize = self.get_var_itemsize(var_name)
        nbytes = self.get_var_nbytes(var_name)
        dim = nbytes // itemsize
    
        if dim != 1:
            raise ValueError(f"{var_name} is not a scalar")
    
        vtype = self.get_var_type(var_name)
        dest = np.empty(1, dtype=vtype)
        x = self.get_value(var_name, dest)
        return x[0]
    
    def set_scalar(self, var_name, value):
        itemsize = self.get_var_itemsize(var_name)
        nbytes = self.get_var_nbytes(var_name)
        dim = nbytes // itemsize
    
        if dim != 1:
            raise ValueError(f"{var_name} is not a scalar")
        
        vtype = self.get_var_type(var_name)
        dest = np.empty(1, dtype=vtype)
        dest[0] = value
        x = self.set_value(var_name, dest)

    def get_selected_output_headings(self):
        # @todo use bmi methods
        return self.GetSelectedOutputHeadings()
    
    def get_selected_output(self):
        # @todo use bmi methods
        return self.GetSelectedOutput()

def my_basic_callback(x1, x2, str, cookie):
    if str == "HYDRAULIC_K":
        phreeqcrm = cookie
        rm_cell_number = int(x1)
        if rm_cell_number >= 0 and rm_cell_number < phreeqcrm.GetChemistryCellCount():
            back = phreeqcrm.GetBackwardMapping() # dict of list
            return phreeqcrm.hydraulic_K[back[rm_cell_number][0]]
    return -999.9


def advect_bmi():

    # Based on PHREEQC Example 11

    # --------------------------------------------------------------------------
    # Create PhreeqcRM
    # --------------------------------------------------------------------------
    nxyz = 40
    bmi = AdvectBMI(nxyz)
    # Demonstrate add to Basic: Set a function for Basic CALLBACK
    bmi.set_basic_callback(my_basic_callback, bmi)
    if phreeqcrm.has_mpi():
        # Put workers in worker loop
        mpi_myself = bmi.GetMpiMyself()
        if (mpi_myself > 0):
            # cannot call bmi.initialize() on MPI workers
            bmi.MpiWorker()
            bmi.finalize()
            print("Worker success: ", mpi_myself)
            return

    # Initialize with YAML file
    file_prefix = "AdvectBmi"
    yaml_file = WriteYAMLFile(file_prefix, nxyz)
    status = bmi.initialize(yaml_file)

#ifdef USE_MPI
    # Optional callback for MPI
    #TODO status = do_something()   # only root is calling do_something here
#endif

    # print(bmi.get_pointable_var_names())

    components = bmi.get_value_ptr("Components")
    ncomps = bmi.get_value_ptr("ComponentCount")[0]
    nxyz = bmi.get_value_ptr("GridCellCount")[0]

    # print(components)
    # print(ncomps)
    # print(nxyz)

    time = bmi.get_value_ptr("Time")
    time_step = bmi.get_value_ptr("TimeStep")

    hydraulic_K = [0.0] * nxyz
    for i in range(nxyz):
        hydraulic_K[i] = i*2.0

    if phreeqcrm.has_openmp():
        nthreads = bmi.GetThreadCount()
        print(f"Number of threads:                                {nthreads}")

    nbytes = bmi.get_var_nbytes("FilePrefix")
    dest = np.full(1, " " * nbytes)
    prefix = bmi.get_value('FilePrefix', dest)[0]

    nchem = bmi.GetChemistryCellCount()
    print(nchem)

    print(f"Number of components for transport:               {ncomps}")

    # Get component information)
    gfw = bmi.get_value_ptr("Gfw")
    for i in range(ncomps):
        print(f"{components[i].rjust(10,' ')}  {gfw[i]}")
    print()

    # Get initial temperatures
    temperature = bmi.get_value_ptr("Temperature")
    # Get initial saturation
    sat = bmi.get_value_ptr("SaturationCalculated")
    # Get initial porosity
    por = bmi.get_value_ptr("Porosity")

    # print(temperature)
    # print(sat)
    # print(por)

    # Get initial volume
    volume = bmi.get_value_ptr("SolutionVolume")
    # Get initial concentrations
    #c = bmi.get_value_ptr("Concentrations")
    dest = np.empty(nxyz*ncomps, dtype=float)
    c = bmi.get_value("Concentrations", dest)

    # print("volume")
    # print(volume)
    # print("c")
    # print(c)

    # Set density, pressure, and temperature
    density = [1.0] * nxyz
    bmi.set_value("DensityUser", density)
    pressure = [2.0] * nxyz
    bmi.set_value("Pressure", pressure)
    temperature = [20.0] * nxyz
    bmi.set_value("Temperature", temperature)

    # --------------------------------------------------------------------------
    # Set boundary condition
    # --------------------------------------------------------------------------
    nbound = 1
    bc1 = [0]           # solution 0 from Initial IPhreeqc instance
    bc2 = [-1]          # no bc2 solution for mixing
    bc_f1 = [1.0]       # mixing fraction for bc1

    bc_conc = bmi.InitialPhreeqc2Concentrations(bc1)

    bmi.SetScreenOn(True)
    time[0] = 0.0
    time_step[0] = 86400.0

    # --------------------------------------------------------------------------
    # Transient loop
    # --------------------------------------------------------------------------
    dummy_int_scalar = np.empty((1,), dtype=int)
    selected_output = bmi.get_value_ptr("SelectedOutputOn")
    nsteps = 10
    for step in range(nsteps):
        print(f"Beginning transport calculation {time*bmi.GetTimeConversion()} days")
        print(f"          Time step             {time_step*bmi.GetTimeConversion()} days")

        advectionbmi(c, bc_conc, ncomps, nxyz, nbound)

        if (step == (nsteps - 1)):
            print_chemistry_on = 1
        else:
            print_chemistry_on = 0
        selected_output[0] = print_chemistry_on

        bmi.SetPrintChemistryOn(print_chemistry_on==1, False, False)  # workers, initial_phreeqc, utility
        time += time_step
        # bmi.set_value("Time", time)

        # Optionally, if values changed during transport
        ##bmi.set_value("Porosity", por)
        bmi.set_value("SaturationUser", sat)
        ##bmi.set_value("Temperature", temperature)
        ##bmi.set_value("Pressure", pressure)
        ##bmi.set_value("TimeStep", time_step)
        
        # Run cells with transported conditions
        print(f"Beginning reaction calculation  {time*bmi.GetTimeConversion()} days")
        bmi.set_value("Concentrations", c)
        bmi.update()
        bmi.get_value("Concentrations", c)

        # Print results at last time step
        if (step == (nsteps - 1)):
            bmi.display_results()

    # Clean up
    if phreeqcrm.has_mpi():
        bmi.MpiWorkerBreak()
    bmi.finalize()
    print("Done.")
    
def advectionbmi(c: np.ndarray, bc_conc: np.ndarray, ncomps: int, nxyz: int, dim: int):
    # Advect
    for i in range(nxyz - 1, 0, -1):
        for j in range(ncomps):
            c[j * nxyz + i] = c[j * nxyz + i - 1]              # component j
    
    # Cell zero gets boundary condition
    for j in range(ncomps):
        c[j * nxyz] = bc_conc[j * dim];                        # component j

def WriteYAMLFile(file_prefix, nxyz=40):
    """
    Write YAML file for this example. Based on PHREEQC Example 11.

    Args:
        file_prefix (str): The prefix for output files.
        nxyz (int, optional): The number of cells. Defaults to 40.

    Returns:
        str: The name of the YAML file.
    """	

    database_file = "phreeqc.dat"
    pqi_file = "advect.pqi"
    # Create YAMLPhreeqcRM document
    yrm = yamlphreeqcrm.YAMLPhreeqcRM()

    # Set GridCellCount
    yrm.YAMLSetGridCellCount(nxyz)
    
    if phreeqcrm.has_openmp():
        # Set ThreadCount
        yrm.YAMLThreadCount(3)
    
    # Set some properties
    yrm.YAMLSetErrorHandlerMode(1)
    yrm.YAMLSetComponentH2O(False)
    yrm.YAMLSetRebalanceFraction(0.5)
    yrm.YAMLSetRebalanceByCell(True)
    yrm.YAMLUseSolutionDensityVolume(False)
    yrm.YAMLSetPartitionUZSolids(False)
    # Open files
    yrm.YAMLSetFilePrefix(file_prefix)
    yrm.YAMLOpenFiles()
    # Set concentration units
    yrm.YAMLSetUnitsSolution(2)           # 1, mg/L; 2, mol/L; 3, kg/kgs
    yrm.YAMLSetUnitsPPassemblage(1)       # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    yrm.YAMLSetUnitsExchange(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    yrm.YAMLSetUnitsSurface(1)            # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    yrm.YAMLSetUnitsGasPhase(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    yrm.YAMLSetUnitsSSassemblage(1)       # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    yrm.YAMLSetUnitsKinetics(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    # Set conversion from seconds to user units (days) Only affects one print statement
    time_conversion = 1.0 / 86400.0
    yrm.YAMLSetTimeConversion(time_conversion)
    # Set representative volume
    rv = [1] * nxyz
    yrm.YAMLSetRepresentativeVolume(rv)
    # Set initial density
    density = [1.0] * nxyz                 # a float list
    yrm.YAMLSetDensityUser(density)
    # Set initial porosity
    por = [0.2] * nxyz                     # an integer list
    yrm.YAMLSetPorosity(por)
    # Set initial saturation
    #sat = [1] * nxyz
    sat = tuple(1.0 for i in range(nxyz))   # a float tuple
    yrm.YAMLSetSaturationUser(sat)   
    # Set cells to print chemistry when print chemistry is turned on
    print_chemistry_mask = [0] * nxyz
    for i in range(nxyz // 2):
        print_chemistry_mask[i] = 1
    yrm.YAMLSetPrintChemistryMask(print_chemistry_mask)  
    # Demonstation of mapping, two equivalent rows by symmetry
    # zero-based indexing
    grid2chem = [-1] * nxyz
    for i in range(nxyz // 2):
        grid2chem[i] = i 
        grid2chem[i + nxyz // 2] = i
    yrm.YAMLCreateMapping(grid2chem)
    # Set printing of chemistry file
    yrm.YAMLSetPrintChemistryOn(False, True, False) # workers, initial_phreeqc, utility
    # Load database
    yrm.YAMLLoadDatabase(database_file)    
    # Run file to define solutions and reactants for initial conditions, selected output
    workers = True             # Worker instances do the reaction calculations for transport
    initial_phreeqc = True     # InitialPhreeqc instance accumulates initial and boundary conditions
    utility = True             # Utility instance is available for processing
    yrm.YAMLRunFile(workers, initial_phreeqc, utility, pqi_file)
    # Clear contents of workers and utility
    initial_phreeqc = False
    input = "DELETE; -all"
    yrm.YAMLRunString(workers, initial_phreeqc, utility, input)
    yrm.YAMLAddOutputVars("AddOutputVars", "true")
    # Determine number of components to transport
    yrm.YAMLFindComponents()
    # set array of initial conditions
    ic1 = [-1]*(nxyz*7)            # an integer list
    #ic2 = [-1]*(nxyz*7)
    ic2 = np.full((nxyz*7,), -1)   # a numpy integer array
    #f1 = [1]*(nxyz*7)
    f1 = np.full((nxyz*7,), 1.0)   # a numpy double array
    for i in range(nxyz): 
        ic1[i] = 1;                # Solution 1
        ic1[nxyz + i] = -1;        # Equilibrium phases none
        ic1[2 * nxyz + i] = 1;     # Exchange 1
        ic1[3 * nxyz + i] = -1;    # Surface none
        ic1[4 * nxyz + i] = -1;    # Gas phase none
        ic1[5 * nxyz + i] = -1;    # Solid solutions none
        ic1[6 * nxyz + i] = -1;    # Kinetics none
    yrm.YAMLInitialPhreeqc2Module_mix(ic1, ic2, f1)    
    # No mixing is defined, so the following is equivalent
    #yrm.YAMLInitialPhreeqc2Module(ic1)

    # alternative for setting initial conditions
    # cell number in first argument (-1 indicates last solution, 40 in this case)
    # in advect.pqi and any reactants with the same number--
    # Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
    # will be written to cells 18 and 19 (0 based)
    module_cells = (18, 19)    # an integer tuple
    yrm.YAMLInitialPhreeqcCell2Module(-1, module_cells)
    # Initial equilibration of cells
    time_step = 0.0    # no kinetics
    yrm.YAMLSetTimeStep(time_step)
    time = 0.0
    yrm.YAMLSetTime(time)
    yrm.YAMLRunCells()
    time_step = 86400.0
    yrm.YAMLSetTimeStep(time_step)    
    # Write YAML file
    file_name = file_prefix + ".yaml"	
    yrm.WriteYAMLDoc(file_name)
    print("Wrote YAML file ", file_name, flush=True)
    return file_name

if __name__ == '__main__':
    advect_bmi()
