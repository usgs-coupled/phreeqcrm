import phreeqcrm

nxyz = 20
nthreads = 3

phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)

# Set properties
status = phreeqc_rm.SetComponentH2O(False)
phreeqc_rm.UseSolutionDensityVolume(False)

# Open files
status = phreeqc_rm.SetFilePrefix("SimpleAdvect_cpp")
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

# Set cells to print chemistry when print chemistry is turned on
print_chemistry_mask = [1] * nxyz
status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask)
nchem = phreeqc_rm.GetChemistryCellCount()

# --------------------------------------------------------------------------
# Set initial conditions
# --------------------------------------------------------------------------

# Set printing of chemistry file
status = phreeqc_rm.SetPrintChemistryOn(False, True, False) # workers, initial_phreeqc, utility

# Load database
status = phreeqc_rm.LoadDatabase("phreeqc.dat")

# Run file to define solutions and reactants for initial conditions, selected output
status = phreeqc_rm.RunFile(True, True, True, "advect.pqi")

# Clear contents of workers and utility
input = "DELETE; -all"
status = phreeqc_rm.RunString(True, False, True, input)

# Determine number of components to transport
ncomps = phreeqc_rm.FindComponents()
# Get component information
##const std::vector<std::string>& components = phreeqc_rm.GetComponents();
components = phreeqc_rm.GetComponents()
##print(type(components))
##print(components)

for comp in components:
    phreeqc_rm.OutputMessage(comp)
phreeqc_rm.OutputMessage("\n")

# Set array of initial conditions
ic1 = [-1] * nxyz * 7
f1 = [1] * nxyz * 7
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

##c = [] * nxyz * len(components)
##status = phreeqc_rm.GetConcentrations(c)
