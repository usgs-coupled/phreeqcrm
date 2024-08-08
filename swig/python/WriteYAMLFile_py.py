import numpy as np
import yamlphreeqcrm

from constants import FilePaths

def WriteYAMLFile_py():
	# Create YAMLPhreeqcRM document
	yrm = yamlphreeqcrm.YAMLPhreeqcRM()
	# Number of cells
	nxyz = 40;
	# Set GridCellCount
	yrm.YAMLSetGridCellCount(nxyz)
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
	yrm.YAMLSetFilePrefix("AdvectBMI_py");
	yrm.YAMLOpenFiles();
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
	yrm.YAMLLoadDatabase("phreeqc.dat")    
    # Run file to define solutions and reactants for initial conditions, selected output
	workers = True             # Worker instances do the reaction calculations for transport
	initial_phreeqc = True     # InitialPhreeqc instance accumulates initial and boundary conditions
	utility = True             # Utility instance is available for processing
	yrm.YAMLRunFile(workers, initial_phreeqc, utility, "advect.pqi")
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
	yrm.WriteYAMLDoc(FilePaths.YAML)
	print("Done.")

if __name__ == '__main__':
    WriteYAMLFile_py()
            

