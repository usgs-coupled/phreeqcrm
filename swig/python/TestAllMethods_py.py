
import numpy as np
import phreeqcrm
#ifdef USE_YAML
import yamlphreeqcrm
    #module mydata
    #  double precision, dimension(:), pointer :: K_ptr
    #  integer                                 :: rm_id
    #end module mydata

def testallmethods_py():

	yrm = yamlphreeqcrm.YAMLPhreeqcRM()
	yrm.YAMLSetGridCellCount(40)
	yrm.YAMLThreadCount(3)
	YAML_filename = "testallmethods_py.yaml"
	yrm.WriteYAMLDoc(YAML_filename)
	#
	# Use all BMIPhreeqcRM methods roughly in order of use
	#
	bmi=phreeqcrm.BMIPhreeqcRM()
	print(f"BMIPhreeqcRM {bmi}")
	#---------
	bmi.initialize(YAML_filename)   # void function
	bmi.InitializeYAML(YAML_filename)
	print(f"Initialize")
	#---------
	nxyz = bmi.GetGridCellCount()
	print(f"GetGridCellCount {type(nxyz)}, {nxyz}")
	dest = np.empty((1,), dtype=int)
	nxyz = bmi.get_value("GridCellCount", dest)
	print(f"get_value('GridCellCount') {type(nxyz)}, {nxyz[0]}")
	dest_ptr = bmi.get_value_ptr("GridCellCount")
	print(f"get_value_ptr('GridCellCount') {type(dest_ptr)}, {dest_ptr}")
	#---------
	x=bmi.GetThreadCount()
	print(f"GetThreadCount {type(x)}, {x}")
	#---------
	# Inactive cells or symmetry
	grid2chem = np.full((nxyz), -1)
	for i in range(nxyz[0]//2):
		grid2chem[i] = i
	x=bmi.CreateMapping(grid2chem)
	print(f"CreateMapping {type(x)}, {x}")
	#--------- 
	x=bmi.LoadDatabase("phreeqc.dat")
	print(f"LoadDatabase {type(x)}, {x}")
	#
	# Set properties
	#
	x=bmi.SetComponentH2O(False)
	print(f"SetComponentH2O {type(x)}, {x}")
	#---------
	x=bmi.SetSpeciesSaveOn(True)
	print(f"SetSpeciesSaveOn {type(x)}, {x}")
	#---------
	x=bmi.SetErrorOn(True)
	print(f"SetErrorOn {type(x)}, {x}")
	#---------
	x=bmi.SetErrorHandlerMode(1)
	print(f"SetErrorHandlerMode {type(x)}, {x}")
	#---------
	x=bmi.SetDumpFileName("TestAllMethods_py.dump")
	print(f"SetDumpFileName {type(x)}, {x}")
	#---------
	x=bmi.SetFilePrefix("TestAllMethods_py")
	print(f"SetFilePrefix {type(x)}, {x}")
	dest = np.full(1, "TestAllMethods_py")
	bmi.set_value("FilePrefix", dest)
	print(f"set_value('FilePrefix')")
	#---------
	x=bmi.OpenFiles()
	print(f"OpenFiles {type(x)}, {x}")
	#---------
	x=bmi.SetPartitionUZSolids(False)	
	print(f"SetPartitionUZSolids {type(x)}, {x}")
	#---------
	x=bmi.SetRebalanceByCell(True)
	print(f"SetRebalanceByCell {type(x)}, {x}")
	#---------
	x=bmi.SetRebalanceFraction(0.5)
	print(f"SetRebalanceFraction {type(x)}, {x}")
	#---------
	x=bmi.SetScreenOn(True)
	print(f"SetScreenOn {type(x)}, {x}")
	#---------
	x=bmi.SetSelectedOutputOn(True)
	print(f"SetSelectedOutputOn {type(x)}, {x}")
	src = np.full(1,1)
	bmi.set_value("SelectedOutputOn", src)
	print(f"set_value('SelectedOutputOn') {src[0]}")
	#---------
	x=bmi.SetUnitsExchange(1)
	print(f"SetUnitsExchange {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsGasPhase(1)
	print(f"SetUnitsGasPhase {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsKinetics(1)
	print(f"SetUnitsKinetics {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsPPassemblage(1)
	print(f"SetUnitsPPassemblage {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsSolution(2)
	print(f"SetUnitsSolution {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsSSassemblage(1)
	print(f"SetUnitsSSassemblage {type(x)}, {x}")
	#---------
	x=bmi.SetUnitsSurface(1)
	print(f"SetUnitsSurface {type(x)}, {x}")
	#---------
	bmi.UseSolutionDensityVolume(False)
	print(f"UseSolutionDensityVolume {type(x)}, {x}")
	#---------
	x=bmi.SetTimeConversion(1.0/86400.0)
	print(f"SetTimeConversion {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1.0)
	x=bmi.SetRepresentativeVolume(v)
	print(f"SetRepresentativeVolume {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.SetPrintChemistryMask(v)
	print(f"SetPrintChemistryMask {type(x)}, {x}")
	#---------
	x=bmi.SetPrintChemistryOn(False,True,False)
	print(f"SetPrintChemistryOn {type(x)}, {x}")
	#
	# Define reactants available for initial 
	# and boundary conditions in this file
	#
	x = bmi.RunFile(True, True, True, "all_reactants.pqi")
	print(f"RunFile {type(x)}, {x}")	
	#---------
	x=bmi.add_output_vars("AddOutputVars", "True")
	x=bmi.add_output_vars("SolutionProperties", "True")
	x=bmi.add_output_vars("SolutionTotalMolalities", "True")
	x=bmi.add_output_vars("ExchangeMolalities", "True")
	x=bmi.add_output_vars("SurfaceMolalities", "True")
	x=bmi.add_output_vars("EquilibriumPhases", "True")
	x=bmi.add_output_vars("Gases", "True")
	x=bmi.add_output_vars("KineticReactants", "True")
	x=bmi.add_output_vars("SolidSolutions", "True")
	x=bmi.add_output_vars("CalculateValues", "True")
	x=bmi.add_output_vars("SolutionActivities", "H+ Ca+2 Na+")
	x=bmi.add_output_vars("SolutionMolalities", "OH- Cl-")
	x=bmi.add_output_vars("SaturationIndices", "Calcite Dolomite")
	print(f"add_output_vars {type(x)}, {x}")		
	#---------
	ncomps = bmi.FindComponents()
	print(f"FindComponents {type(ncomps)}, {ncomps}")
	#
	# Methods up to this point are useful 
	# in a YAML initialization file
	# 
	# Lists of reactants found by FindComponents follow
	# 
	nchem=bmi.GetChemistryCellCount()
	print(f"GetChemistryCellCount {type(nchem)}, {nchem}")
	#---------
	ncomps=bmi.GetComponentCount()
	print(f"GetComponentCount {type(ncomps)}, {ncomps}")
	dest = np.empty((1,), dtype=int)
	ncomps = bmi.get_value("ComponentCount", dest)
	print(f"get_value('ComponentCount') {type(ncomps)}, {ncomps}")
	dest_ptr = bmi.get_value_ptr("ComponentCount")
	print(f"get_value_ptr('ComponentCount') {type(dest_ptr)}, {dest_ptr}")
	#---------
	x=bmi.GetComponents()	
	print(f"GetComponents {type(x)}, {x}")
	itemsize = bmi.get_var_itemsize("Components")
	nbytes = bmi.get_var_nbytes("Components")
	dim = nbytes // itemsize
	s = " " * itemsize
	dest = np.full(dim, s)
	x = bmi.get_value("Components", dest)
	print(f"get_value('Components', dest) {type(dest)}, {dest}")
	print(f"get_value('Components') {type(x)}, {x}")
	# Species info
	x=bmi.GetSpeciesCount()
	print(f"GetSpeciesCount {type(x)}, {x}")
	#---------
	x=bmi.GetSpeciesD25()
	print(f"GetSpeciesD25 {type(x)}, {x[0]}")
	#---------
	x=bmi.GetSpeciesNames()
	print(f"GetSpeciesNames {type(x)}, {x[0]}")
	#---------
	x=bmi.GetSpeciesStoichiometry()
	print(f"GetSpeciesStoichiometry {type(x)}, {x['AlOH+2']}")
	#---------
	x=bmi.GetSpeciesZ()
	print(f"GetSpeciesZ {type(x)}, {x[0]}")
	# Reactant lists
	x=bmi.GetEquilibriumPhases()
	print(f"GetEquilibriumPhases {type(x)}, {x}")
	#---------
	x=bmi.GetEquilibriumPhasesCount()
	print(f"GetEquilibriumPhasesCount {type(x)}, {x}")
	#---------
	x=bmi.GetExchangeNames()
	print(f"GetExchangeNames {type(x)}, {x}")
	#---------
	x=bmi.GetExchangeSpecies()
	print(f"GetExchangeSpecies {type(x)}, {x}")
	#---------
	x=bmi.GetExchangeSpeciesCount()
	print(f"GetExchangeSpeciesCount {type(x)}, {x}")
	#---------
	x=bmi.GetGasComponents()
	print(f"GetGasComponents {type(x)}, {x}")
	#---------
	x=bmi.GetGasComponentsCount()
	print(f"GetGasComponentsCount {type(x)}, {x}")
	#---------
	x=bmi.GetGfw()
	print(f"GetGfw {type(x)}, {x[0]}")
	dest = np.empty(ncomps[0], dtype=float)
	x=bmi.get_value("gfw", dest)
	print(f"get_value('gfw') {type(x)}, {x[0]}")
	x = bmi.get_value_ptr("Gfw")  
	print(f"get_value_ptr('Gfw') {type(x)}, {x[0]}")
	#---------
	x=bmi.GetKineticReactions()
	print(f"GetKineticReactions {type(x)}, {x}")
	#---------
	x=bmi.GetKineticReactionsCount()
	print(f"GetKineticReactionsCount {type(x)}, {x}")
	#---------
	x=bmi.GetSICount()
	print(f"GetSICount {type(x)}, {x}")
	#---------
	x=bmi.GetSINames()
	print(f"GetSINames {type(x)}, {x[0]}")
	#---------
	x=bmi.GetSolidSolutionComponents()
	print(f"GetSolidSolutionComponents {type(x)}, {x}")
	#---------
	x=bmi.GetSolidSolutionComponentsCount()
	print(f"GetSolidSolutionComponentsCount {type(x)}, {x}")
	#---------
	x=bmi.GetSolidSolutionNames()
	print(f"GetSolidSolutionNames {type(x)}, {x}")
	#---------
	x=bmi.GetSurfaceNames()
	print(f"GetSurfaceNames {type(x)}, {x[0]}")
	#---------
	x=bmi.GetSurfaceSpecies()
	print(f"GetSurfaceSpecies {type(x)}, {x[0]}")
	#---------
	x=bmi.GetSurfaceSpeciesCount()
	print(f"GetSurfaceSpeciesCount {type(x)}, {x}")
	#---------
	x=bmi.GetSurfaceTypes()
	print(f"GetSurfaceTypes {type(x)}, {x[0]}")
	#
	# Remove any reactants in workers 
	#
	input = "DELETE; -all"
	bmi.RunString(True, False, False, input)
	print(f"RunString")
	# 
	# Transfer initial conditions
	#
	v = np.full((nxyz), 1)
	x=bmi.InitialEquilibriumPhases2Module(v)
	print(f"InitialEquilibriumPhases2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialExchanges2Module(v)
	print(f"InitialExchanges2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialGasPhases2Module(v)
	print(f"InitialGasPhases2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialKinetics2Module(v)
	print(f"InitialKinetics2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialSolutions2Module(v)
	print(f"InitialSolutions2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialSolidSolutions2Module(v)
	print(f"InitialSolidSolutions2Module {type(x)}, {x}")
	#---------
	v = np.full((nxyz), 1)
	x=bmi.InitialSurfaces2Module(v)	
	print(f"InitialSurfaces2Module {type(x)}, {x}")
	#---------
	# Alternative A. to the previous seven methods
	v = np.full((nxyz*7), 1)
	x=bmi.InitialPhreeqc2Module(v)
	print(f"InitialPhreeqc2Module {type(x)}, {x}")
	#---------
	# Alternative B. to the previous seven methods
	v1 = np.full((nxyz*7), 1)
	v2 = np.full((nxyz*7), -1)
	f1 = np.full((nxyz*7), 1.0)
	x=bmi.InitialPhreeqc2Module_mix(v1, v2, f1)
	print(f"InitialPhreeqc2Module_mix {type(x)}, {x}")
	#---------
	# Alternative C. to the previous seven methods
	cells = np.full((2), 18)
	cells[1] = 19
	x=bmi.InitialPhreeqcCell2Module(1, cells)
	print(f"InitialPhreeqcCell2Module {type(x)}, {x}")
	#
	# Boundary conditions
	# 
	bc1 = np.full((1), 1)
	bc_conc = bmi.InitialPhreeqc2Concentrations(bc1)
	print(f"InitialPhreeqc2Concentrations {type(bc_conc)}, {bc_conc[0]}")
	#--------
	bc1 = np.full((1), 1)
	bc2 = np.full((1), -1)
	f1 = np.full((1), 1.0)
	bc_conc = bmi.InitialPhreeqc2Concentrations_mix(bc1, bc2, f1)
	print(f"InitialPhreeqc2Concentrations_mix {type(bc_conc)}, {bc_conc[0]}")
	#--------
	v = np.full((1), 1)
	species_c=bmi.InitialPhreeqc2SpeciesConcentrations(v)
	print(f"InitialPhreeqc2SpeciesConcentrations {type(species_c)}, {species_c[0]}")
	#----------
	v1 = np.full((1), 1)
	v2 = np.full((1), -1)
	f1 = np.full((1), 1.0)
	species_c=bmi.InitialPhreeqc2SpeciesConcentrations_mix(v1,v2,f1)
	print(f"InitialPhreeqc2SpeciesConcentrations_mix {type(species_c)}, {species_c[0]}")
	#
	# Get/Set methods for time steping
	#
	x=bmi.GetTime()
	print(f"GetTime {type(x)}, {x}")
	x = bmi.get_current_time()
	print(f"get_current_time {type(x)}, {x}")
	x = bmi.get_start_time()
	print(f"get_start_time {type(x)}, {x}")
	x_ptr = bmi.get_value_ptr("Time")
	print(f"get_value_ptr {type(x_ptr)}, {x_ptr}")
	#---------
	x=bmi.SetTime(0.001)
	print(f"SetTime {type(x)}, {x}, {x_ptr}")
	src = np.full(1, 0.002)
	bmi.set_value("Time", src)
	print(f"set_value('Time'), {x_ptr[0]}")
	#---------
	x=bmi.GetTimeStep()
	print(f"GetTimeStep {type(x)}, {x}")
	dest = np.empty(1, dtype=float)
	x= bmi.get_value("TimeStep", dest)
	print(f"get_value('TimeStep') {type(x)}, {x}")
	d_ptr = bmi.get_value_ptr("TimeStep")
	print(f"get_value_ptr {type(d_ptr)}, {d_ptr}")
	#---------
	x=bmi.SetTimeStep(0.001)
	print(f"SetTimeStep {type(x)}, {x}, {d_ptr[0]}")
	src = np.full(1, 0.0)
	bmi.set_value("TimeStep", src)
	print(f"set_value('TimeStep') {d_ptr[0]}")
	#---------
	c=bmi.GetConcentrations()
	print(f"GetConcentrations {type(c)}, {c[0]}")
	dest = np.empty(nxyz[0]*ncomps[0], dtype=float)
	c = bmi.get_value("Concentrations", dest)
	print(f"get_value('Concentrations') {type(c)}, {c[0]}")
	d_ptr = bmi.get_value_ptr("Concentrations")
	print(f"get_value_ptr('Concentrations') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.SetConcentrations(c)
	print(f"SetConcentrations {type(x)}, {x}")
	bmi.set_value("Concentrations", c)
	print(f"set_value('Concentrations') {type(c)}, {c[0]}")
	#---------
	d=bmi.GetDensityCalculated()
	print(f"GetDensityCalculated {type(d)}, {d[0]}") 
	dest = np.empty(nxyz[0], dtype=float)
	d=bmi.get_value("DensityCalculated", dest)
	print(f"get_value('DensityCalculated') {type(d)}, {d[0]}") 
	d_ptr = bmi.get_value_ptr("DensityCalculated")
	print(f"get_value_ptr('DensityCalculated') {type(d_ptr)}, {d_ptr[0]}") 
	#---------
	x=bmi.SetDensityUser(d)
	print(f"SetDensityUser {type(x)}, {x}")
	d = np.full(nxyz[0], 1.0011)
	bmi.set_value("DensityUser", d)
	print(f"set_value('DensityUser') {type(d)}, {d[0]}")
	#---------
	g = bmi.GetGasCompMoles() 
	print(f"GetGasCompMoles {type(g)}, {g[0]}")
	#---------
	x=bmi.SetGasCompMoles(g)
	print(f"SetGasCompMoles {type(x)}, {x}")
	#---------
	v=bmi.GetGasCompPhi()
	print(f"GetGasCompPhi {type(v)}, {v[0]}")
	#---------
	v=bmi.GetGasCompPressures()
	print(f"GetGasCompPressures {type(v)}, {v[0]}")
	#---------
	v=bmi.GetGasPhaseVolume()
	print(f"GetGasPhaseVolume {type(v)}, {v[0]}")
	#---------
	x=bmi.SetGasPhaseVolume(v)
	print(f"SetGasPhaseVolume {type(x)}, {x}")
	#---------
	for i in range(bmi.GetComponentCount()):
		v = bmi.GetIthConcentration(i)
		#-------
		x = bmi.SetIthConcentration(i, v)
	print(f"GetIthConcentration ")
	print(f"SetIthConcentration ")
	#-------
	for i in range(bmi.GetSpeciesCount()):
		v = bmi.GetIthSpeciesConcentration(i)
		#-------
		x = bmi.SetIthSpeciesConcentration(i, v)
	print(f"GetIthSpeciesConcentration ")
	print(f"SetIthSpeciesConcentration ")
	#-------
	x=bmi.GetPorosity()
	print(f"GetPorosity {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x = bmi.get_value("Porosity", dest)
	print(f"get_value('Porosity') {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("Porosity")
	print(f"get_value_ptr('Porosity') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x = np.full(nxyz[0], 1.001)
	x=bmi.SetPorosity(x)
	print(f"SetPorosity {type(x)}, {x}, {d_ptr[0]}")
	x = np.full(nxyz[0], 1.002)
	bmi.set_value("Porosity", x)
	print(f"set_value('Porosity') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.GetPressure()
	print(f"GetPressure {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x = bmi.get_value("Pressure", dest)
	print(f"get_value('Pressure' {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("Pressure");
	print(f"get_value_ptr('Pressure') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.SetPressure(x)
	print(f"SetPressure {type(x)}, {x}")
	x = np.full(nxyz[0], 2.2)
	bmi.set_value("Pressure", x)
	print(f"set_value('Pressure') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.GetSaturationCalculated() 
	print(f"GetSaturationCalculated {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x=bmi.get_value("SaturationCalculated", dest)
	print(f"get_value('SaturationCalculated') {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("SaturationCalculated");
	print(f"get_value_ptr('SaturationCalculated') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.SetSaturationUser(d_ptr)
	print(f"SetSaturationUser {type(x)}, {x}")
	bmi.set_value("SaturationUser", d_ptr)
	print(f"set_value('SaturationUser') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.GetSolutionVolume()
	print(f"GetSolutionVolume {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x=bmi.get_value("SolutionVolume", dest)
	print(f"get_value('SolutionVolume') {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("SolutionVolume")
	print(f"get_value_ptr('SolutionVolume') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.GetSpeciesConcentrations()           
	print(f"GetSpeciesConcentrations {type(x)}, {x[0]}")
	#---------
	x=bmi.SpeciesConcentrations2Module(x)
	print(f"SpeciesConcentrations2Module {type(x)}, {x}")
	#---------
	v=bmi.GetSpeciesLog10Gammas()  
	print(f"GetSpeciesLog10Gammas {type(v)}, {v[0]}")
	#---------
	v=bmi.GetSpeciesLog10Molalities()   
	print(f"GetSpeciesLog10Molalities {type(v)}, {v[0]}")
	#---------
	x=bmi.GetTemperature()
	print(f"GetTemperature {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x=bmi.get_value("Temperature", dest)
	print(f"get_value('Temperature') {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("Temperature")
	print(f"get_value_ptr('Temperature') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	d = np.full(nxyz[0], 26.0)
	x=bmi.SetTemperature(d)	
	print(f"SetTemperature {type(x)}, {x}, {d_ptr[0]}")	
	d = np.full(nxyz[0], 27.0)
	bmi.set_value("Temperature", d)
	print(f"SetTemperature {type(d)}, {d[0]}, {d_ptr[0]}")	
	#-------
	x = bmi.GetViscosity()
	print(f"GetViscosity  {type(x)}, {x[0]}")
	dest = np.empty(nxyz[0], dtype=float)
	x=bmi.get_value("Viscosity", dest)
	print(f"get_value('Viscosity')  {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("Viscosity")
	print(f"get_value_ptr('Viscosity')  {type(d_ptr)}, {d_ptr[0]}")
	#
	# Take a time step
	#
	bmi.update()    # void method
	print(f"update")
	#---------
	x=bmi.RunCells()	
	print(f"RunCells {x}")
	#-------
	bmi.update_until(86400.0) 
	print(f"update_until ")
	#
	# Selected output
	#
	x=bmi.SetNthSelectedOutput(0)
	print(f"SetNthSelectedOutput {type(x)}, {x}")
	src = np.full(1, 0)
	bmi.set_value("NthSelectedOutput", src)
	print(f"set_value('NthSelectedOutput') {type(src)}, {src[0]}")
	#---------
	x=bmi.GetCurrentSelectedOutputUserNumber()
	print(f"GetCurrentSelectedOutputUserNumber {type(x)}, {x}")
	dest = np.empty(1, dtype=int)
	x=bmi.get_value("CurrentSelectedOutputUserNumber", dest)
	print(f"get_value('CurrentSelectedOutputUserNumber') {type(x)}, {x[0]}")
	#---------
	x=bmi.GetNthSelectedOutputUserNumber(0)
	print(f"GetNthSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	v=bmi.GetSelectedOutput()
	print(f"GetSelectedOutput {type(v)}, {v[0]}")
	n = np.empty(1, dtype=int)
	x = bmi.get_value("SelectedOutputColumnCount", n)
	dest = np.empty(n[0]*nxyz[0], dtype=float)
	v=bmi.get_value("SelectedOutput", dest)
	print(f"get_value('SelectedOutput') {type(v)}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutputColumnCount()
	print(f"GetSelectedOutputColumnCount {type(x)}, {x}")
	dest = np.empty(1, dtype=int)
	x = bmi.get_value("SelectedOutputColumnCount", dest)
	print(f"get_value('SelectedOutputColumnCount') {type(x)}, {x[0]}")
	#---------
	x = bmi.GetSelectedOutputCount()
	print(f"GetSelectedOutputCount {type(x)}, {x}")
	dest = np.empty(1, dtype=int)
	x =bmi.get_value("SelectedOutputCount", dest)
	print(f"get_value('SelectedOutputCount') {type(x)}, {x[0]}")
	#---------
	v=bmi.GetSelectedOutputHeadings()
	print(f"GetSelectedOutputHeadings {type(v)}, {v[87]}, {len(v)}")
	itemsize = bmi.get_var_itemsize("SelectedOutputHeadings")
	s = " " * itemsize
	dest = np.empty(1, dtype=int)
	n = bmi.get_value('SelectedOutputColumnCount', dest)
	dest = np.full(n[0], s)
	v=bmi.get_value("SelectedOutputHeadings", dest)
	print(f"get_value('SelectedOutputHeadings') {type(v)}, {v[87]}")
	#---------
	x=bmi.GetSelectedOutputOn()
	print(f"GetSelectedOutputOn {type(x)}, {x}")
	dest = np.empty(1, dtype=int)
	x=bmi.get_value("SelectedOutputOn", dest)
	print(f"get_value('SelectedOutputOn') {type(x)}, {x[0]}")
	d_ptr = bmi.get_value_ptr("SelectedOutputOn")
	print(f"get_value_ptr('SelectedOutputOn') {type(d_ptr)}, {d_ptr[0]}")
	#---------
	x=bmi.GetSelectedOutputRowCount()
	print(f"GetSelectedOutputRowCount {type(x)}, {x}")
	dest = np.empty(1, dtype=int)
	x=bmi.get_value("SelectedOutputRowCount", dest)
	print(f"get_value('SelectedOutputRowCount') {type(x)}, {x[0]}")
	#---------
	x=bmi.SetCurrentSelectedOutputUserNumber(333)
	print(f"SetCurrentSelectedOutputUserNumber {type(x)}, {x}")
	#
	# Getters
	#
	x=bmi.GetBackwardMapping()
	print(f"GetBackwardMapping {type(x)}, {x[5]}")
	#---------
	x=bmi.GetDatabaseFileName()
	print(f"GetDatabaseFileName {type(x)}, {x}")
	#---------
	x=bmi.GetEndCell()
	print(f"GetEndCell {type(x)}, {x}")
	#---------
	x=bmi.GetErrorHandlerMode()
	print(f"GetErrorHandlerMode {type(x)}, {x}")
	#---------
	x=bmi.GetErrorString()
	print(f"GetErrorString {type(x)}, {x}")
	itemsize = bmi.get_var_nbytes("ErrorString")
	s = " " * (itemsize + 1)
	dest = np.full(1,s)
	x=bmi.get_value("ErrorString", dest)
	print(f"get_value('ErrorString') {type(x)}, '{x[0]}'")
	#---------
	x=bmi.GetFilePrefix()
	print(f"GetFilePrefix {type(x)}, {x}")
	n = bmi.get_var_itemsize("FilePrefix")
	s = " " * n
	dest = np.full(1, s)
	x=bmi.get_value("FilePrefix", dest)
	print(f"GetFilePrefix {type(x)}, {x[0]}")
	#---------
	x=bmi.GetForwardMapping()
	print(f"GetForwardMapping {type(x)}, {x[0]}")
	#---------
	#x=bmi.GetIPhreeqcPointer(0)        #Not Implemented
	#print(f"GetIPhreeqcPointer {x}")
	#---------
	x=bmi.GetMpiMyself()
	print(f"GetMpiMyself {type(x)}, {x}")
	#---------
	x=bmi.GetMpiTasks()
	print(f"GetMpiTasks {type(x)}, {x}")
	#---------
	x=bmi.GetPartitionUZSolids()
	print(f"GetPartitionUZSolids {type(x)}, {x}")
	#---------
	x=bmi.GetPrintChemistryMask()
	print(f"GetPrintChemistryMask {type(x)}, {x[0]}")
	#---------
	x=bmi.GetPrintChemistryOn()
	print(f"GetPrintChemistryOn {type(x)}, {x[0]} {x[1]} {x[2]}")
	#---------
	x=bmi.GetRebalanceByCell()
	print(f"GetRebalanceByCell {type(x)}, {x}")
	#---------
	x=bmi.GetRebalanceFraction()
	print(f"GetRebalanceFraction {type(x)}, {x}")
	#---------
	x=bmi.GetSpeciesSaveOn()
	print(f"GetSpeciesSaveOn {type(x)}, {x}")
	#---------
	x=bmi.GetStartCell()
	print(f"GetStartCell {type(x)}, {x}")
	#---------
	x=bmi.GetTimeConversion()
	print(f"GetTimeConversion {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsExchange()
	print(f"GetUnitsExchange {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsGasPhase()
	print(f"GetUnitsGasPhase {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsKinetics()
	print(f"GetUnitsKinetics {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsPPassemblage()
	print(f"GetUnitsPPassemblage {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsSolution()
	print(f"GetUnitsSolution {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsSSassemblage()
	print(f"GetUnitsSSassemblage {type(x)}, {x}")
	#---------
	x=bmi.GetUnitsSurface()
	print(f"GetUnitsSurface {type(x)}, {x}")
	#---------
	#x=bmi.GetWorkers()                         #Not Implemented
	#print(f"GetWorkers {type(x)}, {x}")
	#
	# Utilities
	#
	bmi2=phreeqcrm.BMIPhreeqcRM()
	x=bmi2.CloseFiles()  
	print(f"CloseFiles {type(x)}, {x}")
	#---------
	bmi2.finalize()
	#---------
	ic1 = np.full((1), 1)
	x = bmi.InitialPhreeqc2Concentrations(ic1)
	tc = np.full((1), 30.0)
	p_atm = np.full((1), 1.5)
	x=bmi.Concentrations2Utility(x, tc, p_atm)
	print(f"Concentrations2Utility {type(x)}, {x}")
	#---------
	bmi.DecodeError(-2)	    
	print(f"DecodeError {x}")
	#---------
	x=bmi.DumpModule(True)
	print(f"DumpModule {type(x)}, {x}")
	#---------	
	bmi.ErrorHandler(0, "string") 
	print(f"OK, just a test: ErrorHandler {type(x)}, {x}")
	#---------
	bmi.ErrorMessage("my error")  
	print(f"OK, just a test: ErrorMessage {type(x)}, {x}")
	#---------
	bmi.LogMessage("Log message")  
	print(f"LogMessage")
	#---------	
	bmi.OutputMessage("Output message")  
	print(f"OutputMessage")
	#---------
	bmi.ScreenMessage("Screen message\n")  
	print(f"ScreenMessage")
	#---------
	x=bmi.StateSave(1)
	print(f"StateSave {type(x)}, {x}")
	#---------
	x=bmi.StateApply(1)
	print(f"StateApply {type(x)}, {x}")
	#---------
	x=bmi.StateDelete(1)
	print(f"StateDelete {type(x)}, {x}")
	#---------
	bmi.WarningMessage("Warning message")  
	print(f"WarningMessage")
	#
	# BMI methods, some have been used above
	#
	#x=bmi.add_output_vars("AddOutputVars", "True")
	#(see above)
	print(f"add_output_vars {type(x)}, {x}")
	x=bmi.get_component_name()
	#---------
	print(f"get_component_name {type(x)}, {x}")
	#---------
	x=bmi.get_current_time()
	print(f"get_current_time {type(x)}, {x}")
	#---------
	x=bmi.get_end_time()
	print(f"get_end_time {type(x)}, {x}")
	#---------
	x = bmi.get_grid_rank(0)
	print(f"get_grid_rank {type(x)}, {x}")
	#-------
	n = bmi.get_grid_size(0)
	print(f"get_grid_size {type(x)}, {x}")
	#-------
	x = bmi.get_grid_type(0)
	print(f"get_grid_type {type(x)}, {x}")
	#-------
	x=bmi.get_input_item_count()
	print(f"get_input_item_count {type(x)}, {x}")
	#---------
	x=bmi.get_input_var_names()
	print(f"get_input_var_names {type(x)}, {x[0]}")
	#---------
	x=bmi.get_output_item_count()
	print(f"get_output_item_count {type(x)}, {x}")
	#---------
	x=bmi.get_output_var_names()
	print(f"get_output_var_names {type(x)}, {x[0]}")
	#---------
	x=bmi.get_pointable_item_count()
	print(f"get_pointable_item_count {type(x)}, {x}")
	#---------
	x=bmi.get_pointable_var_names()
	print(f"get_pointable_var_names {type(x)}, {x[0]}")
	#---------
	x=bmi.get_time_step()
	print(f"get_time_step {type(x)}, {x}")
	#---------
	x=bmi.get_time_units()
	print(f"get_time_units {type(x)}, {x}")
	#---------
	dest = np.empty(nxyz[0], dtype=float)
	x=bmi.get_value("solution_saturation_index_Calcite", dest)
	print(f"get_value {type(x)}, {x[0]}")
	#---------
	x=bmi.get_var_itemsize("solution_saturation_index_Calcite")
	print(f"get_var_itemsize {type(x)}, {x}")
	#---------
	x=bmi.get_var_nbytes("solution_saturation_index_Calcite")
	print(f"get_var_nbytes {type(x)}, {x}")
	#---------
	x=bmi.get_var_type("solution_saturation_index_Calcite")
	print(f"get_var_type {type(x)}, {x}")
	#---------
	x=bmi.get_var_units("solution_saturation_index_Calcite")
	print(f"get_var_units {type(x)}, {x}")
	#---------
	print(f"AddOutputVars")
	v = bmi.get_output_var_names()
	for i in range(len(v)):
		itemsize = bmi.get_var_itemsize(v[i])
		nbytes = bmi.get_var_nbytes(v[i])
		vtype = bmi.get_var_type(v[i])
		if itemsize == 0:
			itemsize = 1
		if nbytes == 0:
			nbytes = 1
		dim = nbytes // itemsize
		if (vtype == "float64"):
			dest  = np.empty(dim, dtype=float)
			x = bmi.get_value(v[i], dest)
			print(f"     {v[i]}, {dest[0]}")
		if (vtype == "int32"):
			dest = np.empty(dim, dtype=int)
			x = bmi.get_value(v[i], dest)
			print(f"     {v[i]}, {dest[0]}")
		if vtype.startswith("<U"):
			s = " " * itemsize
			dest = np.full(dim, s)
			x = bmi.get_value(v[i], dest)
			print(f"     {v[i]}, {dest[0]}")
	#x=bmi.initialize("file")
	# (See above)
	print(f"initialize")
	#-------	
	src = np.full(1, 1.0)
	bmi.set_value("Time", src) 
	print(f"set_value")
	#---------
	bmi.update()    
	print(f"update")
	#-------	
	bmi.update_until(864000.0)
	print(f"update_until")
	#---------
	bmi.finalize()    
	print(f"finalize ")
	#---------
	print("Success.")
	return

	#Should be private: x=bmi.ReturnHandler()
	#TODO x=bmi.SetMpiWorkerCallbackC()
	#TODO x=bmi.SetMpiWorkerCallbackCookie()
	#TODO x=bmi.MpiAbort()
	#TODO x=bmi.MpiWorker()
	#TODO x=bmi.MpiWorkerBreak()

if __name__ == '__main__':
	testallmethods_py()

#endif # USE_YAML
