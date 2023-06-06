
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
	print(f"Initialize")
	#---------
	nxyz = bmi.GetGridCellCount()
	nxyz = bmi.get_value("GridCellCount")
	#int *i_ptr = (int*)bmi.GetValuePtr("GridCellCount")  # pointer
	print(f"GetValue('GridCellCount') {type(nxyz)}, {nxyz}")
	#---------
	x=bmi.GetThreadCount()
	print(f"GetThreadCount {type(x)}, {x}")
	#---------
	# Inactive cells or symmetry
	grid2chem = np.full((nxyz), -1)
	for i in range(nxyz//2):
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
	bmi.SetValue("FilePrefix", "TestAllMethods_py")
	print(f"SetFilePrefix {type(x)}, {x}")
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
	bmi.set_value("SelectedOutputOn", True)
	print(f"SetSelectedOutputOn {type(x)}, {x}")
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
	print(f"GridCellCount {type(nchem)}, {nchem}")
	#---------
	ncomps=bmi.GetComponentCount()
	ncomps = bmi.get_value("ComponentCount")
	#i_ptr = (int*) bmi.GetValuePtr("ComponentCount")   # Pointer
	print(f"GetComponentCount {type(ncomps)}, {ncomps}")
	#---------
	x=bmi.GetComponents()	
	x = bmi.get_value("Components")
	print(f"GetComponents {type(x)}, {x}")
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
	x=bmi.get_value("gfw")
	#double *d_ptr = (double*)bmi.GetValuePtr("Gfw")  # Pointer
	print(f"GetGfw {type(x)}, {x[0]}")
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
	d = bmi.get_current_time()
	d = bmi.get_start_time()
	#d_ptr = bmi.GetValuePtr("Time")
	print(f"GetTime {type(x)}, {x}")
	#---------
	x=bmi.SetTime(0.0)
	bmi.set_value("Time", 0.0)
	print(f"SetTime {type(x)}, {x}")
	#---------
	x=bmi.GetTimeStep()
	x= bmi.get_value("TimeStep")
	#d_ptr = (double*)bmi.GetValuePtr("TimeStep")
	print(f"GetTimeStep {type(x)}, {x}")
	#---------
	x=bmi.SetTimeStep(0.0)
	bmi.set_value("TimeStep", 0.0)
	print(f"SetTimeStep {type(x)}, {x}")
	#---------
	c=bmi.GetConcentrations()
	c = bmi.get_value("Concentrations")
	#d_ptr = (double*)bmi.GetValuePtr("Concentrations")
	print(f"GetConcentrations {type(c)}, {c[0]}")
	#---------
	x=bmi.SetConcentrations(c)
	bmi.set_value("Concentrations", c)
	print(f"SetConcentrations {type(x)}, {x}")
	#---------
	d=bmi.GetDensityCalculated()
	d=bmi.get_value("DensityCalculated")
	#d_ptr = (double*)bmi.GetValuePtr("DensityCalculated")
	print(f"GetDensityCalculated {type(d)}, {d[0]}") 
	#---------
	x=bmi.SetDensityUser(d)
	bmi.set_value("DensityUser",d)
	print(f"SetDensityUser {type(x)}, {x}")
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
	x = bmi.get_value("Porosity")
	#d_ptr = (double*)bmi.GetValuePtr("Porosity")
	print(f"GetPorosity {type(x)}, {x[0]}")
	#---------
	bmi.set_value("Porosity", x)
	x=bmi.SetPorosity(x)
	print(f"SetPorosity {type(x)}, {x}")
	#---------
	x=bmi.GetPressure()
	x = bmi.get_value("Pressure")
	#d_ptr = (double*)bmi.GetValuePtr("Pressure");
	print(f"GetPressure {type(x)}, {x[0]}")
	#---------
	bmi.set_value("Pressure", x)
	x=bmi.SetPressure(x)
	print(f"SetPressure {type(x)}, {x}")
	#---------
	x=bmi.get_value("SaturationCalculated")
	x=bmi.GetSaturationCalculated() 
	#d_ptr = (double*)bmi.GetValuePtr("SaturationCalculated");
	print(f"GetSaturationCalculated {type(x)}, {x[0]}")
	#---------
	bmi.set_value("SaturationUser", x)
	x=bmi.SetSaturationUser(x)
	print(f"SetSaturationUser {type(x)}, {x}")
	#---------
	x=bmi.get_value("SolutionVolume")
	x=bmi.GetSolutionVolume()
	#d_ptr = (double*)bmi.GetValuePtr("SolutionVolume");
	print(f"GetSolutionVolume {type(x)}, {x[0]}")
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
	x=bmi.get_value("Temperature")
	print(f"GetTemperature {type(x)}, {x[0]}")
	#---------
	bmi.set_value("Temperature", x)
	x=bmi.SetTemperature(x)	
	#d_ptr = (double*)bmi.GetValuePtr("Temperature");
	print(f"SetTemperature {type(x)}, {x}")
	#-------
	x = bmi.GetViscosity()
	x=bmi.get_value("Viscosity")
	#d_ptr = (double*)bmi.GetValuePtr("Viscosity")
	print(f"GetViscosity  {type(x)}, {x[0]}")
	#
	# Take a time step
	#
	bmi.update()    # void method
	print(f"Update")
	#---------
	x=bmi.RunCells()	
	print(f"RunCells {x}")
	#-------
	bmi.update_until(86400.0) 
	print(f"UpdateUntil ")
	#
	# Selected output
	#
	x=bmi.SetNthSelectedOutput(0)
	bmi.set_value("NthSelectedOutput", 0)
	print(f"SetNthSelectedOutput {type(x)}, {x}")
	#---------
	x=bmi.GetCurrentSelectedOutputUserNumber()
	x=bmi.get_value("CurrentSelectedOutputUserNumber")
	print(f"GetCurrentSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	x=bmi.GetNthSelectedOutputUserNumber(0)
	print(f"GetNthSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	v=bmi.GetSelectedOutput()
	v=bmi.get_value("SelectedOutput")
	print(f"GetSelectedOutput {type(v)}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutputColumnCount()
	x=bmi.get_value("SelectedOutputColumnCount")
	print(f"GetSelectedOutputColumnCount {type(x)}, {x}")
	#---------
	x = bmi.GetSelectedOutputCount()
	x =bmi.get_value("SelectedOutputCount")
	print(f"GetSelectedOutputCount {type(x)}, {x}")
	#---------
	v=bmi.GetSelectedOutputHeadings()
	v=bmi.get_value("SelectedOutputHeadings")
	print(f"GetSelectedOutputHeadings {type(v)}, {v[87]}")
	#---------
	x=bmi.GetSelectedOutputOn()
	x=bmi.get_value("SelectedOutputOn")
	#bool *b_ptr = (bool*)bmi.GetValuePtr("SelectedOutputOn")
	print(f"GetSelectedOutputOn {type(x)}, {x}")
	#---------
	x=bmi.GetSelectedOutputRowCount()
	x=bmi.get_value("SelectedOutputRowCount")
	print(f"GetSelectedOutputRowCount {type(x)}, {x}")
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
	x=bmi.get_value("ErrorString")
	print(f"GetErrorString {type(x)}, {x}")
	#---------
	x=bmi.GetFilePrefix()
	x=bmi.get_value("FilePrefix")
	print(f"GetFilePrefix {type(x)}, {x}")
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
	x=bmi.get_value("solution_saturation_index_Calcite")
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
	#x=bmi.initialize("file")
	# (See above)
	print(f"initialize")
	#-------	
	bmi.set_value("Time", 1.0) 
	print(f"set_value")
	#---------
	bmi.update()    
	print(f"update")
	#-------	
	bmi.update_until(864000.0)
	print(f"UpdateUntil")
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
