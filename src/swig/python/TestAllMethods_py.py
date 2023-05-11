import phreeqcrm
import yamlphreeqcrm
import numpy as np
#ifdef USE_YAML
    #module mydata
    #  double precision, dimension(:), pointer :: K_ptr
    #  integer                                 :: rm_id
    #end module mydata

def testallmethods_py():

	yrm = yamlphreeqcrm.YAMLPhreeqcRM()
	yrm.YAMLSetGridCellCount(40)
	yrm.YAMLThreadCount(3)
	YAML_filename = "testallmethods_py.yaml"
	yrm.WriteYAMLDoc(YAML_filename);
	#
	# Use all BMIPhreeqcRM methods roughly in order of use
	#
	bmi=phreeqcrm.BMIPhreeqcRM()
	print(f"BMIPhreeqcRM {bmi}")
	#---------
	bmi.Initialize(YAML_filename)   # void function
	print(f"Initialize")
	#---------
	nxyz = bmi.GetValue("GridCellCount")
	print(f"GetValue('GridCellCount') {type(nxyz)}, {nxyz}")
	#---------
	nxyz = bmi.GetGridCellCount()
	print(f"GetValue('GetGridCellCount') {type(nxyz)}, {nxyz}")
	#---------
	x=bmi.GetThreadCount()
	print(f"GetThreadCount {type(x)}, {x}")
	#---------
	grid2chem = phreeqcrm.IntVector(nxyz, -1)
	for i in range(nxyz//2):
		grid2chem[i] = i
	x=bmi.CreateMapping(grid2chem)
	print(f"CreateMapping {type(x)}, {x}")
	#--------- LoadDatabase removes definitions in instance	
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
	x=bmi.UseSolutionDensityVolume(False)
	print(f"UseSolutionDensityVolume {type(x)}, {x}")
	#---------
	x=bmi.SetTimeConversion(1.0/86400.0)
	print(f"SetTimeConversion {type(x)}, {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetRepresentativeVolume(v)
	print(f"SetRepresentativeVolume {type(x)}, {x}")
	#---------Chemistry cells may be fewer than GridCellCount
	nchem=bmi.GetChemistryCellCount()
	print(f"GridCellCount {type(nchem)}, {nchem}")
	#
	# Begin initialization
	#
	v = phreeqcrm.IntVector(nxyz, 1)
	x=bmi.SetPrintChemistryMask(v)
	print(f"SetPrintChemistryMask {type(x)}, {x}")
	#---------
	x=bmi.SetPrintChemistryOn(True,True,True)
	print(f"SetPrintChemistryOn {type(x)}, {x}")
	#---------
	x = bmi.RunFile(True, True, True, "all_reactants.pqi")
	print(f"RunFile {type(x)}, {x}")	
	#---------
	#x=bmi.AddOutputVars("SolutionActivities", "True")
	#x=bmi.AddOutputVars("SolutionMolalities", "True")
	#x=bmi.AddOutputVars("SaturationIndices", "True")
	x=bmi.AddOutputVars("SolutionActivities", "H+ Ca+2 Na+")
	x=bmi.AddOutputVars("SolutionMolalities", "OH- Cl-")
	x=bmi.AddOutputVars("SaturationIndices", "Calcite Dolomite")
	print(f"AddOutputVars {type(x)}, {x}")		
	#---------
	ncomps = bmi.FindComponents()
	print(f"FindComponents {type(ncomps)}, {ncomps}")
	#---------
	ncomps=bmi.GetComponentCount()
	print(f"GetComponentCount {type(ncomps)}, {ncomps}")
	#---------
	ncomps=bmi.GetValue("ComponentCount")
	print(f"bmi.GetValue('ComponentCount') {type(ncomps)}, {ncomps}")
	#---------
	x=bmi.GetComponents()
	print(f"GetComponents {type(x)}, {x}")
	#---------
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
	v = phreeqcrm.IntVector(nxyz, 1)
	x=bmi.InitialEquilibriumPhases2Module(v)
	print(f"InitialEquilibriumPhases2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialExchanges2Module(v)
	print(f"InitialExchanges2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialGasPhases2Module(v)
	print(f"InitialGasPhases2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialKinetics2Module(v)
	print(f"InitialKinetics2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialSolutions2Module(v)
	print(f"InitialSolutions2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialSolidSolutions2Module(v)
	print(f"InitialSolidSolutions2Module {type(x)}, {x}")
	#---------
	x=bmi.InitialSurfaces2Module(v)	
	print(f"InitialSurfaces2Module {type(x)}, {x}")
	#---------
	# Alternative A. to the previous seven methods
	v = phreeqcrm.IntVector(nxyz*7, 1)
	x=bmi.InitialPhreeqc2Module(v)
	print(f"InitialPhreeqc2Module {type(x)}, {x}")
	#---------
	# Alternative B. to the previous seven methods
	v1 = phreeqcrm.IntVector(nxyz*7, 1)
	v2 = phreeqcrm.IntVector(nxyz*7, -1)
	f1 = phreeqcrm.DoubleVector(nxyz*7, 1.0)
	x=bmi.InitialPhreeqc2Module(v1, v2, f1)
	print(f"InitialPhreeqc2Module_mix {type(x)}, {x}")
	#---------
	# Alternative C. to the previous seven methods
	cells = [18, 19]
	x=bmi.InitialPhreeqcCell2Module(1, cells)
	print(f"InitialPhreeqcCell2Module {type(x)}, {x}")
	#
	# Boundary conditions
	# 
	bc1 = phreeqcrm.IntVector(1,1)
	bc_conc = bmi.InitialPhreeqc2Concentrations(bc1)
	print(f"InitialPhreeqc2Concentrations {type(bc_conc)}, {bc_conc[0]}")
	#--------
	bc1 = phreeqcrm.IntVector(1,1)
	bc2 = phreeqcrm.IntVector(1,-1)
	f1 = phreeqcrm.DoubleVector(1,1.0)
	bc_conc = bmi.InitialPhreeqc2Concentrations_mix(bc1, bc2, f1)
	print(f"InitialPhreeqc2Concentrations_mix {type(bc_conc)}, {bc_conc[0]}")
	#--------
	v = phreeqcrm.IntVector(1,1)
	species_c=bmi.InitialPhreeqc2SpeciesConcentrations(v)
	print(f"InitialPhreeqc2SpeciesConcentrations {type(species_c)}, {species_c[0]}")
	#----------
	v1 = phreeqcrm.IntVector(1,1)
	v2 = phreeqcrm.IntVector(1,-1)
	f1 = phreeqcrm.DoubleVector(1,1.0)
	species_c=bmi.InitialPhreeqc2SpeciesConcentrations_mix(v1,v2,f1)
	print(f"InitialPhreeqc2SpeciesConcentrations_mix {type(species_c)}, {species_c[0]}")
	#
	# Get/Set methods for time steping
	#
	x=bmi.SetTime(0.0)
	print(f"SetTime {type(x)}, {x}")
	#---------
	x=bmi.SetTimeStep(0.0)
	print(f"SetTimeStep {type(x)}, {x}")
	#---------
	c = phreeqcrm.DoubleVector()
	x=bmi.GetConcentrations(c)
	print(f"GetConcentrations {type(x)}, {x}, {c[0]}.")
		#---------
	x=bmi.SetConcentrations(c)
	print(f"SetConcentrations {type(x)}, {x}")
	#---------
	c = bmi.GetValue("Concentrations")
	print(type(c))
	print(f"GetValue('Concentrations') {c[0]}")
	#---------
	bmi.SetValue("Concentrations", c)
	print(f"SetValue('Concentrations')")
	#---------
	d = phreeqcrm.DoubleVector(nxyz, 0)
	x=bmi.GetDensity(d)
	print(f"GetDensity {type(x)}, {x}, {d[0]}") 
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetDensity(v)
	print(f"SetDensity {x}")
	#---------
	v = phreeqcrm.DoubleVector()
	x = bmi.GetGasCompMoles(v)
	print(f"GetGasCompMoles {type(x)}, {x}")
	#---------
	x=bmi.SetGasCompMoles(v)
	print(f"SetGasCompMoles {type(x)}, {x}")
	#---------
	x=bmi.GetGasCompPhi(v)
	print(f"GetGasCompPhi {type(x)}, {x}, {v[0]}")
	#---------
	x=bmi.GetGasCompPressures(v)
	print(f"GetGasCompPressures {type(x)}, {x}, {v[0]}")
	#---------
	x=bmi.GetGasPhaseVolume(v)
	print(f"GetGasPhaseVolume {type(x)}, {x}, {v[0]}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetGasPhaseVolume(v)
	print(f"SetGasPhaseVolume {type(x)}, {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 0.21)
	x=bmi.SetPorosity(v)
	print(f"SetPorosity {type(x)}, {x}")
	#---------
	x=bmi.GetPressure()
	print(f"GetPressure {type(x)}, {x[0]}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 3.0)
	x=bmi.SetPressure(v)
	print(f"SetPressure {type(x)}, {x}")
	#---------
	x=bmi.SetSaturation(v)
	print(f"SetSaturation {type(x)}, {x}")
	#---------
	x=bmi.GetSolutionVolume()
	print(f"GetSolutionVolume {type(x)}, {x[0]}")
	#---------
	v = phreeqcrm.DoubleVector()
	x=bmi.GetSpeciesConcentrations(v)
	print(f"GetSpeciesConcentrations {type(x)}, {v[0]}")
	#---------
	x=bmi.SpeciesConcentrations2Module(v)
	print(f"SpeciesConcentrations2Module {type(x)}, {x}")
	#---------
	x=bmi.GetSpeciesLog10Gammas(v)
	print(f"GetSpeciesLog10Gammas {type(x)}, {x}, {v[0]}")
	#---------
	x=bmi.GetSpeciesLog10Molalities(v)
	print(f"GetSpeciesLog10Molalities {type(x)}, {x}, {v[0]}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 26.0)
	x=bmi.SetTemperature(v)	
	print(f"SetTemperature {type(x)}, {x}")
	#
	# Take a time step
	#
	bmi.Update()    # void method
	print(f"Update")
	#---------
	x=bmi.RunCells()	
	print(f"RunCells {x}")
	#
	# Selected output
	#
	x=bmi.GetCurrentSelectedOutputUserNumber()
	print(f"GetCurrentSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	x=bmi.SetCurrentSelectedOutputUserNumber(1)
	print(f"SetCurrentSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	x=bmi.GetNthSelectedOutputUserNumber(0)
	print(f"GetNthSelectedOutputUserNumber {type(x)}, {x}")
	#---------
	x=bmi.SetNthSelectedOutput(0)
	print(f"SetNthSelectedOutput {type(x)}, {x}")
	#---------
	x=bmi.GetSelectedOutput(v)
	print(f"GetSelectedOutput {type(x)}, {x}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutputColumnCount()
	print(f"GetSelectedOutputColumnCount {type(x)}, {x}")
	#---------
	v = phreeqcrm.StringVector()
	x=bmi.GetSelectedOutputHeadings(v)
	print(f"GetSelectedOutputHeadings {type(x)}, {x}, {v[87]}")
	#---------
	x=bmi.GetSelectedOutputRowCount()
	print(f"GetSelectedOutputRowCount {type(x)}, {x}")
	#
	# Getters
	#
	x=bmi.GetBackwardMapping()
	print(type(x))
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
	#---------
	x=bmi.GetFilePrefix()
	print(f"GetFilePrefix {type(x)}, {x}")
	#---------
	x=bmi.GetForwardMapping()
	print(f"GetForwardMapping {type(x)}, {x[0]}")
	#---------
	x=bmi.GetGfw()
	print(f"GetGfw {type(x)}, {x[0]}")
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
	x=bmi.GetPorosity()
	print(f"GetPorosity {type(x)}, {x[0]}")
	#---------
	x=bmi.GetPrintChemistryMask()
	print(f"GetPrintChemistryMask {type(x)}, {x[0]}")
	#---------
	x=bmi.GetPrintChemistryOn()
	print(type(x))
	print(f"GetPrintChemistryOn {type(x)}, {x[0]} {x[1]} {x[2]}")
	#---------
	x=bmi.GetRebalanceByCell()
	print(f"GetRebalanceByCell {type(x)}, {x}")
	#---------
	x=bmi.GetRebalanceFraction()
	print(f"GetRebalanceFraction {type(x)}, {x}")
	#---------
	v = phreeqcrm.DoubleVector()
	x=bmi.GetSaturation(v)
	print(f"GetSaturation {type(x)}, {x}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutputOn()
	print(f"GetSelectedOutputOn {type(x)}, {x}")
	#---------
	x=bmi.GetSpeciesSaveOn()
	print(f"GetSpeciesSaveOn {type(x)}, {x}")
	#---------
	x=bmi.GetStartCell()
	print(f"GetStartCell {type(x)}, {x}")
	#---------
	x=bmi.GetTemperature()
	print(f"GetTemperature {type(x)}, {x[0]}")
	#---------
	x=bmi.GetTime()
	print(f"GetTime {type(x)}, {x}")
	#---------
	x=bmi.GetTimeConversion()
	print(f"GetTimeConversion {type(x)}, {x}")
	#---------
	x=bmi.GetTimeStep()
	print(f"GetTimeStep {type(x)}, {x}")
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
	bc1 = phreeqcrm.IntVector(1,1)
	#out = bmi.InitialPhreeqc2Concentrations(bc1)
	#status = out[0]
	#print(type(out[1]))
	#bc_conc = list(out[1])
	conc =bmi.InitialPhreeqc2Concentrations(bc1)
	tc = [30.0]
	p_atm = [1.5]
	bc_double_vect = phreeqcrm.DoubleVector(len(bc_conc))
	for i in range(len(bc_conc)):
		bc_double_vect[i] = bc_conc[i]
	x=bmi.Concentrations2Utility(bc_double_vect, tc, p_atm)
	print(f"Concentrations2Utility {type(x)}, {x}")
	#---------
	bmi.DecodeError(-2)	         # void function
	print(f"DecodeError {x}")
	#---------
	x=bmi.DumpModule(True)
	print(f"DumpModule {type(x)}, {x}")
	#---------	
	bmi.ErrorHandler(0, "string") # void function
	print(f"OK, just a test: ErrorHandler {type(x)}, {x}")
	#---------
	bmi.ErrorMessage("my error")  # void function
	print(f"OK, just a test: ErrorMessage {type(x)}, {x}")
	#---------
	bmi.LogMessage("Log message")  # void method
	print(f"LogMessage")
	#---------	
	bmi.OutputMessage("Output message")  # void method
	print(f"OutputMessage")
	#---------
	bmi.ScreenMessage("Screen message\n")  # void method
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
	bmi.WarningMessage("Warning message")  # void method
	print(f"WarningMessage")
	#
	# BMI methods, some have been used above
	#
	x=bmi.GetComponentName()
	print(f"GetComponentName {type(x)}, {x}")
	#---------
	x=bmi.GetCurrentTime()
	print(f"GetCurrentTime {type(x)}, {x}")
	#---------
	x=bmi.GetEndTime()
	print(f"GetEndTime {type(x)}, {x}")
	#---------
	x=bmi.GetInputItemCount()
	print(f"GetInputItemCount {type(x)}, {x}")
	#---------
	x=bmi.GetInputVarNames()
	print(f"GetInputVarNames {type(x)}, {x[0]}")
	#---------
	x=bmi.GetOutputItemCount()
	print(f"GetOutputItemCount {type(x)}, {x}")
	#---------
	x=bmi.GetOutputVarNames()
	print(f"GetOutputVarNames {type(x)}, {x[0]}")
	#---------
	x=bmi.GetTimeStep()
	print(f"GetTimeStep {type(x)}, {x}")
	#---------
	x=bmi.GetTimeUnits()
	print(f"GetTimeUnits {type(x)}, {x}")
	#---------
	x=bmi.GetValue("solution_saturation_index_Calcite")
	print(f"GetValue {type(x)}, {x[0]}")
	#---------
	x=bmi.GetVarItemsize("solution_saturation_index_Calcite")
	print(f"GetVarItemsize {type(x)}, {x}")
	#---------
	x=bmi.GetVarNbytes("solution_saturation_index_Calcite")
	print(f"GetVarNbytes {type(x)}, {x}")
	#---------
	x=bmi.GetVarType("solution_saturation_index_Calcite")
	print(f"GetVarType {type(x)}, {x}")
	#---------
	x=bmi.GetVarUnits("solution_saturation_index_Calcite")
	print(f"GetVarUnits {type(x)}, {x}")
	#---------
	#x=bmi.Initialize("file")
	# See above
	bmi.Update()    # void method
	print(f"Update")
	#---------
	bmi.SetValue("Time", 1.0)    # void method
	print(f"SetValue")
	#---------
	x=bmi.CloseFiles()  # not a BMI method
	print(f"CloseFiles {type(x)}, {x}")
	#---------
	bmi.Finalize()    # void method
	print(f"Finalize {x}")
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
