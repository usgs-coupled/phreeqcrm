import phreeqcrm
import numpy as np
#ifdef USE_YAML
    #module mydata
    #  double precision, dimension(:), pointer :: K_ptr
    #  integer                                 :: rm_id
    #end module mydata

def testbmi_py():

	#x=GetGridCellCountYAML("AdvectBMI_test_py.yaml") #not defined
	
	bmi=phreeqcrm.BMIPhreeqcRM()
	print(f"BMIPhreeqcRM {bmi}")
	#---------
	bmi.Initialize("AdvectBMI_test_py.yaml")   # void function
	print(f"Initialize")
	#---------
	nxyz = bmi.GetValue("GridCellCount")
	print(f"GetValue('GridCellCount') {nxyz}")
	#---------
	##x=bmi.LoadDatabase("phreeqc.dat")
	print(f"LoadDatabase")
	#---------
	x=bmi.SetSpeciesSaveOn(True)
	print(f"SetSpeciesSaveOn {x}")	
	#---------
	x=bmi.SetPrintChemistryOn(False,True,False)
	print(f"SetPrintChemistryOn {x}")
	#---------
	x = bmi.RunFile(True, True, True, "advectBMI_test.pqi")
	print(f"RunFile {x}")	
	#---------
	x=bmi.SetComponentH2O(False)
	print(f"SetComponentH2O {x}")
	#---------
	x=bmi.SetCurrentSelectedOutputUserNumber(1)
	print(f"SetCurrentSelectedOutputUserNumber {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetDensity(v)
	print(f"SetDensity {x}")
	#---------
	x=bmi.SetDumpFileName("AdvectBMI_test_py.dump")
	print(f"SetDumpFileName {x}")
	#---------
	x=bmi.SetErrorHandlerMode(1)
	print(f"SetErrorHandlerMode {x}")
	#---------
	x=bmi.SetErrorOn(True)
	print(f"SetErrorOn {x}")
	#---------
	x=bmi.SetFilePrefix("AdvectBMI_test_py")
	print(f"SetFilePrefix {x}")
	#---------
	x=bmi.OpenFiles()
	print(f"OpenFiles {x}")
	#---------
	x=bmi.SetPartitionUZSolids(False)	
	print(f"SetPartitionUZSolids {x}")
	#---------
	x=bmi.SetRebalanceByCell(True)
	print(f"SetRebalanceByCell {x}")
	#---------
	x=bmi.SetRebalanceFraction(0.5)
	print(f"SetRebalanceFraction {x}")
	#---------
	x=bmi.SetScreenOn(True)
	print(f"SetScreenOn {x}")
	#---------
	x=bmi.SetSelectedOutputOn(True)
	print(f"SetSelectedOutputOn {x}")
	#---------
	x=bmi.SetTime(0.0)
	print(f"SetTime {x}")	
	#---------
	time_conversion = 1.0 / 86400.0
	bmi.SetTimeConversion(time_conversion)
	print(f"SetTimeConversion {x}")
	#---------
	x=bmi.SetTimeStep(0.0)
	print(f"SetTimeStep {x}")	
	#---------
	x=bmi.SetUnitsExchange(1)
	print(f"SetUnitsExchange {x}")
	#---------
	x=bmi.SetUnitsGasPhase(1)
	print(f"SetUnitsGasPhase {x}")
	#---------
	x=bmi.SetUnitsKinetics(1)
	print(f"SetUnitsKinetics {x}")
	#---------
	x=bmi.SetUnitsPPassemblage(1)
	print(f"SetUnitsPPassemblage {x}")
	#---------
	x=bmi.SetUnitsSolution(2)
	print(f"SetUnitsSolution {x}")
	#---------
	x=bmi.SetUnitsSSassemblage(1)
	print(f"SetUnitsSSassemblage {x}")
	#---------
	x=bmi.SetUnitsSurface(1)
	print(f"SetUnitsSurface {x}")
	#---------
	x=bmi.UseSolutionDensityVolume(False)
	print(f"UseSolutionDensityVolume {x}")
	#---------
	x=bmi.SetSpeciesSaveOn(True)
	print(f"SetSpeciesSaveOn {x}")	
	#---------
	x=bmi.AddOutputVars("SolutionActivities", "True")
	x=bmi.AddOutputVars("SolutionMolalities", "True")
	x=bmi.AddOutputVars("SaturationIndices", "True")
	x=bmi.AddOutputVars("SolutionActivities", "H+ Ca+2 Na+")
	x=bmi.AddOutputVars("SolutionMolalities", "OH- Cl-")
	x=bmi.AddOutputVars("SaturationIndices", "Calcite Dolomite")
	print(f"AddOutputVars {x}")		
	#---------
	x = bmi.FindComponents()
	print(f"FindComponents {x}")
	#---------
	ncomps = bmi.GetValue("ComponentCount")
	print(f"GetValue('ComponentCount') {ncomps}")
	#---------
	input = "DELETE; -all"
	bmi.RunString(True, False, True, input)
	print(f"RunString {x}")
	#---------
	v = phreeqcrm.IntVector(nxyz, 1)
	x=bmi.InitialEquilibriumPhases2Module(v)
	print(f"InitialEquilibriumPhases2Module {x}")
	#---------
	x=bmi.InitialExchanges2Module(v)
	print(f"InitialExchanges2Module {x}")
	#---------
	x=bmi.InitialGasPhases2Module(v)
	print(f"InitialGasPhases2Module {x}")
	#---------
	x=bmi.InitialKinetics2Module(v)
	print(f"InitialKinetics2Module {x}")
	#---------
	x=bmi.InitialSolutions2Module(v)
	print(f"InitialSolutions2Module {x}")
	#---------
	x=bmi.InitialSolidSolutions2Module(v)
	print(f"InitialSolidSolutions2Module {x}")
	#---------
	x=bmi.InitialSurfaces2Module(v)	
	print(f"InitialSurfaces2Module {x}")
	#---------
	# Alternative A. to the previous seven methods
	v = phreeqcrm.IntVector(nxyz*7, 1)
	x=bmi.InitialPhreeqc2Module(v)
	print(f"InitialPhreeqc2Module {x}")
	#---------
	# Alternative B. to the previous seven methods
	v1 = phreeqcrm.IntVector(nxyz*7, 1)
	v2 = phreeqcrm.IntVector(nxyz*7, -1)
	f1 = phreeqcrm.DoubleVector(nxyz*7, 1.0)
	x=bmi.InitialPhreeqc2Module(v1, v2, f1)
	print(f"InitialPhreeqc2Module_mix {x}")
	#---------
	# Alternative C. to the previous seven methods
	cells = [18, 19]
	x=bmi.InitialPhreeqcCell2Module(1, cells)
	print(f"InitialPhreeqcCell2Module {x}")
	#---------
	ngas = bmi.GetGasComponentsCount()
	print(f"GetGasComponentsCount {x}")
	#---------
	v = phreeqcrm.DoubleVector()
	x = bmi.GetGasCompMoles(v)
	print(f"GetGasCompMoles {x}")
	#---------
	x=bmi.SetGasCompMoles(v)
	print(f"SetGasCompMoles {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetGasPhaseVolume(v)
	print(f"SetGasPhaseVolume {x}")
	#---------
	#TODO x=bmi.SetMpiWorkerCallbackC()
	#TODO x=bmi.SetMpiWorkerCallbackCookie()
	x=bmi.SetNthSelectedOutput(0)
	print(f"SetNthSelectedOutput {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 0.21)
	x=bmi.SetPorosity(v)
	print(f"SetPorosity {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 3.0)
	x=bmi.SetPressure(v)
	print(f"SetPressure {x}")
	#---------
	v = phreeqcrm.IntVector(nxyz, 1)
	x=bmi.SetPrintChemistryMask(v)
	print(f"SetPrintChemistryMask {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 1.0)
	x=bmi.SetRepresentativeVolume(v)
	print(f"SetRepresentativeVolume {x}")
	#---------
	x=bmi.SetSaturation(v)
	print(f"SetSaturation {x}")
	#---------
	v = phreeqcrm.DoubleVector(nxyz, 26.0)
	x=bmi.SetTemperature(v)	
	print(f"SetTemperature {x}")
	#---------
	v = phreeqcrm.IntVector(1,1)
	x=bmi.InitialPhreeqc2SpeciesConcentrations(v)
	status = x[0]
	species_c = x[1]
	print(f"InitialPhreeqc2SpeciesConcentrations {status}, {species_c[0]}")
	#---------
	v1 = phreeqcrm.IntVector(1,-1)
	f1 = phreeqcrm.DoubleVector(1,1.0)
	x=bmi.InitialPhreeqc2SpeciesConcentrations(v,v1,f1)
	status = x[0]
	species_c = x[1]
	print(f"InitialPhreeqc2SpeciesConcentrations_mix {status}, {species_c[0]}")
	#---------
	bc1 = phreeqcrm.IntVector(1,1)
	out = bmi.InitialPhreeqc2Concentrations(bc1)
	status = out[0]
	print(type(out[1]))
	bc_conc = list(out[1])
	tc = [30.0]
	p_atm = [1.5]
	bc_double_vect = phreeqcrm.DoubleVector(len(bc_conc))
	for i in range(len(bc_conc)):
		bc_double_vect[i] = bc_conc[i]
	x=bmi.Concentrations2Utility(bc_double_vect, tc, p_atm)
	print(f"Concentrations2Utility {x}")
	#---------
	grid2chem = phreeqcrm.IntVector(nxyz, -1)
	for i in range(nxyz//2):
		grid2chem[i] = i
	x=bmi.CreateMapping(grid2chem)
	print(f"CreateMapping {x}")
	#---------
	bmi.DecodeError(-2)	         # void function
	print(f"DecodeError {x}")
	#---------
	x=bmi.DumpModule(True)
	print(f"DumpModule {x}")
	#---------
	bmi.ErrorHandler(0, "string") # void function
	print(f"OK, just a test: ErrorHandler {x}")
	#---------
	bmi.ErrorMessage("my error")  # void function
	print(f"OK, just a test: ErrorMessage {x}")
	#---------
	bmi.Update()      # void function
	print(f"Update")
	#---------
	x=bmi.GetBackwardMapping()
	print(type(x))
	print("GetBackwardMapping what is it?")
	#---------
	x=bmi.GetChemistryCellCount()
	print(f"GridCellCount {x}")
	#---------
	x=bmi.GetComponentCount()
	print(f"GetComponentCount {x}")
	#---------
	x=bmi.GetComponents()
	print(f"GetComponents {x}")
	#---------
	c = phreeqcrm.DoubleVector()
	x=bmi.GetConcentrations(c)
	print(f"GetConcentrations {x}, {c[0]}.")
	#---------
	x=bmi.GetCurrentSelectedOutputUserNumber()
	print(f"GetCurrentSelectedOutputUserNumber {x}")
	#---------
	x=bmi.GetDatabaseFileName()
	print(f"GetDatabaseFileName {x}")
	#---------
	d = phreeqcrm.DoubleVector(nxyz, 0)
	x=bmi.GetDensity(d)
	print(f"GetDensity {x}, {d[0]}") 
	#---------
	x=bmi.GetEndCell()
	print(f"GetEndCell, FAILS WITH THREADS {x}")
	#---------
	x=bmi.GetEquilibriumPhases()
	print(f"GetEquilibriumPhases {x}")
	#---------
	x=bmi.GetEquilibriumPhasesCount()
	print(f"GetEquilibriumPhasesCount {x}")
	#---------
	x=bmi.GetErrorHandlerMode()
	print(f"GetErrorHandlerMode {x}")
	#---------
	x=bmi.GetErrorString()
	print(f"GetErrorString {x}")
	#---------
	x=bmi.GetExchangeNames()
	print(f"GetExchangeNames {x}")
	#---------
	x=bmi.GetExchangeSpecies()
	print(f"GetExchangeSpecies {x}")
	#---------
	x=bmi.GetExchangeSpeciesCount()
	print(f"GetExchangeSpeciesCount {x}")
	#---------
	x=bmi.GetFilePrefix()
	print(f"GetFilePrefix {x}")
	#---------
	x=bmi.GetForwardMapping()
	print(f"GetForwardMapping {x[0]}")
	#---------
	x=bmi.GetGasComponents()
	print(f"GetGasComponents {x}")
	#---------
	x=bmi.GetGasComponentsCount()
	print(f"GetGasComponentsCount {x}")
	#---------
	v = phreeqcrm.DoubleVector()
	x=bmi.GetGasCompMoles(v)
	print(f"GetGasCompMoles {x}, {v[0]}")
	#---------
	x=bmi.GetGasCompPressures(v)
	print(f"GetGasCompPressures {x}, {v[0]}")
	#---------
	x=bmi.GetGasCompPhi(v)
	print(f"GetGasCompPhi {x}, {v[0]}")
	#---------
	x=bmi.GetGasPhaseVolume(v)
	print(f"GetGasPhaseVolume {x}, {v[0]}")
	#---------
	x=bmi.GetGfw()
	print(f"GetGfw {x[0]}")
	#---------
	x=bmi.GetGridCellCount()
	print(f"GetGridCellCount {x}")
	#---------
	x=bmi.GetIPhreeqcPointer(0)
	print(f"GetIPhreeqcPointer {x}")
	#---------
	x=bmi.GetKineticReactions()
	print(f"GetKineticReactions {x}")
	#---------
	x=bmi.GetKineticReactionsCount()
	print(f"GetKineticReactionsCount {x}")
	#---------
	x=bmi.GetMpiMyself()
	print(f"GetMpiMyself {x}")
	#---------
	x=bmi.GetMpiTasks()
	print(f"GetMpiTasks {x}")
	#---------
	x=bmi.GetNthSelectedOutputUserNumber(0)
	print(f"GetNthSelectedOutputUserNumber {x}")
	#---------
	x=bmi.GetPartitionUZSolids()
	print(f"GetPartitionUZSolids {x}")
	#---------
	x=bmi.GetPorosity()
	print(f"GetPorosity {x[0]}")
	#---------
	x=bmi.GetPressure()
	print(f"GetPressure {x[0]}")
	#---------
	x=bmi.GetPrintChemistryMask()
	print(f"GetPrintChemistryMask {x[0]}")
	#---------
	x=bmi.GetPrintChemistryOn()
	print(type(x))
	print(f"What is it? GetPrintChemistryOn {x}")
	#---------
	x=bmi.GetRebalanceByCell()
	print(f"GetRebalanceByCell {x}")
	#---------
	x=bmi.GetRebalanceFraction()
	print(f"GetRebalanceFraction {x}")
	#---------
	x=bmi.GetSaturation(v)
	print(f"GetSaturation {x}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutput(v)
	print(f"GetSelectedOutput {x}, {v[0]}")
	#---------
	x=bmi.GetSelectedOutputColumnCount()
	print(f"GetSelectedOutputColumnCount {x}")
	#---------
	v = phreeqcrm.StringVector()
	x=bmi.GetSelectedOutputHeadings(v)
	print(f"GetSelectedOutputHeadings {x}, {v[87]}")
	#---------
	x=bmi.GetSelectedOutputOn()
	print(f"GetSelectedOutputOn {x}")
	#---------
	x=bmi.GetSelectedOutputRowCount()
	print(f"GetSelectedOutputRowCount {x}")
	#---------
	x=bmi.GetSICount()
	print(f"GetSICount {x}")
	#---------
	x=bmi.GetSINames()
	print(f"GetSINames {x[0]}")
	#---------
	x=bmi.GetSolidSolutionComponents()
	print(f"GetSolidSolutionComponents {x}")
	#---------
	x=bmi.GetSolidSolutionComponentsCount()
	print(f"GetSolidSolutionComponentsCount {x}")
	#---------
	x=bmi.GetSolidSolutionNames()
	print(f"GetSolidSolutionNames {x}")
	#---------
	x=bmi.GetSolutionVolume()
	print(f"GetSolutionVolume {x[0]}")
	#---------
	v = phreeqcrm.DoubleVector()
	x=bmi.GetSpeciesConcentrations(v)
	print(f"GetSpeciesConcentrations {v[0]}")
	#---------
	x=bmi.GetSpeciesCount()
	print(f"GetSpeciesCount {x}")
	#---------
	x=bmi.GetSpeciesD25()
	print(f"GetSpeciesD25 {x[0]}")
	#---------
	x=bmi.GetSpeciesLog10Gammas(v)
	print(f"GetSpeciesLog10Gammas {x}, {v[0]}")
	#---------
	x=bmi.GetSpeciesLog10Molalities(v)
	print(f"GetSpeciesLog10Molalities {x}, {v[0]}")
	#---------
	x=bmi.GetSpeciesNames()
	print(f"GetSpeciesNames {x[0]}")
	#---------
	x=bmi.GetSpeciesSaveOn()
	print(f"GetSpeciesSaveOn {x}")
	#---------
	x=bmi.GetSpeciesStoichiometry()
	print(f"What is it? GetSpeciesStoichiometry {x}")
	print(type(x))
	#---------
	x=bmi.GetSpeciesZ()
	print(f"GetSpeciesZ {x[0]}")
	#---------
	x=bmi.GetStartCell()
	print(f"GetStartCell, FAILS WITH THREADS {x}")
	#---------
	x=bmi.GetSurfaceNames()
	print(f"GetSurfaceNames {x[0]}")
	#---------
	x=bmi.GetSurfaceSpecies()
	print(f"GetSurfaceSpecies {x[0]}")
	#---------
	x=bmi.GetSurfaceSpeciesCount()
	print(f"GetSurfaceSpeciesCount {x}")
	#---------
	x=bmi.GetSurfaceTypes()
	print(f"GetSurfaceTypes {x[0]}")
	#---------
	x=bmi.GetTemperature()
	print(f"GetTemperature {x[0]}")
	#---------
	x=bmi.GetThreadCount()
	print(f"Wrong: GetThreadCount {x}")
	#---------
	x=bmi.GetTime()
	print(f"GetTime {x}")
	#---------
	x=bmi.GetTimeConversion()
	print(f"GetTimeConversion {x}")
	#---------
	x=bmi.GetTimeStep()
	print(f"GetTimeStep {x}")
	#---------
	x=bmi.GetUnitsExchange()
	print(f"GetUnitsExchange {x}")
	#---------
	x=bmi.GetUnitsGasPhase()
	print(f"GetUnitsGasPhase {x}")
	#---------
	x=bmi.GetUnitsKinetics()
	print(f"GetUnitsKinetics {x}")
	#---------
	x=bmi.GetUnitsPPassemblage()
	print(f"GetUnitsPPassemblage {x}")
	#---------
	x=bmi.GetUnitsSolution()
	print(f"GetUnitsSolution {x}")
	#---------
	x=bmi.GetUnitsSSassemblage()
	print(f"GetUnitsSSassemblage {x}")
	#---------
	x=bmi.GetUnitsSurface()
	print(f"GetUnitsSurface {x}")
	#---------
	x=bmi.GetWorkers()
	print(f"What is it? GetWorkers {x}")
	#---------
	bmi.LogMessage("Log message")  # void method
	print(f"LogMessage")
	#---------
	#TODO x=bmi.MpiAbort()
	#TODO x=bmi.MpiWorker()
	#TODO x=bmi.MpiWorkerBreak()
	bmi.OutputMessage("Output message")  # void method
	print(f"OutputMessage")
	#---------
	#Should be private: x=bmi.ReturnHandler()
	bmi.ScreenMessage("Screen message\n")  # void method
	print(f"ScreenMessage")
	#---------
	v = bmi.GetValue("Concentrations")
	print(type(v))
	x=bmi.SetConcentrations(v)
	print(f"SetConcentrations {x}")
	#---------
	v = phreeqcrm.DoubleVector()
	x=bmi.GetSpeciesConcentrations(v)
	x=bmi.SpeciesConcentrations2Module(v)
	print(f"SpeciesConcentrations2Module {x}")
	#---------
	x=bmi.StateSave(1)
	print(f"StateSave {x}")
	#---------
	x=bmi.StateApply(1)
	print(f"StateApply {x}")
	#---------
	x=bmi.StateDelete(1)
	print(f"StateDelete {x}")
	#---------
	bmi.WarningMessage("Warning message")  # void method
	print(f"WarningMessage")
	#---------
	# BMI methods, some have been used above
	x=bmi.GetComponentName()
	print(f"GetComponentName {x}")
	#---------
	x=bmi.GetCurrentTime()
	print(f"GetCurrentTime {x}")
	#---------
	x=bmi.GetEndTime()
	print(f"GetEndTime {x}")
	#---------
	x=bmi.GetInputItemCount()
	print(f"GetInputItemCount {x}")
	#---------
	x=bmi.GetInputVarNames()
	print(f"GetInputVarNames {x[0]}")
	#---------
	x=bmi.GetOutputItemCount()
	print(f"GetOutputItemCount {x}")
	#---------
	x=bmi.GetOutputVarNames()
	print(f"GetOutputVarNames {x[0]}")
	#---------
	x=bmi.GetTimeStep()
	print(f"GetTimeStep {x}")
	#---------
	x=bmi.GetTimeUnits()
	print(f"GetTimeUnits {x}")
	#---------
	x=bmi.GetValue("solution_saturation_index_Calcite")
	print(f"GetValue {x[0]}")
	#---------
	x=bmi.GetVarItemsize("solution_saturation_index_Calcite")
	print(f"GetVarItemsize {x}")
	#---------
	x=bmi.GetVarNbytes("solution_saturation_index_Calcite")
	print(f"GetVarNbytes {x}")
	#---------
	x=bmi.GetVarType("solution_saturation_index_Calcite")
	print(f"GetVarType {x}")
	#---------
	x=bmi.GetVarUnits("solution_saturation_index_Calcite")
	print(f"GetVarUnits {x}")
	#---------
	#x=bmi.Initialize("AdvectBMI_test_py.yaml")
	# See above
	bmi.SetValue("Time", 1.0)    # void method
	print(f"SetValue")
	#---------
	bmi.Update()    # void method
	print(f"Update")
	#---------
	x=bmi.CloseFiles()
	print(f"CloseFiles {x}")
	#---------
	bmi.Finalize()    # void method
	print(f"Finalize {x}")
	#---------
	print("Success.")
	return

if __name__ == '__main__':
	testbmi_py()

#endif # USE_YAML
