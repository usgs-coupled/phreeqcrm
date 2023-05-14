#ifdef USE_YAML
subroutine TestAllMethods_f90()  BIND(C, NAME='TestAllMethods_f90')
  USE, intrinsic :: ISO_C_BINDING
  USE PhreeqcRM
  USE IPhreeqc
  USE YAMLPreeqcRM
  implicit none
#ifdef USE_MPI
  INCLUDE 'mpif.h'
#endif

  ! Based on PHREEQC Example 11
  integer :: mpi_myself
  integer                      :: i, j
  integer                      :: id, nxyz, nthreads, nchem, ncomps, status
  integer                      :: nbound, isteps, nsteps
  character(100)               :: string
  character(len:), allocatable :: StringVector(:)
  real(kin=8), allocatable     :: IntVector(:)
  integer, allocatable         :: DoubleVector(:), f1(:), temperature(:)
  integer, allocatable         :: ic1(:,:), ic2(:,:), bc_conc(:,:), c(:,:)
  integer, allocatable         :: bc1(:)
  real(kind=8)                 :: time, time_step, t
  ! --------------------------------------------------------------------------
  ! Create PhreeqcRM
  ! --------------------------------------------------------------------------
#ifdef SKIP
    ! Write YAML file
    id = CreateYAMLPhreeqcRM()
    nxyz = 40
    status = YAMLSetGridCellCount(id, nxyz)
    status = YAMLThreadCount(id, 3)
	string = "TestAllMethods_cpp.yaml"
	string = WriteYAMLDoc(id, string)
	status = YAMLClear(id)
    status = DestroyYAMLPhreeqcRM(id)
	!
	! Use all BMIPhreeqcRM methods roughly in order of use
	!
	!
	id = RM_create()
	write(*,*) "RM_create"
	!-------
	status = RM_Initialize(id, string)
	write(*,*) "RM_Initialize"
	!-------
	status = bmif_getvalue(id, "GridCellCount", nxyz)
	write(*,*) "bmif_getvalue('GridCellCount')"
	!-------
	nxyz = RM_GetGridCellCount(id)
	write(*,*) "RM_GetGridCellCount "
	!-------
	int n = RM_GetThreadCount(id)
	write(*,*) "RM_GetThreadCount " << n << ""
	! Inactive cells or symmetry
    deallocate(IntVector)
	allocate(IntVector(nxyz, -1)
	do i = 1, nxyz / 2 
        IntVector(i) = i
    enddo
	status = RM_CreateMapping(id, IntVector)
	write(*,*) "RM_CreateMapping "
	!-------
	status = RM_LoadDatabase(id, "phreeqc.dat")
	write(*,*) "RM_LoadDatabase"
	!
	! Set properties
	!
	!-------
	status = RM_SetComponentH2O(id, .false.)
	write(*,*) "RM_SetComponentH2O "
	!-------
	status = RM_SetSpeciesSaveOn(id, .true.)
	write(*,*) "RM_SetSpeciesSaveOn "
	!-------
	status =RM_SetErrorOn(id, .true.)
	write(*,*) "RM_SetErrorOn "
	!-------
	status =RM_SetErrorHandlerMode(id, 1)
	write(*,*) "RM_SetErrorHandlerMode "
	!-------
	status =RM_SetDumpFileName(id, "TestAllMethods_cpp.dump")
	write(*,*) "RM_SetDumpFileName "
	!-------
	status =RM_SetFilePrefix(id, "TestAllMethods_cpp")
	write(*,*) "RM_SetFilePrefix "
	!-------
	status =RM_OpenFiles(id)
	write(*,*) "RM_OpenFiles "
	!-------
	status =RM_SetPartitionUZSolids(id, .false.)
	write(*,*) "RM_SetPartitionUZSolids "
	!-------
	status =RM_SetRebalanceByCell(id, .true.)
	write(*,*) "RM_SetRebalanceByCell "
	!-------
	status =RM_SetRebalanceFraction(id, 0.5)
	write(*,*) "RM_SetRebalanceFraction "
	!-------
	status =RM_SetScreenOn(id, .true.)
	write(*,*) "RM_SetScreenOn "
	!-------
	status =RM_SetSelectedOutputOn(id, .true.)
	write(*,*) "RM_SetSelectedOutputOn "
	!-------
	status = RM_SetUnitsExchange(id, 1)
	write(*,*) "RM_SetUnitsExchange "
	!-------
	status =RM_SetUnitsGasPhase(id, 1)
	write(*,*) "RM_SetUnitsGasPhase "
	!-------
	status =RM_SetUnitsKinetics(id, 1)
	write(*,*) "RM_SetUnitsKinetics "
	!-------
	status =RM_SetUnitsPPassemblage(id, 1)
	write(*,*) "RM_SetUnitsPPassemblage "
	!-------
	status =RM_SetUnitsSolution(id, 2)
	write(*,*) "RM_SetUnitsSolution "
	!-------
	status =RM_SetUnitsSSassemblage(id, 1)
	write(*,*) "RM_SetUnitsSSassemblage "
	!-------
	status =RM_SetUnitsSurface(id, 1)
	write(*,*) "RM_SetUnitsSurface "
	!-------
	status = RM_UseSolutionDensityVolume(id, .false.)
	write(*,*) "RM_UseSolutionDensityVolume "
	!-------
	t = 1.0 / 86400.0
	status = RM_SetTimeConversion(id, t)
	write(*,*) "RM_SetTimeConversion "
	!-------
	allocate(DoubleVector(nxyz, 1.0))
	status = RM_SetRepresentativeVolume(id, DoubleVector)
	write(*,*) "RM_SetRepresentativeVolume "

	!-------Chemistry cells may be fewer than GridCellCount
	nchem = RM_GetChemistryCellCount(id)
	write(*,*) "RM_GetChemistryCellCount "
	!
	! Begin initialization
	!
	status =RM_SetPrintChemistryOn(id, .false., .true., .false.)
	write(*,*) "RM_SetPrintChemistryOn "
	!-------
	status =RM_RunFile(.true., .true., .true., "all_reactants.pqi")
	write(*,*) "RunFile "
	!-------
	!status = RM_AddOutputVars("SolutionActivities", ".true.")
	!status = RM_AddOutputVars("SolutionMolalities", ".true.")
	!status = RM_AddOutputVars("SaturationIndices", ".true.")
	status = RM_AddOutputVars("SolutionActivities", "H+ Ca+2 Na+")
	status = RM_AddOutputVars("SolutionMolalities", "OH- Cl-")
	status = RM_AddOutputVars("SaturationIndices", "Calcite Dolomite")
	write(*,*) "AddOutputVars "
	!-------
	int ncomps = status = RM_FindComponents()
	write(*,*) "FindComponents "
	! Component names
	ncomps = status = RM_GetComponentCount()
	write(*,*) "GetComponentCount "
	!-------
	status = RM_GetValue("ComponentCount", ncomps)
	write(*,*) "GetValue('ComponentCount')" << ncomps << ""
	!-------
	std::vector<std::string> str_vector = status = RM_GetComponents()
	write(*,*) "GetComponents "
	! Species info
	n = status = RM_GetSpeciesCount()
	write(*,*) "GetSpeciesCount "
	!-------
	str_vector = status = RM_GetSpeciesNames()
	write(*,*) "GetSpeciesNames "
	!-------
	v = status = RM_GetSpeciesD25()
	write(*,*) "GetSpeciesD25 "
	!-------
	v = status = RM_GetSpeciesZ()
	write(*,*) "GetSpeciesZ "
	! Reactant lists
	std::vector<std::string> equiuilibrium_phases = status = RM_GetEquilibriumPhases()
	write(*,*) "GetEquilibriumPhases "
	!-------
	n = status = RM_GetEquilibriumPhasesCount()
	write(*,*) "GetEquilibriumPhasesCount "
	!-------
	str_vector = status = RM_GetExchangeNames()
	write(*,*) "GetExchangeNames "
	!-------
	str_vector = status = RM_GetExchangeSpecies()
	write(*,*) "GetExchangeSpecies "
	!-------
	n = status = RM_GetExchangeSpeciesCount()
	write(*,*) "GetExchangeSpeciesCount "
	!-------
	str_vector = status = RM_GetGasComponents()
	write(*,*) "GetGasComponents "
	!-------
	n = status = RM_GetGasComponentsCount()
	write(*,*) "GetGasComponentsCount "
	!-------
	str_vector = status = RM_GetKineticReactions()
	write(*,*) "GetKineticReactions "
	!-------
	n = status = RM_GetKineticReactionsCount()
	write(*,*) "GetKineticReactionsCount "
	!-------
	n = status = RM_GetSICount()
	write(*,*) "GetSICount "
	!-------
	str_vector = status = RM_GetSINames()
	write(*,*) "GetSINames "
	!-------
	str_vector = status = RM_GetSolidSolutionComponents()
	write(*,*) "GetSolidSolutionComponents "
	!-------
	n = status = RM_GetSolidSolutionComponentsCount()
	write(*,*) "GetSolidSolutionComponentsCount "
	!-------
	str_vector = status = RM_GetSolidSolutionNames()
	write(*,*) "GetSolidSolutionNames "
	!-------
	str_vector = status = RM_GetSurfaceNames()
	write(*,*) "GetSurfaceNames "
	!-------
	str_vector = status = RM_GetSurfaceSpecies()
	write(*,*) "GetSurfaceSpecies "
	!-------
	n = status = RM_GetSurfaceSpeciesCount()
	write(*,*) "GetSurfaceSpeciesCount "
	!-------
	str_vector = status = RM_GetSurfaceTypes()
	write(*,*) "GetSurfaceTypes "
	!
	! Remove any reactants in workers
	!
	std::string input = "DELETE -all"
	status = RM_RunString(.true., .false., .false., input)
	write(*,*) "RunString "
	!-------
	!
	! Transfer initial conditions
	!
	std::vector<int> vi(nxyz, 1)
	status =RM_InitialEquilibriumPhases2Module(vi)
	write(*,*) "InitialEquilibriumPhases2Module "
	!-------
	status =RM_InitialExchanges2Module(vi)
	write(*,*) "InitialExchanges2Module "
	!-------
	status =RM_InitialGasPhases2Module(vi)
	write(*,*) "InitialGasPhases2Module "
	!-------
	status =RM_InitialKinetics2Module(vi)
	write(*,*) "InitialKinetics2Module "
	!-------
	status =RM_InitialSolutions2Module(vi)
	write(*,*) "InitialSolutions2Module "
	!-------
	status =RM_InitialSolidSolutions2Module(vi)
	write(*,*) "InitialSolidSolutions2Module "
	!-------
	status =RM_InitialSurfaces2Module(vi)
	write(*,*) "InitialSurfaces2Module "
	!-------
	! Alternative A.to the previous seven methods
	std::vector<int> ic(nxyz * 7, 1)
	status =RM_InitialPhreeqc2Module(ic)
	write(*,*) "InitialPhreeqc2Module "
	!-------
	! Alternative B.to the previous seven methods, possible mixing
	std::vector<int> v1(nxyz * 7, 1)
	std::vector<int> v2(nxyz * 7, -1)
	std::vector<double> f1(nxyz * 7, 1.0)
	status =RM_InitialPhreeqc2Module(v1, v2, f1)
	write(*,*) "InitialPhreeqc2Modul mix "
	!-------
	! Alternative C.to the previous seven methods, initialize cells 18 and 19
	std::vector<int> cells(2)
	cells(0) = 18
	cells(1) = 19
	status =RM_InitialPhreeqcCell2Module(1, cells)
	write(*,*) "InitialPhreeqcCell2Module "
	!
	! Boundary conditions
	!
	std::vector<double> bc
	vi.resize(1, 1)
	status =RM_InitialPhreeqc2SpeciesConcentrations(bc, vi)
	write(*,*) "InitialPhreeqc2SpeciesConcentrations "
	!-------
	std::vector<double> bc_species
	std::vector<int> vi1(1, 1)
	std::vector<int> vi2(1, -1)
	f1.resize(1, 1.0)
	status =RM_InitialPhreeqc2SpeciesConcentrations(bc_species, vi1, vi2, f1)
	write(*,*) "InitialPhreeqc2SpeciesConcentrations mix "
	!
	! Get/Set methods for time steping
	!
	status =RM_SetTime(0.0)
	write(*,*) "SetTime "
	!-------
	status =RM_SetTimeStep(0.0)
	write(*,*) "SetTimeStep "
	!-------
	status =RM_GetGasCompMoles(v)
	write(*,*) "GetGasCompMoles "
	!-------
	status =RM_SetGasCompMoles(v)
	write(*,*) "SetGasCompMoles "
	!-------
	vi.resize(nxyz, 1)
	status =RM_SetPrintChemistryMask(vi)
	write(*,*) "SetPrintChemistryMask "
	!
	! Get/Set methods for time stepping
	!
	std::vector<double> c
	status =RM_GetConcentrations(c)
	write(*,*) "GetConcentrations "
	!-------
	status =RM_SetConcentrations(c)
	write(*,*) "SetConcentrations "
	!-------
	v.resize(nxyz, 0)
	status =RM_GetDensity(v)
	write(*,*) "GetDensity "
	!-------
	v.resize(nxyz, 1.1)
	status =RM_SetDensity(v)
	write(*,*) "SetDensity "
	!-------
	status =RM_GetGasCompMoles(v)
	write(*,*) "GetGasCompMoles "
	!-------
	status =RM_SetGasCompMoles(v)
	write(*,*) "GetGasCompMoles "
	!-------
	status =RM_GetGasCompPressures(v)
	write(*,*) "GetGasCompPressures "
	!-------
	status =RM_GetGasPhaseVolume(v)
	write(*,*) "GetGasPhaseVolume "
	!-------
	v.resize(nxyz, 1.0)
	status =RM_SetGasPhaseVolume(v)
	write(*,*) "SetGasPhaseVolume "
	!-------
	v.resize(nxyz, 0.21)
	status =RM_SetPorosity(v)
	write(*,*) "SetPorosity "
	!-------
	v = status = RM_GetPressure()
	write(*,*) "GetPressure "
	!-------
	v.resize(nxyz, 3.0)
	status =RM_SetPressure(v)
	write(*,*) "SetPressure "
	!-------
	v.resize(nxyz, 1.0)
	status =RM_SetSaturation(v)
	write(*,*) "SetSaturation "
	!-------
	v = status = RM_GetSolutionVolume()
	write(*,*) "GetSolutionVolume "
	!-------
	status =RM_GetSpeciesConcentrations(v)
	write(*,*) "GetSpeciesConcentrations "
	!-------
	status =RM_SpeciesConcentrations2Module(v)
	write(*,*) "SpeciesConcentrations2Module "
	!-------
	status =RM_GetSpeciesLog10Gammas(v)
	write(*,*) "GetSpeciesLog10Gammas "
	!-------
	status =RM_GetSpeciesLog10Molalities(v)
	write(*,*) "GetSpeciesLog10Molalities "
	!-------
	v.resize(nxyz, 26.0)
	status =RM_SetTemperature(v)
	write(*,*) "SetTemperature "
	!
	! Take a time step
	!
	status = RM_Update()      ! void function
	write(*,*) "Update"
	!-------
	status =RM_RunCells()
	write(*,*) "RunCells"
	!
	! Selected output
	!
	status =RM_SetNthSelectedOutput(0)
	write(*,*) "SetNthSelectedOutput "
	!-------
	int n_user = status = RM_GetCurrentSelectedOutputUserNumber()
	write(*,*) "GetCurrentSelectedOutputUserNumber "
	!-------
	status =RM_SetCurrentSelectedOutputUserNumber(333)
	write(*,*) "SetCurrentSelectedOutputUserNumber "
	!-------
	n = status = RM_GetNthSelectedOutputUserNumber(0)
	write(*,*) "GetNthSelectedOutputUserNumber "
	!-------
	status =RM_GetSelectedOutput(v)
	write(*,*) "GetSelectedOutput "
	!-------
	n = status = RM_GetSelectedOutputColumnCount()
	write(*,*) "GetSelectedOutputColumnCount "
	!-------
	n = status = RM_GetSelectedOutputRowCount()
	write(*,*) "GetSelectedOutputRowCount "
	!
	! Getters
	!
	std::vector< std::vector<int> > back_map = status = RM_GetBackwardMapping()
	write(*,*) "GetBackwardMapping "
	!-------
	std::string db_name = status = RM_GetDatabaseFileName()
	write(*,*) "GetDatabaseFileName "
	!-------
	vi = status = RM_GetEndCell()
	write(*,*) "GetEndCell"
	!-------
	n = status = RM_GetErrorHandlerMode()
	write(*,*) "GetErrorHandlerMode "
	!-------
	std::string str = status = RM_GetErrorString()
	write(*,*) "GetErrorString "
	!-------
	str = status = RM_GetFilePrefix()
	write(*,*) "GetFilePrefix "
	!-------
	vi = status = RM_GetForwardMapping()
	write(*,*) "GetForwardMapping "
	!-------
	status =RM_GetGasCompPhi(v)
	write(*,*) "GetGasCompPhi "
	!-------
	v = status = RM_GetGfw()
	write(*,*) "GetGfw "
	!-------
	IPhreeqc* ipq = status = RM_GetIPhreeqcPointer(0)
	write(*,*) "GetIPhreeqcPointer "
	!-------
	n = status = RM_GetMpiMyself()
	write(*,*) "GetMpiMyself "
	!-------
	n = status = RM_GetMpiTasks()
	write(*,*) "GetMpiTasks "
	!-------
	bool b = status = RM_GetPartitionUZSolids()
	write(*,*) "GetPartitionUZSolids "
	!-------
	v = status = RM_GetPorosity()
	write(*,*) "GetPorosity "
	!-------
	vi = status = RM_GetPrintChemistryMask()
	write(*,*) "GetPrintChemistryMask "
	!-------
	std::vector<bool> vb = status = RM_GetPrintChemistryOn()
	write(*,*) "GetPrintChemistryOn "
	!-------
	b = status = RM_GetRebalanceByCell()
	write(*,*) "GetRebalanceByCell "
	!-------
	double d = status = RM_GetRebalanceFraction()
	write(*,*) "GetRebalanceFraction "
	!-------
	status =RM_GetSaturation(v)
	write(*,*) "GetSaturation "
	!-------
	status =RM_GetSelectedOutputHeadings(str_vector)
	write(*,*) "GetSelectedOutputHeadings "
	!-------
	b = status = RM_GetSelectedOutputOn()
	write(*,*) "GetSelectedOutputOn "
	!-------
	b = status = RM_GetSpeciesSaveOn()
	write(*,*) "GetSpeciesSaveOn "
	!-------
	std::vector< cxxNameDouble > s = status = RM_GetSpeciesStoichiometry()
	write(*,*) "GetSpeciesStoichiometry "
	!-------
	vi = status = RM_GetStartCell()
	write(*,*) "GetStartCell "
	!-------
	v = status = RM_GetTemperature()
	write(*,*) "GetTemperature "
	!-------
	d = status = RM_GetTime()
	write(*,*) "GetTime "
	!-------
	d = status = RM_GetTimeConversion()
	write(*,*) "GetTimeConversion "
	!-------
	d = status = RM_GetTimeStep()
	write(*,*) "GetTimeStep "
	!-------
	n = status = RM_GetUnitsExchange()
	write(*,*) "GetUnitsExchange "
	!-------
	n = status = RM_GetUnitsGasPhase()
	write(*,*) "GetUnitsGasPhase "
	!-------
	n = status = RM_GetUnitsKinetics()
	write(*,*) "GetUnitsKinetics "
	!-------
	n = status = RM_GetUnitsPPassemblage()
	write(*,*) "GetUnitsPPassemblage "
	!-------
	n = status = RM_GetUnitsSolution()
	write(*,*) "GetUnitsSolution "
	n = status = RM_GetUnitsSSassemblage()
	write(*,*) "GetUnitsSSassemblage "
	!-------
	n = status = RM_GetUnitsSurface()
	write(*,*) "GetUnitsSurface "
	!-------
	std::vector<IPhreeqcPhast *> w = status = RM_GetWorkers()
	write(*,*) "GetWorkers "
	!
	! Utilities
	!
	vi.resize(1, 1)
	status =RM_InitialPhreeqc2Concentrations(bc, vi)
	std::vector<double> tc(1, 30.0)
	std::vector<double> p_atm(1, 1.5)
	IPhreeqc* utility_ptr = status = RM_Concentrations2Utility(bc, tc, p_atm)
	write(*,*) "Concentrations2Utility "
	!-------
	status = RM_DecodeError(-2)	         ! void function
	write(*,*) "DecodeError "
	!-------
	status =RM_DumpModule(.true.)
	write(*,*) "DumpModule "
	!-------
	status = RM_ErrorHandler(0, "string") ! void function
	write(*,*) "OK, just a test: ErrorHandler "
	!-------
	status = RM_ErrorMessage("my error")  ! void function
	write(*,*) "OK, just a test: ErrorMessage "
	!-------
	status = RM_LogMessage("Log message")  ! void method
	write(*,*) "LogMessage "
	!-------
	status = RM_OutputMessage("Output message")  ! void method
	write(*,*) "OutputMessage "
	!-------
	status = RM_ScreenMessage("Screen message")  ! void method
	write(*,*) "ScreenMessage "
	!-------
	status =RM_StateSave(1)
	write(*,*) "StateSave "
	!-------
	status =RM_StateApply(1)
	write(*,*) "StateApply "
	!-------
	status =RM_StateDelete(1)
	write(*,*) "StateDelete "
	!-------
	status = RM_WarningMessage("Warning message")  ! void method
	write(*,*) "WarningMessage "
	!-------
	! BMI Methods
	str = status = RM_GetComponentName()
	write(*,*) "GetComponentName "
	!-------
	d = status = RM_GetCurrentTime()
	write(*,*) "GetCurrentTime "
	!-------
	d = status = RM_GetEndTime()
	write(*,*) "GetEndTime "
	!-------
	n = status = RM_GetInputItemCount()
	write(*,*) "GetInputItemCount "
	!-------
	str_vector = status = RM_GetInputVarNames()
	write(*,*) "GetInputVarNames "
	!-------
	n = status = RM_GetOutputItemCount()
	write(*,*) "GetOutputItemCount "
	!-------
	str_vector = status = RM_GetOutputVarNames()
	write(*,*) "GetOutputVarNames "
	!-------
	d = status = RM_GetTimeStep()
	write(*,*) "GetTimeStep "
	!-------
	str = status = RM_GetTimeUnits()
	write(*,*) "GetTimeUnits "
	!-------
	status = RM_GetValue("solution_saturation_index_Calcite", v)
	write(*,*) "GetValue "
	!-------
	n = status = RM_GetVarItemsize("solution_saturation_index_Calcite")
	write(*,*) "GetVarItemsize "
	!-------
	n = status = RM_GetVarNbytes("solution_saturation_index_Calcite")
	write(*,*) "GetVarNbytes "
	!-------
	str = status = RM_GetVarType("solution_saturation_index_Calcite")
	write(*,*) "GetVarType "
	!-------
	str = status = RM_GetVarUnits("solution_saturation_index_Calcite")
	write(*,*) "GetVarUnits "
	!status =RM_Initialize(YAML_filename)
	! See above
	status = RM_SetValue("Time", 1.0)    ! void method
	write(*,*) "SetValue"
	!-------
	status = RM_Update()    ! void method
	write(*,*) "Update"
	!-------
	status =RM_CloseFiles() ! not a BMI method, but needs to be last
	write(*,*) "CloseFiles "
	status = RM_Finalize()    ! void method
	write(*,*) "Finalize "
	!Should be private: status =RM_ReturnHandler()
	!TODO status =RM_MpiAbort()
	!TODO status =RM_MpiWorker()
	!TODO status =RM_MpiWorkerBreak()
	!TODO status =RM_SetMpiWorkerCallbackC()
	!TODO status =RM_SetMpiWorkerCallbackCookie()
	write(*,*) "Success."
	return

#endif
  ! Deallocate
  deallocate(por)
  deallocate(print_chemistry_mask)
  deallocate(components)
  deallocate(ic1)
  deallocate(bc1)
  deallocate(bc_conc)
  deallocate(c)
  deallocate(temperature)
  deallocate(pressure)
  return
end subroutine TestAllMethods_f90

#endif ! YAML
