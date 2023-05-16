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
  integer                      :: mpi_myself
  integer                      :: i, j
  integer                      :: id, nxyz, nthreads, nchem, ncomps, status
  integer                      :: nbound, isteps, nsteps
  character(100)               :: string, yaml_filename
  character(len:), allocatable :: StringVector(:)
  real(kin=8), allocatable     :: IntVector(:)
  integer, allocatable         :: DoubleVector(:), f1(:), temperature(:)
  integer, allocatable         :: ic1(:,:), ic2(:,:), bc_conc(:,:), c(:,:)
  integer, allocatable         :: bc1(:), cells(:), v1(:), v2(:)
  real(kind=8)                 :: time, time_step, t
  ! --------------------------------------------------------------------------
  ! Create PhreeqcRM
  ! --------------------------------------------------------------------------
    ! Write YAML file
    id = CreateYAMLPhreeqcRM()
    nxyz = 40
    status = YAMLSetGridCellCount(id, nxyz)
    status = YAMLThreadCount(id, 3)
	yaml_filename = "TestAllMethods_f90.yaml"
	status = WriteYAMLDoc(id, string)
	status = YAMLClear(id)
    status = DestroyYAMLPhreeqcRM(id)
	!
	! Use all BMIPhreeqcRM methods roughly in order of use
	!
	id = bmif_create()
	write(*,*) "bmif_create"
	!-------
	status = bmif_initialize(id, yaml_filename)
	write(*,*) "bmif_initialize"
	!-------
	status = bmif_get_value(id, "GridCellCount", nxyz)
	nxyz = RM_GetGridCellCount();
	write(*,*) "GetGridCellCount"
	!-------
	int n = RM_GetThreadCount(id)
	write(*,*) "GetThreadCount " << n << ""
	!-------
	! Inactive cells or symmetry
	allocate(IntVector(nxyz, -1)
	do i = 1, nxyz / 2 
        IntVector(i) = i
    enddo
	status = RM_CreateMapping(id, IntVector)
	write(*,*) "CreateMapping "
	!-------
	status = RM_LoadDatabase(id, "phreeqc.dat")
	write(*,*) "LoadDatabase"
	!
	! Set properties
	!
	!-------
	status = RM_SetComponentH2O(id, .false.)
	write(*,*) "SetComponentH2O "
	!-------
	status = RM_SetSpeciesSaveOn(id, .true.)
	write(*,*) "SetSpeciesSaveOn "
	!-------
	status = RM_SetErrorOn(id, .true.)
	write(*,*) "SetErrorOn "
	!-------
	status = RM_SetErrorHandlerMode(id, 1)
	write(*,*) "SetErrorHandlerMode "
	!-------
	status = RM_SetDumpFileName(id, "TestAllMethods_py.dump")
	write(*,*) "SetDumpFileName "
	!-------
	status = bmif_set_value(id, TestAllMethods_py)
	status = RM_SetFilePrefix(id, "TestAllMethods_py")
	write(*,*) "SetFilePrefix "
	!-------
	status =RM_OpenFiles(id)
	write(*,*) "OpenFiles "
	!-------
	status =RM_SetPartitionUZSolids(id, .false.)
	write(*,*) "SetPartitionUZSolids "
	!-------
	status =RM_SetRebalanceByCell(id, .true.)
	write(*,*) "SetRebalanceByCell "
	!-------
	status =RM_SetRebalanceFraction(id, 0.5)
	write(*,*) "SetRebalanceFraction "
	!-------
	status =RM_SetScreenOn(id, .true.)
	write(*,*) "SetScreenOn "
	!-------
	status = RM_SetSelectedOutputOn(id, .true.)
	status = bmif_set_value(id, "SelectedOutputOn", .true.)
	write(*,*) "SetSelectedOutputOn "
	!-------
	status = RM_SetUnitsExchange(id, 1)
	write(*,*) "SetUnitsExchange "
	!-------
	status =RM_SetUnitsGasPhase(id, 1)
	write(*,*) "SetUnitsGasPhase "
	!-------
	status =RM_SetUnitsKinetics(id, 1)
	write(*,*) "SetUnitsKinetics "
	!-------
	status =RM_SetUnitsPPassemblage(id, 1)
	write(*,*) "SetUnitsPPassemblage "
	!-------
	status =RM_SetUnitsSolution(id, 2)
	write(*,*) "SetUnitsSolution "
	!-------
	status =RM_SetUnitsSSassemblage(id, 1)
	write(*,*) "SetUnitsSSassemblage "
	!-------
	status =RM_SetUnitsSurface(id, 1)
	write(*,*) "SetUnitsSurface "
	!-------
	status = RM_UseSolutionDensityVolume(id, .false.)
	write(*,*) "UseSolutionDensityVolume "
	!-------
	t = 1.0 / 86400.0
	status = RM_SetTimeConversion(id, t)
	write(*,*) "SetTimeConversion "
	!-------
	allocate(DoubleVector(nxyz, 1.0))
	status = RM_SetRepresentativeVolume(id, DoubleVector)
	write(*,*) "SetRepresentativeVolume "

	!-------Chemistry cells may be fewer than GridCellCount
	deallocate(IntVector)
	allocate(IntVector(nxyz))
	IntVector = 1
	status =RM_SetPrintChemistryMask(id, vi)
	write(*,*) "SetPrintChemistryMask "
	!-------
	status =RM_SetPrintChemistryOn(id, .false., .true., .false.)
	write(*,*) "RM_SetPrintChemistryOn "
	!
	! Define reactants available for initial 
	! and boundary conditions in this file
	!
	status =RM_RunFile(id, .true., .true., .true., "all_reactants.pqi")
	write(*,*) "RunFile "
	!-------
	!status = RM_AddOutputVars(id, "SolutionActivities", ".true.")
	!status = RM_AddOutputVars(id, "SolutionMolalities", ".true.")
	!status = RM_AddOutputVars(id, "SaturationIndices", ".true.")
	status = RM_AddOutputVars(id, "SolutionActivities", "H+ Ca+2 Na+")
	status = RM_AddOutputVars(id, "SolutionMolalities", "OH- Cl-")
	status = RM_AddOutputVars(id, "SaturationIndices", "Calcite Dolomite")
	write(*,*) "AddOutputVars "
	!-------
	int ncomps = RM_FindComponents(id)
	write(*,*) "FindComponents "
	!
	! Methods up to this point are useful 
	! in a YAML initialization file
	! 
	! Lists of reactants found by FindComponents follow
	! 
	nchem = RM_GetChemistryCellCount(id)
	write(*,*) "GetChemistryCellCount "
	!-------
	ncomps = status = RM_GetComponentCount(id)
	status = bmif_get_value(id, "ComponentCount", ncomps)
	i_ptr = bmi.get_value_ptr(id, "ComponentCount");
	write(*,*) "GetComponentCount)" 
	!-------
	status = RM_GetComponents(id, StringVector)
	status = bmif_get_value(id, "Components", StringVector)
	write(*,*) "GetComponents)" 
	! Species info
	status = RM_GetSpeciesCount(id, n)
	write(*,*) "GetSpeciesCount "
	!-------
	status = RM_GetSpeciesName(id, 1, string)
	write(*,*) "GetSpeciesName "
	!-------
	status = RM_GetSpeciesD25(id, DoubleVector)
	write(*,*) "GetSpeciesD25 "
	!-------
	status = RM_GetSpeciesZ(id, DoubleVector)
	write(*,*) "GetSpeciesZ "
	! Reactant lists
	status = RM_GetEquilibriumPhasesName(id, 1, string)
	write(*,*) "GetEquilibriumPhasesName "
	!-------
	n = RM_GetEquilibriumPhasesCount(id)
	write(*,*) "GetEquilibriumPhasesCount "
	!-------
	status = RM_GetExchangeName(id, 1, string)
	write(*,*) "GetExchangeName "
	!-------
	status = RM_GetExchangeSpeciesName(id, 1, string)
	write(*,*) "GetExchangeSpeciesName "
	!-------
	n = RM_GetExchangeSpeciesCount(id)
	write(*,*) "GetExchangeSpeciesCount "
	!-------
	status = RM_GetGasComponentsName(id, 1, string)
	write(*,*) "GetGasComponentsName "
	!-------
	n = RM_GetGasComponentsCount(id)
	write(*,*) "GetGasComponentsCount "
	!-------
	status = RM_GetGfw(id, DoubleVector)
	status = bmif_get_value(id, DoubleVector)
	d_ptr = bmif_get_value_ptr(id, "Gfw", DoubleVector);
	write(*,*) "GetGfw "
	!-------
	n = RM_GetKineticReactionsCount(id)
	write(*,*) "GetKineticReactionsCount "
	!-------
	status = RM_GetKineticReactionsName(id, 1, string)
	write(*,*) "GetKineticReactionsName "
	!-------
	n = RM_GetSICount(id)
	write(*,*) "GetSICount "
	!-------
	str_vector = status = RM_GetSIName(d, 1, string)
	write(*,*) "GetSIName "
	!-------
	n = RM_GetSolidSolutionComponentsCount(id)
	write(*,*) "GetSolidSolutionComponentsCount "
	!-------
	status = RM_GetSolidSolutionComponentsName(id, 1, string)
	write(*,*) "GetSolidSolutionComponentsName "
	!-------
	status = RM_GetSolidSolutionName(id, 1, string)
	write(*,*) "GetSolidSolutionName "
	!-------
	status = RM_GetSurfaceName(id, 1, string)
	write(*,*) "GetSurfaceName "
	!-------
	n = RM_GetSurfaceSpeciesCount(id)
	write(*,*) "GetSurfaceSpeciesCount "
	!-------
	status = RM_GetSurfaceSpeciesName(id, 1, string)
	write(*,*) "GetSurfaceSpeciesName "
	!-------
	status = RM_GetSurfaceType(id, 1, string)
	write(*,*) "GetSurfaceType "
	!
	! Remove any reactants in workers 
	! before populating cells with reactants
	!
	std::string input = "DELETE -all"
	status = RM_RunString(id, .true., .false., .false., input)
	write(*,*) "RunString "
	!-------
	!
	! Transfer initial conditions
	!
	deallocate(IntVector)
	allocate(IntVector(nxyz)
	IntVector = 1
	status = RM_InitialEquilibriumPhases2Module(id, IntVector)
	write(*,*) "InitialEquilibriumPhases2Module "
	!-------
	status =RM_InitialExchanges2Module(id, IntVector)
	write(*,*) "InitialExchanges2Module "
	!-------
	status =RM_InitialGasPhases2Module(id, IntVector)
	write(*,*) "InitialGasPhases2Module "
	!-------
	status =RM_InitialKinetics2Module(id, IntVector)
	write(*,*) "InitialKinetics2Module "
	!-------
	status =RM_InitialSolutions2Module(id, IntVector)
	write(*,*) "InitialSolutions2Module "
	!-------
	status =RM_InitialSolidSolutions2Module(id, IntVector)
	write(*,*) "InitialSolidSolutions2Module "
	!-------
	status = RM_InitialSurfaces2Module(id, IntVector)
	write(*,*) "InitialSurfaces2Module "
	!-------
	! Alternative A.to the previous seven methods
	deallocate(IntVector)
	allocate(IntVector(nxyz * 7))
	IntVector = 1
	status = RM_InitialPhreeqc2Module(id, IntVector)
	write(*,*) "InitialPhreeqc2Module "
	!-------
	! Alternative B.to the previous seven methods, possible mixing
	allocate(v1(nxyz * 7)
	v1 = 1
	allocate(v2(nxyz * 7)
	v2 = -1
	deallocate(f1)
	allocate(f1(nxyz*7))
	f1 = 1.0
	status = RM_InitialPhreeqc2Module(ic, v1, v2, f1)
	write(*,*) "InitialPhreeqc2Modul mix "
	!-------
	! Alternative C.to the previous seven methods, initialize cells 18 and 19
	allocate(cells(2))
	cells(1) = 18
	cells(2) = 19
	status = RM_InitialPhreeqcCell2Module(id, 1, cells, size(cells))
	write(*,*) "InitialPhreeqcCell2Module "
	!
	! Boundary conditions
	!
	deallocate(v1)
	allocate(v1(1))
	v1 = 1
	deallocate(v2)
	allocate(v2(1))
	v2 = -1
	deallocate(f1)
	allocate(f1(1))
	f1 = 1
	status = RM_InitialPhreeqc2Concentrations(id, bc, size(v1), v1, v2, f1)
	write(*,*) "InitialPhreeqc2Concentrations mix "
	!-------
	deallocate(v1)
	allocate(v1(1))
	v1 = 1
	deallocate(v2)
	allocate(v2(1))
	v2 = -1
	deallocate(f1)
	allocate(f1(1))
	f1 = 1
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_species, size(v1), v1, v2, f1)
	write(*,*) "InitialPhreeqc2SpeciesConcentrations mix "
	!
	! Get/Set methods for time steping
	!
	d = status = RM_GetTime(id)
	status = bmif_get_value(id, "Time", d)
	d = bmif_get_current_time(id)
	d = bmif_get_start_time();
	d_ptr = bmif_get_value_ptr(id, "Time")
	write(*,*) "GetTime "
	!-------
	status = RM_SetTime(id, 0.0)
	status = bmif_set_value(id, "Time", 0.0)
	write(*,*) "SetTime "
	!-------
	d = status = RM_GetTimeStep()
	status = bmif_get_value(id, "TimeStep", d)
	d_ptr = bmif_get_value_ptr(id, "TimeStep")
	write(*,*) "GetTimeStep "
	!-------
	status = RM_SetTimeStep(id, 0.0)
	status = bmif_set_value(id, "TimeStep", 0.0)
	write(*,*) "SetTimeStep "
	!-------
	status =RM_GetConcentrations(id, c)
	status = bmif_get_value(id, "Concentrations", c)
	d_ptr = bmif_get_value_ptr(id, "Concentrations")
	write(*,*) "GetConcentrations "
	!-------
	status =RM_SetConcentrations(id, c)
	status = bmif_set_value(id, "Concentrations", c)
	write(*,*) "SetConcentrations "
	!-------
	status = RM_GetDensity(id, v)
	status = bmif_get_value(id, "Density", v)
	d_ptr = bmif_get_value_ptr(id, "Density")
	write(*,*) "GetDensity "
	!-------
	status = RM_SetDensity(id, v)
	status = bmif_set_value(id, "Density", v)
	write(*,*) "SetDensity "
	!-------
	status = RM_GetGasCompMoles(id, v)
	write(*,*) "GetGasCompMoles "
	!-------
	status = RM_SetGasCompMoles(id, v)
	write(*,*) "SetGasCompMoles "
	!-------
	status = RM_GetGasCompPhi(id, v)
	write(*,*) "GetGasCompPhi "
	!-------
	status = RM_GetGasCompPressures(id, v)
	write(*,*) "GetGasCompPressures "
	!-------
	status = RM_GetGasPhaseVolume(id, v)
	write(*,*) "GetGasPhaseVolume "
	!-------
	status =RM_SetGasPhaseVolume(id, v)
	write(*,*) "SetGasPhaseVolume "
	!-------
	status = RM_GetIthConcentration(id, 1, v);
	std::cerr << "GetIthConcentration \n";
	//-------
	status = RM_GetIthSpeciesConcentration(id, 1, v);
	std::cerr << "GetIthSpeciesConcentration \n";
	//-------
	status = RM_GetPorosity(id, v)
	status = bmif_get_value(id, "Porosity", v)
	d_ptr = bmif_get_value_ptr(id, "Porosity")
	write(*,*) "GetPorosity "
	!-------
	status =RM_SetPorosity(id, v)
	status = bmif_set_value(id, "Porosity", v)
	write(*,*) "SetPorosity "
	!-------
	status = RM_GetPressure(id, v)
	status = bmif_get_value(id, "Pressure", v)
	d_ptr = bmif_get_value_ptr(id, "Pressure")
	write(*,*) "GetPressure "
	!-------
	status = RM_SetPressure(id, v)
	status = bmif_set_value(id, "Pressure", v)
	write(*,*) "SetPressure "
	!-------
	status = RM_GetSaturation(id, v)
	status = bmif_get_value(id, "Saturation", v)
	d_ptr = bmif_get_value_ptr(id, "Saturation")
	write(*,*) "GetSaturation "
	!-------
	status = RM_SetSaturation(id, v)
	status = bmif_set_value(id, "Saturation", v)
	write(*,*) "SetSaturation "
	!-------
	status = RM_GetSolutionVolume(id, v)
	status = bmif_get_value(id, "SolutionVolume", v)
	d_ptr = bmif_get_value_ptr(id, "SolutionVolume")
	write(*,*) "GetSolutionVolume "
	!-------
	status = RM_GetSpeciesConcentrations(id, v)
	write(*,*) "GetSpeciesConcentrations "
	!-------
	status = RM_SpeciesConcentrations2Module(id, v)
	write(*,*) "SpeciesConcentrations2Module "
	!-------
	status = RM_GetSpeciesLog10Gammas(id, v)
	write(*,*) "GetSpeciesLog10Gammas "
	!-------
	status = RM_GetSpeciesLog10Molalities(id, v)
	write(*,*) "GetSpeciesLog10Molalities "
	!-------
	status = RM_GetTemperature(id, v)
	status = bmif_get_value(id, "Temperature", v)
	d_ptr = bmif_get_value_ptr(id, "Temperature", v)
	write(*,*) "GetTemperature "
	!-------
	status = RM_SetTemperature(id, v)
	status = bmif_set_value(id, "Temperature", v)
	write(*,*) "SetTemperature "
	!-------
	status = RM_GetViscosity(id, v);
	status = bmif_get_value(id, "Viscosity", v);
	d_ptr = bmif_get_value_ptr(id, "Viscosity");	
	std::cerr << "GetViscosity \n";
	!
	! Take a time step
	!
	status = bmif_update(id)
	write(*,*) "Update"
	!-------
	status =RM_RunCells(id)
	write(*,*) "RunCells"
	!-------
	status = bmif_update_until(id, 86400.)
	write(*,*) "UpdateUntil"
	!
	! Selected output
	!
	status = RM_SetNthSelectedOutput(id, 1)
	status = bmif_set_value(id, "NthSelectedOutput", 1)
	write(*,*) "SetNthSelectedOutput "
	!-------
	n_user = RM_GetCurrentSelectedOutputUserNumber(id)
	status = bmif_get_value(id, "CurrentSelectedOutputUserNumber", n_user)
	write(*,*) "GetCurrentSelectedOutputUserNumber "
	!-------
	n = RM_GetNthSelectedOutputUserNumber(id, 1)
	write(*,*) "GetNthSelectedOutputUserNumber "
	!-------
	status = RM_GetSelectedOutput(id, v)
	status = bmif_get_value(id, "SelectedOutput", v)
	write(*,*) "GetSelectedOutput "
	!-------
	n = RM_GetSelectedOutputColumnCount(id)
	status = bmif_get_value(id, "SelectedOutputColumnCount", n)
	write(*,*) "GetSelectedOutputColumnCount "
	!-------
	n = RM_GetSelectedOutputCount(id)
	status = bmif_get_value(id, "SelectedOutputCount", n)
	write(*,*) "GetSelectedOutputCount "
	!-------
	status = RM_GetSelectedOutputHeadings(id, str_vector)
	status = bmif_get_value(id, "SelectedOutputHeadings", str_vector)
	write(*,*) "GetSelectedOutputHeadings "
	!-------
	b = RM_GetSelectedOutputOn(id)
	status = bmif_get_value(id, "SelectedOutputOn", b)
	b_ptr = bmif_get_value_ptr(id, "SelectedOutputOn");	
	write(*,*) "GetSelectedOutputOn "
	!-------
	n = RM_GetSelectedOutputRowCount(id)
	status = bmif_get_value(id, "SelectedOutputOn", n)
	write(*,*) "GetSelectedOutputRowCount "
	!-------
	status = RM_SetCurrentSelectedOutputUserNumber(id, 333)
	write(*,*) "SetCurrentSelectedOutputUserNumber "
	!
	! Getters
	!
#ifdef SKIP
	status = RM_GetBackwardMapping(id, 1, IntVector, n)
	write(*,*) "GetBackwardMapping "
	!-------
	status = RM_GetDatabaseFileName(id, string)
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
	b = status = RM_GetSpeciesSaveOn()
	write(*,*) "GetSpeciesSaveOn "
	!-------
	std::vector< cxxNameDouble > s = status = RM_GetSpeciesStoichiometry()
	write(*,*) "GetSpeciesStoichiometry "
	!-------
	vi = status = RM_GetStartCell()
	write(*,*) "GetStartCell "
	!-------
	d = status = RM_GetTimeConversion()
	write(*,*) "GetTimeConversion "
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
	deallocate(IntVector)
	allocate(IntVector(1))
	IntVector = 1
	status =RM_InitialPhreeqc2Concentrations(bc, IntVector)
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
	!status = bmif_nitialize(YAML_filename)
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
