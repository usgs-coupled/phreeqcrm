#ifdef USE_YAML
subroutine TestAllMethods_f90()  BIND(C, NAME='TestAllMethods_f90')
  USE, intrinsic :: ISO_C_BINDING
  USE BMIPhreeqcRM
  USE IPhreeqc
  USE YAMLPhreeqcRM
  
  implicit none
#ifdef USE_MPI
  INCLUDE 'mpif.h'
#endif

  ! Based on PHREEQC Example 11
  integer                      :: mpi_myself
  integer                      :: i, j, n
  integer                      :: id, yid, id1, nxyz, nthreads, nchem, ncomps, nspecies
  integer                      :: ngas
  integer                      :: nbound, isteps, nsteps, n_user, status
  character(100)               :: string, yaml_filename, vtype
  character(len=:), allocatable :: Names(:), StringVector(:), AllocString
  integer, allocatable         :: IntVector(:), IntVector2(:,:)
  real(kind=8), allocatable    :: p_atm(:), tc(:)
  real(kind=8), allocatable    :: DoubleVector(:), f1(:), DoubleVector2(:,:), f2(:,:)
  real(kind=8)                 :: d
  integer, allocatable         :: ic1(:,:), ic2(:,:) 
  real(kind=8), allocatable    :: bc1(:), c(:,:), bc_species(:,:), bc2(:,:)
  
  integer, allocatable         :: cells(:), v1(:,:), v2(:,:)
  integer, allocatable         :: u1(:), u2(:) , ic(:)
  real(kind=8)                 :: time, time_step, t
  real(kind=8), pointer        :: d_ptr
  integer, pointer             :: i_ptr
  logical(kind=1), pointer     :: b_ptr
  logical                      :: l
  integer                      :: itemsize, nbytes, dim
  type(bmi) :: bmif, bmif1
  ! --------------------------------------------------------------------------
  ! Create PhreeqcRM
  ! --------------------------------------------------------------------------
    nxyz = 40
#ifdef USE_MPI
    ! MPI
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
    if (status .ne. MPI_SUCCESS) then
        stop "Failed to get mpi_myself"
    endif  
    id = bmif%bmif_create(nxyz, MPI_COMM_WORLD)
    if (mpi_myself > 0) then
        status = RM_MpiWorker(id)
        status = bmif%bmif_finalize(id)
        return
    endif
#else
    ! OpenMP
    !id = bmif%bmif_create()
	id = bmif%bmif_initialize()
#endif
    ! Write YAML file
    yid = CreateYAMLPhreeqcRM()
    status = YAMLSetGridCellCount(yid, nxyz)
    status = YAMLThreadCount(yid, 3)
	yaml_filename = "TestAllMethods_f90.yaml"
	status = WriteYAMLDoc(yid, yaml_filename)
	status = YAMLClear(yid)
    status = DestroyYAMLPhreeqcRM(yid)  
  
	write(*,*) "bmif%bmif_create"
	!-------
	status = bmif%bmif_initialize(yaml_filename)
	status = RM_InitializeYAML(id, yaml_filename)
	write(*,*) "bmif%bmif_initialize"
	!-------
	nxyz = RM_GetGridCellCount(id)
	status = bmif%bmif_get_value("GridCellCount", nxyz)
	write(*,*) "GetGridCellCount"
	!-------
	n = RM_GetThreadCount(id)
	write(*,*) "GetThreadCount " 
	!-------
	! Inactive cells or symmetry
	allocate(IntVector(nxyz))
    IntVector = -1
	do i = 1, nxyz / 2 
        IntVector(i) = i
    enddo
	status = RM_LoadDatabase(id, "phreeqc.dat")
	write(*,*) "LoadDatabase"
	!
	! Set properties
	!
	!-------
	status = RM_SetComponentH2O(id, 0)
	write(*,*) "SetComponentH2O "
	!-------
	status = RM_SetSpeciesSaveOn(id, 1)
	write(*,*) "SetSpeciesSaveOn "
	!-------
	status = RM_SetErrorOn(id, 1)
	write(*,*) "SetErrorOn "
	!-------
	status = RM_SetErrorHandlerMode(id, 1)
	write(*,*) "SetErrorHandlerMode "
	!-------
	status = RM_SetDumpFileName(id, "TestAllMethods_py.dump")
	write(*,*) "SetDumpFileName "
	!-------
	status = bmif%bmif_set_value("FilePrefix", "TestAllMethods_py")
	status = RM_SetFilePrefix(id, "TestAllMethods_py")
	write(*,*) "SetFilePrefix "
	!-------
	status =RM_OpenFiles(id)
	write(*,*) "OpenFiles "
	!-------
	status =RM_SetPartitionUZSolids(id, 0)
	write(*,*) "SetPartitionUZSolids "
	!-------
	status =RM_SetRebalanceByCell(id, 1)
	write(*,*) "SetRebalanceByCell "
	!-------
	status =RM_SetRebalanceFraction(id, 0.5d0)
	write(*,*) "SetRebalanceFraction "
	!-------
	status =RM_SetScreenOn(id, 1)
	write(*,*) "SetScreenOn "
	!-------
	status = RM_SetSelectedOutputOn(id, 1)
	status = bmif%bmif_set_value("SelectedOutputOn", .true.)
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
	status = RM_UseSolutionDensityVolume(id, 0)
	write(*,*) "UseSolutionDensityVolume "
	!-------
	t = 1.0 / 86400.0
	status = RM_SetTimeConversion(id, t)
	write(*,*) "SetTimeConversion "
	!-------
	allocate(DoubleVector(nxyz))
    DoubleVector = 1.0d0
	status = RM_SetRepresentativeVolume(id, DoubleVector)
	write(*,*) "SetRepresentativeVolume "

	!-------Chemistry cells may be fewer than GridCellCount
	deallocate(IntVector)
	allocate(IntVector(nxyz))
	IntVector = 1
	status = RM_SetPrintChemistryMask(id, IntVector)
	write(*,*) "SetPrintChemistryMask "
	!-------
	status = RM_SetPrintChemistryOn(id, 0, 1, 0)
	write(*,*) "RM_SetPrintChemistryOn "
	!
	! Define reactants available for initial 
	! and boundary conditions in this file
	!
	status = RM_RunFile(id, 1, 1, 1, "all_reactants.pqi")
	write(*,*) "RunFile "
	!-------
	status = bmif%bmif_add_output_vars("AddOutputVars", "True")
	status = bmif%bmif_add_output_vars("SolutionProperties", "True")
	status = bmif%bmif_add_output_vars("SolutionTotalMolalities", "True")
	status = bmif%bmif_add_output_vars("ExchangeMolalities", "True")
	status = bmif%bmif_add_output_vars("SurfaceMolalities", "True")
	status = bmif%bmif_add_output_vars("EquilibriumPhases", "True")
	status = bmif%bmif_add_output_vars("Gases", "True")
	status = bmif%bmif_add_output_vars("KineticReactants","True")
	status = bmif%bmif_add_output_vars("SolidSolutions", "True")
	status = bmif%bmif_add_output_vars("CalculateValues", "True")
	status = bmif%bmif_add_output_vars("SolutionActivities", "H+ Ca+2 Na+")
	status = bmif%bmif_add_output_vars("SolutionMolalities", "OH- Cl-")
	status = bmif%bmif_add_output_vars("SaturationIndices", "Calcite Dolomite")
	write(*,*) "AddOutputVars "
	!-------
	ncomps = RM_FindComponents(id)
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
	ncomps = RM_GetComponentCount(id)
	status = bmif%bmif_get_value("ComponentCount", ncomps)
	status = bmif%bmif_get_value_ptr("ComponentCount", i_ptr)
	write(*,*) "GetComponentCount)" 
	!-------
	status = bmif%bmif_get_value("Components", StringVector)
	status = RM_GetComponents(id, StringVector)
	write(*,*) "GetComponents)" 
	! Species info
	nspecies = RM_GetSpeciesCount(id)
	write(*,*) "GetSpeciesCount "
	!-------
	status = RM_GetSpeciesNames(id, StringVector)
	write(*,*) "GetSpeciesNames "
	!-------
    if(allocated(DoubleVector)) deallocate(DoubleVector)
    allocate(DoubleVector(nxyz*nspecies))
	status = RM_GetSpeciesD25(id, DoubleVector)
	write(*,*) "GetSpeciesD25 "
	!-------
	status = RM_GetSpeciesZ(id, DoubleVector)
	write(*,*) "GetSpeciesZ "
	! Reactant lists
	!status = RM_GetEquilibriumPhasesName(id, 1, string)
    !write(*,*) "GetEquilibriumPhasesName "
	status = RM_GetEquilibriumPhasesNames(id, StringVector)
    write(*,*) "GetEquilibriumPhasesNames "
	!-------
	n = RM_GetEquilibriumPhasesCount(id)
	write(*,*) "GetEquilibriumPhasesCount "
	!-------
	!status = RM_GetExchangeName(id, 1, string)
	!write(*,*) "GetExchangeName "
	status = RM_GetExchangeNames(id, StringVector)
	write(*,*) "GetExchangeNames "
	!-------
	!status = RM_GetExchangeSpeciesName(id, 1, string)
	!write(*,*) "GetExchangeSpeciesName "
	status = RM_GetExchangeSpeciesNames(id, StringVector)
	write(*,*) "GetExchangeSpeciesNames "
	!-------
	n = RM_GetExchangeSpeciesCount(id)
	write(*,*) "GetExchangeSpeciesCount "
	!-------
	!status = RM_GetGasComponentsName(id, 1, string)
	!write(*,*) "GetGasComponentsName "
	status = RM_GetGasComponentsNames(id, StringVector)
	write(*,*) "GetGasComponentsNames "
	!-------
	ngas = RM_GetGasComponentsCount(id)
	write(*,*) "GetGasComponentsCount "
	!-------
	status = bmif%bmif_get_value("Gfw", DoubleVector)
	status = RM_GetGfw(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("Gfw", d_ptr)
	write(*,*) "GetGfw "
	!-------
	n = RM_GetKineticReactionsCount(id)
	write(*,*) "GetKineticReactionsCount "
	!-------
	!status = RM_GetKineticReactionsName(id, 1, string)
	!write(*,*) "GetKineticReactionsName "
	status = RM_GetKineticReactionsNames(id, StringVector)
	write(*,*) "GetKineticReactionsNames "
	!-------
	n = RM_GetSICount(id)
	write(*,*) "GetSICount "
	!-------
	!status = RM_GetSIName(id, 1, string)
	!write(*,*) "GetSIName "
	status = RM_GetSINames(id, StringVector)
	write(*,*) "GetSINames "
	!-------
	n = RM_GetSolidSolutionComponentsCount(id)
	write(*,*) "GetSolidSolutionComponentsCount "
	!-------
	!status = RM_GetSolidSolutionComponentsName(id, 1, string)
	!write(*,*) "GetSolidSolutionComponentsName "
	status = RM_GetSolidSolutionComponentsNames(id, StringVector)
	write(*,*) "GetSolidSolutionComponentsNames "
	!-------
	!status = RM_GetSolidSolutionName(id, 1, string)
	!write(*,*) "GetSolidSolutionName "
	status = RM_GetSolidSolutionNames(id, StringVector)
	write(*,*) "GetSolidSolutionNames "
	!-------
	!status = RM_GetSurfaceName(id, 1, string)
	!write(*,*) "GetSurfaceName "
	status = RM_GetSurfaceNames(id, StringVector)
	write(*,*) "GetSurfaceNames "
	!-------
	n = RM_GetSurfaceSpeciesCount(id)
	write(*,*) "GetSurfaceSpeciesCount "
	!-------
	!status = RM_GetSurfaceSpeciesName(id, 1, string)
	!write(*,*) "GetSurfaceSpeciesName "
	status = RM_GetSurfaceSpeciesNames(id, StringVector)
	write(*,*) "GetSurfaceSpeciesNames "
	!-------
	!status = RM_GetSurfaceType(id, 1, string)
	!write(*,*) "GetSurfaceType "
	status = RM_GetSurfaceTypes(id, StringVector)
	write(*,*) "GetSurfaceTypes "
	!
	! Remove any reactants in workers 
	! before populating cells with reactants
	!
	string = "DELETE -all"
	status = RM_RunString(id, 1, 0, 0, string)
	write(*,*) "RunString "
	!-------
	!
	! Transfer initial conditions
	!
	deallocate(IntVector)
	allocate(IntVector(nxyz))
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
	allocate(IntVector2(nxyz, 7))
	IntVector2 = 1
	status = RM_InitialPhreeqc2Module(id, IntVector2)
	write(*,*) "InitialPhreeqc2Module "
	!-------
	! Alternative B.to the previous seven methods, possible mixing
	allocate(v1(nxyz, 7))
	v1 = 1
	allocate(v2(nxyz, 7))
	v2 = -1
	allocate(f2(nxyz,7))
	f2 = 1.0
	status = RM_InitialPhreeqc2Module(id, v1, v2, f2)
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
	allocate(u1(1), u2(1), f1(1))
	u1 = 1
	u2 = -1
	f1 = 1
    allocate(bc2(size(u1), ncomps))
	status = RM_InitialPhreeqc2Concentrations(id, bc2, size(u1), u1, u2, f1)
	write(*,*) "InitialPhreeqc2Concentrations mix "
	!-------
	deallocate(u1, u2, f1)
	allocate(u1(1), u2(1), f1(1))
	u1 = 1
	u2 = -1
	f1 = 1.0d0
    allocate(bc_species(size(u1), nspecies))
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_species, size(u1), u1, u2, f1)
	write(*,*) "InitialPhreeqc2SpeciesConcentrations mix "
	!
	! Get/Set methods for time steping
	!
	d = RM_GetTime(id)
	status = bmif%bmif_get_value("Time", d)
	status = bmif%bmif_get_current_time(d)
	status = bmif%bmif_get_start_time(d)
	status = bmif%bmif_get_value_ptr("Time", d_ptr)
	write(*,*) "GetTime "
	!-------
	status = RM_SetTime(id, 0.0d0)
	status = bmif%bmif_set_value("Time", 0.0d0)
	write(*,*) "SetTime "
	!-------
	d = RM_GetTimeStep(id)
	status = bmif%bmif_get_value("TimeStep", d)
	status = bmif%bmif_get_value_ptr("TimeStep", d_ptr)
	write(*,*) "GetTimeStep "
	!-------
	status = RM_SetTimeStep(id, 0.0d0)
	status = bmif%bmif_set_value("TimeStep", 0.0d0)
	write(*,*) "SetTimeStep "
	!-------
	status = bmif%bmif_get_value("Concentrations", DoubleVector)
    allocate(DoubleVector2(nxyz, ncomps))
	status = RM_GetConcentrations(id, DoubleVector2)
	status = bmif%bmif_get_value_ptr("Concentrations", d_ptr)
	write(*,*) "GetConcentrations "
	!-------
	status =RM_SetConcentrations(id, DoubleVector2)
	status = bmif%bmif_set_value("Concentrations", DoubleVector)
	write(*,*) "SetConcentrations "
	!-------
	status = bmif%bmif_get_value("DensityCalculated", DoubleVector)
	status = RM_GetDensityCalculated(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("DensityCalculated", d_ptr)
	write(*,*) "GetDensityCalculated "
	!-------
	status = RM_SetDensityUser(id, DoubleVector)
	status = bmif%bmif_set_value("DensityUser", DoubleVector)
	write(*,*) "SetDensityUser "
	!-------
    deallocate(DoubleVector2)
    allocate(DoubleVector2(nxyz, ngas))
	status = RM_GetGasCompMoles(id, DoubleVector2)
	write(*,*) "GetGasCompMoles "
	!-------
	status = RM_SetGasCompMoles(id, DoubleVector2)
	write(*,*) "SetGasCompMoles "
	!-------
	status = RM_GetGasCompPhi(id, DoubleVector2)
	write(*,*) "GetGasCompPhi "
	!-------
	status = RM_GetGasCompPressures(id, DoubleVector2)
	write(*,*) "GetGasCompPressures "
	!-------
    deallocate(DoubleVector)
    allocate(DoubleVector(nxyz))
	status = RM_GetGasPhaseVolume(id, DoubleVector)
	write(*,*) "GetGasPhaseVolume "
	!-------
	status =RM_SetGasPhaseVolume(id, DoubleVector)
	write(*,*) "SetGasPhaseVolume "
	!-------
    do i = 1, RM_GetComponentCount(id)
		status = RM_GetIthConcentration(id, i, DoubleVector)
		!-------
		status = RM_SetIthConcentration(id, i, DoubleVector)
    enddo
	write(*,*)  "GetIthConcentration "
	write(*,*)  "SetIthConcentration "
	!-------
    do i = 1, RM_GetSpeciesCount(id)
		status = RM_GetIthSpeciesConcentration(id, i, DoubleVector)
		!-------
		status = RM_SetIthSpeciesConcentration(id, i, DoubleVector)
    enddo
    write(*,*) "GetIthSpeciesConcentration "
    write(*,*) "SetIthSpeciesConcentration "
	!-------
	status = bmif%bmif_get_value("Porosity", DoubleVector)
	status = RM_GetPorosity(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("Porosity", d_ptr)
	write(*,*) "GetPorosity "
	!-------
	status = bmif%bmif_set_value("Porosity", DoubleVector)
	status = RM_SetPorosity(id, DoubleVector)
	write(*,*) "SetPorosity "
	!-------
	status = bmif%bmif_get_value("Pressure", DoubleVector)
	status = RM_GetPressure(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("Pressure", d_ptr)
	write(*,*) "GetPressure "
	!-------
	status = bmif%bmif_set_value("Pressure", DoubleVector)
	status = RM_SetPressure(id, DoubleVector)
	write(*,*) "SetPressure "
	!-------
	status = bmif%bmif_get_value("SaturationCalculated", DoubleVector)
	status = RM_GetSaturationCalculated(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("SaturationCalculated", d_ptr)
	write(*,*) "GetSaturationCalculated "
	!-------
	status = RM_SetSaturationUser(id, DoubleVector)
	status = bmif%bmif_set_value("SaturationUser", DoubleVector)
	write(*,*) "SetSaturationUser "
	!-------
	status = bmif%bmif_get_value("SolutionVolume", DoubleVector)
	status = RM_GetSolutionVolume(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("SolutionVolume", d_ptr)
	write(*,*) "GetSolutionVolume "
	!-------
    deallocate(DoubleVector2)
    allocate(DoubleVector2(nxyz, nspecies))
	status = RM_GetSpeciesConcentrations(id, DoubleVector2)
	write(*,*) "GetSpeciesConcentrations "
	!-------
	status = RM_SpeciesConcentrations2Module(id, DoubleVector2)
	write(*,*) "SpeciesConcentrations2Module "
	!-------
	status = RM_GetSpeciesLog10Gammas(id, DoubleVector2)
	write(*,*) "GetSpeciesLog10Gammas "
	!-------
	status = RM_GetSpeciesLog10Molalities(id, DoubleVector2)
	write(*,*) "GetSpeciesLog10Molalities "
	!-------
	status = bmif%bmif_get_value("Temperature", DoubleVector)
	status = RM_GetTemperature(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("Temperature", d_ptr)
	write(*,*) "GetTemperature "
	!-------
	status = RM_SetTemperature(id, DoubleVector)
	status = bmif%bmif_set_value("Temperature", DoubleVector)
	write(*,*) "SetTemperature "
	!-------
	status = bmif%bmif_get_value("Viscosity", DoubleVector)
	status = RM_GetViscosity(id, DoubleVector)
	status = bmif%bmif_get_value_ptr("Viscosity", d_ptr)	
	write(*,*) "GetViscosity "
	!
	! Take a time step
	!
	status = bmif%bmif_update()
	write(*,*) "Update"
	!-------
	status =RM_RunCells(id)
	write(*,*) "RunCells"
	!-------
	status = bmif%bmif_update_until(86400.0d0)
	write(*,*) "UpdateUntil"
	!
	! Selected output
	!
	status = RM_SetNthSelectedOutput(id, 1)
	status = bmif%bmif_set_value("NthSelectedOutput", 1)
	write(*,*) "SetNthSelectedOutput "
	!-------
	n_user = RM_GetCurrentSelectedOutputUserNumber(id)
	status = bmif%bmif_get_value("CurrentSelectedOutputUserNumber", n_user)
	write(*,*) "GetCurrentSelectedOutputUserNumber "
	!-------
	n = RM_GetNthSelectedOutputUserNumber(id, 1)
	write(*,*) "GetNthSelectedOutputUserNumber "
	!-------
	status = bmif%bmif_get_value("SelectedOutput", DoubleVector)
    if(allocated(DoubleVector2)) deallocate(DoubleVector2)
    n = RM_GetSelectedOutputColumnCount(id)
    allocate(DoubleVector2(nxyz, n))
	status = RM_GetSelectedOutput(id, DoubleVector2)
	write(*,*) "GetSelectedOutput "
	!-------
	n = RM_GetSelectedOutputColumnCount(id)
	status = bmif%bmif_get_value("SelectedOutputColumnCount", n)
	write(*,*) "GetSelectedOutputColumnCount "
	!-------
	n = RM_GetSelectedOutputCount(id)
	status = bmif%bmif_get_value("SelectedOutputCount", n)
	write(*,*) "GetSelectedOutputCount "
	!-------
	status = bmif%bmif_get_value("SelectedOutputHeadings", StringVector)
	status = RM_GetSelectedOutputHeadings(id, StringVector)
	write(*,*) "GetSelectedOutputHeadings "
	!-------
	!b = RM_GetSelectedOutputOn(id)
	status = bmif%bmif_get_value("SelectedOutputOn", l)
	status = bmif%bmif_get_value_ptr("SelectedOutputOn", b_ptr)	
	write(*,*) "GetSelectedOutputOn "
	!-------
	n = RM_GetSelectedOutputRowCount(id)
	status = bmif%bmif_get_value("SelectedOutputRowCount", n)
	write(*,*) "GetSelectedOutputRowCount "
	!-------
	status = RM_SetCurrentSelectedOutputUserNumber(id, 333)
	write(*,*) "SetCurrentSelectedOutputUserNumber "
	!
	! Getters
	!
	status = RM_GetBackwardMapping(id, 1, IntVector)
	write(*,*) "GetBackwardMapping "
	!-------
	!status = RM_GetDatabaseFileName(id, string)
	!write(*,*) "GetDatabaseFileName "
	!-------
	status = RM_GetEndCell(id, IntVector)
	write(*,*) "GetEndCell"
	!-------
	!n = RM_GetErrorHandlerMode()
	!write(*,*) "GetErrorHandlerMode "
	!-------
	status = RM_GetErrorString(id, AllocString)
	status = bmif%bmif_get_value("ErrorString", AllocString)
	write(*,*) "GetErrorString "
	!-------
	status = RM_GetFilePrefix(id, AllocString)
	status = bmif%bmif_get_value("FilePrefix", AllocString)
	write(*,*) "GetFilePrefix "
	!-------
	!status = RM_GetForwardMapping()  ! not implemented
	!write(*,*) "GetForwardMapping "
	!-------
	n = RM_GetIPhreeqcID(id, 0)
	write(*,*) "GetIPhreeqcID "
	!-------
	n = RM_GetMpiMyself(id)
	write(*,*) "GetMpiMyself "
	!-------
	n = RM_GetMpiTasks(id)
	write(*,*) "GetMpiTasks "
	!-------
	!status = RM_GetPartitionUZSolids(id, n)  ! Not implemented
	!write(*,*) "GetPartitionUZSolids "
	!-------
	!status = RM_GetPrintChemistryMask(id, n)   ! Not implemented
	!write(*,*) "GetPrintChemistryMask "
	!-------
	!status = RM_GetPrintChemistryOn(id, n)  ! Not implemented
	!write(*,*) "GetPrintChemistryOn "
	!-------
	!status = RM_GetRebalanceByCell(id, n)   ! Not implemented
	!write(*,*) "GetRebalanceByCell "
	!-------
	!status = RM_GetRebalanceFraction(id, d)  ! Not implemented
	!write(*,*) "GetRebalanceFraction "
	!-------
	!status = RM_GetSpeciesSaveOn(id, n)       ! Not implemented
	!write(*,*) "GetSpeciesSaveOn "
	!-------
	!status = RM_GetSpeciesStoichiometry(id, i, IntVector)   ! Not implemented
	!write(*,*) "GetSpeciesStoichiometry "
	!-------
	status = RM_GetStartCell(id, IntVector)
	write(*,*) "GetStartCell "
	!-------
	d = RM_GetTimeConversion(id)
	write(*,*) "GetTimeConversion "
	!-------
	!status = RM_GetUnitsExchange(id, n)   ! Not implemented
	!write(*,*) "GetUnitsExchange "
	!-------
	!status = RM_GetUnitsGasPhase(id, n)   ! Not implemented
	!write(*,*) "GetUnitsGasPhase "
	!-------
	!status = RM_GetUnitsKinetics(id, n)   ! Not implemented
	!write(*,*) "GetUnitsKinetics "
	!-------
	!status = RM_GetUnitsPPassemblage(id, n)   ! Not implemented
	!write(*,*) "GetUnitsPPassemblage "
	!-------
	!status = RM_GetUnitsSolution(id, n)   ! Not implemented
	!write(*,*) "GetUnitsSolution "
	!-------
	!status = RM_GetUnitsSSassemblage(id, n)   ! Not implemented
	!write(*,*) "GetUnitsSSassemblage "
	!-------
	!status = RM_GetUnitsSurface(id, n)   ! Not implemented
	!write(*,*) "GetUnitsSurface "
	!-------
	!std::vector<IPhreeqcPhast *> w = status = RM_GetWorkers()   ! Not implemented
	!write(*,*) "GetWorkers "
	!
	! Utilities
	!
#ifndef USE_MPI    
	!id1 = bmif%bmif_create(10, 1)  ! make another bmiphreeqcrm
	id1 = bmif1%bmif_initialize(10, 1)
	status = RM_CloseFiles(id1) 
	status = bmif1%bmif_finalize()   ! destroy the new bmiphreeqcrm
	write(*,*) "CloseFiles "
#endif    
	!-------
	deallocate(IntVector)
	allocate(IntVector(1))
	IntVector = 1
	n = size(IntVector)
    deallocate(bc2)
	allocate(bc2(n,ncomps))
	status = RM_InitialPhreeqc2Concentrations(id, bc2, size(IntVector), IntVector)
	n = size(IntVector)
	allocate(tc(n), p_atm(n))
	tc = 30
	p_atm = 1.5
	status = RM_Concentrations2Utility(id, bc2, n, tc, p_atm)
	write(*,*) "Concentrations2Utility "
	!-------
	status = RM_DecodeError(id, -2)	   
	write(*,*) "DecodeError "
	!-------
	status =RM_DumpModule(id, 1, 0)
	write(*,*) "DumpModule "
	!-------
	!status = RM_ErrorHandler(id, 0, "string") 
	!write(*,*) "OK, just a test: ErrorHandler "
	!-------
	status = RM_ErrorMessage(id, "my error")  
	write(*,*) "OK, just a test: ErrorMessage "
	!-------
	status = RM_LogMessage(id, "Log message") 
	write(*,*) "LogMessage "
	!-------
	status = RM_OutputMessage(id, "Output message") 
	write(*,*) "OutputMessage "
	!-------
	status = RM_ScreenMessage(id, "Screen message") 
	write(*,*) "ScreenMessage "
	!-------
	status =RM_StateSave(id, 1)
	write(*,*) "StateSave "
	!-------
	status =RM_StateApply(id, 1)
	write(*,*) "StateApply "
	!-------
	status =RM_StateDelete(id, 1)
	write(*,*) "StateDelete "
	!-------
	status = RM_WarningMessage(id, "Warning message")  
	write(*,*) "WarningMessage "
	!
	! BMI Methods
	!
	status = bmif%bmif_get_component_name(string)
	write(*,*) "bmif%bmif_get_component_name "
	!-------
	status = bmif%bmif_get_current_time(d)
	write(*,*) "bmif%bmif_get_current_time "
	!-------
	status = bmif%bmif_get_end_time(d)
	write(*,*) "bmif%bmif_get_end_time "
	!-------
	status = bmif%bmif_grid_rank(0, n)
	write(*,*) "bmif%bmif_grid_rank "
	!-------
	status = bmif%bmif_grid_size(0, n)
	write(*,*) "bmif%bmif_grid_size "
	!-------
	status = bmif%bmif_grid_type(0, string)
	write(*,*) "bmif%bmif_grid_type "
	!-------
	status = bmif%bmif_get_input_item_count(n)
	write(*,*) "bmif%bmif_get_input_item_count "
	!-------
	status = bmif%bmif_get_input_var_names(StringVector)
	write(*,*) "bmif%bmif_get_input_var_names "
	!-------
	status = bmif%bmif_get_output_item_count(n)
	write(*,*) "bmif%bmif_get_output_item_count "
	!-------
	status = bmif%bmif_get_output_var_names(StringVector)
	write(*,*) "bmif%bmif_get_output_var_names "
	!-------
	status = bmif%bmif_get_pointable_item_count(n)
	write(*,*) "bmif%bmif_get_pointable_item_count "
	!-------
	status = bmif%bmif_get_pointable_var_names(StringVector)
	write(*,*) "bmif%bmif_get_pointable_var_names "
	!-------
	status = bmif%bmif_get_time_step(d)
	write(*,*) "bmif%bmif_get_time_step "
	!-------
	status = bmif%bmif_get_time_units(string)
	write(*,*) "bmif%bmif_get_time_units "
	!-------
	status = bmif%bmif_get_value("solution_saturation_index_Calcite", DoubleVector)
	write(*,*) "bmif%bmif_get_value "
	!-------
	status = bmif%bmif_get_var_itemsize("solution_saturation_index_Calcite", n)
	write(*,*) "bmif%bmif_get_var_itemsize "
	!-------
	status = bmif%bmif_get_var_nbytes("solution_saturation_index_Calcite", n)
	write(*,*) "bmif%bmif_get_var_nbytes "
	!-------
	status = bmif%bmif_get_var_type("solution_saturation_index_Calcite", string)
	write(*,*) "bmif%bmif_get_var_type "
	!-------
	status = bmif%bmif_get_var_units("solution_saturation_index_Calcite", string)
	write(*,*) "bmif%bmif_get_var_units "
	!status = bmif%bmif_initialize(YAML_filename)
	! See above
	status = bmif%bmif_set_value("Time", 1.0d0) 
	write(*,*) "bmif%bmif_set_value"
	!-------
	status = bmif%bmif_update() 
	write(*,*) "bmif%bmif_update"
	!-------
	status = bmif%bmif_update_until(864000.0d0) 
	write(*,*) "bmif%bmif_update_until"
	!-------
    
	write(*,*) "AddOutputVars"
	status = bmif%bmif_get_output_var_names(Names)
	do i = 1, size(StringVector)
		status = bmif%bmif_get_var_itemsize(Names(i), itemsize)
		status = bmif%bmif_get_var_nbytes(Names(i), nbytes)
		status = bmif%bmif_get_var_type(Names(i), string)
		if (itemsize .eq. 0) itemsize=1
		if (nbytes .eq. 0) nbytes=1
		dim = nbytes / itemsize
		status = bmif%bmif_get_var_type(Names(i), vtype)
		if (vtype .eq. "real(kind=8)") then
			if (dim .eq. 1) then
				status = bmif%bmif_get_value(Names(i), d)
                write(*,*) "     ", Names(i), "  ", d
			else
				status = bmif%bmif_get_value(Names(i), DoubleVector)
                write(*,*) "     ", Names(i), "  ", DoubleVector(1)
			endif
		else if (vtype .eq. "integer") then
			if (dim .eq. 1) then
				status = bmif%bmif_get_value(Names(i), j)
                write(*,*) "     ", Names(i), "  ", j
			else
				status = bmif%bmif_get_value(Names(i), IntVector)
                write(*,*) "     ", Names(i), "  ", IntVector(1)
			endif
		else if (vtype .eq. "logical") then
			if (dim == 1) then
				status = bmif%bmif_get_value(Names(i), l)
                write(*,*) "     ", Names(i), "  ", l
			endif
		else if (vtype(1:9) .eq. "character") then
			if (dim == 1) then
				status = bmif%bmif_get_value(Names(i), string)
                write(*,*) "     ", Names(i), "  ", trim(string)
			else
				status = bmif%bmif_get_value(Names(i), StringVector)
                write(*,*) "     ", Names(i), "  ", trim(StringVector(1))
			endif
		endif
	enddo    
 
#ifdef USE_MPI
	status = RM_MpiWorkerBreak(id)
#endif    
	status = bmif%bmif_finalize()    ! void method
	write(*,*) "bmif%bmif_finalize "

	!Should be private: status =RM_ReturnHandler()
	!TODO status =RM_MpiAbort()
	!TODO status =RM_MpiWorker()
	!TODO status =RM_MpiWorkerBreak()
	!TODO status =RM_SetMpiWorkerCallbackC()
	!TODO status =RM_SetMpiWorkerCallbackCookie()
	write(*,*) "Success."
	return

  ! Deallocate
  deallocate(IntVector, IntVector2)
  deallocate(DoubleVector, DoubleVector2)
  deallocate(StringVector, AllocString)
  deallocate(IntVector)
  deallocate(Names)
  deallocate(v1, v2, f1)
  deallocate(ic1, ic2)
  deallocate(u1, u2, ic)
  deallocate(cells, v1, v2)
  deallocate(bc1, c, bc_species, bc2)
  deallocate(tc)
  deallocate(p_atm)
  return
end subroutine TestAllMethods_f90

! YAML
#endif
