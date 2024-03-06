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
  real(kind=4), pointer        :: float_ptr
  real(kind=4)                 :: float
  real(kind=4), allocatable    :: FloatVector(:)
  logical                      :: l
  integer                      :: itemsize, nbytes, dim
  type(bmi) :: brm, brm1
  type(YAML_PhreeqcRM) :: yrm
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
    id = brm%bmif_create(nxyz, MPI_COMM_WORLD)
    if (mpi_myself > 0) then
        !status = RM_MpiWorker(id)     ! Deprecated
        status = brm%MpiWorker()
        status = brm%bmif_finalize()
        return
    endif
#else
    ! OpenMP
    id = brm%bmif_create_default()
#endif
    ! Write YAML file
    yid = yrm%CreateYAMLPhreeqcRM()
    status = yrm%YAMLSetGridCellCount(nxyz)
    status = yrm%YAMLThreadCount(3)
	yaml_filename = "TestAllMethods_f90.yaml"
	status = yrm%WriteYAMLDoc(yaml_filename)
	status = yrm%YAMLClear()
    status = yrm%DestroyYAMLPhreeqcRM()  
  
	write(*,*) "brm%bmif_create"
	!-------
	status = brm%bmif_initialize(yaml_filename)
	status = RM_InitializeYAML(id, yaml_filename)     ! Deprecated
	write(*,*) "brm%bmif_initialize"
	!-------
	nxyz = RM_GetGridCellCount(id)     ! Deprecated
	nxyz = brm%GetGridCellCount()
	status = brm%bmif_get_value("GridCellCount", nxyz)
	write(*,*) "GetGridCellCount"
	!-------
	n = RM_GetThreadCount(id)     ! Deprecated
	n = brm%GetThreadCount()
	write(*,*) "GetThreadCount " 
	!-------
	! Inactive cells or symmetry
	allocate(IntVector(nxyz))
    IntVector = -1
	do i = 1, nxyz / 2 
        IntVector(i) = i
    enddo
	status = RM_LoadDatabase(id, "phreeqc.dat")     ! Deprecated
	status = brm%LoadDatabase("phreeqc.dat")
	write(*,*) "LoadDatabase"
	!
	! Set properties
	!
	!-------
	status = RM_SetComponentH2O(id, 0)     ! Deprecated
	status = brm%SetComponentH2O(0)
	write(*,*) "SetComponentH2O "
	!-------
	status = RM_SetSpeciesSaveOn(id, 1)     ! Deprecated
	status = brm%SetSpeciesSaveOn(1)
	write(*,*) "SetSpeciesSaveOn "
	!-------
	status = RM_SetErrorOn(id, 1)     ! Deprecated
	status = brm%SetErrorOn(1)
	write(*,*) "SetErrorOn "
	!-------
	status = RM_SetErrorHandlerMode(id, 1)     ! Deprecated
	status = brm%SetErrorHandlerMode(1)
	write(*,*) "SetErrorHandlerMode "
	!-------
	status = RM_SetDumpFileName(id, "TestAllMethods_py.dump")     ! Deprecated
	status = brm%SetDumpFileName("TestAllMethods_py.dump")
	write(*,*) "SetDumpFileName "
	!-------
	status = brm%bmif_set_value("FilePrefix", "TestAllMethods_py")
	status = RM_SetFilePrefix(id, "TestAllMethods_py")     ! Deprecated
	status = brm%SetFilePrefix("TestAllMethods_py")
	write(*,*) "SetFilePrefix "
	!-------
	status =RM_OpenFiles(id)     ! Deprecated
	status =brm%OpenFiles()
	write(*,*) "OpenFiles "
	!-------
	status =RM_SetPartitionUZSolids(id, 0)     ! Deprecated
	status =brm%SetPartitionUZSolids(0)
	write(*,*) "SetPartitionUZSolids "
	!-------
	status =RM_SetRebalanceByCell(id, 1)     ! Deprecated
	status =brm%SetRebalanceByCell(1)
	write(*,*) "SetRebalanceByCell "
	!-------
	status =RM_SetRebalanceFraction(id, 0.5d0)     ! Deprecated
	status =brm%SetRebalanceFraction(0.5d0)
	write(*,*) "SetRebalanceFraction "
	!-------
	status =RM_SetScreenOn(id, 1)     ! Deprecated
	status =brm%SetScreenOn(1)
	write(*,*) "SetScreenOn "
	!-------
	status = RM_SetSelectedOutputOn(id, 1)     ! Deprecated
	status = brm%SetSelectedOutputOn(1)
	status = brm%bmif_set_value("SelectedOutputOn", .true.)
	write(*,*) "SetSelectedOutputOn "
	!-------
	status = RM_SetUnitsExchange(id, 1)     ! Deprecated
	status = brm%SetUnitsExchange(1)
	write(*,*) "SetUnitsExchange "
	!-------
	status =RM_SetUnitsGasPhase(id, 1)     ! Deprecated
	status =brm%SetUnitsGasPhase(1)
	write(*,*) "SetUnitsGasPhase "
	!-------
	status =RM_SetUnitsKinetics(id, 1)     ! Deprecated
	status =brm%SetUnitsKinetics(1)
	write(*,*) "SetUnitsKinetics "
	!-------
	status =RM_SetUnitsPPassemblage(id, 1)     ! Deprecated
	status =brm%SetUnitsPPassemblage(1)
	write(*,*) "SetUnitsPPassemblage "
	!-------
	status =RM_SetUnitsSolution(id, 2)     ! Deprecated
	status =brm%SetUnitsSolution(2)
	write(*,*) "SetUnitsSolution "
	!-------
	status =RM_SetUnitsSSassemblage(id, 1)     ! Deprecated
	status =brm%SetUnitsSSassemblage(1)
	write(*,*) "SetUnitsSSassemblage "
	!-------
	status =RM_SetUnitsSurface(id, 1)     ! Deprecated
	status =brm%SetUnitsSurface(1)
	write(*,*) "SetUnitsSurface "
	!-------
	status = RM_UseSolutionDensityVolume(id, 0)     ! Deprecated
	status = brm%UseSolutionDensityVolume(0)
	write(*,*) "UseSolutionDensityVolume "
	!-------
	t = 1.0 / 86400.0
	status = RM_SetTimeConversion(id, t)     ! Deprecated
	status = brm%SetTimeConversion(t)
	write(*,*) "SetTimeConversion "
	!-------
	allocate(DoubleVector(nxyz))
    DoubleVector = 1.0d0
	status = RM_SetRepresentativeVolume(id, DoubleVector)     ! Deprecated
	status = brm%SetRepresentativeVolume(DoubleVector)
	write(*,*) "SetRepresentativeVolume "

	!-------Chemistry cells may be fewer than GridCellCount
	deallocate(IntVector)
	allocate(IntVector(nxyz))
	IntVector = 1
	status = RM_SetPrintChemistryMask(id, IntVector)     ! Deprecated
	status = brm%SetPrintChemistryMask(IntVector)
	write(*,*) "SetPrintChemistryMask "
	!-------
	status = RM_SetPrintChemistryOn(id, 0, 1, 0)     ! Deprecated
	status = brm%SetPrintChemistryOn(0, 1, 0)
	write(*,*) "RM_SetPrintChemistryOn "     ! Deprecated
	write(*,*) "brm%SetPrintChemistryOn "
	!
	! Define reactants available for initial 
	! and boundary conditions in this file
	!
	status = RM_RunFile(id, 1, 1, 1, "all_reactants.pqi")     ! Deprecated
	status = brm%RunFile(1, 1, 1, "all_reactants.pqi")
	write(*,*) "RunFile "
	!-------
	status = brm%bmif_add_output_vars("AddOutputVars", "True")
	status = brm%bmif_add_output_vars("SolutionProperties", "True")
	status = brm%bmif_add_output_vars("SolutionTotalMolalities", "True")
	status = brm%bmif_add_output_vars("ExchangeMolalities", "True")
	status = brm%bmif_add_output_vars("SurfaceMolalities", "True")
	status = brm%bmif_add_output_vars("EquilibriumPhases", "True")
	status = brm%bmif_add_output_vars("Gases", "True")
	status = brm%bmif_add_output_vars("KineticReactants","True")
	status = brm%bmif_add_output_vars("SolidSolutions", "True")
	status = brm%bmif_add_output_vars("CalculateValues", "True")
	status = brm%bmif_add_output_vars("SolutionActivities", "H+ Ca+2 Na+")
	status = brm%bmif_add_output_vars("SolutionMolalities", "OH- Cl-")
	status = brm%bmif_add_output_vars("SaturationIndices", "Calcite Dolomite")
	write(*,*) "AddOutputVars "
	!-------
	ncomps = RM_FindComponents(id)     ! Deprecated
	ncomps = brm%FindComponents()
	write(*,*) "FindComponents "
	!
	! Methods up to this point are useful 
	! in a YAML initialization file
	! 
	! Lists of reactants found by FindComponents follow
	! 
	nchem = RM_GetChemistryCellCount(id)     ! Deprecated
	nchem = brm%GetChemistryCellCount()
	write(*,*) "GetChemistryCellCount "
	!-------
	ncomps = RM_GetComponentCount(id)     ! Deprecated
	ncomps = brm%GetComponentCount()
	status = brm%bmif_get_value("ComponentCount", ncomps)
	status = brm%bmif_get_value_ptr("ComponentCount", i_ptr)
	write(*,*) "GetComponentCount)" 
	!-------
	status = brm%bmif_get_value("Components", StringVector)
	status = RM_GetComponents(id, StringVector)     ! Deprecated
	status = brm%GetComponents(StringVector)
	write(*,*) "GetComponents)" 
	! Species info
	nspecies = RM_GetSpeciesCount(id)     ! Deprecated
	nspecies = brm%GetSpeciesCount()
	write(*,*) "GetSpeciesCount "
	!-------
	status = RM_GetSpeciesNames(id, StringVector)     ! Deprecated
	status = brm%GetSpeciesNames(StringVector)
	write(*,*) "GetSpeciesNames "
	!-------
    if(allocated(DoubleVector)) deallocate(DoubleVector)
    allocate(DoubleVector(nspecies))
	status = RM_GetSpeciesD25(id, DoubleVector)     ! Deprecated
	status = brm%GetSpeciesD25(DoubleVector)
	write(*,*) "GetSpeciesD25 "
	!-------
	status = RM_GetSpeciesZ(id, DoubleVector)     ! Deprecated
	status = brm%GetSpeciesZ(DoubleVector)
	write(*,*) "GetSpeciesZ "
	! Reactant lists
	!status = RM_GetEquilibriumPhasesName(id, 1, string)     ! Deprecated
	!status = brm%GetEquilibriumPhasesName(1, string)
    !write(*,*) "GetEquilibriumPhasesName "
	status = RM_GetEquilibriumPhasesNames(id, StringVector)     ! Deprecated
	status = brm%GetEquilibriumPhasesNames(StringVector)
    write(*,*) "GetEquilibriumPhasesNames "
	!-------
	n = RM_GetEquilibriumPhasesCount(id)     ! Deprecated
	n = brm%GetEquilibriumPhasesCount()
	write(*,*) "GetEquilibriumPhasesCount "
	!-------
	!status = RM_GetExchangeName(id, 1, string)     ! Deprecated
	!status = brm%GetExchangeName(1, string)
	!write(*,*) "GetExchangeName "
	status = RM_GetExchangeNames(id, StringVector)     ! Deprecated
	status = brm%GetExchangeNames(StringVector)
	write(*,*) "GetExchangeNames "
	!-------
	!status = RM_GetExchangeSpeciesName(id, 1, string)     ! Deprecated
	!status = brm%GetExchangeSpeciesName(1, string)
	!write(*,*) "GetExchangeSpeciesName "
	status = RM_GetExchangeSpeciesNames(id, StringVector)     ! Deprecated
	status = brm%GetExchangeSpeciesNames(StringVector)
	write(*,*) "GetExchangeSpeciesNames "
	!-------
	n = RM_GetExchangeSpeciesCount(id)     ! Deprecated
	n = brm%GetExchangeSpeciesCount()
	write(*,*) "GetExchangeSpeciesCount "
	!-------
	!status = RM_GetGasComponentsName(id, 1, string)     ! Deprecated
	!status = brm%GetGasComponentsName(1, string)
	!write(*,*) "GetGasComponentsName "
	status = RM_GetGasComponentsNames(id, StringVector)     ! Deprecated
	status = brm%GetGasComponentsNames(StringVector)
	write(*,*) "GetGasComponentsNames "
	!-------
	ngas = RM_GetGasComponentsCount(id)     ! Deprecated
	ngas = brm%GetGasComponentsCount()
	write(*,*) "GetGasComponentsCount "
	!-------
	status = brm%bmif_get_value("Gfw", DoubleVector)
	status = RM_GetGfw(id, DoubleVector)     ! Deprecated
	status = brm%GetGfw(DoubleVector)
	status = brm%bmif_get_value_ptr("Gfw", d_ptr)
	write(*,*) "GetGfw "
	!-------
	n = RM_GetKineticReactionsCount(id)     ! Deprecated
	n = brm%GetKineticReactionsCount()
	write(*,*) "GetKineticReactionsCount "
	!-------
	!status = RM_GetKineticReactionsName(id, 1, string)     ! Deprecated
	!status = brm%GetKineticReactionsName(1, string)
	!write(*,*) "GetKineticReactionsName "
	status = RM_GetKineticReactionsNames(id, StringVector)     ! Deprecated
	status = brm%GetKineticReactionsNames(StringVector)
	write(*,*) "GetKineticReactionsNames "
	!-------
	n = RM_GetSICount(id)     ! Deprecated
	n = brm%GetSICount()
	write(*,*) "GetSICount "
	!-------
	!status = RM_GetSIName(id, 1, string)     ! Deprecated
	!status = brm%GetSIName(1, string)
	!write(*,*) "GetSIName "
	status = RM_GetSINames(id, StringVector)     ! Deprecated
	status = brm%GetSINames(StringVector)
	write(*,*) "GetSINames "
	!-------
	n = RM_GetSolidSolutionComponentsCount(id)     ! Deprecated
	n = brm%GetSolidSolutionComponentsCount()
	write(*,*) "GetSolidSolutionComponentsCount "
	!-------
	!status = RM_GetSolidSolutionComponentsName(id, 1, string)     ! Deprecated
	!status = brm%GetSolidSolutionComponentsName(1, string)
	!write(*,*) "GetSolidSolutionComponentsName "
	status = RM_GetSolidSolutionComponentsNames(id, StringVector)     ! Deprecated
	status = brm%GetSolidSolutionComponentsNames(StringVector)
	write(*,*) "GetSolidSolutionComponentsNames "
	!-------
	!status = RM_GetSolidSolutionName(id, 1, string)     ! Deprecated
	!status = brm%GetSolidSolutionName(1, string)
	!write(*,*) "GetSolidSolutionName "
	status = RM_GetSolidSolutionNames(id, StringVector)     ! Deprecated
	status = brm%GetSolidSolutionNames(StringVector)
	write(*,*) "GetSolidSolutionNames "
	!-------
	!status = RM_GetSurfaceName(id, 1, string)     ! Deprecated
	!status = brm%GetSurfaceName(1, string)
	!write(*,*) "GetSurfaceName "
	status = RM_GetSurfaceNames(id, StringVector)     ! Deprecated
	status = brm%GetSurfaceNames(StringVector)
	write(*,*) "GetSurfaceNames "
	!-------
	n = RM_GetSurfaceSpeciesCount(id)     ! Deprecated
	n = brm%GetSurfaceSpeciesCount()
	write(*,*) "GetSurfaceSpeciesCount "
	!-------
	!status = RM_GetSurfaceSpeciesName(id, 1, string)     ! Deprecated
	!status = brm%GetSurfaceSpeciesName(1, string)
	!write(*,*) "GetSurfaceSpeciesName "
	status = RM_GetSurfaceSpeciesNames(id, StringVector)     ! Deprecated
	status = brm%GetSurfaceSpeciesNames(StringVector)
	write(*,*) "GetSurfaceSpeciesNames "
	!-------
	!status = RM_GetSurfaceType(id, 1, string)     ! Deprecated
	!status = brm%GetSurfaceType(1, string)
	!write(*,*) "GetSurfaceType "
	status = RM_GetSurfaceTypes(id, StringVector)     ! Deprecated
	status = brm%GetSurfaceTypes(StringVector)
	write(*,*) "GetSurfaceTypes "
	!
	! Remove any reactants in workers 
	! before populating cells with reactants
	!
	string = "DELETE -all"
	status = RM_RunString(id, 1, 0, 0, string)     ! Deprecated
	status = brm%RunString(1, 0, 0, string)
	write(*,*) "RunString "
	!-------
	!
	! Transfer initial conditions
	!
	deallocate(IntVector)
	allocate(IntVector(nxyz))
	IntVector = 1
	status = RM_InitialEquilibriumPhases2Module(id, IntVector)     ! Deprecated
	status = brm%InitialEquilibriumPhases2Module(IntVector)
	write(*,*) "InitialEquilibriumPhases2Module "
	!-------
	status =RM_InitialExchanges2Module(id, IntVector)     ! Deprecated
	status =brm%InitialExchanges2Module(IntVector)
	write(*,*) "InitialExchanges2Module "
	!-------
	status =RM_InitialGasPhases2Module(id, IntVector)     ! Deprecated
	status =brm%InitialGasPhases2Module(IntVector)
	write(*,*) "InitialGasPhases2Module "
	!-------
	status =RM_InitialKinetics2Module(id, IntVector)     ! Deprecated
	status =brm%InitialKinetics2Module(IntVector)
	write(*,*) "InitialKinetics2Module "
	!-------
	status =RM_InitialSolutions2Module(id, IntVector)     ! Deprecated
	status =brm%InitialSolutions2Module(IntVector)
	write(*,*) "InitialSolutions2Module "
	!-------
	status =RM_InitialSolidSolutions2Module(id, IntVector)     ! Deprecated
	status =brm%InitialSolidSolutions2Module(IntVector)
	write(*,*) "InitialSolidSolutions2Module "
	!-------
	status = RM_InitialSurfaces2Module(id, IntVector)     ! Deprecated
	status = brm%InitialSurfaces2Module(IntVector)
	write(*,*) "InitialSurfaces2Module "
	!-------
	! Alternative A.to the previous seven methods
	allocate(IntVector2(nxyz, 7))
	IntVector2 = 1
	status = RM_InitialPhreeqc2Module(id, IntVector2)     ! Deprecated
	status = brm%InitialPhreeqc2Module(IntVector2)
	write(*,*) "InitialPhreeqc2Module "
	!-------
	! Alternative B.to the previous seven methods, possible mixing
	allocate(v1(nxyz, 7))
	v1 = 1
	allocate(v2(nxyz, 7))
	v2 = -1
	allocate(f2(nxyz,7))
	f2 = 1.0
	status = RM_InitialPhreeqc2Module(id, v1, v2, f2)     ! Deprecated
	status = brm%InitialPhreeqc2Module(v1, v2, f2)
	write(*,*) "InitialPhreeqc2Modul mix "
	!-------
	! Alternative C.to the previous seven methods, initialize cells 18 and 19
	allocate(cells(2))
	cells(1) = 18
	cells(2) = 19
	status = RM_InitialPhreeqcCell2Module(id, 1, cells, size(cells))     ! Deprecated
	status = brm%InitialPhreeqcCell2Module(1, cells, size(cells))
	write(*,*) "InitialPhreeqcCell2Module "
	!
	! Boundary conditions
	!
	allocate(u1(1), u2(1), f1(1))
	u1 = 1
	u2 = -1
	f1 = 1
    allocate(bc2(size(u1), ncomps))
	status = RM_InitialPhreeqc2Concentrations(id, bc2, size(u1), u1, u2, f1)     ! Deprecated
	status = brm%InitialPhreeqc2Concentrations(bc2, size(u1), u1, u2, f1)
	write(*,*) "InitialPhreeqc2Concentrations mix "
	!-------
	deallocate(u1, u2, f1)
	allocate(u1(1), u2(1), f1(1))
	u1 = 1
	u2 = -1
	f1 = 1.0d0
    allocate(bc_species(size(u1), nspecies))
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_species, size(u1), u1, u2, f1)     ! Deprecated
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_species, size(u1), u1, u2, f1)     ! Deprecated
	write(*,*) "InitialPhreeqc2SpeciesConcentrations mix "
	!
	! Get/Set methods for time steping
	!
	d = RM_GetTime(id)     ! Deprecated
	d = brm%GetTime()
	status = brm%bmif_get_value("Time", d)
	status = brm%bmif_get_current_time(d)
	status = brm%bmif_get_start_time(d)
	status = brm%bmif_get_value_ptr("Time", d_ptr)
	write(*,*) "GetTime "
	!-------
	status = RM_SetTime(id, 0.0d0)     ! Deprecated
	status = brm%SetTime(0.0d0)
	status = brm%bmif_set_value("Time", 0.0d0)
	write(*,*) "SetTime "
	!-------
	d = RM_GetTimeStep(id)     ! Deprecated
	d = brm%GetTimeStep()
	status = brm%bmif_get_value("TimeStep", d)
	status = brm%bmif_get_value_ptr("TimeStep", d_ptr)
	write(*,*) "GetTimeStep "
	!-------
	status = RM_SetTimeStep(id, 0.0d0)     ! Deprecated
	status = brm%SetTimeStep(0.0d0)
	status = brm%bmif_set_value("TimeStep", 0.0d0)
	write(*,*) "SetTimeStep "
	!-------
	status = brm%bmif_get_value("Concentrations", DoubleVector)
    allocate(DoubleVector2(nxyz, ncomps))
	status = RM_GetConcentrations(id, DoubleVector2)     ! Deprecated
	status = brm%GetConcentrations(DoubleVector2)
	status = brm%bmif_get_value_ptr("Concentrations", d_ptr)
	write(*,*) "GetConcentrations "
	!-------
	status =RM_SetConcentrations(id, DoubleVector2)     ! Deprecated
	status =brm%SetConcentrations(DoubleVector2)
	status = brm%bmif_set_value("Concentrations", DoubleVector)
	write(*,*) "SetConcentrations "
	!-------
	status = brm%bmif_get_value("DensityCalculated", DoubleVector)
	status = RM_GetDensityCalculated(id, DoubleVector)     ! Deprecated
	status = brm%GetDensityCalculated(DoubleVector)
	status = brm%bmif_get_value_ptr("DensityCalculated", d_ptr)
	write(*,*) "GetDensityCalculated "
	!-------
	status = RM_SetDensityUser(id, DoubleVector)     ! Deprecated
	status = brm%SetDensityUser(DoubleVector)
	status = brm%bmif_set_value("DensityUser", DoubleVector)
	write(*,*) "SetDensityUser "
	!-------
    deallocate(DoubleVector2)
    allocate(DoubleVector2(nxyz, ngas))
	status = RM_GetGasCompMoles(id, DoubleVector2)     ! Deprecated
	status = brm%GetGasCompMoles(DoubleVector2)
	write(*,*) "GetGasCompMoles "
	!-------
	status = RM_SetGasCompMoles(id, DoubleVector2)     ! Deprecated
	status = brm%SetGasCompMoles(DoubleVector2)
	write(*,*) "SetGasCompMoles "
	!-------
	status = RM_GetGasCompPhi(id, DoubleVector2)     ! Deprecated
	status = brm%GetGasCompPhi(DoubleVector2)
	write(*,*) "GetGasCompPhi "
	!-------
	status = RM_GetGasCompPressures(id, DoubleVector2)     ! Deprecated
	status = brm%GetGasCompPressures(DoubleVector2)
	write(*,*) "GetGasCompPressures "
	!-------
    deallocate(DoubleVector)
    allocate(DoubleVector(nxyz))
	status = RM_GetGasPhaseVolume(id, DoubleVector)     ! Deprecated
	status = brm%GetGasPhaseVolume(DoubleVector)
	write(*,*) "GetGasPhaseVolume "
	!-------
	status =RM_SetGasPhaseVolume(id, DoubleVector)     ! Deprecated
	status =brm%SetGasPhaseVolume(DoubleVector)
	write(*,*) "SetGasPhaseVolume "
	!-------
    !do i = 1, RM_GetComponentCount(id)     ! Deprecated
    do i = 1, brm%GetComponentCount()
		status = RM_GetIthConcentration(id, i, DoubleVector)     ! Deprecated
		status = brm%GetIthConcentration(i, DoubleVector)
		!-------
		status = RM_SetIthConcentration(id, i, DoubleVector)     ! Deprecated
		status = brm%SetIthConcentration(i, DoubleVector)
    enddo
	write(*,*)  "GetIthConcentration "
	write(*,*)  "SetIthConcentration "
	!-------
    !do i = 1, RM_GetSpeciesCount(id)     ! Deprecated
    do i = 1, brm%GetSpeciesCount()
		status = RM_GetIthSpeciesConcentration(id, i, DoubleVector)     ! Deprecated
		status = brm%GetIthSpeciesConcentration(i, DoubleVector)
		!-------
		status = RM_SetIthSpeciesConcentration(id, i, DoubleVector)     ! Deprecated
		status = brm%SetIthSpeciesConcentration(i, DoubleVector)
    enddo
    write(*,*) "GetIthSpeciesConcentration "
    write(*,*) "SetIthSpeciesConcentration "
	!-------
	status = brm%bmif_get_value("Porosity", DoubleVector)
	status = RM_GetPorosity(id, DoubleVector)     ! Deprecated
	status = brm%GetPorosity(DoubleVector)
	status = brm%bmif_get_value_ptr("Porosity", d_ptr)
	write(*,*) "GetPorosity "
	!-------
	status = brm%bmif_set_value("Porosity", DoubleVector)
	status = RM_SetPorosity(id, DoubleVector)     ! Deprecated
	status = brm%SetPorosity(DoubleVector)
	write(*,*) "SetPorosity "
	!-------
	status = brm%bmif_get_value("Pressure", DoubleVector)
	status = RM_GetPressure(id, DoubleVector)     ! Deprecated
	status = brm%GetPressure(DoubleVector)
	status = brm%bmif_get_value_ptr("Pressure", d_ptr)
	write(*,*) "GetPressure "
	!-------
	status = brm%bmif_set_value("Pressure", DoubleVector)
	status = RM_SetPressure(id, DoubleVector)     ! Deprecated
	status = brm%SetPressure(DoubleVector)
	write(*,*) "SetPressure "
	!-------
	status = brm%bmif_get_value("SaturationCalculated", DoubleVector)
	status = RM_GetSaturationCalculated(id, DoubleVector)     ! Deprecated
	status = brm%GetSaturationCalculated(DoubleVector)
	status = brm%bmif_get_value_ptr("SaturationCalculated", d_ptr)
	write(*,*) "GetSaturationCalculated "
	!-------
	status = RM_SetSaturationUser(id, DoubleVector)     ! Deprecated
	status = brm%SetSaturationUser(DoubleVector)
	status = brm%bmif_set_value("SaturationUser", DoubleVector)
	write(*,*) "SetSaturationUser "
	!-------
	status = brm%bmif_get_value("SolutionVolume", DoubleVector)
	status = RM_GetSolutionVolume(id, DoubleVector)     ! Deprecated
	status = brm%GetSolutionVolume(DoubleVector)
	status = brm%bmif_get_value_ptr("SolutionVolume", d_ptr)
	write(*,*) "GetSolutionVolume "
	!-------
    deallocate(DoubleVector2)
    allocate(DoubleVector2(nxyz, nspecies))
	status = RM_GetSpeciesConcentrations(id, DoubleVector2)     ! Deprecated
	status = brm%GetSpeciesConcentrations(DoubleVector2)
	write(*,*) "GetSpeciesConcentrations "
	!-------
	status = RM_SpeciesConcentrations2Module(id, DoubleVector2)     ! Deprecated
	status = brm%SpeciesConcentrations2Module(DoubleVector2)
	write(*,*) "SpeciesConcentrations2Module "
	!-------
	status = RM_GetSpeciesLog10Gammas(id, DoubleVector2)     ! Deprecated
	status = brm%GetSpeciesLog10Gammas(DoubleVector2)
	write(*,*) "GetSpeciesLog10Gammas "
	!-------
	status = RM_GetSpeciesLog10Molalities(id, DoubleVector2)     ! Deprecated
	status = brm%GetSpeciesLog10Molalities(DoubleVector2)
	write(*,*) "GetSpeciesLog10Molalities "
	!-------
	status = brm%bmif_get_value("Temperature", DoubleVector)
	status = RM_GetTemperature(id, DoubleVector)     ! Deprecated
	status = brm%GetTemperature(DoubleVector)
	status = brm%bmif_get_value_ptr("Temperature", d_ptr)
	write(*,*) "GetTemperature "
	!-------
	status = RM_SetTemperature(id, DoubleVector)     ! Deprecated
	status = brm%SetTemperature(DoubleVector)
	status = brm%bmif_set_value("Temperature", DoubleVector)
	write(*,*) "SetTemperature "
	!-------
	status = brm%bmif_get_value("Viscosity", DoubleVector)
	status = RM_GetViscosity(id, DoubleVector)     ! Deprecated
	status = brm%GetViscosity(DoubleVector)
	status = brm%bmif_get_value_ptr("Viscosity", d_ptr)	
	write(*,*) "GetViscosity "
	!
	! Take a time step
	!
	status = brm%bmif_update()
	write(*,*) "Update"
	!-------
	status =RM_RunCells(id)     ! Deprecated
	status =brm%RunCells()
	write(*,*) "RunCells"
	!-------
	status = brm%bmif_update_until(86400.0d0)
	write(*,*) "UpdateUntil"
	!
	! Selected output
	!
	status = RM_SetNthSelectedOutput(id, 1)     ! Deprecated
	status = brm%SetNthSelectedOutput(1)
	status = brm%bmif_set_value("NthSelectedOutput", 1)
	write(*,*) "SetNthSelectedOutput "
	!-------
	n_user = RM_GetCurrentSelectedOutputUserNumber(id)     ! Deprecated
	n_user = brm%GetCurrentSelectedOutputUserNumber()
	status = brm%bmif_get_value("CurrentSelectedOutputUserNumber", n_user)
	write(*,*) "GetCurrentSelectedOutputUserNumber "
	!-------
	n = RM_GetNthSelectedOutputUserNumber(id, 1)     ! Deprecated
	n = brm%GetNthSelectedOutputUserNumber(1)
	write(*,*) "GetNthSelectedOutputUserNumber "
	!-------
	status = brm%bmif_get_value("SelectedOutput", DoubleVector)
    if(allocated(DoubleVector2)) deallocate(DoubleVector2)
    n = RM_GetSelectedOutputColumnCount(id)     ! Deprecated
    n = brm%GetSelectedOutputColumnCount()
    allocate(DoubleVector2(nxyz, n))
	status = RM_GetSelectedOutput(id, DoubleVector2)     ! Deprecated
	status = brm%GetSelectedOutput(DoubleVector2)
	write(*,*) "GetSelectedOutput "
	!-------
	n = RM_GetSelectedOutputColumnCount(id)     ! Deprecated
	n = brm%GetSelectedOutputColumnCount()
	status = brm%bmif_get_value("SelectedOutputColumnCount", n)
	write(*,*) "GetSelectedOutputColumnCount "
	!-------
	n = RM_GetSelectedOutputCount(id)     ! Deprecated
	n = brm%GetSelectedOutputCount()
	status = brm%bmif_get_value("SelectedOutputCount", n)
	write(*,*) "GetSelectedOutputCount "
	!-------
	status = brm%bmif_get_value("SelectedOutputHeadings", StringVector)
	status = RM_GetSelectedOutputHeadings(id, StringVector)     ! Deprecated
	status = brm%GetSelectedOutputHeadings(StringVector)
	write(*,*) "GetSelectedOutputHeadings "
	!-------
	!b = RM_GetSelectedOutputOn(id)     ! Deprecated
	!b = brm%GetSelectedOutputOn()
	status = brm%bmif_get_value("SelectedOutputOn", l)
	status = brm%bmif_get_value_ptr("SelectedOutputOn", b_ptr)	
	write(*,*) "GetSelectedOutputOn "
	!-------
	n = RM_GetSelectedOutputRowCount(id)     ! Deprecated
	n = brm%GetSelectedOutputRowCount()
	status = brm%bmif_get_value("SelectedOutputRowCount", n)
	write(*,*) "GetSelectedOutputRowCount "
	!-------
	status = RM_SetCurrentSelectedOutputUserNumber(id, 333)     ! Deprecated
	status = brm%SetCurrentSelectedOutputUserNumber(333)
	write(*,*) "SetCurrentSelectedOutputUserNumber "
	!
	! Getters
	!
	status = RM_GetBackwardMapping(id, 1, IntVector)     ! Deprecated
	status = brm%GetBackwardMapping(1, IntVector)
	write(*,*) "GetBackwardMapping "
	!-------
	!status = RM_GetDatabaseFileName(id, string)     ! Deprecated
	!status = brm%GetDatabaseFileName(string)
	!write(*,*) "GetDatabaseFileName "
	!-------
	status = RM_GetEndCell(id, IntVector)     ! Deprecated
	status = brm%GetEndCell(IntVector)
	write(*,*) "GetEndCell"
	!-------
	!n = RM_GetErrorHandlerMode()     ! Deprecated
	!n = brm%GetErrorHandlerMode()
	!write(*,*) "GetErrorHandlerMode "
	!-------
	status = RM_GetErrorString(id, AllocString)     ! Deprecated
	status = brm%bmif_get_value("ErrorString", AllocString)
	write(*,*) "GetErrorString "
	!-------
	status = RM_GetFilePrefix(id, AllocString)     ! Deprecated
	status = brm%GetFilePrefix(AllocString)
	status = brm%bmif_get_value("FilePrefix", AllocString)
	write(*,*) "GetFilePrefix "
	!-------
	!status = RM_GetForwardMapping()  ! not implemented
	!status = brm%GetForwardMapping()  ! not implemented
	!write(*,*) "GetForwardMapping "
	!-------
	n = RM_GetIPhreeqcID(id, 0)     ! Deprecated
    n = brm%GetIPhreeqcID(0)
	write(*,*) "GetIPhreeqcID "
	!-------
	n = RM_GetMpiMyself(id)     ! Deprecated
	n = brm%GetMpiMyself()
	write(*,*) "GetMpiMyself "
	!-------
	n = RM_GetMpiTasks(id)     ! Deprecated
	n = brm%GetMpiTasks()
	write(*,*) "GetMpiTasks "
	!-------
	!status = RM_GetPartitionUZSolids(id, n)  ! not implemented
	!status = brm%GetPartitionUZSolids(n)  ! Not implemented
	!write(*,*) "GetPartitionUZSolids "
	!-------
	!status = RM_GetPrintChemistryMask(id, n)   ! not implemented
	!status = brm%GetPrintChemistryMask(n)   ! Not implemented
	!write(*,*) "GetPrintChemistryMask "
	!-------
	!status = RM_GetPrintChemistryOn(id, n)  ! not implemented
	!status = brm%GetPrintChemistryOn(n)  ! Not implemented
	!write(*,*) "GetPrintChemistryOn "
	!-------
	!status = RM_GetRebalanceByCell(id, n)   ! not implemented
	!status = brm%GetRebalanceByCell(n)   ! Not implemented
	!write(*,*) "GetRebalanceByCell "
	!-------
	!status = RM_GetRebalanceFraction(id, d)  ! not implemented
	!status = brm%GetRebalanceFraction(d)  ! Not implemented
	!write(*,*) "GetRebalanceFraction "
	!-------
	!status = RM_GetSpeciesSaveOn(id, n)       ! not implemented
	!status = brm%GetSpeciesSaveOn(n)       ! Not implemented
	!write(*,*) "GetSpeciesSaveOn "
	!-------
	!status = RM_GetSpeciesStoichiometry(id, i, IntVector)   ! not implemented
	!status = RM_GetSpeciesStoichiometry(id, i, IntVector)   ! not implemented
	!write(*,*) "GetSpeciesStoichiometry "
	!-------
	status = RM_GetStartCell(id, IntVector)     ! Deprecated
	status = brm%GetStartCell(IntVector)
	write(*,*) "GetStartCell "
	!-------
	d = RM_GetTimeConversion(id)     ! Deprecated
	d = brm%GetTimeConversion()
	write(*,*) "GetTimeConversion "
	!-------
	!status = RM_GetUnitsExchange(id, n)   ! not implemented
	!status = brm%GetUnitsExchange(n)   ! Not implemented
	!write(*,*) "GetUnitsExchange "
	!-------
	!status = RM_GetUnitsGasPhase(id, n)   ! not implemented
	!status = brm%GetUnitsGasPhase(n)   ! Not implemented
	!write(*,*) "GetUnitsGasPhase "
	!-------
	!status = RM_GetUnitsKinetics(id, n)   ! not implemented
	!status = brm%GetUnitsKinetics(n)   ! Not implemented
	!write(*,*) "GetUnitsKinetics "
	!-------
	!status = RM_GetUnitsPPassemblage(id, n)   ! not implemented
	!status = brm%GetUnitsPPassemblage(n)   ! Not implemented
	!write(*,*) "GetUnitsPPassemblage "
	!-------
	!status = RM_GetUnitsSolution(id, n)   ! not implemented
	!status = brm%GetUnitsSolution(n)   ! Not implemented
	!write(*,*) "GetUnitsSolution "
	!-------
	!status = RM_GetUnitsSSassemblage(id, n)   ! not implemented
	!status = brm%GetUnitsSSassemblage(n)   ! Not implemented
	!write(*,*) "GetUnitsSSassemblage "
	!-------
	!status = RM_GetUnitsSurface(id, n)   ! not implemented
	!status = brm%GetUnitsSurface(n)   ! Not implemented
	!write(*,*) "GetUnitsSurface "
	!-------
	!std::vector<IPhreeqcPhast *> w = status = RM_GetWorkers()   ! not implemented
	!std::vector<IPhreeqcPhast *> w = status = RM_GetWorkers()   ! not implemented
	!write(*,*) "GetWorkers "
	!
	! Utilities
	!
#ifndef USE_MPI    
	!id1 = brm%bmif_create(10, 1)  ! make another bmiphreeqcrm
	id1 = brm1%bmif_create(10,1)
	status = brm1%CloseFiles() 
	status = brm1%bmif_finalize()   ! destroy the new bmiphreeqcrm
	write(*,*) "CloseFiles "
#endif    
	!-------
	deallocate(IntVector)
	allocate(IntVector(1))
	IntVector = 1
	n = size(IntVector)
    deallocate(bc2)
	allocate(bc2(n,ncomps))
	status = RM_InitialPhreeqc2Concentrations(id, bc2, size(IntVector), IntVector)     ! Deprecated
	status = RM_InitialPhreeqc2Concentrations(id, bc2, size(IntVector), IntVector)     ! Deprecated
	n = size(IntVector)
	allocate(tc(n), p_atm(n))
	tc = 30
	p_atm = 1.5
	status = RM_Concentrations2Utility(id, bc2, n, tc, p_atm)     ! Deprecated
	status = brm%Concentrations2Utility(bc2, n, tc, p_atm)
	write(*,*) "Concentrations2Utility "
	!-------
	status = RM_DecodeError(id, -2)	        ! Deprecated
	status = brm%DecodeError(-2)	   
	write(*,*) "DecodeError "
	!-------
	status =RM_DumpModule(id, 1, 0)     ! Deprecated
	status =brm%DumpModule(1, 0)
	write(*,*) "DumpModule "
	!-------
	!status = RM_ErrorHandler(id, 0, "string")      ! Deprecated
	!status = brm%ErrorHandler(0, "string") 
	!write(*,*) "OK, just a test: ErrorHandler "
	!-------
	status = RM_ErrorMessage(id, "my error")       ! Deprecated
	status = brm%ErrorMessage("my error")  
	write(*,*) "OK, just a test: ErrorMessage "
	!-------
	status = RM_LogMessage(id, "Log message")      ! Deprecated
	status = brm%LogMessage("Log message") 
	write(*,*) "LogMessage "
	!-------
	status = RM_OutputMessage(id, "Output message")      ! Deprecated
	status = brm%OutputMessage("Output message") 
	write(*,*) "OutputMessage "
	!-------
	status = RM_ScreenMessage(id, "Screen message")      ! Deprecated
	status = brm%ScreenMessage("Screen message") 
	write(*,*) "ScreenMessage "
	!-------
	status =RM_StateSave(id, 1)     ! Deprecated
	status =brm%StateSave(1)
	write(*,*) "StateSave "
	!-------
	status =RM_StateApply(id, 1)     ! Deprecated
	status =brm%StateApply(1)
	write(*,*) "StateApply "
	!-------
	status =RM_StateDelete(id, 1)     ! Deprecated
	status =brm%StateDelete(1)
	write(*,*) "StateDelete "
	!-------
	status = RM_WarningMessage(id, "Warning message")       ! Deprecated
	status = brm%WarningMessage("Warning message")  
	write(*,*) "WarningMessage "
	!
	! BMI Methods
	!
	status = brm%bmif_get_component_name(string)
	write(*,*) "brm%bmif_get_component_name "
	!-------
	status = brm%bmif_get_current_time(d)
	write(*,*) "brm%bmif_get_current_time "
	!-------
	status = brm%bmif_get_end_time(d)
	write(*,*) "brm%bmif_get_end_time "
	!-------
	status = brm%bmif_grid_rank(0, n)
	write(*,*) "brm%bmif_grid_rank "
	!-------
	status = brm%bmif_grid_size(0, n)
	write(*,*) "brm%bmif_grid_size "
	!-------
	status = brm%bmif_grid_type(0, string)
	write(*,*) "brm%bmif_grid_type "
	!-------
	status = brm%bmif_get_input_item_count(n)
	write(*,*) "brm%bmif_get_input_item_count "
	!-------
	status = brm%bmif_get_input_var_names(StringVector)
	write(*,*) "brm%bmif_get_input_var_names "
	!-------
	status = brm%bmif_get_output_item_count(n)
	write(*,*) "brm%bmif_get_output_item_count "
	!-------
	status = brm%bmif_get_output_var_names(StringVector)
	write(*,*) "brm%bmif_get_output_var_names "
	!-------
	status = brm%bmif_get_pointable_item_count(n)
	write(*,*) "brm%bmif_get_pointable_item_count "
	!-------
	status = brm%bmif_get_pointable_var_names(StringVector)
	write(*,*) "brm%bmif_get_pointable_var_names "
	!-------
	status = brm%bmif_get_time_step(d)
	write(*,*) "brm%bmif_get_time_step "
	!-------
	status = brm%bmif_get_time_units(string)
	write(*,*) "brm%bmif_get_time_units "
	!-------
	status = brm%bmif_get_value("solution_saturation_index_Calcite", DoubleVector)
	write(*,*) "brm%bmif_get_value "
	!-------
	status = brm%bmif_get_var_itemsize("solution_saturation_index_Calcite", n)
	write(*,*) "brm%bmif_get_var_itemsize "
	!-------
	status = brm%bmif_get_var_nbytes("solution_saturation_index_Calcite", n)
	write(*,*) "brm%bmif_get_var_nbytes "
	!-------
	status = brm%bmif_get_var_type("solution_saturation_index_Calcite", string)
	write(*,*) "brm%bmif_get_var_type "
	!-------
	status = brm%bmif_get_var_units("solution_saturation_index_Calcite", string)
	write(*,*) "brm%bmif_get_var_units "
	!status = brm%bmif_initialize(YAML_filename)
	! See above
	status = brm%bmif_set_value("Time", 1.0d0) 
	write(*,*) "brm%bmif_set_value"
	!-------
	status = brm%bmif_update() 
	write(*,*) "brm%bmif_update"
	!-------
	status = brm%bmif_update_until(864000.0d0) 
	write(*,*) "brm%bmif_update_until"
	!-------
    
	write(*,*) "AddOutputVars"
	status = brm%bmif_get_output_var_names(Names)
	do i = 1, size(StringVector)
		status = brm%bmif_get_var_itemsize(Names(i), itemsize)
		status = brm%bmif_get_var_nbytes(Names(i), nbytes)
		status = brm%bmif_get_var_type(Names(i), string)
		if (itemsize .eq. 0) itemsize=1
		if (nbytes .eq. 0) nbytes=1
		dim = nbytes / itemsize
		status = brm%bmif_get_var_type(Names(i), vtype)
		if (vtype .eq. "real(kind=8)") then
			if (dim .eq. 1) then
				status = brm%bmif_get_value(Names(i), d)
                write(*,*) "     ", Names(i), "  ", d
			else
				status = brm%bmif_get_value(Names(i), DoubleVector)
                write(*,*) "     ", Names(i), "  ", DoubleVector(1)
			endif
		else if (vtype .eq. "integer") then
			if (dim .eq. 1) then
				status = brm%bmif_get_value(Names(i), j)
                write(*,*) "     ", Names(i), "  ", j
			else
				status = brm%bmif_get_value(Names(i), IntVector)
                write(*,*) "     ", Names(i), "  ", IntVector(1)
			endif
		else if (vtype .eq. "logical") then
			if (dim == 1) then
				status = brm%bmif_get_value(Names(i), l)
                write(*,*) "     ", Names(i), "  ", l
			endif
		else if (vtype(1:9) .eq. "character") then
			if (dim == 1) then
				status = brm%bmif_get_value(Names(i), string)
                write(*,*) "     ", Names(i), "  ", trim(string)
			else
				status = brm%bmif_get_value(Names(i), StringVector)
                write(*,*) "     ", Names(i), "  ", trim(StringVector(1))
			endif
		endif
    enddo    
    
    ! Not implemented
	status = brm%bmif_get_var_location(string, string)
	status = brm%bmif_get_value_float(string, float)
	status = brm%bmif_get_value_ptr_float(string, float_ptr)
	status = brm%bmif_get_value_at_indices_double(string, DoubleVector, IntVector)
	status = brm%bmif_get_value_at_indices_float(string, FloatVector, IntVector)
	status = brm%bmif_get_value_at_indices_int(string, IntVector, IntVector)
	status = brm%bmif_set_value_at_indices_double(string, DoubleVector, IntVector)
	status = brm%bmif_set_value_at_indices_float(string, FloatVector, IntVector)
	status = brm%bmif_set_value_at_indices_int(string, IntVector, IntVector)
	status = brm%bmif_set_value_float(string, float)
	status = brm%bmif_get_grid_shape(n, IntVector)
	status = brm%bmif_get_grid_spacing(n, DoubleVector)
	status = brm%bmif_get_grid_origin(n, DoubleVector)
	status = brm%bmif_get_grid_x(n, DoubleVector)
	status = brm%bmif_get_grid_y(n, DoubleVector)
	status = brm%bmif_get_grid_z(n, DoubleVector)
	status = brm%bmif_get_grid_node_count(n, n)
	status = brm%bmif_get_grid_edge_count(n, n)
	status = brm%bmif_get_grid_face_count(n, n)
	status = brm%bmif_get_grid_edge_nodes(n, IntVector)
	status = brm%bmif_get_grid_face_edges(n, IntVector)
	status = brm%bmif_get_grid_face_nodes(n, IntVector)
	status = brm%bmif_get_grid_nodes_per_face(n, IntVector)
#ifdef USE_MPI
!	status = RM_MpiWorkerBreak(id)     ! Deprecated
	status = brm%MpiWorkerBreak()
#endif    
	status = brm%bmif_finalize()    ! void method
	write(*,*) "brm%bmif_finalize "

	!TODO status =RM_MpiAbort()
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
