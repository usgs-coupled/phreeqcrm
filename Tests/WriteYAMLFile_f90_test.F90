#ifdef USE_YAML
    subroutine WriteYAMLFile_f90_test()  BIND(C, NAME='WriteYAMLFile_f90_test')
    USE, intrinsic :: ISO_C_BINDING
    USE YAMLPhreeqcRM
    implicit none
    integer :: id
    integer :: status, i
    integer :: nxyz
    character(len=100) :: input, YAML_filename
    real(kind=8) :: time_conversion
    logical :: workers, initial_phreeqc, utility
    real(kind=8) :: time, time_step
    real(kind=8), allocatable, dimension(:)   :: density
    real(kind=8), allocatable, dimension(:)   :: rv
    real(kind=8), allocatable, dimension(:)   :: por
    real(kind=8), allocatable, dimension(:)   :: sat
    integer,          allocatable, dimension(:)   :: print_chemistry_mask
    integer,          allocatable, dimension(:)   :: grid2chem
    integer,          allocatable, dimension(:)   :: ic
    ! Create YAMLPhreeqcRM document
    id = CreateYAMLPhreeqcRM()   
    ! Number of cells
    nxyz = 40;
	! Set GridCellCount
	status = YAMLSetGridCellCount(id, nxyz)
    status = YAMLThreadCount(id, 4)
	! Set some properties
	status = YAMLSetErrorHandlerMode(id, 1)
	status = YAMLSetComponentH2O(id, .false.)
	status = YAMLSetRebalanceFraction(id, 0.5d0)
	status = YAMLSetRebalanceByCell(id, .true.)
	status = YAMLUseSolutionDensityVolume(id, .false.)
	status = YAMLSetPartitionUZSolids(id, .false.)
    status = YAMLSetFilePrefix(id, "AdvectBMI_f90_test")
    status = YAMLOpenFiles(id)
    ! Set concentration units
	status = YAMLSetUnitsSolution(id, 2)           ! 1, mg/L; 2, mol/L; 3, kg/kgs
	status = YAMLSetUnitsPPassemblage(id, 1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = YAMLSetUnitsExchange(id, 1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = YAMLSetUnitsSurface(id, 1)            ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = YAMLSetUnitsGasPhase(id, 1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = YAMLSetUnitsSSassemblage(id, 1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = YAMLSetUnitsKinetics(id, 1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock

	! Set conversion from seconds to user units (days) Only affects one print statement
	time_conversion = 1.0d0 / 86400.0d0
	status = YAMLSetTimeConversion(id, time_conversion)
    
	! Set representative volume
    allocate(rv(nxyz))
    rv = 1.0d0
	status = YAMLSetRepresentativeVolume(id, rv)
	! Set initial density
    allocate(density(nxyz))
    density = 1.0d0
	status = YAMLSetDensityUser(id, density)
    ! Set initial porosity
    allocate(por(nxyz))
    por = 0.2d0
	status = YAMLSetPorosity(id, por)
	! Set initial saturation
    allocate(sat(nxyz))
    sat = 1.0d0
	status = YAMLSetSaturationUser(id, sat)   
	! Set cells to print chemistry when print chemistry is turned on
    allocate(print_chemistry_mask(nxyz))
    print_chemistry_mask = 0
	do i = 1, nxyz / 2 
		print_chemistry_mask(i) = 1
	enddo
	status = YAMLSetPrintChemistryMask(id, print_chemistry_mask)  
	! Demonstation of mapping, two equivalent rows by symmetry
    ! zero-based indexing
    allocate(grid2chem(nxyz))
    grid2chem = -1
	do i = 1, nxyz / 2 
		grid2chem(i) = i - 1
		grid2chem(i + nxyz / 2) = i - 1
	enddo
	status = YAMLCreateMapping(id, grid2chem)
	! Set printing of chemistry file
	status = YAMLSetPrintChemistryOn(id, .false., .true., .false.) ! workers, initial_phreeqc, utility
	! Load database
	status = YAMLLoadDatabase(id, "phreeqc.dat")    
    ! Run file to define solutions and reactants for initial conditions, selected output
	workers = .true.             ! Worker instances do the reaction calculations for transport
	initial_phreeqc = .true.     ! InitialPhreeqc instance accumulates initial and boundary conditions
	utility = .true.             ! Utility instance is available for processing
	status = YAMLRunFile(id, workers, initial_phreeqc, utility, "all_reactants.pqi")
	! Clear contents of workers and utility
	initial_phreeqc = .false.
	input = "DELETE; -all"
	status = YAMLRunString(id, workers, initial_phreeqc, utility, input)
    
	status = YAMLAddOutputVars(id, "AddOutputVars", "true")
	status = YAMLAddOutputVars(id, "SolutionProperties", "true")
	status = YAMLAddOutputVars(id, "SolutionTotalMolalities", "true")
	status = YAMLAddOutputVars(id, "ExchangeMolalities", "true")
	status = YAMLAddOutputVars(id, "SurfaceMolalities", "true")
	status = YAMLAddOutputVars(id, "EquilibriumPhases", "true")
	status = YAMLAddOutputVars(id, "Gases", "true")
	status = YAMLAddOutputVars(id, "KineticReactants", "true")
	status = YAMLAddOutputVars(id, "SolidSolutions", "true")
	status = YAMLAddOutputVars(id, "CalculateValues", "true")
	status = YAMLAddOutputVars(id, "SolutionActivities", "false")
	status = YAMLAddOutputVars(id, "SolutionMolalities", "false")
	status = YAMLAddOutputVars(id, "SaturationIndices", "false")
	! Determine number of components to transport
	status = YAMLFindComponents(id)
	! set array of initial conditions
    allocate(ic(nxyz))
    ic = 1
    status = YAMLInitialSolutions2Module(id, ic)
    status = YAMLInitialEquilibriumPhases2Module(id, ic)
    status = YAMLInitialExchanges2Module(id, ic)
    status = YAMLInitialGasPhases2Module(id, ic)
    status = YAMLInitialKinetics2Module(id, ic)
    status = YAMLInitialSolidSolutions2Module(id, ic)
    status = YAMLInitialSurfaces2Module(id, ic)
	! Initial equilibration of cells
	time_step = 0.0d0  ! no kinetics
	status = YAMLSetTimeStep(id, time_step)
	time = 0.0d0
	status = YAMLSetTime(id, time)
	status = YAMLRunCells(id)
    time_step = 86400.0d0 
	status = YAMLSetTimeStep(id, time_step)    

	! Write YAML file
    YAML_filename = "AdvectBMI_f90_test.yaml"
	status = WriteYAMLDoc(id, YAML_filename)
	status = YAMLClear(id)  
    status = DestroyYAMLPhreeqcRM(id)
    status = 0
    end subroutine WriteYAMLFile_f90_test
#endif     