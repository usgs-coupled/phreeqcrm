#ifdef USE_YAML
    subroutine WriteYAMLFile_f90()  BIND(C, NAME='WriteYAMLFile_f90')
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
    integer,          allocatable, dimension(:,:) :: ic1
    integer,          allocatable, dimension(:,:) :: ic2
    real(kind=8), allocatable, dimension(:,:) :: f1
    integer,          allocatable, dimension(:)   :: module_cells
    ! Create YAMLPhreeqcRM document
    id = CreateYAMLPhreeqcRM()   
    ! Number of cells
    nxyz = 40;
	! Set GridCellCount
	status = YAMLSetGridCellCount(id, nxyz)
	status = YAMLThreadCount(id, 3)
	! Set some properties
	status = YAMLSetErrorHandlerMode(id, 1)
	status = YAMLSetComponentH2O(id, .false.)
	status = YAMLSetRebalanceFraction(id, 0.5d0)
	status = YAMLSetRebalanceByCell(id, .true.)
	status = YAMLUseSolutionDensityVolume(id, .false.)
	status = YAMLSetPartitionUZSolids(id, .false.)
    status = YAMLSetFilePrefix(id, "AdvectBMI_f90")
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
	status = YAMLRunFile(id, workers, initial_phreeqc, utility, "advect.pqi")
	! Clear contents of workers and utility
	initial_phreeqc = .false.
	input = "DELETE; -all"
	status = YAMLRunString(id, workers, initial_phreeqc, utility, input)
	! Determine number of components to transport
	status = YAMLFindComponents(id)
	! set array of initial conditions
    allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
    ic1 = -1
    ic2 = -1
    f1 = 1.0d0
	do i = 1, nxyz
		ic1(i,1) = 1     ! Solution 1
		ic1(i,2) = -1    ! Equilibrium phases none
		ic1(i,3) = 1     ! Exchange 1
		ic1(i,4) = -1    ! Surface none
		ic1(i,5) = -1    ! Gas phase none
		ic1(i,6) = -1    ! Solid solutions none
		ic1(i,7) = -1    ! Kinetics none
	enddo
	status = YAMLInitialPhreeqc2Module_mix(id, ic1, ic2, f1)    
	! No mixing is defined, so the following is equivalent
	!status = YAMLInitialPhreeqc2Module(id, id, ic1)

	! alternative for setting initial conditions
	! cell number in first argument (id, -1 indicates last solution, 40 in this case)
	! in advect.pqi and any reactants with the same number--
	! Equilibrium phases, exchange, surface, gas phase, solid solution, and (id, or) kinetics--
	! will be written to cells 18 and 19 (id, 0 based)
    allocate(module_cells(2))
	module_cells(1) = 18
	module_cells(2) = 19
	status = YAMLInitialPhreeqcCell2Module(id, -1, module_cells)
	! Initial equilibration of cells
	time_step = 0.0d0  ! no kinetics
	status = YAMLSetTimeStep(id, time_step)
	time = 0.0d0
	status = YAMLSetTime(id, time)
	status = YAMLRunCells(id)
    time_step = 86400.0d0 
	status = YAMLSetTimeStep(id, time_step)    

	! Write YAML file
    YAML_filename = "AdvectBMI_f90.yaml"
	status = WriteYAMLDoc(id, YAML_filename)
	status = YAMLClear(id)  
    status = DestroyYAMLPhreeqcRM(id)
    status = 0
    end subroutine WriteYAMLFile_f90
#endif     