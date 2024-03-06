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
    type(YAML_PhreeqcRM) :: yrm
    ! Create YAMLPhreeqcRM document
    id = yrm%CreateYAMLPhreeqcRM()   
    ! Number of cells
    nxyz = 40;
	! Set GridCellCount
	status = yrm%YAMLSetGridCellCount(nxyz)
	status = yrm%YAMLThreadCount(3)
	! Set some properties
	status = yrm%YAMLSetErrorHandlerMode(1)
	status = yrm%YAMLSetComponentH2O(.false.)
	status = yrm%YAMLSetRebalanceFraction(0.5d0)
	status = yrm%YAMLSetRebalanceByCell(.true.)
	status = yrm%YAMLUseSolutionDensityVolume(.false.)
	status = yrm%YAMLSetPartitionUZSolids(.false.)
    status = yrm%YAMLSetFilePrefix("AdvectBMI_f90")
    status = yrm%YAMLOpenFiles()
    ! Set concentration units
	status = yrm%YAMLSetUnitsSolution(2)           ! 1, mg/L; 2, mol/L; 3, kg/kgs
	status = yrm%YAMLSetUnitsPPassemblage(1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = yrm%YAMLSetUnitsExchange(1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = yrm%YAMLSetUnitsSurface(1)            ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = yrm%YAMLSetUnitsGasPhase(1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = yrm%YAMLSetUnitsSSassemblage(1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = yrm%YAMLSetUnitsKinetics(1)           ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock

	! Set conversion from seconds to user units (days) Only affects one print statement
	time_conversion = 1.0d0 / 86400.0d0
	status = yrm%YAMLSetTimeConversion(time_conversion)
    
	! Set representative volume
    allocate(rv(nxyz))
    rv = 1.0d0
	status = yrm%YAMLSetRepresentativeVolume(rv)
	! Set initial density
    allocate(density(nxyz))
    density = 1.0d0
	status = yrm%YAMLSetDensityUser(density)
    ! Set initial porosity
    allocate(por(nxyz))
    por = 0.2d0
	status = yrm%YAMLSetPorosity(por)
	! Set initial saturation
    allocate(sat(nxyz))
    sat = 1.0d0
	status = yrm%YAMLSetSaturationUser(sat)   
	! Set cells to print chemistry when print chemistry is turned on
    allocate(print_chemistry_mask(nxyz))
    print_chemistry_mask = 0
	do i = 1, nxyz / 2 
		print_chemistry_mask(i) = 1
	enddo
	status = yrm%YAMLSetPrintChemistryMask(print_chemistry_mask)  
	! Demonstation of mapping, two equivalent rows by symmetry
    ! zero-based indexing
    allocate(grid2chem(nxyz))
    grid2chem = -1
	do i = 1, nxyz / 2 
		grid2chem(i) = i - 1
		grid2chem(i + nxyz / 2) = i - 1
	enddo
	status = yrm%YAMLCreateMapping(grid2chem)
	! Set printing of chemistry file
	status = yrm%YAMLSetPrintChemistryOn(.false., .true., .false.) ! workers, initial_phreeqc, utility
	! Load database
	status = yrm%YAMLLoadDatabase("phreeqc.dat")    
    ! Run file to define solutions and reactants for initial conditions, selected output
	workers = .true.             ! Worker instances do the reaction calculations for transport
	initial_phreeqc = .true.     ! InitialPhreeqc instance accumulates initial and boundary conditions
	utility = .true.             ! Utility instance is available for processing
	status = yrm%YAMLRunFile(workers, initial_phreeqc, utility, "advect.pqi")
	! Clear contents of workers and utility
	initial_phreeqc = .false.
	input = "DELETE; -all"
	status = yrm%YAMLRunString(workers, initial_phreeqc, utility, input)
	! Determine number of components to transport
	status = yrm%YAMLFindComponents()
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
	status = yrm%YAMLInitialPhreeqc2Module_mix(ic1, ic2, f1)    
	! No mixing is defined, so the following is equivalent
	!status = yrm%YAMLInitialPhreeqc2Module(id, ic1)

	! alternative for setting initial conditions
	! cell number in first argument (-1 indicates last solution, 40 in this case)
	! in advect.pqi and any reactants with the same number--
	! Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
	! will be written to cells 18 and 19 (0 based)
    allocate(module_cells(2))
	module_cells(1) = 18
	module_cells(2) = 19
	status = yrm%YAMLInitialPhreeqcCell2Module(-1, module_cells)
	! Initial equilibration of cells
	time_step = 0.0d0  ! no kinetics
	status = yrm%YAMLSetTimeStep(time_step)
	time = 0.0d0
	status = yrm%YAMLSetTime(time)
	status = yrm%YAMLRunCells()
    time_step = 86400.0d0 
	status = yrm%YAMLSetTimeStep(time_step)    

	! Write YAML file
    YAML_filename = "AdvectBMI_f90.yaml"
	status = yrm%WriteYAMLDoc(YAML_filename)
	status = yrm%YAMLClear()  
    status = Yrm%DestroyYAMLPhreeqcRM()
    status = 0
    end subroutine WriteYAMLFile_f90
#endif     