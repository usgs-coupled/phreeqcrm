

module mydata
  double precision, dimension(:), pointer :: K_ptr
  integer                                 :: rm_id
end module mydata
    
subroutine advection_f90()  BIND(C)
  USE, intrinsic :: ISO_C_BINDING
  USE PhreeqcRM
  USE IPhreeqc
  USE mydata
  implicit none
#ifdef USE_MPI    
  INCLUDE 'mpif.h'
#endif
    interface
        subroutine advect_f90(c, bc_conc, ncomps, nxyz)
            implicit none
            double precision, dimension(:,:), allocatable, intent(inout) :: c 
            double precision, dimension(:,:), allocatable, intent(in) :: bc_conc
            integer, intent(in)                                       :: ncomps, nxyz
        end subroutine advect_f90
        integer function do_something()
        end function do_something
        integer(kind=C_INT) function worker_tasks_f(method_number) BIND(C, NAME='worker_tasks_f')
            USE ISO_C_BINDING
            implicit none
            integer(kind=c_int), intent(in) :: method_number
        end function worker_tasks_f
        SUBROUTINE register_basic_callback_fortran()
            implicit none
        END SUBROUTINE register_basic_callback_fortran
    end interface

  ! Based on PHREEQC Example 11
  
  integer :: mpi_myself
  integer :: i, j
  integer :: nxyz
  integer :: nthreads
  integer :: id
  integer :: status
  double precision, dimension(:), allocatable, target :: hydraulic_K
  double precision, dimension(:), allocatable   :: rv
  double precision, dimension(:), allocatable   :: por
  double precision, dimension(:), allocatable   :: sat
  integer,          dimension(:), allocatable   :: print_chemistry_mask
  integer,          dimension(:), allocatable   :: grid2chem
  integer                                       :: nchem
  character(100)                                :: string
  character(200)                                :: string1
  integer                                       :: ncomps, ncomps1
  character(100),   dimension(:), allocatable   :: components
  double precision, dimension(:), allocatable   :: gfw
  integer,          dimension(:,:), allocatable :: ic1, ic2
  double precision, dimension(:,:), allocatable :: f1
  integer                                       :: nbound
  integer,          dimension(:), allocatable   :: bc1, bc2
  double precision, dimension(:), allocatable   :: bc_f1
  integer,          dimension(:), allocatable   :: module_cells
  double precision, dimension(:,:), allocatable :: bc_conc
  double precision, dimension(:,:), allocatable :: c
  double precision                              :: time, time_step
  double precision, dimension(:), allocatable   :: density
  double precision, dimension(:), allocatable   :: sat_calc
  double precision, dimension(:), allocatable   :: volume
  double precision, dimension(:), allocatable   :: temperature
  double precision, dimension(:), allocatable   :: pressure
  integer                                       :: isteps, nsteps
  double precision, dimension(:,:), allocatable :: selected_out
  integer                                       :: col, isel, n_user
  character(100)                                :: heading
  double precision, dimension(:,:), allocatable :: c_well
  double precision, dimension(:), allocatable   :: tc, p_atm
  integer                                       :: vtype
  double precision                              :: pH
  character(100)                                :: svalue
  integer                                       :: iphreeqc_id, iphreeqc_id1
  integer                                       :: dump_on, append
  !character(LEN=1), dimension(:), allocatable   :: errstr
#ifdef FORTRAN_2003
  character(LEN=:), allocatable                 :: errstr
#else
  character(LEN=10000)                          :: errstr
#endif
  integer                                       :: l, n
  integer, dimension(:), allocatable            :: sc, ec
  ! --------------------------------------------------------------------------
  ! Create PhreeqcRM
  ! --------------------------------------------------------------------------

  nxyz = 40
  ! Bogus conductivity field for Basic callback demonstration
  allocate(hydraulic_K(nxyz))
  do i = 1, nxyz
      hydraulic_K(i) = i * 2.0
  enddo
  K_ptr => hydraulic_K
#ifdef USE_MPI
  ! MPI
  id = RM_Create(nxyz, MPI_COMM_WORLD)
  rm_id = id
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
  if (status .ne. MPI_SUCCESS) then
     stop "Failed to get mpi_myself"
  endif
  if (mpi_myself > 0) then
      status = RM_SetMpiWorkerCallback(id, worker_tasks_f)
      status = RM_MpiWorker(id)
      status = RM_Destroy(id)
      return
  endif
#else
  ! OpenMP
  nthreads = 3
  id = RM_Create(nxyz, nthreads)
  rm_id = id
#endif
  ! Set properties
  status = RM_SetErrorHandlerMode(id, 2)  ! exit on error
  status = RM_SetComponentH2O(id, 0)
  status = RM_SetRebalanceFraction(id, 0.5d0)
  status = RM_SetRebalanceByCell(id, 1)
  status = RM_UseSolutionDensityVolume(id, 0)
  status = RM_SetPartitionUZSolids(id, 0)
  ! Open files
  status = RM_SetFilePrefix(id, "Advect_f90")
  status = RM_OpenFiles(id)
#ifdef USE_MPI
  ! Optional callback for MPI
  status = do_something()   ! only root is calling do_something here
#endif
  ! Set concentration units
  status = RM_SetUnitsSolution(id, 2)      ! 1, mg/L; 2, mol/L; 3, kg/kgs
  status = RM_SetUnitsPPassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  status = RM_SetUnitsExchange(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  status = RM_SetUnitsSurface(id, 1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  status = RM_SetUnitsGasPhase(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  status = RM_SetUnitsSSassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  status = RM_SetUnitsKinetics(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  ! Set conversion from seconds to user units (days)
  status = RM_SetTimeConversion(id, dble(1.0 / 86400.0))
  ! Set representative volume
  allocate(rv(nxyz))
  rv = 1.0
  status = RM_SetRepresentativeVolume(id, rv)
  ! Set initial porosity
  allocate(por(nxyz))
  por = 0.2
  status = RM_SetPorosity(id, por)
  ! Set initial saturation
  allocate(sat(nxyz))
  sat = 1.0
  status = RM_SetSaturation(id, sat)
  ! Set cells to print chemistry when print chemistry is turned on
  allocate(print_chemistry_mask(nxyz))
  do i = 1, nxyz/2
     print_chemistry_mask(i) = 1
     print_chemistry_mask(i+nxyz/2) = 0
  enddo
  status = RM_SetPrintChemistryMask(id, print_chemistry_mask)
  ! Demonstation of mapping, two equivalent rows by symmetry
  allocate(grid2chem(nxyz))
  do i = 1, nxyz/2
     grid2chem(i) = i - 1
     grid2chem(i+nxyz/2) = i - 1
  enddo
  status = RM_CreateMapping(id, grid2chem)  
  if (status < 0) status = RM_DecodeError(id, status) 
  nchem = RM_GetChemistryCellCount(id)

  ! --------------------------------------------------------------------------
  ! Set initial conditions
  ! --------------------------------------------------------------------------

  ! Set printing of chemistry file to false
  status = RM_SetPrintChemistryOn(id, 0, 1, 0)  ! workers, initial_phreeqc, utility
  ! Load database
  status = RM_LoadDatabase(id, "phreeqc.dat") 
  
  ! Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
  CALL register_basic_callback_fortran()
  
  ! Demonstration of error handling if ErrorHandlerMode is 0
  if (status .ne. 0) then
     l = RM_GetErrorStringLength(id)
#ifdef FORTRAN_2003
     allocate (character(len=l) :: errstr)
#endif
     write(*,*) "Start of error string: "
     status = RM_GetErrorString(id, errstr)
     write(*,"(A)") trim(errstr)
     write(*,*) "End of error string."
#ifdef FORTRAN_2003
     deallocate(errstr)
#endif
     status = RM_Destroy(id)
     stop
  endif
  ! Run file to define solutions and reactants for initial conditions, selected output
  ! There are three types of IPhreeqc instances in PhreeqcRM
  ! Argument 1 refers to the worker IPhreeqcs for doing reaction calculations for transport
  ! Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
  ! Argument 3 refers to the Utility instance available for processing
  status = RM_RunFile(id, 1, 1, 1, "advect.pqi")
  ! Clear contents of workers and utility
  string = "DELETE; -all"
  status = RM_RunString(id, 1, 0, 1, string)  ! workers, initial_phreeqc, utility
  ! Determine number of components to transport
  ncomps = RM_FindComponents(id)
  ! Print some of the reaction module information		
  write(string1, "(A,I10)") "Number of threads:                                ", RM_GetThreadCount(id)
  status = RM_OutputMessage(id, string1)
  write(string1, "(A,I10)") "Number of MPI processes:                          ", RM_GetMpiTasks(id)
  status = RM_OutputMessage(id, string1)
  write(string1, "(A,I10)") "MPI task number:                                  ", RM_GetMpiMyself(id)
  status = RM_OutputMessage(id, string1)
  status = RM_GetFilePrefix(id, string)
  write(string1, "(A,A)") "File prefix:                                      ", string
  status = RM_OutputMessage(id, trim(string1))
  write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", RM_GetGridCellCount(id)
  status = RM_OutputMessage(id, trim(string1))
  write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", RM_GetChemistryCellCount(id)
  status = RM_OutputMessage(id, trim(string1))
  write(string1, "(A,I10)") "Number of components for transport:               ", RM_GetComponentCount(id)
  status = RM_OutputMessage(id, trim(string1))
  ! Get component information
  allocate(components(ncomps))
  allocate(gfw(ncomps))
  status = RM_GetGfw(id, gfw)
  do i = 1, ncomps
     status = RM_GetComponent(id, i, components(i))
     write(string,"(A10, F15.4)") components(i), gfw(i)
     status = RM_OutputMessage(id, string)
  enddo
  status = RM_OutputMessage(id, " ")
  ! Set array of initial conditions
  allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
  ic1 = -1
  ic2 = -1
  f1 = 1.0
  do i = 1, nxyz
     ic1(i,1) = 1       ! Solution 1
     ic1(i,2) = -1      ! Equilibrium phases none
     ic1(i,3) = 1       ! Exchange 1
     ic1(i,4) = -1      ! Surface none
     ic1(i,5) = -1      ! Gas phase none
     ic1(i,6) = -1      ! Solid solutions none
     ic1(i,7) = -1      ! Kinetics none
  enddo
  status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)
  ! No mixing is defined, so the following is equivalent
  ! status = RM_InitialPhreeqc2Module(id, ic1)

  ! alternative for setting initial conditions
  ! cell number in second argument (-1 indicates last solution, 40 in this case)
  ! in advect.pqi and any reactants with the same number--
  ! Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
  ! will be written to cells 18 and 19 (0 based)
  allocate (module_cells(2))
  module_cells(1) = 18
  module_cells(2) = 19
  status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2)
  ! Initial equilibration of cells
  time = 0.0
  time_step = 0.0
  allocate(c(nxyz, ncomps))
  status = RM_SetTime(id, time)
  status = RM_SetTimeStep(id, time_step)
  status = RM_RunCells(id) 
  status = RM_GetConcentrations(id, c)
  
  ! --------------------------------------------------------------------------
  ! Set boundary condition
  ! --------------------------------------------------------------------------

  nbound = 1
  allocate(bc1(nbound), bc2(nbound), bc_f1(nbound))
  allocate(bc_conc(nbound, ncomps))  
  bc1 = 0           ! solution 0 from Initial IPhreeqc instance
  bc2 = -1          ! no bc2 solution for mixing
  bc_f1 = 1.0       ! mixing fraction for bc1 
  status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1)

  ! --------------------------------------------------------------------------
  ! Transient loop
  ! --------------------------------------------------------------------------

  nsteps = 10
  allocate(density(nxyz), pressure(nxyz), temperature(nxyz), volume(nxyz), sat_calc(nxyz))
  volume = 1.0
  density = 1.0
  pressure = 2.0
  temperature = 20.0
  status = RM_SetDensity(id, density)
  status = RM_SetTemperature(id, temperature)
  status = RM_SetPressure(id, pressure)
  time_step = 86400.0
  status = RM_SetTimeStep(id, time_step)
  do isteps = 1, nsteps
     ! Transport calculation here
     write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
          RM_GetTime(id) * RM_GetTimeConversion(id), " days"
     status = RM_LogMessage(id, string)
     status = RM_SetScreenOn(id, 1)
     status = RM_ScreenMessage(id, string)
     write(string, "(A32,F15.1,A)") "          Time step             ", &
          RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days"
     status = RM_LogMessage(id, string)
     status = RM_ScreenMessage(id, string)        
     call advect_f90(c, bc_conc, ncomps, nxyz)
     ! print at last time step
     if (isteps == nsteps) then
        status = RM_SetSelectedOutputOn(id, 1)         ! enable selected output
        status = RM_SetPrintChemistryOn(id, 1, 0, 0)   ! workers, initial_phreeqc, utility
     else
        status = RM_SetSelectedOutputOn(id, 0)         ! disable selected output
        status = RM_SetPrintChemistryOn(id, 0, 0, 0)   ! workers, initial_phreeqc, utility
     endif
     ! Transfer data to PhreeqcRM for reactions
     status = RM_SetPorosity(id, por)                ! If pore volume changes 
     status = RM_SetSaturation(id, sat)              ! If saturation changes
     status = RM_SetTemperature(id, temperature)     ! If temperature changes
     status = RM_SetPressure(id, pressure)           ! If pressure changes
     status = RM_SetConcentrations(id, c)            ! Transported concentrations
     status = RM_SetTimeStep(id, time_step)          ! Time step for kinetic reactions
     time = time + time_step
     status = RM_SetTime(id, time)                   ! Current time
     ! Run cells with transported conditions
     write(string, "(A32,F15.1,A)") "Beginning reaction calculation  ", &
          time * RM_GetTimeConversion(id), " days"
     status = RM_LogMessage(id, string) 
     status = RM_ScreenMessage(id, string) 
     status = RM_RunCells(id)  
     ! Transfer data from PhreeqcRM for transport
     status = RM_GetConcentrations(id, c)            ! Concentrations after reaction
     status = RM_GetDensity(id, density)             ! Density after reaction
     status = RM_GetSolutionVolume(id, volume)       ! Solution volume after reaction
     status = RM_GetSaturation(id, sat_calc)         ! Saturation after reaction
     ! Print results at last time step
     if (isteps == nsteps) then
		write(*,*) "Current distribution of cells for workers"
		write(*,*) "Worker      First cell        Last Cell"
		n = RM_GetThreadCount(id) * RM_GetMpiTasks(id)
        allocate(sc(n), ec(n))
        status = RM_GetStartCell(id, sc)
        status = RM_GetEndCell(id, ec)
        do i = 1, n
			write(*,*) i,"           ", sc(i),"                 ",ec(i)
        enddo
        
        ! Loop through possible multiple selected output definitions
        do isel = 1, RM_GetSelectedOutputCount(id)
           n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
           status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
           write(*,*) "Selected output sequence number: ", isel
           write(*,*) "Selected output user number:     ", n_user
           ! Get double array of selected output values
           col = RM_GetSelectedOutputColumnCount(id)
           allocate(selected_out(nxyz,col))
           status = RM_GetSelectedOutput(id, selected_out)
           ! Print results
           do i = 1, RM_GetSelectedOutputRowCount(id)/2
              write(*,*) "Cell number ", i
              write(*,*) "     Density: ", density(i)
              write(*,*) "     Volume:  ", volume(i)
              write(*,*) "     Components: "
              do j = 1, ncomps
                 write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ",trim(components(j)), ": ", c(i,j)
              enddo
              write(*,*) "     Selected output: "
              do j = 1, col
                 status = RM_GetSelectedOutputHeading(id, j, heading)    
                 write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ", trim(heading),": ", selected_out(i,j)
              enddo
           enddo
           deallocate(selected_out)
        enddo
     endif
  enddo

  ! --------------------------------------------------------------------------
  ! Additional features and finalize
  ! --------------------------------------------------------------------------

  ! Use utility instance of PhreeqcRM to calculate pH of a mixture
  allocate (c_well(1,ncomps))
  do i = 1, ncomps
     c_well(1,i) = 0.5 * c(1,i) + 0.5 * c(10,i)
  enddo
  allocate(tc(1), p_atm(1))
  tc(1) = 15.0
  p_atm(1) = 3.0
  iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm)
  string = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1"
  ! Alternatively, utility pointer is worker number nthreads + 1 
  iphreeqc_id1 = RM_GetIPhreeqcId(id, RM_GetThreadCount(id) + 1)
  status = SetOutputFileName(iphreeqc_id, "utility_f90.txt")
  status = SetOutputFileOn(iphreeqc_id, .true.)
  status = RunString(iphreeqc_id, string)
  if (status .ne. 0) status = RM_Abort(id, status, "IPhreeqc RunString failed") 
  status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5) 
  status = GetSelectedOutputValue(iphreeqc_id, 1, 1, vtype, pH, svalue)
  ! Dump results   
  status = RM_SetDumpFileName(id, "advection_f90.dmp")  
  dump_on = 1
  append = 0  
  status = RM_DumpModule(id, dump_on, append)    
  ! Clean up
  status = RM_CloseFiles(id)
  status = RM_MpiWorkerBreak(id)
  status = RM_Destroy(id)
  ! Deallocate
  deallocate(rv) 
  deallocate(por) 
  deallocate(sat) 
  deallocate(print_chemistry_mask) 
  deallocate(grid2chem) 
  deallocate(components) 
  deallocate(ic1) 
  deallocate(ic2) 
  deallocate(f1) 
  deallocate(bc1) 
  deallocate(bc2) 
  deallocate(bc_f1) 
  deallocate(bc_conc) 
  deallocate(c) 
  deallocate(density) 
  deallocate(temperature) 
  deallocate(c_well) 
  deallocate(pressure) 
  deallocate(tc) 
  deallocate(p_atm) 
  return 
end subroutine advection_f90

SUBROUTINE ADVECT_F90(c, bc_conc, ncomps, nxyz)
  implicit none
  double precision, dimension(:,:), allocatable, intent(inout) :: c 
  double precision, dimension(:,:), allocatable, intent(in)    :: bc_conc
  integer, intent(in)                                          :: ncomps, nxyz
  integer                                                      :: i, j
  ! Advect
  do i = nxyz/2, 2, -1
     do j = 1, ncomps
        c(i,j) = c(i-1,j)
     enddo
  enddo
  ! Cell 1 gets boundary condition
  do j = 1, ncomps
     c(1,j) = bc_conc(1,j)
  enddo
END SUBROUTINE ADVECT_F90

#ifdef USE_MPI
integer(kind=C_INT) function worker_tasks_f(method_number) BIND(C, NAME='worker_tasks_f')
    USE ISO_C_BINDING
    implicit none
    interface
        integer function do_something()
        end function do_something

        SUBROUTINE register_basic_callback_fortran()
            implicit none
        END SUBROUTINE register_basic_callback_fortran        
    end interface
    integer(kind=c_int), intent(in) :: method_number
    integer :: status
	if (method_number .eq. 1000) then
		status = do_something()
	else if (method_number .eq. 1001) then
		call register_basic_callback_fortran()
    endif
    worker_tasks_f = 0
end function worker_tasks_f
    
integer function do_something()
    implicit none  
	INCLUDE 'mpif.h'
	integer status
	integer i, method_number, mpi_myself, mpi_task, mpi_tasks, worker_number 
    method_number = 1000
	call MPI_Comm_size(MPI_COMM_WORLD, mpi_tasks, status)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
	if (mpi_myself .eq. 0) then     
		CALL MPI_Bcast(method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status)
		write(*,*) "I am root."        
		do i = 1, mpi_tasks-1
			CALL MPI_Recv(worker_number, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, status)
			write(*,*) "Recieved data from worker number ", worker_number, "."
		enddo
	else
		CALL MPI_Send(mpi_myself, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status)
    endif
	do_something = 0
end function do_something
#endif
SUBROUTINE register_basic_callback_fortran()
    USE IPhreeqc
    USE mydata
    USE PhreeqcRM
    USE ISO_C_BINDING
    implicit none
    INTERFACE
        REAL(kind=C_DOUBLE) FUNCTION my_basic_fortran_callback(x1, x2, str, l) BIND(C)
            USE ISO_C_BINDING
            IMPLICIT none
            REAL(kind=C_DOUBLE), INTENT(in)           :: x1, x2
            CHARACTER(kind=C_CHAR), INTENT(in)        :: str(*)
            INTEGER(kind=C_INT),    INTENT(in), value :: l
        END FUNCTION my_basic_fortran_callback   
    END INTERFACE    
#ifdef USE_MPI    
  INCLUDE 'mpif.h'
#endif
	integer status, method_number, mpi_tasks, mpi_myself
    integer i, j
    
    method_number = 1001

#ifdef USE_MPI
	CALL MPI_Comm_size(MPI_COMM_WORLD, mpi_tasks, status)
	CALL MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
	if (mpi_myself == 0) then
		CALL MPI_Bcast(method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status)
	endif
#endif

    do i = 1, RM_GetThreadCount(rm_id) + 2
		j = RM_GetIPhreeqcId(rm_id, i-1)
		j = SetBasicFortranCallback(j, my_basic_fortran_callback)
    enddo
END SUBROUTINE register_basic_callback_fortran
    
REAL(kind=C_DOUBLE) FUNCTION my_basic_fortran_callback(x1, x2, str, l) BIND(C, name='my_basic_fortran_callback')
    USE ISO_C_BINDING
    USE PhreeqcRM
    USE mydata
    IMPLICIT none
    REAL(kind=C_DOUBLE),    INTENT(in)        :: x1, x2
    CHARACTER(kind=C_CHAR), INTENT(in)        :: str(*)
    INTEGER(kind=C_INT),    INTENT(in), value :: l
    character(100) fstr

    INTEGER :: list(4), i
    INTEGER :: size=4, rm_cell_number
    
    do i = 1, l
        fstr(i:i) = str(i)
    enddo
	rm_cell_number = DINT(x1)
    my_basic_fortran_callback = -999.9
	if (rm_cell_number .ge. 0 .and. rm_cell_number < RM_GetChemistryCellCount(rm_id)) then
		if (RM_GetBackwardMapping(rm_id, rm_cell_number, list, size) .eq. 0) then
			if (fstr(1:l) .eq. "HYDRAULIC_K") then
				my_basic_fortran_callback = K_ptr(list(1)+1)
            endif
        endif
    endif
    END FUNCTION my_basic_fortran_callback
    
INTEGER FUNCTION example_selected_output(id)
    USE ISO_C_BINDING
    USE PhreeqcRM
    implicit none
    integer, intent(in) :: id
    character*100 :: line, line1, line2
    character*100, allocatable :: lines(:)
    integer iline, i, nlines, status
    nlines = 50
    nlines = nlines + RM_GetComponentCount(id) 
    nlines = nlines + RM_GetSpeciesCount(id)
    nlines = nlines + RM_GetExchangeSpeciesCount(id)
    nlines = nlines + RM_GetSurfaceSpeciesCount(id)
    nlines = nlines + RM_GetEquilibriumPhasesCount(id)
    nlines = nlines + RM_GetGasComponentsCount(id)
    nlines = nlines + RM_GetKineticReactionsCount(id)
    nlines = nlines + RM_GetSolidSolutionComponentsCount(id)
    nlines = nlines + RM_GetSICount(id)
    allocate (lines(nlines))
    
    iline = 1
    lines(iline) = "SELECTED_OUTPUT 2"
    iline = iline + 1
    ! totals
    lines(iline) = "  -totals"
    iline = iline + 1
    do i = 1, RM_GetComponentCount(id)
        status = RM_GetComponent(id, i, line)
        if (trim(line) .ne. "H" .and. &
            trim(line) .ne. "O" .and. &
            trim(line) .ne. "charge" .and. &
            trim(line) .ne. "Charge" .and. &
            trim(line) .ne. "H2O" ) then
            lines(iline) = "    " // trim(line)
            iline = iline + 1
         endif 
    enddo
    ! aqueous species
    lines(iline) = "  -molalities"
    iline = iline + 1
    do i = 1, RM_GetSpeciesCount(id)
        status = RM_GetSpeciesName(id, i, line)
        lines(iline) = trim(line)
        iline = iline + 1
    enddo 
    ! exchange species
    do i = 1, RM_GetExchangeSpeciesCount(id)
        status = RM_GetExchangeSpeciesName(id, i, line)
        status = RM_GetExchangeName(id, i, line1)
        lines(iline) = "    " // trim(line) // " # " // trim(line1)
        iline = iline + 1
    enddo  
    ! surface species species
    do i = 1, RM_GetSurfaceSpeciesCount(id)
        status = RM_GetSurfaceSpeciesName(id, i, line)
        status = RM_GetSurfaceType(id, i, line1)
        status = RM_GetSurfaceName(id, i, line2)
        lines(iline) = "    " // trim(line) // " # " // trim(line1) // "  " // trim(line2)
        iline = iline + 1
    enddo   
    ! equilibrium phases
    lines(iline) = "  -equilibrium_phases "
    iline = iline + 1
    do i = 1, RM_GetEquilibriumPhasesCount(id)
        status = RM_GetEquilibriumPhasesName(id, i, line)
        lines(iline) = "    " // trim(line)
        iline = iline + 1
    enddo
    ! gas components
    lines(iline) = "  -gas "
    iline = iline + 1
    do i = 1, RM_GetGasComponentsCount(id)
        status = RM_GetGasComponentsName(id, i, line)
        lines(iline) = "    " // trim(line)
        iline = iline + 1
    enddo
    ! kinetic reactions
    lines(iline) = "  -kinetics "
    iline = iline + 1
    do i = 1, RM_GetKineticReactionsCount(id)
        status = RM_GetKineticReactionsName(id, i, line)
        lines(iline) = "    " // trim(line)
        iline = iline + 1
    enddo
    ! solid solution components
    lines(iline) = "  -solid_solutions "
    iline = iline + 1
    do i = 1, RM_GetSolidSolutionComponentsCount(id)
        status = RM_GetSolidSolutionComponentsName(id, i, line)
        status = RM_GetSolidSolutionName(id, i, line1)
        lines(iline) = "    " // trim(line) // " # " // trim(line1)
        iline = iline + 1
    enddo
    ! saturation indices
    lines(iline) = "  -si "
    iline = iline + 1
    do i = 1, RM_GetSICount(id)
        status = RM_GetSIName(id, i, line)
        lines(iline) = "    " // trim(line)
        iline = iline + 1
    enddo    
#ifdef SKIP

	oss << "  -solid_solutions " << "\n";
	{
		// solid solutions 
		const std::vector<std::string> &ss_comps = phreeqc_rm.GetSolidSolutionComponents();
		const std::vector<std::string> &ss_names = phreeqc_rm.GetSolidSolutionNames();
		for (size_t i = 0; i < phreeqc_rm.GetSolidSolutionComponentsCount(); i++)
		{
			
			oss << "    ";
			oss.width(15);
			oss  << std::left << ss_comps[i];
			oss << " # " << ss_names[i] << "\n";
		}
	}
	oss << "  -saturation_indices " << "\n";
	{
		// molalities of aqueous species 
		const std::vector<std::string> &si = phreeqc_rm.GetSINames();
		for (size_t i = 0; i < phreeqc_rm.GetSICount(); i++)
		{
			oss << "    " << si[i] << "\n";
		}
	}

	// Uncommenting the following line would define SELECTED_OUTPUT 2 with all species, reactants, and SIs
	// int status = phreeqc_rm.RunString(true, false, false, oss.str().c_str());
	std::cerr << oss.str();

	return(0);
}
#endif    

    example_selected_output = 0
END FUNCTION example_selected_output