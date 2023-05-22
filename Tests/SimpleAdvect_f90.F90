    
subroutine SimpleAdvect_f90()  BIND(C, NAME='SimpleAdvect_f90')
  USE, intrinsic :: ISO_C_BINDING
  USE PhreeqcRM
  USE IPhreeqc
  USE mydata
  implicit none
#ifdef USE_MPI    
  INCLUDE 'mpif.h'
#endif
    interface
        subroutine simpleadvection_f90(c, bc_conc, ncomps, nxyz)
            implicit none
            real(kind=8), dimension(:,:), allocatable, intent(inout) :: c 
            real(kind=8), dimension(:,:), allocatable, intent(in) :: bc_conc
            integer, intent(in)                                       :: ncomps, nxyz
        end subroutine simpleadvection_f90
    end interface

  ! Based on PHREEQC Example 11
  integer :: mpi_myself
  integer :: i, j
  integer :: nxyz
  integer :: nthreads
  integer :: id
  integer :: status
  real(kind=8), dimension(:), allocatable   :: por
  integer,          dimension(:), allocatable   :: print_chemistry_mask
  integer                                       :: nchem
  character(100)                                :: string
  integer                                       :: ncomps
  character(len=:), dimension(:), allocatable   :: components
  integer,          dimension(:,:), allocatable :: ic1
  integer                                       :: nbound
  integer,          dimension(:), allocatable   :: bc1
  real(kind=8), dimension(:,:), allocatable :: bc_conc
  real(kind=8), dimension(:,:), allocatable :: c
  real(kind=8)                              :: time, time_step
  real(kind=8), dimension(:), allocatable   :: temperature
  real(kind=8), dimension(:), allocatable   :: pressure
  integer                                       :: isteps, nsteps
  ! --------------------------------------------------------------------------
  ! Create PhreeqcRM
  ! --------------------------------------------------------------------------
  nxyz = 20
#ifdef USE_MPI
  ! MPI
  id = RM_Create(nxyz, MPI_COMM_WORLD)
  rm_id = id
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
  if (status .ne. MPI_SUCCESS) then
     stop "Failed to get mpi_myself"
  endif
  if (mpi_myself > 0) then
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
  status = RM_SetComponentH2O(id, 0)
  status = RM_UseSolutionDensityVolume(id, 0)
  ! Open files
  status = RM_SetFilePrefix(id, "SimpleAdvect_f90")
  status = RM_OpenFiles(id)
  ! Set concentration units
  status = RM_SetUnitsSolution(id, 2)      ! 1, mg/L; 2, mol/L; 3, kg/kgs
  status = RM_SetUnitsExchange(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  ! Set conversion from seconds to user units (days)
  status = RM_SetTimeConversion(id, dble(1.0 / 86400.0))
  ! Set initial porosity
  allocate(por(nxyz))
  por = 0.2
  status = RM_SetPorosity(id, por)
  ! Set cells to print chemistry when print chemistry is turned on
  allocate(print_chemistry_mask(nxyz))
  print_chemistry_mask = 1
  status = RM_SetPrintChemistryMask(id, print_chemistry_mask)
  nchem = RM_GetChemistryCellCount(id)
  ! --------------------------------------------------------------------------
  ! Set initial conditions
  ! --------------------------------------------------------------------------
  ! Set printing of chemistry file to false
  status = RM_SetPrintChemistryOn(id, 0, 1, 0)  ! workers, initial_phreeqc, utility
  ! Load database
  status = RM_LoadDatabase(id, "phreeqc.dat") 
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
  ! Get component information
  !allocate(components(ncomps))
  status = RM_GetComponents(id, components)
  do i = 1, ncomps
     !status = RM_GetComponent(id, i, components(i))
     write(string,"(A10, F15.4)") components(i)
     status = RM_OutputMessage(id, string)
  enddo
  ! Set array of initial conditions
  allocate(ic1(nxyz,7))
  ic1 = -1
  do i = 1, nxyz
     ic1(i,1) = 1       ! Solution 1
     ic1(i,2) = -1      ! Equilibrium phases none
     ic1(i,3) = 1       ! Exchange 1
     ic1(i,4) = -1      ! Surface none
     ic1(i,5) = -1      ! Gas phase none
     ic1(i,6) = -1      ! Solid solutions none
     ic1(i,7) = -1      ! Kinetics none
  enddo
  status = RM_InitialPhreeqc2Module(id, ic1)
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
  allocate(bc1(nbound))
  allocate(bc_conc(nbound, ncomps))  
  bc1 = 0           ! solution 0 from Initial IPhreeqc instance
  status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1)
  ! --------------------------------------------------------------------------
  ! Transient loop
  ! --------------------------------------------------------------------------
  nsteps = 10
  allocate(pressure(nxyz), temperature(nxyz))
  pressure = 2.0
  temperature = 20.0
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
     ! advect one step
     call simpleadvection_f90(c, bc_conc, ncomps, nxyz)
     ! print at last time step
     if (isteps == nsteps) then
        status = RM_SetSelectedOutputOn(id, 1)         ! enable selected output
        status = RM_SetPrintChemistryOn(id, 1, 0, 0)   ! workers, initial_phreeqc, utility
     else
        status = RM_SetSelectedOutputOn(id, 0)         ! disable selected output
        status = RM_SetPrintChemistryOn(id, 0, 0, 0)   ! workers, initial_phreeqc, utility
     endif
     ! Transfer data to PhreeqcRM for reactions
     status = RM_SetConcentrations(id, c)            ! Transported concentrations
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
  enddo
  status = RM_CloseFiles(id)
  status = RM_MpiWorkerBreak(id)
  status = RM_Destroy(id)
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
end subroutine SimpleAdvect_F90

SUBROUTINE simpleadvection_f90(c, bc_conc, ncomps, nxyz)
  implicit none
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: c 
  real(kind=8), dimension(:,:), allocatable, intent(in)    :: bc_conc
  integer, intent(in)                                          :: ncomps, nxyz
  integer                                                      :: i, j
  ! Advect
  do i = nxyz, 2, -1
     do j = 1, ncomps
        c(i,j) = c(i-1,j)
     enddo
  enddo
  ! Cell 1 gets boundary condition
  do j = 1, ncomps
     c(1,j) = bc_conc(1,j)
  enddo
END SUBROUTINE simpleadvection_f90

