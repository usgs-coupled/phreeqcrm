subroutine Gas_f90()  BIND(C, NAME='Gas_f90')
  USE, intrinsic :: ISO_C_BINDING
  USE PhreeqcRM
  USE IPhreeqc
  USE mydata
  implicit none
#ifdef USE_MPI    
  INCLUDE 'mpif.h'
#endif
interface
    subroutine PrintCells(gas_comps, gas_moles, gas_p, gas_phi, str)
      character(len=:), dimension(:), allocatable, intent(in) :: gas_comps
      real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_moles
      real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_p
      real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_phi
      character(*), intent(in) :: str
    end subroutine PrintCells
end interface

  ! Based on PHREEQC Example 11

  integer :: mpi_myself
  integer :: i, j
  integer :: nxyz
  integer :: nthreads
  integer :: id
  integer :: status
  integer                                       :: ncomps
  integer                                       :: ngas
  real(kind=8),     dimension(:), allocatable   :: por
  real(kind=8),     dimension(:), allocatable   :: sat
  character(len=:), dimension(:), allocatable   :: gas_comps
  character(len=:), dimension(:), allocatable   :: components
  integer,          dimension(:,:), allocatable :: ic1, ic2
  real(kind=8), dimension(:,:), allocatable :: f1
  real(kind=8), dimension(:,:), allocatable :: gas_moles, gas_p, gas_phi
  real(kind=8), dimension(:), allocatable   :: gas_volume
  !character(LEN=1), dimension(:), allocatable   :: errstr
#ifdef FORTRAN_2003
  character(LEN=:), allocatable                 :: errstr
#else
  character(LEN=10000)                          :: errstr
#endif
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
  ! Open files
  status = RM_SetFilePrefix(id, "Gas_f90")
  status = RM_OpenFiles(id)

  ! Set concentration units
  status = RM_SetUnitsSolution(id, 2)      ! 1, mg/L; 2, mol/L; 3, kg/kgs
  status = RM_SetUnitsGasPhase(id, 0)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock

  ! Set initial porosity
  allocate(por(nxyz))
  por = 0.2
  status = RM_SetPorosity(id, por)
  ! Set initial saturation
  allocate(sat(nxyz))
  sat = 0.5
  status = RM_SetSaturationUser(id, sat)

  ! --------------------------------------------------------------------------
  ! Set initial conditions
  ! --------------------------------------------------------------------------

  ! Set printing of chemistry file to false
  status = RM_SetPrintChemistryOn(id, 0, 1, 0)  ! workers, initial_phreeqc, utility
  ! Load database
  status = RM_LoadDatabase(id, "phreeqc.dat")
  if (status < 0) then
      status = RM_OutputMessage(id, "Unable to open database.")
  endif

  status = RM_RunFile(id, 0, 1, 0, "gas.pqi")
  if (status < 0) then
      status = RM_OutputMessage(id, "Unable to run input file.")
  endif

  ! Determine number of components to transport
  ncomps = RM_FindComponents(id)
  ngas = RM_GetGasComponentsCount(id)

  ! Determine number of components and gas components
  !allocate(gas_comps(ngas))
  status = RM_GetGasComponentsNames(id, gas_comps)
  !do i = 1, ngas
  !    status = RM_GetGasComponentsName(id, i, gas_comps(i))
  !enddo
  status = RM_GetComponents(id, components)
  !allocate(components(ncomps))
  !do i = 1, ncomps
  !    status = RM_GetComponent(id, i, components(i))
  !enddo

  ! Set array of initial conditions
  allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
  ic1 = -1
  ic2 = -1
  f1 = 1.0
  do i = 1, nxyz
      ic1(i,1) = 1                 ! Solution 1
      ic1(i,5) = MOD(i-1,3) + 1    ! Gas phase none
  enddo
  status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)

  ! Get gases
  allocate(gas_moles(nxyz, ngas))
  allocate(gas_p(nxyz, ngas))
  allocate(gas_phi(nxyz, ngas))
  status = RM_GetGasCompMoles(id, gas_moles)
  status = RM_GetGasCompPressures(id, gas_p)
  status = RM_GetGasCompPhi(id, gas_phi)
  call PrintCells(gas_comps, gas_moles, gas_p, gas_phi, "Initial condition")
  
  ! multiply by 2
  do i = 1,nxyz
      do j = 1,ngas
          gas_moles(i,j) = gas_moles(i,j) * 2.0
      enddo
  enddo
  status = RM_SetGasCompMoles(id, gas_moles)
  status = RM_GetGasCompMoles(id, gas_moles)
  status = RM_GetGasCompPressures(id, gas_p)
  status = RM_GetGasCompPhi(id, gas_phi)
  call PrintCells(gas_comps, gas_moles, gas_p, gas_phi, "Initial condition times 2")

  ! Eliminate CH4(g) from cell 1, all gases from cell 2
  gas_moles(1,1) = -1.0
  gas_moles(2,1) = -1.0
  gas_moles(2,2) = -1.0
  gas_moles(2,3) = -1.0
  status = RM_SetGasCompMoles(id, gas_moles)
  status = RM_RunCells(id)
  status = RM_GetGasCompMoles(id, gas_moles)
  status = RM_GetGasCompPressures(id, gas_p)
  status = RM_GetGasCompPhi(id, gas_phi)
  call PrintCells(gas_comps, gas_moles, gas_p, gas_phi, "Remove some components")
  
    ! add CH4 in cell 0
    gas_moles(1,1) = 0.02
    ! Gas phase is added to cell 1; fixed pressure by default
    gas_moles(2,1) = 0.01
	gas_moles(2,2) = 0.02
	gas_moles(2,3) = 0.03
    status = RM_SetGasCompMoles(id, gas_moles)
	! Set volume for cell 1 and convert to fixed pressure gas phase
    allocate(gas_volume(nxyz))
    gas_volume = -1.0
	gas_volume(2) = 12.25
	status = RM_SetGasPhaseVolume(id, gas_volume)
	status = RM_GetGasPhaseVolume(id, gas_volume)
	status = RM_RunCells(id)
    status = RM_GetGasCompMoles(id, gas_moles)
    status = RM_GetGasCompPressures(id, gas_p)
    status = RM_GetGasCompPhi(id, gas_phi)
    call PrintCells(gas_comps, gas_moles, gas_p, gas_phi, "Add components back")
      
  ! Clean up
  status = RM_CloseFiles(id)
#ifdef USE_MPI
  status = RM_MpiWorkerBreak(id)
#endif
  status = RM_Destroy(id)
  ! Deallocate
  deallocate(por)
  deallocate(sat)
  deallocate(gas_comps)
  deallocate(ic1)
  deallocate(ic2)
  deallocate(f1)
  deallocate(gas_moles)

  return
end subroutine gas_f90
   
subroutine PrintCells(gas_comps, gas_moles, gas_p, gas_phi, str)
  implicit none
  character(len=:), dimension(:), allocatable, intent(in) :: gas_comps
  real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_moles
  real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_p
  real(kind=8), dimension(:,:), allocatable, intent(in) :: gas_phi
  integer :: n
  character(*), intent(in) :: str
  integer :: i,j
  write(*, *)
  write(*,"(A)") str
  do i = 1,3  ! cells
    write(*,"(A6,I2)") "Cell: ",i
    write(*,"(8x,3A10)") "Moles", "P", "Phi"
	do j = 1, size(gas_comps) ! gas components
	  write(*,"(2x,A6,3f10.6)") gas_comps(j), gas_moles(i,j), gas_p(i,j), gas_phi(i,j)
	enddo
  enddo    
  return
end subroutine PrintCells
    
