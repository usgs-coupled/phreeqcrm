    !module mydata
    !  double precision, dimension(:), pointer :: K_ptr
    !  integer                                 :: rm_id
    !end module mydata

    subroutine AdvectBMI_f90()  BIND(C, NAME='AdvectBMI_f90')
    USE, intrinsic :: ISO_C_BINDING
    USE PhreeqcRM
    USE IPhreeqc
    USE mydata
    USE BMI_PhreeqcRM
    implicit none
#ifdef USE_MPI    
    INCLUDE 'mpif.h'
#endif
    interface
        subroutine advectionbmi_f90(c, bc_conc, ncomps, nxyz)
            implicit none
            double precision, dimension(:,:), allocatable, intent(inout) :: c
            double precision, dimension(:,:), allocatable, intent(in) :: bc_conc
            integer, intent(in)                                       :: ncomps, nxyz
        end subroutine advectionbmi_f90
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
    double precision, pointer :: d1_ptr(:)
    character(100) :: yaml_file
    integer :: mpi_myself
    integer :: i, j
    integer :: nxyz
    integer :: nthreads
    integer :: id
    integer :: status
    integer :: bytes, nbytes
    double precision, dimension(:), allocatable, target :: hydraulic_K
    double precision, dimension(:), allocatable   :: por
    double precision, dimension(:), allocatable   :: sat
    integer                                       :: nchem
    character(100)                                :: string
    character(200)                                :: string1
    integer                                       :: ncomps, ncomps1
    character(LEN=:), dimension(:), allocatable   :: components
    double precision, dimension(:), allocatable   :: gfw
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
    integer                                       :: col, isel, n_user, rows
    character(LEN=:), dimension(:), allocatable   :: headings
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
    yaml_file = "AdvectBMI_cpp.yaml"
    ! RM_GetGridCellCountYAML must be called BEFORE
    ! the PhreeqcRM instance is created. The
    ! return value can be used to create the
    ! PhreeqcRM instance.
    !
    ! If the YAML file does not contain
    ! a node "SetGridCellCount:" (usually written
    ! using the YAMLPhreeqcRM class and the method
    ! YAMLSetGridCellCount), the return
    ! value is zero.
    nxyz = RM_GetGridCellCountYAML(yaml_file)
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
    ! Open files
    status = RM_SetFilePrefix(id, "AdvectBMI_f90")
    status = RM_OpenFiles(id)
    ! Initialize with YAML file
    status = RM_BMI_Initialize(id, yaml_file)

    ! Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
    CALL register_basic_callback_fortran()
#ifdef USE_MPI
    ! Optional callback for MPI
    status = do_something()   ! only root is calling do_something here
#endif

    nchem = RM_GetChemistryCellCount(id)
    status = RM_BMI_GetValue(id, "ComponentCount", ncomps)
    
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
    write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", nxyz
    status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", nchem
    status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of components for transport:               ", ncomps
    status = RM_OutputMessage(id, trim(string1))
    ! Get component information
    bytes = RM_BMI_GetVarItemsize(id, "Components")
    allocate(character(len=bytes) :: components(ncomps))    
    status = RM_BMI_GetValue(id, "Components", components)
    allocate(gfw(ncomps))
    status = RM_BMI_GetValue(id, "Gfw", gfw)
    do i = 1, ncomps
        write(string,"(A10, F15.4)") trim(components(i)), gfw(i)
        status = RM_OutputMessage(id, string)
    enddo
    status = RM_OutputMessage(id, " ")	
    ! Get initial temperatures
    allocate(temperature(nxyz))
    status = RM_BMI_GetValue(id, "Temperature", temperature)
    ! Get initial saturation
    allocate(sat(nxyz))
    status = RM_BMI_GetValue(id, "Saturation", sat)
    ! Get initial porosity
    allocate(por(nxyz))
    status = RM_BMI_GetValue(id, "Porosity", por)
    allocate(volume(nxyz))
    status = RM_BMI_GetValue(id, "SolutionVolume", volume)
    ! Get initial concentrations
    allocate(c(nxyz, ncomps))
    status = RM_BMI_GetValue(id, "Concentrations", c)
    ! Set density, pressure, and temperature
    allocate(density(nxyz), pressure(nxyz))
    density = 1.0
    status = RM_BMI_SetValue(id, "Density", density)
    pressure = 2.0
    status = RM_BMI_SetValue(id, "Pressure", pressure)  
    temperature = 20.0
    status = RM_BMI_SetValue(id, "Temperature", temperature)  
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
    time = 0.0
    status = RM_BMI_SetValue(id, "Time", time)
    time_step = 86400.0
    status = RM_BMI_SetValue(id, "TimeStep", time_step)   
    do isteps = 1, nsteps
        write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
            time/86400., " days"
        status = RM_LogMessage(id, string)
        status = RM_SetScreenOn(id, 1)
        status = RM_ScreenMessage(id, string)
        write(string, "(A32,F15.1,A)") "          Time step             ", &
            time_step / 86400., " days"
        status = RM_LogMessage(id, string)
        status = RM_ScreenMessage(id, string)    
        ! Transport calculation here, changes c
        call advectionbmi_f90(c, bc_conc, ncomps, nxyz)
    
        ! print at last time step
        if (isteps == nsteps) then
            status = RM_SetSelectedOutputOn(id, 1)         ! enable selected output
            status = RM_SetPrintChemistryOn(id, 1, 0, 0)   ! workers, initial_phreeqc, utility
        else
            status = RM_SetSelectedOutputOn(id, 0)         ! disable selected output
            status = RM_SetPrintChemistryOn(id, 0, 0, 0)   ! workers, initial_phreeqc, utility
        endif
        ! Transfer data to PhreeqcRM after transport      
        status = RM_BMI_SetValue(id, "Concentrations", c)  ! Transported concentrations
        ! Optionally, if values changed during transport
        status = RM_BMI_SetValue(id, "Porosity", por)              
        status = RM_BMI_SetValue(id, "Saturation", sat)            
        status = RM_BMI_SetValue(id, "Temperature", temperature) 
        status = RM_BMI_SetValue(id, "Pressure", pressure)          
        status = RM_BMI_SetValue(id, "TimeStep", time_step) 
        ! Set new time
        time = time + time_step              
        status = RM_BMI_SetValue(id, "Time", time)  ! Current time
        ! Run cells with transported conditions
        write(string, "(A32,F15.1,A)") "Beginning reaction calculation  ", &
            time / 86400., " days"
        status = RM_LogMessage(id, string)
        status = RM_ScreenMessage(id, string)
        ! Demonstration of state 
        status = RM_StateSave(id, 1)
        status = RM_StateApply(id, 1)
        status = RM_StateDelete(id, 1)
        ! Run chemistry
        status = RM_BMI_Update(id)
        ! Get new data calculated by PhreeqcRM for transport
        status = RM_BMI_GetValue(id, "Concentrations", c)
        status = RM_BMI_GetValue(id, "Density", density)
        status = RM_BMI_GetValue(id, "SolutionVolume", volume)   
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
            status = RM_BMI_GetValue(id, "SelectedOutputCount", n)
            do isel = 0, n-1
                i = isel
                status = RM_BMI_SetValue(id, "NthSelectedOutput", i)
                status = RM_BMI_GetValue(id, "CurrentSelectedOutputUserNumber", n_user)
                write(*,*) "Selected output sequence number: ", isel
                write(*,*) "Selected output user number:     ", n_user
                ! Get 2D array of selected output values
                status = RM_BMI_GetValue(id, "SelectedOutputColumnCount", col)
                status = RM_BMI_GetValue(id, "SelectedOutputRowCount", rows)
                allocate(selected_out(rows,col))
                ! Get headings
                bytes = RM_BMI_GetVarItemsize(id, "SelectedOutputHeadings")
                allocate(character(len=bytes) :: headings(col))    
                status = RM_BMI_GetValue(id, "SelectedOutputHeadings", headings)
                ! Get selected output
                status = RM_BMI_GetValue(id, "SelectedOutput", selected_out)
                ! Print results
                do i = 1, rows/2
                    write(*,*) "Cell number ", i
                    write(*,*) "     Density: ", density(i)
                    write(*,*) "     Volume:  ", volume(i)
                    write(*,*) "     Components: "
                    do j = 1, ncomps
                        write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ",trim(components(j)), ": ", c(i,j)
                    enddo
                    write(*,*) "     Selected output: "
                    do j = 1, col
                        write(*,'(10x,i3,A2,A20,A2,1pe15.4)') j, " ", trim(headings(j)),": ", selected_out(i,j)
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
    status = SetOutputFileName(iphreeqc_id, "Advect_f90_utility.txt")
    status = SetOutputFileOn(iphreeqc_id, .true.)
    status = RunString(iphreeqc_id, string)
    if (status .ne. 0) status = RM_Abort(id, status, "IPhreeqc RunString failed")
    status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5)
    status = GetSelectedOutputValue(iphreeqc_id, 1, 1, vtype, pH, svalue)
    ! Dump results
    status = RM_SetDumpFileName(id, "AdvectBMI_f90.dmp")
    dump_on = 1
    append = 0
    status = RM_DumpModule(id, dump_on, append)
    ! Clean up
    status = RM_CloseFiles(id)
    status = RM_MpiWorkerBreak(id)
    status = RM_BMI_Finalize(id)
    ! Deallocate
    deallocate(por)
    deallocate(sat)
    deallocate(components)
    deallocate(headings)
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
    end subroutine AdvectBMI_f90

    SUBROUTINE advectionbmi_f90(c, bc_conc, ncomps, nxyz)
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
    END SUBROUTINE advectionbmi_f90

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

