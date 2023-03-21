#ifdef USE_YAML
    !module mydata
    !  double precision, dimension(:), pointer :: K_ptr
    !  integer                                 :: rm_id
    !end module mydata

    subroutine AdvectBMI_f90()  BIND(C, NAME='AdvectBMI_f90')
    USE, intrinsic :: ISO_C_BINDING
    USE PhreeqcRM
    USE IPhreeqc
    USE mydata
    !USE BMI_PhreeqcRM
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
        subroutine BMI_testing(id)
            implicit none
            integer, intent(in) :: id
        end subroutine BMI_testing
    end interface

    ! Based on PHREEQC Example 11
    double precision, pointer :: d1_ptr(:)
    character(100) :: yaml_file
    integer :: mpi_myself
    integer :: i, j
    logical :: tf
    integer :: nxyz
    integer :: nthreads
    integer :: id
    integer :: status
    integer :: bytes, nbytes
    double precision, dimension(:), allocatable, target :: hydraulic_K
    double precision, dimension(:), allocatable   :: por
    double precision, dimension(:), allocatable   :: sat
    integer                                       :: nchem
    character(len=:), allocatable                 :: prefix
    character(len=:), allocatable                 :: alloc_string
    character(100)                                :: string
    character(200)                                :: string1
    integer                                       :: ncomps, ncomps1
    character(len=:), dimension(:), allocatable          :: components
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
    character(len=:), dimension(:), allocatable   :: headings
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
    yaml_file = "AdvectBMI_f90.yaml"
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
    nxyz = GetGridCellCountYAML(yaml_file)

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
    !status = RM_BMI_SetValue(id, "FilePrefix", "AdvectBMI_f90")
    !status = RM_OpenFiles(id)
    ! Initialize with YAML file
    status = RM_BMI_Initialize(id, yaml_file)

    ! Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
    CALL register_basic_callback_fortran()
#ifdef USE_MPI
    ! Optional callback for MPI
    status = do_something()   ! only root is calling do_something here
#endif
    status = RM_BMI_GetValue(id, "ComponentCount", ncomps)
    
    ! Print some of the reaction module information
    write(string1, "(A,I10)") "Number of threads:                                ", RM_GetThreadCount(id)
    status = RM_OutputMessage(id, string1)
    write(string1, "(A,I10)") "Number of MPI processes:                          ", RM_GetMpiTasks(id)
    status = RM_OutputMessage(id, string1)
    write(string1, "(A,I10)") "MPI task number:                                  ", RM_GetMpiMyself(id)
    status = RM_OutputMessage(id, string1)
    status = RM_BMI_GetValue(id, "FilePrefix", prefix)
    write(string1, "(A,A)") "File prefix:                                        ", prefix
    status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", nxyz
    status = RM_OutputMessage(id, trim(string1))
    !write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", nchem
    !status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of components for transport:               ", ncomps
    status = RM_OutputMessage(id, trim(string1))
    ! Get component information
    status = RM_BMI_GetValue(id, "Components", components)
    status = RM_BMI_GetValue(id, "Gfw", gfw)
    do i = 1, ncomps
        write(string,"(A10, F15.4)") trim(components(i)), gfw(i)
        status = RM_OutputMessage(id, string)
    enddo
    status = RM_OutputMessage(id, " ")	
    ! Get initial temperatures
    status = RM_BMI_GetValue(id, "Temperature", temperature)
    ! Get initial temperature
    status = RM_BMI_GetValue(id, "Saturation", sat)
    ! Get initial porosity
    status = RM_BMI_GetValue(id, "Porosity", por)
    ! Get initial temperature
    status = RM_BMI_GetValue(id, "SolutionVolume", volume)
    ! Get initial concentrations
    status = RM_BMI_GetValue(id, "Concentrations", c)
    ! Set density, pressure, and temperature (previously allocated)
    allocate(density(nxyz))
    allocate(pressure(nxyz))
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
            status = RM_BMI_SetValue(id, "SelectedOutputOn", .true.)    ! enable selected output
            status = RM_SetPrintChemistryOn(id, 1, 0, 0)                ! workers, initial_phreeqc, utility
        else        
            status = RM_BMI_SetValue(id, "SelectedOutputOn", .false.)   ! disable selected output
            status = RM_SetPrintChemistryOn(id, 0, 0, 0)                ! workers, initial_phreeqc, utility
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
            do isel = 1, n  ! one based
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
    call BMI_testing(id)
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

subroutine BMI_testing(id)
USE, intrinsic :: ISO_C_BINDING
    USE PhreeqcRM
    USE IPhreeqc
    !USE BMI_PhreeqcRM
    implicit none
    interface
        integer function assert(tf)
        logical, intent(in) :: tf
        end function assert
    end interface
    integer, intent(in) :: id
    character(100) :: string
    integer :: status, dim, isel
    integer :: i, j, n, bytes, nbytes
    integer                                              :: nxyz, rm_nxyz
    integer                                              :: ncomps, rm_ncomps
    character(len=:), dimension(:), allocatable          :: components, rm_components
    character(len=:), allocatable                        :: component
    character(len=:), dimension(:), allocatable          :: inputvars
    character(len=:), dimension(:), allocatable          :: outputvars 
    character(len=:), allocatable                        :: prefix, rm_prefix
    double precision, dimension(:), allocatable          :: gfw, rm_gfw
    double precision, dimension(:,:), allocatable        :: c, c_rm
    integer                                              :: so_count, rm_so_count
    integer                                              :: nuser, rm_nuser
    integer                                              :: col_count, rm_col_count
    integer                                              :: row_count, rm_row_count
    integer                                              :: itemsize
    double precision, dimension(:,:), allocatable        :: so, rm_so
    character(LEN=:), dimension(:), allocatable          :: headings, rm_headings
    character(LEN=:), allocatable                        :: heading
    double precision                                     :: time, rm_time
    double precision                                     :: time_step, rm_time_step
    double precision, dimension(:), allocatable          :: density, rm_density
    double precision, dimension(:), allocatable          :: porosity, rm_porosity
    double precision, dimension(:), allocatable          :: pressure, rm_pressure
    double precision, dimension(:), allocatable          :: saturation, rm_saturation
    double precision, dimension(:), allocatable          :: temperature, rm_temperature
    double precision, dimension(:), allocatable          :: volume, rm_volume
    logical :: tf
    character(LEN=:), allocatable :: alloc_string

    status = RM_BMI_GetComponentName(id, string)
    write(*,*) trim(string)
    write(*,*) RM_BMI_GetCurrentTime(id)
    write(*,*) RM_BMI_GetEndTime(id)
    status = RM_BMI_GetTimeUnits(id, string)
    write(*,*) string
    write(*,*) RM_BMI_GetTimeStep(id)
    ! Time
    rm_time = 3600
    status = RM_BMI_SetValue(id, "Time", rm_time)
    status = RM_BMI_GetValue(id, "Time", time)
    status = assert(rm_time .eq. time)
    rm_time = RM_GetTime(id)
    status = assert(time .eq. rm_time)
    time = RM_BMI_GetCurrentTime(id)
    status = assert(time .eq. rm_time)
    ! TimeStep
    rm_time_step = 60
    status = RM_BMI_SetValue(id, "TimeStep", rm_time_step)
    status = RM_BMI_GetValue(id, "TimeStep", time_step)
    status = assert(time_step .eq. rm_time_step)
    rm_time_step = RM_GetTimeStep(id)
    status = assert(time_step .eq. rm_time_step)
    time_step = RM_BMI_GetTimeStep(id)
    status = assert(time_step .eq. rm_time_step)
    ! InputVarNames
    status = RM_BMI_GetInputVarNames(id, inputvars)
    status = RM_BMI_GetValue(id, "InputVarNames", inputvars)
    write(*,*) "Input variables (setters)"
    do i = 1, size(inputvars)
        write(*,"(1x, I4, A40)") i, trim(inputvars(i))
        status = RM_BMI_GetVarUnits(id, inputvars(i), string)
        write(*,"(5x, A15)") trim(string)
        status = RM_BMI_GetVarType(id, inputvars(i), string)
        write(*,"(5x, A15)") trim(string)
        write(*, "(5x, I15)") RM_BMI_GetVarItemsize(id, inputvars(i))
        write(*, "(5x, I15)") RM_BMI_GetVarNbytes(id, inputvars(i))
    enddo
    ! OutputVarNames
    status = RM_BMI_GetOutputVarNames(id, outputvars)
    status = RM_BMI_GetValue(id, "OutputVarNames", outputvars)
    write(*,*) "Output variables (getters)"
    do i = 1, size(outputvars)
        write(*,"(1x, I4, A40)") i, trim(outputvars(i))
        status = RM_BMI_GetVarUnits(id, outputvars(i), string)
        write(*,"(5x, A15)") trim(string)
        status = RM_BMI_GetVarType(id, outputvars(i), string)
        write(*,"(5x, A15)") trim(string)
        write(*, "(5x, I15)") RM_BMI_GetVarItemsize(id, outputvars(i))
        write(*, "(5x, I15)") RM_BMI_GetVarNbytes(id, outputvars(i))
    enddo

    ! ComponentCount
    status = RM_BMI_GetValue(id, "ComponentCount", ncomps)
    rm_ncomps = RM_GetComponentCount(id)
    status = assert(ncomps .eq. rm_ncomps)
    ! Components
    status = RM_BMI_GetValue(id, "Components", components)
    itemsize = RM_BMI_GetVarItemsize(id, "Components")
    nbytes = RM_BMI_GetVarNBytes(id, "Components")
    dim = nbytes / itemsize
    status = assert(dim .eq. size(components))
    allocate(character(len=itemsize) :: component)
    do i = 1, dim
        j = i
        status = RM_GetComponent(id, j, component)
        status = assert(component .eq. components(i))
    enddo
    status = RM_GetComponents(id, rm_components)
    do i = 1, dim
        status = assert(components(i) .eq. rm_components(i))
    enddo

    status = RM_BMI_GetValue(id, "Gfw", gfw)
    do i = 1, size(components)
        write(*,"(A10, F15.4)") trim(components(i)), gfw(i)
    enddo
    write(*,*)	   
    
    ! Concentrations
    status = RM_BMI_GetValue(id, "GridCellCount", nxyz)
    status = assert(nxyz .eq. RM_GetGridCellCount(id))
    status = RM_BMI_GetValue(id, "Concentrations", c)
    allocate(c_rm(nxyz, ncomps))
    status = RM_GetConcentrations(id, c_rm)
    do j = 1, ncomps
        do i = 1, nxyz
            if (c(i,j) .ne. c_rm(i,j)) then
                status = assert(.false.)
                exit
            endif
        enddo
    enddo
   
	! GetValue("Density")
    ! RM_GetDensity and RM_BMI_GetValue("Density) always return 
    ! the calculated solution density
    itemsize = RM_BMI_GetVarItemsize(id, "Density")
    nbytes = RM_BMI_GetVarNBytes(id, "Density")
    dim = nbytes / itemsize
    allocate(rm_density(dim))
    status = RM_GetDensity(id, rm_density)
    status = RM_BMI_GetValue(id, "Density", density)
    do i = 1, nxyz
        if (density(i) .ne. rm_density(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
    ! FilePrefix
    string = "NewPrefix"
    status = RM_BMI_SetValue(id, "FilePrefix", string)
    itemsize = RM_BMI_GetVarItemsize(id, "FilePrefix")
    nbytes = RM_BMI_GetVarNBytes(id, "FilePrefix")
    allocate(character(len=itemsize) :: rm_prefix)
    status = assert(itemsize .eq. nbytes)
    status = RM_GetFilePrefix(id, rm_prefix)
    status = assert(string .eq. rm_prefix)
    allocate(character(len=itemsize) :: prefix)
    status = RM_BMI_GetValue(id, "FilePrefix", prefix)
    status = assert(prefix .eq. rm_prefix)
         
	! GetValue("Gfw")
    status = RM_BMI_GetValue(id, "Gfw", gfw)
    itemsize = RM_BMI_GetVarItemsize(id, "Gfw")
    nbytes = RM_BMI_GetVarNBytes(id, "Gfw")
    dim = nbytes / itemsize
    allocate(rm_gfw(dim))
	status = RM_GetGfw(id, rm_gfw);
    do i = 1, ncomps
        if (gfw(i) .ne. rm_gfw(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
    ! GridCellCount
    status = RM_BMI_GetValue(id, "GridCellCount", nxyz)
    rm_nxyz = RM_GetGridCellCount(id)
    status = assert(nxyz .eq. rm_nxyz)
    
	! GetValue("Porosity")
    itemsize = RM_BMI_GetVarItemsize(id, "Porosity")
    nbytes = RM_BMI_GetVarNBytes(id, "Porosity")
    dim = nbytes / itemsize
    allocate(rm_porosity(dim))
    rm_porosity = 0.25
    status = RM_BMI_SetValue(id, "Porosity", rm_porosity)
    status = RM_BMI_GetValue(id, "Porosity", porosity)
    do i = 1, nxyz
        if (porosity(i) .ne. rm_porosity(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
	status = RM_GetPorosity(id, rm_porosity);
        do i = 1, nxyz
        if (porosity(i) .ne. rm_porosity(i)) then
            status = assert(.false.)
            exit
        endif
    enddo

	! GetValue("Pressure")
    itemsize = RM_BMI_GetVarItemsize(id, "Pressure")
    nbytes = RM_BMI_GetVarNBytes(id, "Pressure")
    dim = nbytes / itemsize
    allocate(rm_pressure(dim))
    rm_pressure = 10.
    status = RM_BMI_SetValue(id, "Pressure", rm_pressure)
    status = RM_BMI_GetValue(id, "Pressure", pressure)
    do i = 1, nxyz
        if (pressure(i) .ne. rm_pressure(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    status = RM_GetPressure(id, rm_pressure);
    do i = 1, nxyz
        if (pressure(i) .ne. rm_pressure(i)) then
            status = assert(.false.)
            exit
        endif
    enddo

	! GetValue("Saturation")
    ! Always returns solution_volume/(rv * porosity) for each cell
    itemsize = RM_BMI_GetVarItemsize(id, "Saturation")
    nbytes = RM_BMI_GetVarNBytes(id, "Saturation")
    dim = nbytes / itemsize
    allocate(rm_saturation(dim))
    status = RM_BMI_GetValue(id, "Saturation", saturation)
	status = RM_GetSaturation(id, rm_saturation);
    do i = 1, nxyz
        if (saturation(i) .ne. rm_saturation(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
    ! GetValue("SolutionVolume")
    itemsize = RM_BMI_GetVarItemsize(id, "SolutionVolume")
    nbytes = RM_BMI_GetVarNBytes(id, "SolutionVolume")
    dim = nbytes / itemsize
    allocate(rm_volume(dim))
    status = RM_BMI_GetValue(id, "SolutionVolume", volume)
	status = RM_GetSolutionVolume(id, rm_volume);
    do i = 1, nxyz
        if (volume(i) .ne. rm_volume(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
	! GetValue("Temperature")
    itemsize = RM_BMI_GetVarItemsize(id, "Temperature")
    nbytes = RM_BMI_GetVarNBytes(id, "Temperature")
    dim = nbytes / itemsize
    allocate(rm_temperature(dim))
    rm_temperature = 11.
    status = RM_BMI_SetValue(id, "Temperature", rm_temperature)
    status = RM_BMI_GetValue(id, "Temperature", temperature)
    do i = 1, nxyz
        if (temperature(i) .ne. rm_temperature(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    status = RM_GetTemperature(id, rm_temperature);
    do i = 1, nxyz
        if (temperature(i) .ne. rm_temperature(i)) then
            status = assert(.false.)
            exit
        endif
    enddo

	! GetValue("SelectedOutput")
    status = RM_BMI_GetValue(id, "SelectedOutputCount", so_count);
	rm_so_count = RM_GetSelectedOutputCount(id);
    status = assert(so_count .eq. rm_so_count)
    do isel = 1, so_count ! one based
        i = isel
		status = RM_BMI_SetValue(id, "NthSelectedOutput", i)
        status = RM_BMI_GetValue(id, "CurrentSelectedOutputUserNumber", nuser)
        rm_nuser = RM_GetCurrentSelectedOutputUserNumber(id)
        status = assert(nuser .eq. rm_nuser)
        i = isel        
		rm_nuser = RM_GetNthSelectedOutputUserNumber(id, i)
        status = assert(nuser .eq. rm_nuser)

        status = RM_BMI_GetValue(id, "SelectedOutputColumnCount", col_count)
        rm_col_count =RM_GetSelectedOutputColumnCount(id)
        status = assert(col_count .eq. rm_col_count)

        status = RM_BMI_GetValue(id, "SelectedOutputRowCount", row_count)
        rm_row_count = RM_GetSelectedOutputRowCount(id)
        status = assert(row_count .eq. rm_row_count)

        nbytes = RM_BMI_GetVarNbytes(id, "SelectedOutput");
        itemsize = RM_BMI_GetVarItemsize(id, "SelectedOutput");
        dim = nbytes / itemsize;
        status = assert(dim .eq. rm_row_count*rm_col_count)
        status = RM_BMI_GetValue(id, "SelectedOutput", so)
        allocate(rm_so(rm_row_count, rm_col_count))
        status = RM_GetSelectedOutput(id, rm_so);
        do j = 1, col_count
            do i = 1, row_count
                if (so(i,j) .ne. rm_so(i,j)) then
                    status = assert(.false.)
                    exit
                endif
            enddo
        enddo
        ! check headings
        nbytes = RM_BMI_GetVarNbytes(id, "SelectedOutputHeadings")
        itemsize = RM_BMI_GetVarItemsize(id, "SelectedOutputHeadings")
        dim = nbytes / itemsize
        status = assert(dim .eq. RM_GetSelectedOutputColumnCount(id))
        status = RM_BMI_GetValue(id, "SelectedOutputHeadings", headings)
        !allocate(character(len=itemsize) :: rm_headings(dim))
        status = RM_GetSelectedOutputHeadings(id, rm_headings)
        do j = 1, col_count
            if (headings(j) .ne. rm_headings(j)) then
                status = assert(.false.)
            endif
        enddo
        itemsize = RM_BMI_GetVarItemsize(id, "SelectedOutputHeadings")
        allocate(character(len=itemsize) :: heading)
        do j = 1, col_count
            status = RM_GetSelectedOutputHeading(id, j, heading)
            if (heading .ne. rm_headings(j)) then
                status = assert(.false.)
            endif
        enddo
    enddo

END subroutine BMI_testing
integer function assert(tf)
logical, intent(in) :: tf
if(tf) then
    assert = 0
    return
endif
stop "Assert failed"
end function assert

#endif 
