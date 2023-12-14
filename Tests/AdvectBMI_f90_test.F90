#ifdef USE_YAML
    subroutine AdvectBMI_f90_test()  BIND(C, NAME='AdvectBMI_f90_test')
    USE, intrinsic :: ISO_C_BINDING
    USE BMIPhreeqcRM
    implicit none
#ifdef USE_MPI    
    INCLUDE 'mpif.h'
#endif
    interface
        subroutine advectionbmi_f90_test(c, bc_conc, ncomps, nxyz)
            implicit none
            real(kind=8), dimension(:,:), allocatable, intent(inout) :: c
            real(kind=8), dimension(:,:), allocatable, intent(in)    :: bc_conc
            integer, intent(in)                                      :: ncomps, nxyz
        end subroutine advectionbmi_f90_test
        integer function do_something()
        end function do_something
        integer(kind=C_INT) function bmi_worker_tasks_f(method_number) BIND(C, NAME='worker_tasks_f')
            USE ISO_C_BINDING
            implicit none
            integer(kind=c_int), intent(in) :: method_number
        end function bmi_worker_tasks_f
        SUBROUTINE register_basic_callback_fortran()
            implicit none
        END SUBROUTINE register_basic_callback_fortran
        !subroutine BMI_testing(self)
        !    implicit none
        !    class(bmi), intent(inout) :: self
        !end subroutine BMI_testing
    end interface
    ! Based on PHREEQC Example 11
    real(kind=8), pointer :: d1_ptr(:)
    character(100) :: yaml_file
    integer :: mpi_myself
    integer :: i, j
    logical :: tf
    integer :: nthreads
    integer :: status
    integer :: bytes, nbytes
    real(kind=8), dimension(:), allocatable, target :: hydraulic_K
    real(kind=8), dimension(:), allocatable         :: por
    real(kind=8), dimension(:), allocatable         :: sat
    integer                                         :: nchem
    character(len=:), allocatable                   :: prefix
    character(len=2) :: shortprefix
    character(len=:), allocatable                   :: alloc_string
    character(100)                                  :: string
    character(200)                                  :: string1
    character(len=:), dimension(:), allocatable     :: components
    real(kind=8), dimension(:), allocatable         :: gfw
    integer                                         :: nbound
    integer,          dimension(:), allocatable     :: bc1, bc2
    real(kind=8), dimension(:), allocatable         :: bc_f1
    integer,          dimension(:), allocatable     :: module_cells
    real(kind=8), dimension(:,:), allocatable       :: bc_conc
    real(kind=8), dimension(:,:), allocatable       :: c
    real(kind=8), dimension(:), allocatable         :: c1
    real(kind=8)                                    :: time, time_step
    real(kind=8), dimension(:), allocatable         :: density
    real(kind=8), dimension(:), allocatable         :: sat_calc
    real(kind=8), dimension(:), allocatable         :: volume
    real(kind=8), dimension(:), allocatable         :: temperature
    real(kind=8), dimension(:), allocatable         :: pressure
    integer                                         :: isteps, nsteps
    real(kind=8), dimension(:,:), allocatable       :: selected_out
    integer                                         :: col, isel, n_user, rows
    character(len=:), dimension(:), allocatable     :: headings
    real(kind=8), dimension(:,:), allocatable       :: c_well
    real(kind=8), dimension(:), allocatable         :: tc, p_atm
    integer                                         :: vtype
    real(kind=8)                                    :: pH
    character(100)                                  :: svalue
    integer                                         :: iphreeqc_id, iphreeqc_id1
    integer                                         :: dump_on, append
    integer                                         :: dim
    character(len=:), dimension(:), allocatable     :: outputvars 
	real(kind=8), dimension(:), allocatable         :: CaX2, KX, NaX, pH_vector, SAR
    common /i_ptrs/ ComponentCount_ptr, GridCellCount_ptr, SelectedOutputOn_ptr, id
    common /r_ptrs/ Concentrations_ptr, Density_calculated_ptr, Gfw_ptr, &
	    Saturation_ptr, SolutionVolume_ptr, Time_ptr, TimeStep_ptr, &
        Porosity_ptr, Pressure_ptr, Temperature_ptr
    integer :: id
	integer, pointer :: ComponentCount_ptr
	integer, pointer :: GridCellCount_ptr
    integer :: nxyz, ncomps
	logical(kind=1), pointer :: SelectedOutputOn_ptr
	real(kind=8), pointer :: Concentrations_ptr(:)
	real(kind=8), pointer :: Density_calculated_ptr(:)
	real(kind=8), pointer :: Gfw_ptr(:)
	real(kind=8), pointer :: Saturation_ptr(:)
	real(kind=8), pointer :: SolutionVolume_ptr(:)
	real(kind=8), pointer :: Time_ptr
	real(kind=8), pointer :: TimeStep_ptr
	real(kind=8), pointer :: Porosity_ptr(:)
	real(kind=8), pointer :: Pressure_ptr(:)
	real(kind=8), pointer :: Temperature_ptr(:)
    type(bmi) :: bmif
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
    yaml_file = "AdvectBMI_f90_test.yaml"
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
    !nxyz = GetGridCellCountYAML(yaml_file)

#ifdef USE_MPI
    ! MPI
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
    if (status .ne. MPI_SUCCESS) then
        stop "Failed to get mpi_myself"
    endif
    if (mpi_myself == 0) then
        nxyz = GetGridCellCountYAML(yaml_file)
    endif
    CALL MPI_Bcast(nxyz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status)
    id = bmif%bmif_create(nxyz, MPI_COMM_WORLD)
    if (mpi_myself > 0) then
        status = RM_MpiWorker(id)
        status = bmif%bmif_finalize()
        return
    endif
#else
    ! OpenMP
    id = bmif%bmif_create_default()
#endif
    ! Initialize with YAML file
    status = bmif%bmif_initialize(yaml_file)
    id = bmif%bmif_get_id()
    status = bmif%bmif_get_value("GridCellCount", nxyz)
    status = bmif%bmif_get_value("ComponentCount", ncomps)
    ! OutputVarNames
    status = bmif%bmif_get_output_var_names(outputvars)
    write(*,*) "Output variables (getters)"
    do i = 1, size(outputvars)
        status = bmif%bmif_get_var_units(outputvars(i), string)
        write(*,"(1x, I4, A60, 2x, A15)") i, trim(outputvars(i)), string
    enddo

	status = bmif%bmif_get_value_ptr("ComponentCount", ComponentCount_ptr)
	status = bmif%bmif_get_value_ptr("GridCellCount", GridCellCount_ptr)
	status = bmif%bmif_get_value_ptr("SelectedOutputOn", SelectedOutputOn_ptr)
	status = bmif%bmif_get_value_ptr("Concentrations", Concentrations_ptr)
	status = bmif%bmif_get_value_ptr("DensityCalculated", Density_calculated_ptr)
	status = bmif%bmif_get_value_ptr("Gfw", Gfw_ptr)
	status = bmif%bmif_get_value_ptr("SaturationCalculated", Saturation_ptr)
	status = bmif%bmif_get_value_ptr("SolutionVolume", SolutionVolume_ptr)
	status = bmif%bmif_get_value_ptr("Time", Time_ptr)
	status = bmif%bmif_get_value_ptr("TimeStep", TimeStep_ptr)
	status = bmif%bmif_get_value_ptr("Porosity", Porosity_ptr)
	status = bmif%bmif_get_value_ptr("Pressure", Pressure_ptr)
	status = bmif%bmif_get_value_ptr("Temperature", Temperature_ptr)
	status = bmif%bmif_get_var_nbytes("FilePrefix", n)
    allocate(character(len=n) :: prefix)
    status = bmif%bmif_get_value("FilePrefix", prefix)
    write(string1, "(A,A)") "File prefix:                                        ", prefix
    status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", GridCellCount_ptr
    status = RM_OutputMessage(id, trim(string1))
    write(string1, "(A,I10)") "Number of components for transport:               ", ComponentCount_ptr
    status = RM_OutputMessage(id, trim(string1))
    ! Get component information
    status = bmif%bmif_get_value("Components", components)
    status = bmif%bmif_get_value("Gfw", gfw)
    do i = 1, ncomps
        write(string,"(A10, F15.4)") trim(components(i)), gfw(i)
        status = RM_OutputMessage(id, string)
    enddo
    status = RM_OutputMessage(id, " ")	
    ! Get initial temperatures
    status = bmif%bmif_get_value("Temperature", temperature)
    ! Get initial temperature
    status = bmif%bmif_get_value("SaturationCalculated", sat)
    ! Get initial porosity
    status = bmif%bmif_get_value("Porosity", por)
    ! Get initial temperature
    status = bmif%bmif_get_value("SolutionVolume", volume)
    ! Get initial concentrations
    ! flattened version
    status = bmif%bmif_get_value("Concentrations", c1)
    c = reshape(c1, (/nxyz, ncomps/))
    ! non-flattened version
    status = bmif%bmif_get_value("Concentrations", c)
    ! Set density, pressure, and temperature (previously allocated)
    allocate(density(nxyz))
    density = 1.0
    status = bmif%bmif_set_value("DensityUser", density)
    allocate(pressure(nxyz))
    pressure = 2.0
    status = bmif%bmif_set_value("Pressure", pressure)  
    temperature = 20.0
    status = bmif%bmif_set_value("Temperature", temperature)  
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
    call compare_ptrs(bmif)
    ! --------------------------------------------------------------------------
    ! Transient loop
    ! --------------------------------------------------------------------------
    nsteps = 10
    time = 0.0
    status = bmif%bmif_set_value("Time", time)
    time_step = 86400.0
    status = bmif%bmif_set_value("TimeStep", time_step)   
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
        call advectionbmi_f90_test(c, bc_conc, ncomps, nxyz)
    
        ! print at last time step
        if (isteps == nsteps) then     
            status = bmif%bmif_set_value("SelectedOutputOn", .true.)    ! enable selected output
            status = RM_SetPrintChemistryOn(id, 1, 0, 0)                ! workers, initial_phreeqc, utility
        else        
            status = bmif%bmif_set_value("SelectedOutputOn", .false.)   ! disable selected output
            status = RM_SetPrintChemistryOn(id, 0, 0, 0)                ! workers, initial_phreeqc, utility
        endif
        ! Transfer data to PhreeqcRM after transport      
        status = bmif%bmif_set_value("Concentrations", c)  ! Transported concentrations
        ! Optionally, if values changed during transport
        status = bmif%bmif_set_value("Porosity", por)              
        status = bmif%bmif_set_value("SaturationUser", sat)            
        status = bmif%bmif_set_value("Temperature", temperature) 
        status = bmif%bmif_set_value("Pressure", pressure)          
        status = bmif%bmif_set_value("TimeStep", time_step) 
        ! Set new time
        time = time + time_step              
        status = bmif%bmif_set_value("Time", time)  ! Current time
        ! Run cells with transported conditions
        write(string, "(A32,F15.1,A)") "Beginning reaction calculation  ", &
            time / 86400., " days"
        status = RM_LogMessage(id, string)
        status = RM_ScreenMessage(id, string)
        ! Demonstration of state 
        status = RM_StateSave(id, 1)
        status = RM_StateApply(id, 1)
        status = RM_StateDelete(id, 1)
        call compare_ptrs(bmif)
        ! Run chemistry
        status = bmif%bmif_update()
        call compare_ptrs(bmif)
        ! Get new data calculated by PhreeqcRM for transport
        status = bmif%bmif_get_value("Concentrations", c)
        status = bmif%bmif_get_value("DensityCalculated", density)
        status = bmif%bmif_get_value("SolutionVolume", volume)   
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
            status = bmif%bmif_get_value("SelectedOutputCount", n)
            do isel = 1, n  ! one based
                i = isel
                status = bmif%bmif_set_value("NthSelectedOutput", i)
                status = bmif%bmif_get_value("CurrentSelectedOutputUserNumber", n_user)
                write(*,*) "Selected output sequence number: ", isel
                write(*,*) "Selected output user number:     ", n_user
                ! Get 2D array of selected output values
                status = bmif%bmif_get_value("SelectedOutputColumnCount", col)
                status = bmif%bmif_get_value("SelectedOutputRowCount", rows)
                allocate(selected_out(rows,col))
                ! Get headings
                status = bmif%bmif_get_var_itemsize("SelectedOutputHeadings", bytes)
                if (allocated(headings)) deallocate(headings)
                allocate(character(len=bytes) :: headings(col))    
                status = bmif%bmif_get_value("SelectedOutputHeadings", headings)
                ! Get selected output
                status = bmif%bmif_get_value("SelectedOutput", selected_out)
                ! Print results
                do i = rows/2, rows/2
                    write(*,*) "Cell number ", i
                    write(*,*) "     Calculated Density: ", density(i)
                    write(*,*) "     Volume:             ", volume(i)
                    write(*,*) "     Components: "
                    do j = 1, ncomps
                        write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ",trim(components(j)), ": ", c(i,j)
                    enddo
                    write(*,*) "     Selected output: "
                    do j = 1, col
                        write(*,'(10x,i3,A2,A50,A2,1pe15.4)') j, " ", trim(headings(j)),": ", selected_out(i,j)
                    enddo
                enddo
                deallocate(selected_out)
            enddo     
            ! Use GetValue to extract exchange composition and pH
			! YAMLAddOutputVars was called in YAML
			! to select additional OutputVarNames variables
			status = bmif%bmif_get_value("solution_ph", pH_vector)
			status = bmif%bmif_get_value("exchange_X_species_log_molality_CaX2", CaX2)
			status = bmif%bmif_get_value("exchange_X_species_log_molality_KX", KX)
			status = bmif%bmif_get_value("exchange_X_species_log_molality_NaX", NaX)
			status = bmif%bmif_get_value("calculate_value_sar", SAR)
            write(string1, "(A)") "      pH      CaX2     KX       NaX       SAR"
            status = RM_OutputMessage(id, string1)
			do i = 1, nxyz
                write(string1,"(5f10.5)") pH_vector(i), 10.0d0**CaX2(i), \
					10.0d0**KX(i), 10.0d0**NaX(i), SAR(i)
            enddo
        endif
    enddo 
    call BMI_testing(bmif)
    ! Clean up
#ifdef USE_MPI    
    status = RM_MpiWorkerBreak(id)
#endif    
    status = bmif%bmif_finalize()
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
    deallocate(pressure)
    
    return
    end subroutine AdvectBMI_f90_test

    SUBROUTINE advectionbmi_f90_test(c, bc_conc, ncomps, nxyz)
    implicit none
    real(kind=8), dimension(:,:), allocatable, intent(inout) :: c
    real(kind=8), dimension(:,:), allocatable, intent(in)    :: bc_conc
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
    END SUBROUTINE advectionbmi_f90_test

#ifdef USE_MPI
    integer(kind=C_INT) function bmi_worker_tasks_f(method_number) BIND(C, NAME='bmi_worker_tasks_f')
    USE ISO_C_BINDING
    implicit none
    interface
    integer function bmi_do_something()
    end function bmi_do_something

    SUBROUTINE register_basic_callback_fortran()
    implicit none
    END SUBROUTINE register_basic_callback_fortran
    end interface
    integer(kind=c_int), intent(in) :: method_number
    integer :: status
    if (method_number .eq. 1000) then
        status = bmi_do_something()
    else if (method_number .eq. 1001) then
        call register_basic_callback_fortran()
    endif
    bmi_worker_tasks_f = 0
    end function bmi_worker_tasks_f

    integer function bmi_do_something()
    USE BMIPhreeqcRM
    implicit none
    INCLUDE 'mpif.h'
    integer status
    integer i, method_number, mpi_myself, mpi_task, mpi_tasks, worker_number
    type(bmi) :: bmif
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
    bmi_do_something = 0
    end function bmi_do_something
#endif

subroutine BMI_testing(bmif)
USE, intrinsic :: ISO_C_BINDING
    USE BMIPhreeqcRM
    USE IPhreeqc
    implicit none
    interface
        integer function assert(tf)
        logical, intent(in) :: tf
        end function assert
    end interface
    !integer, intent(in) :: id
    character(100) :: string
    integer :: status, dim, isel
    integer :: i, j, n, bytes, nbytes
    integer                                              :: nxyz, rm_nxyz
    integer                                              :: ncomps, rm_ncomps
    character(len=:), dimension(:), allocatable          :: components, rm_components
    character(len=:), allocatable                        :: component
    character(len=:), dimension(:), allocatable          :: inputvars
    character(len=:), dimension(:), allocatable          :: outputvars 
    character(len=:), dimension(:), allocatable          :: pointablevars 
    character(len=:), allocatable                        :: prefix, rm_prefix 
    real(kind=8), dimension(:), allocatable          :: gfw, rm_gfw
    real(kind=8), dimension(:,:), allocatable        :: c, c_rm
    integer                                              :: so_count, rm_so_count
    integer                                              :: nuser, rm_nuser
    integer                                              :: col_count, rm_col_count
    integer                                              :: row_count, rm_row_count
    integer                                              :: itemsize
    real(kind=8), dimension(:,:), allocatable        :: so, rm_so
    character(LEN=:), dimension(:), allocatable          :: headings, rm_headings
    character(LEN=:), allocatable                        :: heading
    real(kind=8)                                     :: time, rm_time
    real(kind=8)                                     :: time_step, rm_time_step
    real(kind=8), dimension(:), allocatable          :: density_calculated, rm_density_calculated
    real(kind=8), dimension(:), allocatable          :: porosity, rm_porosity
    real(kind=8), dimension(:), allocatable          :: pressure, rm_pressure
    real(kind=8), dimension(:), allocatable          :: saturation, rm_saturation
    real(kind=8), dimension(:), allocatable          :: temperature, rm_temperature
    real(kind=8), dimension(:), allocatable          :: volume, rm_volume
    logical :: tf
    character(LEN=:), allocatable :: alloc_string
    real(kind=8), pointer :: real_ptr
    real(kind=c_double), pointer :: real_dim_ptr(:)
    integer, pointer :: integer_ptr
    logical, pointer :: logical_ptr
    type(bmi), intent(inout) :: bmif
    integer :: id
    
    id = bmif%bmif_get_id()
    status = bmif%bmif_get_component_name(string)
    write(*,*) trim(string)
    time = bmif%bmif_get_current_time(time)
    write(*,*) time
    time = bmif%bmif_get_end_time(time)
    write(*,*) time
    status = bmif%bmif_get_time_units(string)
    write(*,*) string
    time = bmif%bmif_get_time_step(time)
    write(*,*) time
    ! Time
    rm_time = 3600
    status = bmif%bmif_set_value("Time", rm_time)
    status = bmif%bmif_get_value("Time", time)
    status = assert(rm_time .eq. time)
    rm_time = RM_GetTime(id)
    status = bmif%bmif_get_value_ptr("Time", real_ptr)
    status = assert(time .eq. rm_time)
    status = bmif%bmif_get_current_time(time)
    status = assert(time .eq. rm_time)
    ! TimeStep
    rm_time_step = 60
    status = bmif%bmif_set_value("TimeStep", rm_time_step)
    status = bmif%bmif_get_value("TimeStep", time_step)
    status = assert(time_step .eq. rm_time_step)
    rm_time_step = RM_GetTimeStep(id)
    status = assert(time_step .eq. rm_time_step)
    status = bmif%bmif_get_time_step(time_step)
    status = assert(time_step .eq. rm_time_step)
    ! InputVarNames
    status = bmif%bmif_get_input_var_names(inputvars)
    write(*,*) "Input variables (setters)"
    do i = 1, size(inputvars)
        write(*,"(1x, I4, A60)") i, trim(inputvars(i))
        !status = bmif%bmif_get_var_units(inputvars(i), string)
        !write(*,"(5x, A60)") trim(string)
        !status = bmif%bmif_get_var_type(inputvars(i), string)
        !write(*,"(5x, A60)") trim(string)
        !status = bmif%bmif_get_var_itemsize(inputvars(i), itemsize)
        !write(*, "(5x, I60)") itemsize
        !status = bmif%bmif_get_var_nbytes(inputvars(i), nbytes)
        !write(*, "(5x, I60)") nbytes
    enddo
    ! OutputVarNames
    status = bmif%bmif_get_output_var_names(outputvars)
    write(*,*) "Output variables (getters)"
    do i = 1, size(outputvars)
        write(*,"(1x, I4, A60)") i, trim(outputvars(i))
        !status = bmif%bmif_get_var_units(outputvars(i), string)
        !write(*,"(5x, A60)") trim(string)
        !status = bmif%bmif_get_var_type(outputvars(i), string)
        !write(*,"(5x, A60)") trim(string)
        !status = bmif%bmif_get_var_itemsize(outputvars(i), itemsize)
        !write(*, "(5x, I60)") itemsize
        !status = bmif%bmif_get_var_nbytes(outputvars(i), nbytes)
        !write(*, "(5x, I60)") nbytes
    enddo
    ! PointableVarNames
    status = bmif%bmif_get_pointable_var_names(pointablevars)
    write(*,*) "Pointable variables (GetValuePtr)"
    do i = 1, size(pointablevars)
        write(*,"(1x, I4, A60)") i, trim(pointablevars(i))
    enddo
    ! ComponentCount
    status = bmif%bmif_get_value("ComponentCount", ncomps)
    rm_ncomps = RM_GetComponentCount(id)
    status = assert(ncomps .eq. rm_ncomps)
    ! Components
    status = bmif%bmif_get_value("Components", components)
    status = bmif%bmif_get_var_itemsize("Components", itemsize)
    status = bmif%bmif_get_var_nbytes("Components", nbytes)
    dim = nbytes / itemsize
    status = assert(dim .eq. size(components))
    !allocate(character(len=itemsize) :: component)
    !do i = 1, dim
    !    j = i
    !    status = RM_GetComponent(j, component)
    !    status = assert(component .eq. components(i))
    !enddo
    status = RM_GetComponents(id, rm_components)
    do i = 1, dim
        status = assert(components(i) .eq. rm_components(i))
    enddo

    status = bmif%bmif_get_value("Gfw", gfw)
    do i = 1, size(components)
        write(*,"(A10, F15.4)") trim(components(i)), gfw(i)
    enddo
    write(*,*)	   
    
    ! Concentrations
    status = bmif%bmif_get_value("GridCellCount", nxyz)
    status = assert(nxyz .eq. RM_GetGridCellCount(id))
    status = bmif%bmif_get_value("Concentrations", c)
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
   
	! GetValue("DensityCalculated")
    ! RM_GetDensityCalculated and bmif%bmif_get_value("DensityCalculated) always return 
    ! the calculated solution density
    status = bmif%bmif_get_var_itemsize("DensityCalculated", itemsize)
    status = bmif%bmif_get_var_nbytes("DensityCalculated", nbytes)
    dim = nbytes / itemsize
    allocate(rm_density_calculated(dim))
    status = RM_GetDensityCalculated(id, rm_density_calculated)
    status = bmif%bmif_get_value("DensityCalculated", density_calculated)
    do i = 1, nxyz
        if (density_calculated(i) .ne. rm_density_calculated(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
    ! FilePrefix
    string = "NewPrefix"
    status = bmif%bmif_set_value("FilePrefix", string)
    status = bmif%bmif_get_var_itemsize("FilePrefix", itemsize)
    status = bmif%bmif_get_var_nbytes("FilePrefix", nbytes)
    allocate(character(len=itemsize) :: rm_prefix)
    status = assert(itemsize .eq. nbytes)
    status = RM_GetFilePrefix(id, rm_prefix)
    status = assert(string .eq. rm_prefix)
    allocate(character(len=itemsize) :: prefix)
    status = bmif%bmif_get_value("FilePrefix", prefix)
    status = assert(prefix .eq. rm_prefix)
         
	! GetValue("Gfw")
    status = bmif%bmif_get_value("Gfw", gfw)
    status = bmif%bmif_get_var_itemsize("Gfw", itemsize)
    status = bmif%bmif_get_var_nbytes("Gfw", nbytes)
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
    status = bmif%bmif_get_value("GridCellCount", nxyz)
    rm_nxyz = RM_GetGridCellCount(id)
    status = assert(nxyz .eq. rm_nxyz)
    
	! GetValue("Porosity")
    status = bmif%bmif_get_value_ptr("Porosity", real_dim_ptr)
    status = bmif%bmif_get_var_itemsize("Porosity", itemsize)
    status = bmif%bmif_get_var_nbytes("Porosity", nbytes)
    dim = nbytes / itemsize
    allocate(rm_porosity(dim))
    rm_porosity = 0.25
    status = bmif%bmif_set_value("Porosity", rm_porosity)
    status = bmif%bmif_get_value("Porosity", porosity)
    do i = 1, nxyz
        status = assert(porosity(i) .eq. rm_porosity(i)) 
        status = assert(porosity(i) .eq. real_dim_ptr(i)) 
    enddo
	status = RM_GetPorosity(id, rm_porosity);
        do i = 1, nxyz
        if (porosity(i) .ne. rm_porosity(i)) then
            status = assert(.false.)
            exit
        endif
    enddo

	! GetValue("Pressure")
    status = bmif%bmif_get_var_itemsize("Pressure", itemsize)
    status = bmif%bmif_get_var_nbytes("Pressure", nbytes)
    dim = nbytes / itemsize
    allocate(rm_pressure(dim))
    rm_pressure = 10.
    status = bmif%bmif_set_value("Pressure", rm_pressure)
    status = bmif%bmif_get_value("Pressure", pressure)
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

	! GetValue("SaturationCalculated")
    ! Always returns solution_volume/(rv * porosity) for each cell
    status = bmif%bmif_get_var_itemsize("SaturationCalculated", itemsize)
    status = bmif%bmif_get_var_nbytes("SaturationCalculated", nbytes)
    dim = nbytes / itemsize
    allocate(rm_saturation(dim))
    status = bmif%bmif_get_value("SaturationCalculated", saturation)
	status = RM_GetSaturationCalculated(id, rm_saturation);
    do i = 1, nxyz
        if (saturation(i) .ne. rm_saturation(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
    ! GetValue("SolutionVolume")
    status = bmif%bmif_get_var_itemsize("SolutionVolume", itemsize)
    status = bmif%bmif_get_var_nbytes("SolutionVolume", nbytes)
    dim = nbytes / itemsize
    allocate(rm_volume(dim))
    status = bmif%bmif_get_value("SolutionVolume", volume)
	status = RM_GetSolutionVolume(id, rm_volume);
    do i = 1, nxyz
        if (volume(i) .ne. rm_volume(i)) then
            status = assert(.false.)
            exit
        endif
    enddo
    
	! GetValue("Temperature")
    status = bmif%bmif_get_var_itemsize("Temperature", itemsize)
    status = bmif%bmif_get_var_nbytes("Temperature", nbytes)
    dim = nbytes / itemsize
    allocate(rm_temperature(dim))
    rm_temperature = 11.
    status = bmif%bmif_set_value("Temperature", rm_temperature)
    status = bmif%bmif_get_value("Temperature", temperature)
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
    status = bmif%bmif_get_value("SelectedOutputCount", so_count);
	rm_so_count = RM_GetSelectedOutputCount(id);
    status = assert(so_count .eq. rm_so_count)
 
    do isel = 1, so_count ! one based
        i = isel
		status = bmif%bmif_set_value("NthSelectedOutput", i)
        status = bmif%bmif_get_value("CurrentSelectedOutputUserNumber", nuser)
        rm_nuser = RM_GetCurrentSelectedOutputUserNumber(id)
        status = assert(nuser .eq. rm_nuser)
        i = isel        
		rm_nuser = RM_GetNthSelectedOutputUserNumber(id, i)
        status = assert(nuser .eq. rm_nuser)

        status = bmif%bmif_get_value("SelectedOutputColumnCount", col_count)
        rm_col_count =RM_GetSelectedOutputColumnCount(id)
        status = assert(col_count .eq. rm_col_count)

        status = bmif%bmif_get_value("SelectedOutputRowCount", row_count)
        rm_row_count = RM_GetSelectedOutputRowCount(id)
        status = assert(row_count .eq. rm_row_count)

        status = bmif%bmif_get_var_nbytes("SelectedOutput", nbytes);
        status = bmif%bmif_get_var_itemsize("SelectedOutput", itemsize);
        dim = nbytes / itemsize;
        status = assert(dim .eq. rm_row_count*rm_col_count)
        status = bmif%bmif_get_value("SelectedOutput", so)
        if (allocated(rm_so)) deallocate(rm_so)
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
        status = bmif%bmif_get_value("SelectedOutputHeadings", headings)
        status = RM_GetSelectedOutputHeadings(id, rm_headings)
        do j = 1, col_count
            if (headings(j) .ne. rm_headings(j)) then
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
write(*,*) "Assert failed"
call exit(-1)
end function assert

subroutine compare_ptrs(bmif)
USE, intrinsic :: ISO_C_BINDING
USE BMIPhreeqcRM
USE IPhreeqc
implicit none
    interface
        integer function assert(tf)
        logical, intent(in) :: tf
        end function assert
    end interface
    common /i_ptrs/ ncomps, nxyz, SelectedOutputOn_ptr, id
    common /r_ptrs/ Concentrations_ptr, Density_calculated_ptr, Gfw_ptr, &
	Saturation_ptr, SolutionVolume_ptr, Time_ptr, TimeStep_ptr, &
    Porosity_ptr, Pressure_ptr, Temperature_ptr
    integer :: id
	integer, pointer :: ncomps
	integer, pointer :: nxyz
	logical(kind=1), pointer :: SelectedOutputOn_ptr
	real(kind=8), pointer :: Concentrations_ptr(:)
	real(kind=8), pointer :: Density_calculated_ptr(:)
	real(kind=8), pointer :: Gfw_ptr(:)
	real(kind=8), pointer :: Saturation_ptr(:)
	real(kind=8), pointer :: SolutionVolume_ptr(:)
	real(kind=8), pointer :: Time_ptr
	real(kind=8), pointer :: TimeStep_ptr
	real(kind=8), pointer :: Porosity_ptr(:)
    real(kind=8), pointer :: Pressure_ptr(:)
    real(kind=8), pointer :: Temperature_ptr(:)

    integer :: status, i, dim
    integer :: componentcount, gridcellcount
    real(kind=8), allocatable, dimension(:) :: gfw, density_calculated, saturation, &
        SolutionVolume, Porosity, Pressure, Temperature, concentrations
    real(kind=8) :: time, timestep
    logical :: selectedoutputon, test_logical
    type(bmi), intent(inout) :: bmif
    ! ComponentCount
    id = bmif%bmif_get_id()
    status = bmif%bmif_get_value("ComponentCount", componentcount)
    status = assert(ncomps .eq. componentcount)
    componentcount = RM_GetComponentCount(id)
    status = assert(ncomps .eq. componentcount)
    ! Gfw
    status = bmif%bmif_get_value("Gfw", gfw)
	do i = 1, ncomps
		status = assert(Gfw_ptr(i) .eq. gfw(i))
    enddo
    status = RM_GetGfw(id, gfw)
	do i = 1, ncomps
		status = assert(Gfw_ptr(i) .eq. gfw(i))
    enddo
    ! GridCellCount
    status = bmif%bmif_get_value("GridCellCount", gridcellcount)
	status = assert(nxyz .eq. gridcellcount)
    ! Density, Saturation, SolutionVolume, Porosity, Pressure, Temperature
	status = bmif%bmif_get_value("DensityCalculated", density_calculated)
	status = bmif%bmif_get_value("SaturationCalculated", saturation)
	status = bmif%bmif_get_value("solutionvolume", solutionvolume)
	status = bmif%bmif_get_value("Porosity", Porosity)
	status = bmif%bmif_get_value("Pressure", Pressure)
	status = bmif%bmif_get_value("Temperature", Temperature)

	do i = 1, nxyz
		status = assert(Density_calculated_ptr(i) .eq. density_calculated(i))
		status = assert(Saturation_ptr(i) .eq. saturation(i))
		status = assert(SolutionVolume_ptr(i) .eq. SolutionVolume(i))
		status = assert(Porosity_ptr(i) .eq. Porosity(i))
		status = assert(Pressure_ptr(i) .eq. Pressure(i))
		status = assert(Temperature_ptr(i) .eq. Temperature(i))
    enddo   
	status = RM_GetDensityCalculated(id, density_calculated)
	status = RM_GetSaturationCalculated(id, saturation)
	status = RM_GetSolutionVolume(id, solutionvolume)
	status = RM_GetPorosity(id, Porosity)
	status = RM_GetPressure(id, Pressure)
	status = RM_GetTemperature(id, Temperature)
	do i = 1, nxyz
		status = assert(Density_calculated_ptr(i) .eq. density_calculated(i))
		status = assert(Saturation_ptr(i) .eq. saturation(i))
		status = assert(SolutionVolume_ptr(i) .eq. SolutionVolume(i))
		status = assert(Porosity_ptr(i) .eq. Porosity(i))
		status = assert(Pressure_ptr(i) .eq. Pressure(i))
		status = assert(Temperature_ptr(i) .eq. Temperature(i))
    enddo   
    ! Concentrations
	status = bmif%bmif_get_value("Concentrations", Concentrations)
	dim = ncomps * nxyz
	do i = 1, dim
		status = assert(Concentrations_ptr(i) .eq. Concentrations(i))
    enddo
    ! Time
    status = bmif%bmif_get_value("Time", time)
	status = assert(Time_ptr .eq. time)
    status = bmif%bmif_get_value("TimeStep", timestep)
	status = assert(TimeStep_ptr .eq. timestep)
    ! SelectedOutputOn
    status = bmif%bmif_get_value("SelectedOutputOn", selectedoutputon)
    test_logical = SelectedOutputOn_ptr
	status = assert(test_logical .eqv. selectedoutputon)
    end subroutine compare_ptrs  
#endif 
