#ifdef USE_YAML
    subroutine AdvectBMI_f90()  BIND(C, NAME='AdvectBMI_f90')
    USE, intrinsic :: ISO_C_BINDING
    USE BMIPhreeqcRM
    implicit none
#ifdef USE_MPI    
    INCLUDE 'mpif.h'
#endif
    interface
        subroutine advectionbmi_f90(c, bc_conc, ncomps, nxyz)
            implicit none
            real(kind=8), dimension(:,:), allocatable, intent(inout) :: c
            real(kind=8), dimension(:,:), allocatable, intent(in) :: bc_conc
            integer, intent(in)                                       :: ncomps, nxyz
        end subroutine advectionbmi_f90
    end interface

    ! Based on PHREEQC Example 11
    real(kind=8), pointer :: d1_ptr(:)
    character(100) :: yaml_file
    integer :: mpi_myself
    integer :: i, j
    logical :: tf
    integer :: nxyz, ncomps
    integer :: nthreads
    integer :: status
    integer :: bytes, nbytes
    real(kind=8), dimension(:), allocatable         :: por
    real(kind=8), dimension(:), allocatable         :: sat
    character(len=:), allocatable                   :: prefix
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
    integer                                         :: vtype
    real(kind=8)                                    :: pH
    character(100)                                  :: svalue
    integer                                         :: iphreeqc_id, iphreeqc_id1
    integer                                         :: dump_on, append
    integer                                         :: dim
	real(kind=8), dimension(:), allocatable         :: CaX2, KX, NaX, pH_vector
    integer :: id
	integer, pointer :: ComponentCount_ptr
	integer, pointer :: GridCellCount_ptr
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
    yaml_file = "AdvectBMI_f90.yaml"
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
    nxyz = GetGridCellCountYAML(yaml_file)
    nthreads = 3
    
    id =bmif%bmif_create(nxyz, nthreads)
#endif
    status = bmif%bmif_initialize(yaml_file)
    
    ! Open files
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
    status = bmif%bmif_get_value("ComponentCount", ncomps)
    ! Print some of the reaction module information
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
    do i = 1, ComponentCount_ptr
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
        call advectionbmi_f90(c, bc_conc, ncomps, nxyz)
    
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
        ! Run chemistry
        status = bmif%bmif_update()
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
            write(string1, "(A)") "      pH      CaX2     KX       NaX"
            status = RM_OutputMessage(id, string1)
			do i = 1, nxyz
                write(string1,"(5f10.5)") pH_vector(i), 10.0d0**CaX2(i), \
					10.0d0**KX(i), 10.0d0**NaX(i)
            enddo
        endif
    enddo 
    ! Clean up
    !status = RM_CloseFiles(id)
    status = RM_MpiWorkerBreak(id)
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
    end subroutine AdvectBMI_f90

    SUBROUTINE advectionbmi_f90(c, bc_conc, ncomps, nxyz)
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
    END SUBROUTINE advectionbmi_f90
#endif 
