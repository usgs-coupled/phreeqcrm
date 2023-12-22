    !> @file BMI_interface.F90
    !> @brief BMIPhreeqcRM module definition 
    !>
    !
    !*MODULE BMIPhreeqcRM PHREEQC Reaction Module for Transport Codes
    !> @brief Fortran documentation for the geochemical reaction module BMIPhreeqcRM.
    !> @par ""
    !> @n "USE BMIPhreeqcRM" defines a module with both Basic Model Interface methods and 
    !> native PhreeqcRM methods for Fortran programs. 
    !> @n For Windows, 
    !> include the files BMI_interface.F90 and RM_interface.F90 in your project.
    !> @n For Linux, configure, compile, and install the PhreeqcRM library and module file.
    !> You will need installed include directory (-I) added to the project to 
    !> reference the module file.
    !> You will need to link to the library to produce the executable for your code.
    !>
    MODULE BMIPhreeqcRM
    USE PhreeqcRM
    IMPLICIT NONE
    PRIVATE :: Lower
    PRIVATE :: success

    integer, parameter :: BMI_MAX_COMPONENT_NAME = 2048
    integer, parameter :: BMI_MAX_VAR_NAME = 2048
    integer, parameter :: BMI_MAX_TYPE_NAME = 2048
    integer, parameter :: BMI_MAX_UNITS_NAME = 2048

    integer, parameter :: BMI_FAILURE = 1
    integer, parameter :: BMI_SUCCESS = 0
    !> @a bmif_create creates a BMIPhreeqcRM instance. 
    !>
    !> Two constructors are available to create a BMIPhreeqcRM
    !> instance. The default method, with no arguments, can be
    !> used if the instance is to be initialized with a YAML file.
    !> The YAML file must provide the number of cells in the user's 
    !> model through the YAMLSetGridCellCount method, and, optionally,
    !> the number of threads with YAMLThreadCount. The other
    !> constructor requires the number of cells in the user model
    !> as the first argument, and the number of threads as the
    !> second argument.  

  !  INTERFACE bmif_create
		!module procedure bmif_create_default 
		!module procedure bmif_create  
  !  END INTERFACE bmif_create


    public bmi
!!!!!!!!!!
    type :: bmi 
!!!!!!!!!!
        INTEGER :: bmiphreeqcrm_id = -1
    contains
		procedure :: bmif_create_default 
		procedure :: bmif_create  
        procedure :: bmif_get_id
        procedure :: bmif_add_output_vars
        procedure :: bmif_finalize
        procedure :: bmif_get_component_name
        procedure :: bmif_get_current_time
        procedure :: bmif_get_end_time
        procedure :: bmif_get_input_item_count
        procedure :: bmif_get_input_var_names
        procedure :: bmif_get_output_item_count
        procedure :: bmif_get_output_var_names
        procedure :: bmif_get_pointable_item_count
        procedure :: bmif_get_pointable_var_names
        procedure :: bmif_grid_rank
        procedure :: bmif_grid_size
        procedure :: bmif_grid_type
        procedure :: bmif_get_start_time
        procedure :: bmif_get_time_step
        procedure :: bmif_get_time_units
        procedure :: bmif_get_value_at_indices_double, &
            bmif_get_value_at_indices_float, bmif_get_value_at_indices_int
        generic :: bmif_get_value_at_indicies => bmif_get_value_at_indices_double, &
            bmif_get_value_at_indices_float, bmif_get_value_at_indices_int
        procedure :: bmif_get_value_logical, bmif_get_value_char, bmif_get_value_char1, &
            bmif_get_value_double, bmif_get_value_double1, bmif_get_value_double2, bmif_get_value_float, &! not implemented
            bmif_get_value_int, bmif_get_value_int1, bmif_get_value_int2
        generic :: bmif_get_value => bmif_get_value_logical, bmif_get_value_char, bmif_get_value_char1, &
            bmif_get_value_double, bmif_get_value_double1, bmif_get_value_double2, bmif_get_value_float, &
            bmif_get_value_int, bmif_get_value_int1, bmif_get_value_int2
        procedure :: bmif_get_value_ptr_logical, bmif_get_value_ptr_integer, &
             bmif_get_value_ptr_double, bmif_get_value_ptr_double1, bmif_get_value_ptr_float
        generic :: bmif_get_value_ptr => bmif_get_value_ptr_logical, bmif_get_value_ptr_integer, &
             bmif_get_value_ptr_double, bmif_get_value_ptr_double1, bmif_get_value_ptr_float
        procedure :: bmif_get_var_location
        
        procedure :: bmif_set_value_at_indices_double, bmif_set_value_at_indices_float, &
            bmif_set_value_at_indices_int
        generic :: bmif_set_value_at_indices => bmif_set_value_at_indices_double, &
            bmif_set_value_at_indices_float, bmif_set_value_at_indices_int
              
        procedure :: bmif_get_var_itemsize
        procedure :: bmif_get_var_nbytes
        procedure :: bmif_get_var_type
        procedure :: bmif_get_var_units
        procedure :: bmif_initialize
        !generic :: bmif_initialize => bmif_initialize, bmif_initialize_default ! procedure declaration

        procedure :: bmif_set_value_b, bmif_set_value_c, bmif_set_value_double, bmif_set_value_double1, &
            bmif_set_value_double2, bmif_set_value_float, &! not implemented
            bmif_set_value_int, bmif_set_value_int1, bmif_set_value_int2
        generic :: bmif_set_value => bmif_set_value_b, bmif_set_value_c, bmif_set_value_double, bmif_set_value_double1, &
            bmif_set_value_double2, bmif_set_value_float, &! not implemented
            bmif_set_value_int, bmif_set_value_int1, bmif_set_value_int2
        
        procedure :: bmif_update
        procedure :: bmif_update_until
        procedure :: bmif_get_grid_shape
        procedure :: bmif_get_grid_spacing
        procedure :: bmif_get_grid_origin
        procedure :: bmif_get_grid_x
        procedure :: bmif_get_grid_y
        procedure :: bmif_get_grid_z
        procedure :: bmif_get_grid_node_count
        procedure :: bmif_get_grid_edge_count
        procedure :: bmif_get_grid_face_count
        procedure :: bmif_get_grid_edge_nodes
        procedure :: bmif_get_grid_face_edges
        procedure :: bmif_get_grid_face_nodes
        procedure :: bmif_get_grid_nodes_per_face
        
        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#define EXTEND_BMIPHREEQCRM
#ifdef EXTEND_BMIPHREEQCRM        

	!procedure :: getgridcellcountyaml
	procedure :: abort
	procedure :: closefiles
	procedure :: concentrations2utility
	procedure :: createmapping
	procedure :: decodeerror
	procedure :: destroy
	procedure :: dumpmodule
	procedure :: errormessage
	procedure :: findcomponents
	procedure :: getbackwardmapping
	procedure :: getchemistrycellcount
	procedure :: getcomponents                      ! bmif_get_var
	procedure :: getcomponentcount                  ! bmif_get_var
	procedure :: getconcentrations                  ! bmif_get_var
	procedure :: getcurrentselectedoutputusernumber ! bmif_get_var
	procedure :: getdensitycalculated               ! bmif_get_var
	procedure :: GetDensity                         ! bmif_get_var
	procedure :: getendcell
	procedure :: getequilibriumphasescount
	procedure :: getequilibriumphasesnames
	procedure :: geterrorstring                     ! bmif_get_var
	procedure :: getexchangenames
	procedure :: getexchangespeciescount
	procedure :: getexchangespeciesnames
	procedure :: getfileprefix                      ! bmif_get_var
	procedure :: getgascomponentscount
	procedure :: getgascomponentsnames
	procedure :: getgascompmoles
	procedure :: getgascomppressures
	procedure :: getgascompphi
	procedure :: getgasphasevolume
	procedure :: getgfw                             ! bmif_get_var
	procedure :: getgridcellcount                   ! bmif_get_var
	procedure :: getiphreeqcid
	procedure :: getithconcentration
	procedure :: getithspeciesconcentration
	procedure :: getkineticreactionscount
	procedure :: getkineticreactionsnames
	procedure :: getmpimyself
	procedure :: getmpitasks
	procedure :: getnthselectedoutputusernumber
	procedure :: getporosity                        ! bmif_get_var
	procedure :: getpressure                        ! bmif_get_var
	procedure :: getsaturationcalculated
	procedure :: GetSaturation
	procedure :: getselectedoutput                  ! bmif_get_var
	procedure :: getselectedoutputcolumncount       ! bmif_get_var
	procedure :: getselectedoutputcount             ! bmif_get_var
	procedure :: getselectedoutputheadings          ! bmif_get_var
    
	procedure :: getselectedoutputrowcount          ! bmif_get_var
	procedure :: getsicount
	procedure :: getsinames
	procedure :: getsolidsolutioncomponentscount
	procedure :: getsolidsolutioncomponentsnames
	procedure :: getsolidsolutionnames
	procedure :: getsolutionvolume                  ! bmif_get_var
	procedure :: getspeciesconcentrations
	procedure :: getspeciescount
	procedure :: getspeciesd25
	procedure :: getspecieslog10gammas
	procedure :: getspecieslog10molalities
	procedure :: getspeciesnames
	procedure :: getspeciessaveon
	procedure :: getspeciesz
	procedure :: getstartcell
	procedure :: getsurfacenames
	procedure :: getsurfacespeciescount
	procedure :: getsurfacespeciesnames
	procedure :: getsurfacetypes
	procedure :: gettemperature                     ! bmif_get_var
	procedure :: getthreadcount
	procedure :: gettime                            ! bmif_get_var
	procedure :: gettimeconversion
	procedure :: gettimestep                        ! bmif_get_var
	procedure :: getviscosity                       ! bmif_get_var
	procedure :: initializeyaml
	procedure :: initialphreeqc2concentrations
	procedure :: initialphreeqc2module
	procedure :: initialsolutions2module
	procedure :: initialequilibriumphases2module
	procedure :: initialexchanges2module
	procedure :: initialgasphases2module
	procedure :: initialsolidsolutions2module
	procedure :: initialsurfaces2module
	procedure :: initialkinetics2module
	procedure :: initialphreeqc2speciesconcentrations
	procedure :: initialphreeqccell2module
	procedure :: loaddatabase
	procedure :: logmessage
	procedure :: mpiworker
	procedure :: mpiworkerbreak
	procedure :: openfiles
	procedure :: outputmessage
	procedure :: runcells
	procedure :: runfile
	procedure :: runstring
	procedure :: screenmessage
	procedure :: setcomponenth2o
	procedure :: setconcentrations                  ! bmif_set_var
	procedure :: setcurrentselectedoutputusernumber
	procedure :: setdensityuser                     ! bmif_set_var
	procedure :: SetDensity                         ! bmif_set_var
	procedure :: setdumpfilename
	procedure :: seterrorhandlermode
	procedure :: seterroron
	procedure :: setfileprefix                      ! bmif_set_var
	procedure :: setgascompmoles
	procedure :: setgasphasevolume
	procedure :: setithconcentration
	procedure :: setithspeciesconcentration
	procedure :: SetmpiworkerCallback
	procedure :: setnthselectedoutput               ! bmif_set_var
	procedure :: setpartitionuzsolids
	procedure :: setporosity                        ! bmif_set_var
	procedure :: setpressure                        ! bmif_set_var
	procedure :: setprintchemistrymask
	procedure :: setprintchemistryon
	procedure :: setrebalancebycell
	procedure :: setrebalancefraction
	procedure :: setrepresentativevolume
	procedure :: setsaturationuser                  ! bmif_set_var
	procedure :: SetSaturation                      ! bmif_set_var
	procedure :: setscreenon
	procedure :: setselectedoutputon                ! bmif_set_var
	procedure :: setspeciessaveon
	procedure :: settemperature                     ! bmif_set_var
	procedure :: settime                            ! bmif_set_var
	procedure :: settimeconversion
	procedure :: settimestep                        ! bmif_set_var
	procedure :: setunitsexchange
	procedure :: setunitsgasphase
	procedure :: setunitskinetics
	procedure :: setunitsppassemblage
	procedure :: setunitssolution
	procedure :: setunitsssassemblage
	procedure :: setunitssurface
	procedure :: speciesconcentrations2module
	procedure :: statesave
	procedure :: stateapply
	procedure :: statedelete
	procedure :: usesolutiondensityvolume
	procedure :: warningmessage
#endif
!End EXTEND_BMIPHREEQCRM
    
    end type         
    CONTAINS

    ! ====================================================
    ! Initialize, run, finalize (IRF) 
    ! ====================================================
    !> @a bmif_create creates a BMIPhreeqcRM instance. 
    !> The default method, with no arguments, can be
    !> used if the instance is to be initialized with a YAML file.
    !> The YAML file must provide the number of cells in the user's 
    !> model through the YAMLSetGridCellCount method, and, optionally,
    !> the number of threads with YAMLThreadCount. The default
    !> method cannot be used with MPI.
    INTEGER FUNCTION bmif_create_default(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RM_BMI_create_default() &
        BIND(C, NAME='RM_BMI_create_default')
    USE ISO_C_BINDING
    IMPLICIT NONE
    END FUNCTION RM_BMI_create_default
    END INTERFACE
    class(bmi), intent(inout) :: self
#if defined(USE_MPI)
    bmif_create_default = -1
    STOP "ERROR: You must use bmif_create(nxyz, COMM) when using MPI."
#endif
    bmif_create_default = RM_BMI_create_default() 
    self%bmiphreeqcrm_id = bmif_create_default
    return
    END FUNCTION bmif_create_default 
    
    !> @a bmif_create creates a reaction module. If the code is compiled with
    !> the preprocessor directive USE_OPENMP, the reaction module is multithreaded.
    !> If the code is compiled with the preprocessor directive USE_MPI, the reaction
    !> module will use MPI and multiple processes. If neither preprocessor directive is used,
    !> the reaction module will be serial (unparallelized).
    !> @param nxyz                         The number of grid cells in the user's model.
    !> @param nthreads (or @a comm, MPI)   When using OPENMP, the argument (@a nthreads) 
    !> is the number of worker threads to be used.
    !> If @a nthreads <= 0, the number of threads is set equal to the number of 
    !> processors of the computer.
    !> When using MPI, the argument (@a comm) is the MPI communicator to use within 
    !> the reaction module.
    !> @retval Id of the BMIPhreeqcRM instance, negative is failure.
    !> @see
    !> @ref bmif_finalize.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> nxyz = 40
    !> nthreads = 3
    !> id = brm%bmif_create(nxyz, nthreads)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.
    INTEGER FUNCTION bmif_create(self, nxyz, nthreads)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RM_BMI_create(nxyz, nthreads) &
        BIND(C, NAME='RM_BMI_create')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: nxyz
    INTEGER(KIND=C_INT), INTENT(in) :: nthreads
    END FUNCTION RM_BMI_create
    END INTERFACE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: nxyz
    INTEGER, INTENT(in) :: nthreads
    bmif_create = RM_BMI_create(nxyz, nthreads)
    self%bmiphreeqcrm_id = bmif_create
    return
    END FUNCTION bmif_create   

    !> @a bmif_initialize uses a YAML file to initialize an instance of BMIPhreeqcRM. Same as
    !> @ref InitializeYAML.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param config_file   String containing the YAML file name.
    !> @retval              0 is success, 1 is failure.
    !> @par
    !> The file contains a YAML map of PhreeqcRM methods
    !> and the arguments corresponding to the methods.
    !> Note that the PhreeqcRM methods do not have the "RM_" prefix
    !> and the id argument is not included.
    !> For example,
    !> @par
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !>  key: loaddatabase
    !>   database: phreeqc.dat
    !> - key: runfile
    !> workers: true
    !> initial_phreeqc: true
    !>   utility: true
    !>  chemistry_name: advect.pqi
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par
    !> @ref bmif_initialize will read the YAML file and execute the specified methods with
    !> the specified arguments. Using YAML
    !> terminology, the argument(s) for a method may be a scalar, a sequence, or a map,
    !> depending if the argument is
    !> a single item, a single vector, or there are multiple arguments.
    !> In the case of a map, the name associated
    !> with each argument (for example "chemistry_name" above) is arbitrary.
    !> The names of the map keys for map
    !> arguments are not used in parsing the YAML file; only the order of
    !> the arguments is important.
    !> @par
    !> The following list gives the PhreeqcRM methods that can be specified in a YAML file
    !> and the arguments that are required. The arguments are described with C++ formats, which
    !> are sufficient to identify which arguments are YAML scalars (single bool/logical,
    !> int, double, string/character argument),
    !> sequences (single vector argument), or maps (multiple arguments).
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> closefiles();
    !> createmapping(std::vector< int >& grid2chem);
    !> dumpmodule();
    !> findcomponents();
    !> initialequilibriumphases2module(std::vector< int > equilibrium_phases);
    !> initialexchanges2module(std::vector< int > exchanges);
    !> initialgasphases2module(std::vector< int > gas_phases);
    !> InitialKineticss2Module(std::vector< int > kinetics);
    !> initialsolidsolutions2module(std::vector< int > solid_solutions);
    !> initialsolutions2module(std::vector< int > solutions);
    !> initialsurfaces2module(std::vector< int > surfaces);
    !> initialphreeqc2module(std::vector< int > initial_conditions1);
    !> initialphreeqc2module(std::vector< int > initial_conditions1, 
    !> std::vector< int > initial_conditions2, std::vector< double > fraction1);
    !> initialphreeqccell2module(int n, std::vector< int > cell_numbers);
    !> loaddatabase(std::string database);
    !> openfiles();
    !> outputmessage(std::string str);
    !> runcells();
    !> runfile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
    !> runstring(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    !> screenmessage(std::string str);
    !> setcomponenth2o(bool tf);
    !> setconcentrations(std::vector< double > c);
    !> setcurrentselectedoutputusernumber(int n_user);
    !> setdensityuser(std::vector< double > density);
    !> setdumpfilename(std::string dump_name);
    !> seterrorhandlermode(int mode);
    !> seterroron(bool tf);
    !> setfileprefix(std::string prefix);
    !> setgascompmoles(std::vector< double > gas_moles);
    !> setgasphasevolume(std::vector< double > gas_volume);
    !> setpartitionuzsolids(bool tf);
    !> setporosity(std::vector< double > por);
    !> setpressure(std::vector< double > p);
    !> setprintchemistrymask(std::vector< int > cell_mask);
    !> setprintchemistryon(bool workers, bool initial_phreeqc, bool utility);
    !> setrebalancebycell(bool tf);
    !> setrebalancefraction(double f);
    !> setrepresentativevolume(std::vector< double > rv);
    !> setsaturationuser(std::vector< double > sat);
    !> setscreenon(bool tf);
    !> setselectedoutputon(bool tf);
    !> setspeciessaveon(bool save_on);
    !> settemperature(std::vector< double > t);
    !> settime(double time);
    !> settimeconversion(double conv_factor);
    !> settimestep(double time_step);
    !> setunitsexchange(int option);
    !> setunitsgasphase(int option);
    !> setunitskinetics(int option);
    !> setunitsppassemblage(int option);
    !> setunitssolution(int option);
    !> setunitsssassemblage(int option);
    !> setunitssurface(int option);
    !> speciesconcentrations2module(std::vector< double > species_conc);
    !> statesave(int istate);
    !> stateapply(int istate);
    !> statedelete(int istate);
    !> usesolutiondensityvolume(bool tf);
    !> warningmessage(std::string warnstr);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !>
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> id = brm%bmif_create(nxyz, nthreads)
    !> status = brm%bmif_initializeYAML("myfile.yaml")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.

#ifdef USE_YAML
    INTEGER FUNCTION bmif_initialize(self, config_file)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_Initialize(id, config_file) &
        BIND(C, NAME='RMF_BMI_Initialize')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: config_file(*)
    END FUNCTION RMF_BMI_Initialize
    END INTERFACE
    class(bmi), intent(inout) :: self
    integer :: return_value
    CHARACTER(len=*), INTENT(in) :: config_file   
    if (self%bmiphreeqcrm_id .lt. 0) then
        STOP "ERROR: You must create the bmif instance before you call bmif_initialize."
    endif
#if defined(USE_MPI)
    if (RM_getmpimyself(self%bmiphreeqcrm_id) .gt. 0) then
        STOP "bmif_initialize with YAML can only be called by root MPI process."
    endif
#endif    
!    if (self%bmiphreeqcrm_id .lt. 0) then
!        self%bmiphreeqcrm_id = bmif_create()
!    endif
    bmif_initialize = success(RMF_BMI_Initialize(self%bmiphreeqcrm_id, trim(config_file)//C_NULL_CHAR))
    return
    END FUNCTION bmif_initialize
#endif
 
    
    INTEGER FUNCTION bmif_get_id(self)
    class(bmi), intent(inout) :: self
    bmif_get_id = self%bmiphreeqcrm_id
    return 
    END FUNCTION bmif_get_id    
    
    !> @a bmif_update runs a reaction step for all of the cells in the reaction module.
    !> Same as @ref RunCells.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval     0 is success, 1 is failure.
    !> Tranport concentrations are transferred to the reaction cells
    !> (@ref bmif_set_value "Concentrations" before
    !> reaction calculations are run. The length of time over which kinetic
    !> reactions are integrated is set
    !> by @ref bmif_set_value "TimeStep". Other properties that may need to be updated
    !> as a result of the transport
    !> calculations include porosity (@ref bmif_set_value "Porosity"),
    !> pressure (@ref bmif_set_value "Pressure"),
    !> saturation (@ref bmif_set_value "SaturationUser"),
    !> temperature (@ref bmif_set_value "Temperature").
    !> @see
    !> @ref bmif_get_value,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_set_value("Porosity", por)                ! If pore volume changes
    !> status = brm%bmif_set_value("SaturationUser", sat)              ! If saturation changes
    !> status = brm%bmif_set_value("Temperature", temperature)     ! If temperature changes
    !> status = brm%bmif_set_value("Pressure", pressure)           ! If pressure changes
    !> status = brm%bmif_set_value("Concentrations", c)            ! Transported concentrations
    !> status = brm%bmif_set_value("TimeStep", time_step)          ! Time step for kinetic reactions
    !> status = brm%bmif_update()
    !> status = brm%bmif_get_value("Concentrations", c)            ! Concentrations after reaction
    !> status = brm%bmif_get_value("DensityCalculated", density)   ! Density after reaction
    !> status = brm%bmif_get_value("SolutionVolume", volume)       ! Solution volume after reaction
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.
    INTEGER FUNCTION bmif_update(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_Update(id) &
        BIND(C, NAME='RMF_BMI_Update')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_Update
    END INTERFACE
    class(bmi), intent(inout) :: self
    bmif_update = success(RMF_BMI_Update(self%bmiphreeqcrm_id))
    return
    END FUNCTION bmif_update
     
    !> @a bmif_update_until is the same as @ref bmif_update, except the time step is calculated
    !> from the argument @a end_time. The time step is calculated to be @a end_time minus 
    !> the current time (@ref bmif_get_current_time).
    !> @param self Fortran-supplied BMIPhreeqcRM instance..
    !> @param end_time Time at the end of the time step. 
    !> @see
    !> @ref bmif_initialize,
    !> @ref bmif_update.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_set_value("Time", time)
    !> status = brm%bmif_set_value("Concentrations", c)
    !> status = brm%bmif_update_until(time + 86400.0)
    !> status = brm%bmif_get_value("Concentrations", c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.    

    INTEGER FUNCTION bmif_update_until(self, end_time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_UpdateUntil(id, time) &
        BIND(C, NAME='RMF_BMI_UpdateUntil')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    real(kind=C_DOUBLE), INTENT(in) :: time
    END FUNCTION RMF_BMI_UpdateUntil
    END INTERFACE
    class(bmi), intent(in) :: self
    real(kind=8), INTENT(in) :: end_time
    bmif_update_until = success(RMF_BMI_UpdateUntil(self%bmiphreeqcrm_id, end_time))
    return
    END FUNCTION bmif_update_until

    !> @a bmif_finalize destroys a reaction module, same as @ref Destroy.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval    0 is success, 1 is failure.
    !> @see
    !> @ref bmif_create.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_finalize()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.
    INTEGER FUNCTION bmif_finalize(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_destroy(id) &
        BIND(C, NAME='RMF_BMI_destroy')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_destroy
    END INTERFACE
    class(bmi), intent(inout) :: self
    INTEGER :: status
#ifdef USE_MPI
    !status = brm%RM_mpiworkerbreak()
#endif    
    bmif_finalize = success(RMF_BMI_destroy(self%bmiphreeqcrm_id))
    return
    END FUNCTION bmif_finalize
    
    ! ====================================================
    ! Exchange items
    ! ====================================================

    !> @a bmif_get_component_name returns the component name--BMIPhreeqcRM. 
    !>
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param component_name Filled with "BMIPhreeqcRM", the name of the component.
    !> @retval               0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_component_name(component_name)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_component_name(self, component_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(inout) :: component_name
    component_name = "BMI_PhreeqcRM"
    bmif_get_component_name = BMI_SUCCESS
    return
    END FUNCTION bmif_get_component_name

    !> @a bmif_get_input_item_count returns count of variables that
    !> can be set with @ref bmif_set_value.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param count  Number of input variables that can be set with @ref bmif_set_value.
    !> @retval       0 is success, 1 is failure.
    !> 
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_input_item_count(count)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    !>
    INTEGER FUNCTION bmif_get_input_item_count(self, count)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetInputItemCount(id) &
        BIND(C, NAME='RMF_BMI_GetInputItemCount')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_GetInputItemCount
    END INTERFACE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetInputItemCount(self%bmiphreeqcrm_id)
    bmif_get_input_item_count = success(count)
    END FUNCTION bmif_get_input_item_count

    !> @a bmif_get_output_item_count returns count of output variables that can be
    !> retrieved with @ref bmif_get_value.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param count  Number of output variables that can be retrieved with @ref bmif_get_value.
    !> @retval       0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_value,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_output_item_count(outputvarcount)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_output_item_count(self, count)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetOutputItemCount(id) &
        BIND(C, NAME='RMF_BMI_GetOutputItemCount')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_GetOutputItemCount
    END INTERFACE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetOutputItemCount(self%bmiphreeqcrm_id)
    bmif_get_output_item_count = success(count)
    END FUNCTION bmif_get_output_item_count
    
    !> @a bmif_get_pointable_item_count returns count of variables for which
    !> pointers can be retrieved with @ref bmif_get_value_ptr.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param count  Number of variables for which pointers can be retrieved with 
    !> @ref bmif_get_value_ptr.
    !> @retval       0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_pointable_item_count(pointablevarcount)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    !>
    INTEGER FUNCTION bmif_get_pointable_item_count(self, count)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetPointableItemCount(id) &
            BIND(C, NAME='RMF_BMI_GetPointableItemCount')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetPointableItemCount
        END INTERFACE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetPointableItemCount(self%bmiphreeqcrm_id)
    bmif_get_pointable_item_count = success(count)
    END FUNCTION bmif_get_pointable_item_count

    !> Basic Model Interface method that returns a list of the variable names that can be set
    !> with @ref bmif_set_value.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var_names   Character array of variable names.
    !> @retval            0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_input_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), dimension(:), allocatable :: inputvars
    !> status = brm%bmif_get_input_var_names(inputvars)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_input_var_names(self, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNames(id, type, dest) &
            BIND(C, NAME='RMF_BMI_GetNames')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: type
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetNames
        END INTERFACE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNamesSize(id, type, dest) &
            BIND(C, NAME='RMF_BMI_GetNamesSize')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: type
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetNamesSize
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "inputvarnames"
    status = RMF_BMI_GetNamesSize(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_input_item_count(self, dim)
    nbytes = dim * itemsize
    dim = nbytes / itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_input_var_names = RMF_BMI_GetNames(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, var_names(1))
    return
    END FUNCTION bmif_get_input_var_names

    !> @a bmif_get_output_var_names returns a list of the variable names that can be
    !> retrieved with @ref bmif_get_value.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var_names   Character array of variable names.
    !> @retval            0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable, dimension(:) :: var_names
    !> status = brm%bmif_get_output_var_names(var_names)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_output_var_names(self, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNames(id, type, dest) &
            BIND(C, NAME='RMF_BMI_GetNames')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: type
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetNames
        END INTERFACE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNamesSize(id, type, dest) &
            BIND(C, NAME='RMF_BMI_GetNamesSize')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: type
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetNamesSize
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "outputvarnames"
    status = RMF_BMI_GetNamesSize(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_output_item_count(self, dim)
    nbytes = dim * itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_output_var_names = RMF_BMI_GetNames(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, var_names(1))
    return
    END FUNCTION bmif_get_output_var_names

    !> @a bmif_get_pointable_var_names returns a list of the variable names for which a pointer can be
    !> retrieved with @ref bmif_get_value_ptr.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var_names   Character array of variable names.
    !> @retval            0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_pointable_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable, dimension(:) :: var_names
    !> status = brm%bmif_get_pointable_var_names(var_names)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_pointable_var_names(self, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNames(id, type, dest) &
        BIND(C, NAME='RMF_BMI_GetNames')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: type
    CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
    END FUNCTION RMF_BMI_GetNames
    END INTERFACE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetNamesSize(id, type, dest) &
        BIND(C, NAME='RMF_BMI_GetNamesSize')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: type
    INTEGER(KIND=C_INT), INTENT(inout) :: dest
    END FUNCTION RMF_BMI_GetNamesSize
    END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "pointablevarnames"
    status = RMF_BMI_GetNamesSize(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_pointable_item_count(self, dim)
    nbytes = dim * itemsize
    dim = nbytes / itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_pointable_var_names = RMF_BMI_GetNames(self%bmiphreeqcrm_id, trim(vartype)//C_NULL_CHAR, var_names(1))
    return
    END FUNCTION bmif_get_pointable_var_names

    ! ====================================================
    ! Variable information
    ! ====================================================

    !> @a bmif_get_var_grid returns a value of 1, indicating points.
    !> BMIPhreeqcRM does not have a grid of its own. The cells
    !> of BMIPhreeqcRM are associated with the user's model grid,
    !> and all spatial characterists are assigned by the user's
    !> model.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var   Varaiable name. (Return value is the same regardless of @a name.)
    !> @param grid  1 (points). BMIPhreeqcRM cells derive meaning from the user's model. 
    !> @retval      0 is success, 1 is failure.
    integer FUNCTION bmif_get_var_grid(self, var, grid)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    INTEGER, INTENT(out) :: grid
    grid = 1
    bmif_get_var_grid = BMI_SUCCESS
    END FUNCTION bmif_get_var_grid

    !> @a bmif_get_var_type retrieves the type of a variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> Types are "character", "real(kind=8)", "integer", or "logical",
    !> or an allocatable array of these types.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var   Name of the variable to retrieve the type.
    !> @param vtype Type of the variable.
    !> @retval      0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_input_item_count,
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_pointable_var_names,
    !> @ref bmif_get_pointable_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_value_ptr,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do i = 1, size(inputvars)
    !>     write(*,"(1x, I4, A40)") i, trim(inputvars(i))
    !>     status = brm%bmif_get_var_units(inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = brm%bmif_get_var_type(inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = brm%bmif_get_var_itemsize(inputvars(i), itemsize)
    !>     write(*, "(5x, I15)") itemsize
    !>     status = brm%bmif_get_var_nbytes(inputvars(i), nbytes)
    !>     write(*, "(5x, I15)") nbytes
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_type(self, var, vtype)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarType(id, var, vtype, l) &
        BIND(C, NAME='RMF_BMI_GetVarType')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id, l
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    CHARACTER(KIND=C_CHAR), INTENT(inout) :: vtype(*)
    END FUNCTION RMF_BMI_GetVarType
    END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: vtype
    bmif_get_var_type = success(RMF_BMI_GetVarType(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, vtype, len(vtype)))
    return
    END FUNCTION bmif_get_var_type
    
    !> @a bmif_get_var_units retrieves the units of a
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var     Name of the variable to retrieve units.
    !> @param units   Units of the variable.
    !> @retval        0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_input_item_count,
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_pointable_var_names,
    !> @ref bmif_get_pointable_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_value_ptr,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do i = 1, size(inputvars)
    !>     write(*,"(1x, I4, A40)") i, trim(inputvars(i))
    !>     status = brm%bmif_get_var_units(inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = brm%bmif_get_var_type(inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = brm%bmif_get_var_itemsize(inputvars(i), itemsize)
    !>     write(*, "(5x, I15)") itemsize
    !>     status = brm%bmif_get_var_nbytes(inputvars(i), nbytes)
    !>     write(*, "(5x, I15)") nbytes
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_units(self, var, units)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarUnits(id, var, units, l) &
        BIND(C, NAME='RMF_BMI_GetVarUnits')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id, l
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    CHARACTER(KIND=C_CHAR), INTENT(inout) :: units(*)
    END FUNCTION RMF_BMI_GetVarUnits
    END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: units
    bmif_get_var_units = success(RMF_BMI_GetVarUnits(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, units, len(units)))
    return
    END FUNCTION bmif_get_var_units

    !> @a bmif_get_var_itemsize retrieves the size, in bytes, of a
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> Sizes may be the size of an integer, real(kind=8), 
    !> or a character length for string variables.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var        Name of the variable to retrieve the item size.
    !> @param itemsize   Size, in bytes, of one element of the variable.
    !> @retval           0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_input_item_count,
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_pointable_var_names,
    !> @ref bmif_get_pointable_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_value_ptr,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer nbytes, item_size, dim
    !> real(kind=8), allocatable, dimension(:) :: bmi_temperature
    !> status = brm%bmif_get_var_nbytes("Temperature", nbytes)
    !> status = brm%bmif_get_var_itemsize("Temperature", item_size)
    !> dim    = nbytes/item_size
    !> allocate(bmi_temperature(dim))
    !> bmi_temperature = 25.0
    !> status = brm%bmif_set_value("Temperature", bmi_temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_itemsize(self, var, itemsize)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarItemsize(id, var) &
        BIND(C, NAME='RMF_BMI_GetVarItemsize')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    END FUNCTION RMF_BMI_GetVarItemsize
    END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    INTEGER, INTENT(out) :: itemsize
    itemsize = RMF_BMI_GetVarItemsize(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR)
    bmif_get_var_itemsize = success(itemsize)
    END FUNCTION bmif_get_var_itemsize

    !> @a bmif_get_var_nbytes retrieves the total number of bytes needed for a 
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var     Name of the variable to retrieve the number of bytes needed to
    !> retrieve or store the variable.
    !> @param nbytes  Total number of bytes needed for the variable.
    !> @retval        0 is success, 1 is failure.
    !>
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_input_item_count,
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_pointable_var_names,
    !> @ref bmif_get_pointable_item_count,
    !> @ref bmif_get_value,
    !> @ref bmif_get_value_ptr,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer nbytes, item_size, dim
    !> real(kind=8), allocatable, dimension(:) :: bmi_temperature
    !> status = brm%bmif_get_var_nbytes("Temperature", nbytes)
    !> status = brm%bmif_get_var_itemsize("Temperature", item_size)
    !> dim    = nbytes/item_size
    !> allocate(bmi_temperature(dim))
    !> bmi_temperature = 25.0
    !> status = brm%bmif_set_value("Temperature", bmi_temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_nbytes(self, var, nbytes)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarNbytes(id, var) &
        BIND(C, NAME='RMF_BMI_GetVarNbytes')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    END FUNCTION RMF_BMI_GetVarNbytes
    END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    INTEGER, INTENT(out) :: nbytes
    nbytes  = RMF_BMI_GetVarNbytes(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR)
    bmif_get_var_nbytes = success(nbytes)
    END FUNCTION bmif_get_var_nbytes

    ! ====================================================
    ! Time information
    ! ====================================================

    !> @a bmif_get_current_time returns the current simulation time, in seconds. (Same as @ref GetTime.)
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param time    The current simulation time, in seconds.
    !> @retval        0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_end_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref GetTime,
    !> @ref GetTimeStep,
    !> @ref SetTime,
    !> @ref SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_current_time(time)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_current_time(self, time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetCurrentTime(id) &
        BIND(C, NAME='RMF_BMI_GetCurrentTime')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_GetCurrentTime
    END INTERFACE
    class(bmi), intent(inout) :: self
    real(kind=8), intent(inout) :: time
    time = RMF_BMI_GetCurrentTime(self%bmiphreeqcrm_id)
    bmif_get_current_time = BMI_SUCCESS
    END FUNCTION bmif_get_current_time

    !> @a bmif_get_start_time returns the current simulation time, in seconds. 
    !> (Same as @ref bmif_get_current_time and @ref GetTime.)
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param start_time  The current simulation time, in seconds.
    !> @retval        0 is success, 1 is failure.
    INTEGER FUNCTION bmif_get_start_time(self, start_time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout) :: start_time
    bmif_get_start_time = bmif_get_current_time(self, start_time)
    END FUNCTION bmif_get_start_time

    !> @a bmif_get_end_time returns @ref bmif_get_current_time plus
    !> @ref bmif_get_time_step, in seconds.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param end_time  The end of the time step, in seconds.
    !> @retval          0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref GetTime,
    !> @ref GetTimeStep,
    !> @ref SetTime,
    !> @ref SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_end_time(end_time)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_end_time(self, end_time)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetEndTime(id) &
            BIND(C, NAME='RMF_BMI_GetEndTime')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetEndTime
        END INTERFACE
    class(bmi), intent(inout) :: self
    real(kind=8), intent(inout) :: end_time
    end_time = RMF_BMI_GetEndTime(self%bmiphreeqcrm_id)
    bmif_get_end_time = BMI_SUCCESS
    END FUNCTION bmif_get_end_time

    !> @a bmif_get_time_units returns the time units of PhreeqcRM.
    !> All time units are seconds for PhreeqcRM.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param time_units    Returns the string "seconds".
    !> @retval              0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_end_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref GetTime,
    !> @ref GetTimeStep,
    !> @ref SetTime,
    !> @ref SetTimeStep,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(20) time_units
    !> status = brm%bmif_get_time_units(time_units)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_time_units(self, time_units)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_gettimeUnits(id, time_units, l) &
            BIND(C, NAME='RMF_BMI_gettimeUnits')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id, l
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: time_units(*)
        END FUNCTION RMF_BMI_gettimeUnits
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(inout) :: time_units
    bmif_get_time_units = success(RMF_BMI_gettimeUnits(self%bmiphreeqcrm_id, time_units, len(time_units)))
    return
    END FUNCTION bmif_get_time_units


    !> @a bmif_get_time_step returns the current simulation time step,
    !> in seconds. (Same as @ref GetTimeStep.)
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param time_step     The current simulation time step, in seconds.
    !> @retval              0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_end_time,
    !> @ref bmif_set_value,
    !> @ref GetTime,
    !> @ref GetTimeStep,
    !> @ref SetTime,
    !> @ref SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_get_time_step(time_step)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_time_step(self, time_step)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_gettimestep(id) &
            BIND(C, NAME='RMF_BMI_gettimestep')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_gettimestep
        END INTERFACE
    class(bmi), intent(inout) :: self
    real(kind=8), intent(inout) :: time_step
    time_step = RMF_BMI_gettimestep(self%bmiphreeqcrm_id)
    bmif_get_time_step = BMI_SUCCESS
    END FUNCTION bmif_get_time_step
        
    ! ====================================================
    ! Getters, by type
    ! ====================================================

    !> @a bmif_get_value retrieves model variables. Only variables in the list
    !> provided by @ref bmif_get_output_var_names can be retrieved. 
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var    Name of the variable to retrieve.
    !> @param dest   Variable in which to place results.
    !> @retval       0 is success, 1 is failure.
    !>
    !> Variable names for the second argument (@a name) and variable type of the
    !> third argument (@a dest).
    !> @n "ComponentCount", @a dest: integer;
    !> @n "Components", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "Concentrations", @a dest: real(kind=8), allocatable, dimension(:,:);
    !> @n "CurrentSelectedOutputUserNumber", @a dest: integer;
    !> @n "DensityCalculated", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "ErrorString", @a dest: character;
    !> @n "FilePrefix", @a dest: character;
    !> @n "Gfw", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "GridCellCount", @a dest: integer;
    !> @n "InputVarNames", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "OutputVarNames", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "Porosity", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "Pressure", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "SaturationCalculated", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "SelectedOutput", @a dest: real(kind=8), allocatable, dimension(:,:);
    !> @n "SelectedOutputColumnCount", @a dest: integer;
    !> @n "SelectedOutputCount", @a dest: integer;
    !> @n "SelectedOutputHeadings", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "SelectedOutputOn", @a dest: logical;
    !> @n "SelectedOutputRowCount", @a dest: integer;
    !> @n "SolutionVolume", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "Temperature", @a dest: real(kind=8), allocatable, dimension(:);
    !> @n "Time", @a dest: real(kind=8);
    !> @n "TimeStep", @a dest: real(kind=8);
    !> @n "Viscosity", @a dest: real(kind=8), allocatable, dimension(:).
    !>
    !> @see
    !> @ref bmif_get_output_var_names,
    !> @ref bmif_get_output_item_count,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable, dimension(:) :: bmi_density
    !> character(len=:), allocatable, dimension(:) :: bmi_comps
    !> status = brm%bmif_get_value("DensityCalculated", bmi_density)
    !> status = brm%bmif_get_value("Components", bmi_comps)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.
!Add NEW_VARIABLE to bmif_get_value Documentation
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_logical(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL, INTENT(inout) :: dest
    character(100) :: vartype
    integer :: status, dest_int
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "logical") then
        write(*,*) vartype, " logical"
        stop "Variable type error."
    endif
    bmif_get_value_logical = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest_int)
    dest = .true.
    if (dest_int .eq. 0) dest = .false.
    return
    END FUNCTION bmif_get_value_logical

    !> \overload
    INTEGER FUNCTION bmif_get_value_char(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: dest
    character(BMI_MAX_TYPE_NAME) :: vartype
    integer :: itemsize, status
    CHARACTER(len=:), allocatable :: temp
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "character") then
        write(*,*) vartype, " character"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(self, var, itemsize)
    allocate(character(len=itemsize) :: temp)
    status = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, temp)
    if (len(dest) .gt. 0) then
        dest = temp
    else
        status = RM_errormessage(self%bmiphreeqcrm_id, "Variable length is zero")   
        status = -1
    endif
    bmif_get_value_char = success(status)
    return
    END FUNCTION bmif_get_value_char
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_char_alloc(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=:), allocatable, INTENT(inout) :: dest
    character(BMI_MAX_TYPE_NAME) :: vartype
    integer :: itemsize, status
    CHARACTER(len=:), allocatable :: temp
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "character") then
        write(*,*) vartype, " character"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(self, var, itemsize)
    allocate(character(len=itemsize) :: temp)
    status = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, temp)
    dest = temp
    bmif_get_value_char_alloc = success(status)
    return
    END FUNCTION bmif_get_value_char_alloc

    !> \overload
    INTEGER FUNCTION bmif_get_value_char1(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "character(len=:),allocatable,dimension(:)") then
        write(*,*) vartype, " character(len=:),allocatable,dimension(:)"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(self, var, itemsize)
    status = bmif_get_var_nbytes(self, var, nbytes)
    dim = nbytes / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
        if (dim2 > 0) then
            dim1 = len(dest(1))
        endif
    endif
    if (dim1 .ne. itemsize .or. dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(character(len=itemsize) :: dest(dim))
    endif
    bmif_get_value_char1 = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_char1

    !> \overload
    INTEGER FUNCTION bmif_get_value_double(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8"
        stop "Variable type error."
    endif
    bmif_get_value_double = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest)
    return
    END FUNCTION bmif_get_value_double

    !> \overload
    INTEGER FUNCTION bmif_get_value_double1(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), allocatable, dimension(:), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8 1d"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(self, var, itemsize)
    status = bmif_get_var_nbytes(self, var, nbytes)
    dim = nbytes / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
    endif
    bmif_get_value_double1 = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_double1

    !> \overload
    INTEGER FUNCTION bmif_get_value_double2(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), allocatable, INTENT(inout) :: dest(:,:)
    character(100) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8 2d"
        stop "Variable type error."
    endif
    varname = Lower(var)
    need_alloc = .true.
    if (varname .eq. "concentrations") then
        dim1 = RM_getgridcellcount(self%bmiphreeqcrm_id)
        dim2 = RM_getcomponentcount(self%bmiphreeqcrm_id)
    else if (varname .eq. "selectedoutput") then
        dim1 = RM_getselectedoutputrowcount(self%bmiphreeqcrm_id)
        dim2 = RM_getselectedoutputcolumncount(self%bmiphreeqcrm_id)
    else
        stop "Unknown 2d variable"
    endif
    if (allocated(dest)) then
        if ((size(dest,1) .eq. dim1) .and. &
            (size(dest,2) .eq. dim2)) then
            need_alloc = .false.
        else
            deallocate(dest)
        endif
    endif
    if (need_alloc) then
        allocate(dest(dim1, dim2))
    endif
    bmif_get_value_double2 = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION bmif_get_value_double2

    !> \overload
    INTEGER FUNCTION bmif_get_value_int(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer") then
        write(*,*) vartype, " integer"
        stop "Variable type error."
    endif
    bmif_get_value_int = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest)
    return
    END FUNCTION bmif_get_value_int

    !> \overload
    INTEGER FUNCTION bmif_get_value_int1(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, allocatable, INTENT(inout) :: dest(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:)") then
        write(*,*) vartype, " integer 1d"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(self, var, itemsize)
    status = bmif_get_var_nbytes(self, var, nbytes)  
    dim = nbytes / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
    endif
    bmif_get_value_int1 = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_int1

    !> \overload
    INTEGER FUNCTION bmif_get_value_int2(self, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, allocatable, INTENT(inout) :: dest(:,:)
    character(100) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer") then
        write(*,*) vartype, " integer 2d"
        stop "Variable type error."
    endif
    varname = Lower(varname)
    need_alloc = .true.
    stop "Unknown 2d variable"
    bmif_get_value_int2 = RMF_BMI_GetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION bmif_get_value_int2
   
    !> @a bmif_get_value_ptr retrieves pointers to model variables. Only variables in the list
    !> provided by @ref bmif_get_pointable_var_names can be pointed to. 
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var    Name of the variable to retrieve.
    !> @param ptr    Pointer to the variable's data.
    !> @retval       0 is success, 1 is failure.
    !> The following list gives the name in the second argument (@a var) and the
    !> data type the pointer (@a ptr):
    !> @n "ComponentCount": integer;
    !> @n "Concentrations": real(kind=8) (:);
    !> @n "DensityCalculated": real(kind=8) (:);
    !> @n "Gfw": real(kind=8) (:);
    !> @n "GridCellCount": integer;
    !> @n "Porosity": real(kind=8) (:);
    !> @n "Pressure": real(kind=8) (:);
    !> @n "SaturationCalculated": real(kind=8) (:);
    !> @n "SelectedOutputOn": logical(kind=1);
    !> @n "SolutionVolume": real(kind=8) (:);
    !> @n "Temperature": real(kind=8) (:);
    !> @n "Time": real(kind=8);
    !> @n "TimeStep": real(kind=8);
    !> @n "Viscosity": real(kind=8) (:);
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.     
!Add NEW_VARIABLE to bmif_get_value_ptr Documentation  
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_double(self, var, ptr)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValuePtr(id, var, src) &
            BIND(C, NAME='RMF_BMI_GetValuePtr')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        type (c_ptr), INTENT(inout) :: src
        END FUNCTION RMF_BMI_GetValuePtr
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=c_double), pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer :: status
    status = RMF_BMI_GetValuePtr(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src)
    call C_F_POINTER(src, ptr)
    bmif_get_value_ptr_double = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_double
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_double1(self, var, ptr)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValuePtr(id, var, src) &
            BIND(C, NAME='RMF_BMI_GetValuePtr')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        type (c_ptr), INTENT(inout) :: src
        END FUNCTION RMF_BMI_GetValuePtr
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=c_double), pointer, INTENT(inout) :: ptr(:)
    type (c_ptr) :: src
    integer nbytes, itemsize, dim, status
    status = bmif_get_var_nbytes(self, var, nbytes)
    status = bmif_get_var_itemsize(self, var, itemsize)
    dim = nbytes/itemsize
    status = RMF_BMI_GetValuePtr(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr, [dim])
    bmif_get_value_ptr_double1 = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_double1
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_integer(self, var, ptr)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValuePtr(id, var, src) &
            BIND(C, NAME='RMF_BMI_GetValuePtr')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        type (c_ptr), INTENT(inout) :: src
        END FUNCTION RMF_BMI_GetValuePtr
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer status
    status = RMF_BMI_GetValuePtr(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr)
    bmif_get_value_ptr_integer = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_integer

    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_logical(self, var, ptr)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValuePtr(id, var, src) &
            BIND(C, NAME='RMF_BMI_GetValuePtr')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        type (c_ptr), INTENT(inout) :: src
        END FUNCTION RMF_BMI_GetValuePtr
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    logical(kind=1), pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer status
    status = RMF_BMI_GetValuePtr(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr)
    bmif_get_value_ptr_logical = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_logical    

    ! ====================================================
    ! Setters, by type
    ! ====================================================

    !> @a bmif_set_value sets model variables. Only variables in the list
    !> provided by @ref bmif_get_input_var_names can be set. 
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param var    Name of variable to set.
    !> @param src    Data to use to set the variable.
    !> @retval       0 is success, 1 is failure.
    !>
    !> Variable names for the second argument (@a var)
    !> and required variable type for the third argument (@a src):
    !> @n "Concentrations", real(kind=8), allocatable, dimension(:,:);
    !> @n "DensityUser", real(kind=8), allocatable, dimension(:);
    !> @n "FilePrefix", character;
    !> @n "NthSelectedOutput", integer;
    !> @n "Porosity", real(kind=8), allocatable, dimension(:);
    !> @n "Pressure", real(kind=8), allocatable, dimension(:);
    !> @n "SaturationUser", real(kind=8), allocatable, dimension(:);
    !> @n "SelectedOutputOn", logical;
    !> @n "Temperature", real(kind=8), allocatable, dimension(:);
    !> @n "Time", real(kind=8);
    !> @n "TimeStep", real(kind=8).
    !>
    !> @see
    !> @ref bmif_get_input_var_names,
    !> @ref bmif_get_input_item_count,,
    !> @ref bmif_get_value,
    !> @ref bmif_get_var_itemsize,
    !> @ref bmif_get_var_nbytes,
    !> @ref bmif_get_var_type,
    !> @ref bmif_get_var_units.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable, dimension(:) :: tc
    !> allocate(tc(nxyz))
    !> tc = 28.0d0
    !> status = brm%bmif_set_value("Temperature", tc)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.
    !> \overload
    INTEGER FUNCTION bmif_set_value_b(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL, INTENT(in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, src_int
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "logical") then
        stop "Variable type error."
    endif
    src_int = 1
    if(.not. src) src_int = 0
    bmif_set_value_b = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src_int)
    return
    END FUNCTION bmif_set_value_b

    !> \overload
    INTEGER FUNCTION bmif_set_value_c(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        CHARACTER(KIND=C_CHAR), INTENT(in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "character") then
        stop "Variable type error."
    endif
    bmif_set_value_c = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, trim(src)//C_NULL_CHAR)
    return
    END FUNCTION bmif_set_value_c

    !> \overload
    INTEGER FUNCTION bmif_set_value_int(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src
    integer :: src_copy
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer") then
        stop "Variable type error."
    endif
    src_copy = src
    if (var .eq. "NthSelectedOutput") src_copy = src - 1
    bmif_set_value_int = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src_copy)
    return
    END FUNCTION bmif_set_value_int

    !> \overload
    INTEGER FUNCTION bmif_set_value_int1(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT( in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(self, var, nbytes)
    status =  bmif_get_var_itemsize(self, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_int1 = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src(1))
    return
    END FUNCTION bmif_set_value_int1

    !> \overload
    INTEGER FUNCTION bmif_set_value_int2(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src(:,:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:,:)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(self, var, nbytes)
    status = bmif_get_var_itemsize(self, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src,1)*size(src,2)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_int2 = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src(1,1))
    return
    END FUNCTION bmif_set_value_int2

    !> \overload
    INTEGER FUNCTION bmif_set_value_double(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT( in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    bmif_set_value_double = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src)
    return
    END FUNCTION bmif_set_value_double

    !> \overload
    INTEGER FUNCTION bmif_set_value_double1(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT( in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(self, var, nbytes)
    status = bmif_get_var_itemsize(self, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_double1 = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src(1))
    return
    END FUNCTION bmif_set_value_double1

    !> \overload
    INTEGER FUNCTION bmif_set_value_double2(self, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, src) &
            BIND(C, NAME='RMF_BMI_SetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT( in) :: src
        END FUNCTION RMF_BMI_SetValue
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src(:,:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(self, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(self, var, nbytes)
    status = bmif_get_var_itemsize(self, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src,1)*size(src,2)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_double2 = RMF_BMI_SetValue(self%bmiphreeqcrm_id, trim(var)//C_NULL_CHAR, src(1,1))
    return
    END FUNCTION bmif_set_value_double2

    ! ====================================================
    ! Grid information
    ! ====================================================

    !> @a bmif_get_grid_rank returns a rank of 1 for grid 0. 
    !> BMIPhreeqcRM has a 1D series of
    !> cells; any grid or spatial information must
    !> be found in the user's model.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param grid   Grid number, only grid 0 is considered.
    !> @param rank   Rank of 1 is returned for grid 0; 0 for
    !> all other values of @a grid.
    !> @retval       0 is success, 1 is failure.    
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_grid_rank(grid, rank)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_rank(self, grid, rank)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    integer, intent(in) :: grid
    INTEGER, INTENT(inout) :: rank
    if (grid .eq. 0) then
        rank = 1
        bmif_grid_rank = BMI_SUCCESS
    else
        rank = 0
        bmif_grid_rank = BMI_FAILURE
    endif
    END FUNCTION bmif_grid_rank

    !> @ref bmif_grid_size returns the number of cells specified
    !> at creation of the BMIPhreeqcRM instance. 
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param grid  Grid number, only grid 0 is considered.
    !> @param ngrid Same value as @ref GetGridCellCount 
    !> or @ref bmif_get_value "GridCellCount" is returned for grid 0; 
    !> 0 for all other values of @a grid.
    !> @retval       0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_grid_size(grid, ngrid)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_size(self, grid, ngrid)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    integer, intent(in) :: grid
    INTEGER, INTENT(inout) :: ngrid
    if (grid .eq. 0) then
        bmif_grid_size = success(self%bmif_get_value("GridCellCount", ngrid))
    else
        ngrid = 0
        bmif_grid_size = BMI_FAILURE
    endif
    END FUNCTION bmif_grid_size

    !> @a bmif_grid_type defines the grid to be points. No grid
    !> information is available in BMIPhreeqcRM; all grid 
    !> information must be found in the user's model.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param grid  Grid number, only grid 0 is considered.
    !> @param str   "Points" is returned for grid 0;
    !> "Undefined grid identifier" is returned for all other 
    !> values of @a grid.
    !> @retval      0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%bmif_grid_type(grid, str)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_type(self, grid, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    integer, intent(in) :: grid
    CHARACTER(len=*), INTENT(inout) :: str
    if (grid .eq. 1) then
        str = "points"
        bmif_grid_type = BMI_SUCCESS
    else
        str = ""
        bmif_grid_type = BMI_FAILURE
    endif
    END FUNCTION bmif_grid_type

!> @a bmif_add_output_vars allows selection of sets of variables that can be retieved
!> by the @ref bmif_get_value method. Sets of variables can be included or excluded with
!> multiple calls to this method. All calls must precede the final call to
!> the PhreeqcRM method findcomponents. findcomponents generates SELECTED_OUTPUT 333 and
!> USER_PUNCH 333 data blocks that make the variables accessible. Variables will
!> only be accessible if the system includes the given reactant; for example, no
!> gas variables will be created if there are no GAS_PHASEs in the model. 
!>
!> @param self Fortran-supplied BMIPhreeqcRM instance.
!> @param option A string value, among those listed below, that includes or
!> excludes variables from @ref bmif_get_output_var_names, @ref bmif_get_value,
!> and other BMI methods.
!> @param def A string value that can be "false", "true", or a list of items to be included as
!> accessible variables. A value of "false", excludes all variables of the given type; a 
!> value of "true" includes all variables of the given type for the current system; a list
!> specifies a subset of items of the given type. 
!>
!> Values for the the parameter @a option:
!> @n AddOutputVars: False excludes all variables; True causes the settings for each variable group
!> to determine the variables that will be defined. Default True;
!> @n SolutionProperties: False excludes all solution property variables; True includes variables pH, pe,
!> alkalinity, ionic strength, water mass, charge balance, percent error, and specific conductance.
!> Default True.
!> @n SolutionTotalMolalities: False excludes all total element and element redox state variables;
!> True includes all elements and element redox state variables for the system defined for the 
!> calculation; list restricts variables to the specified elements and redox states.
!> Default True.
!> @n ExchangeMolalities: False excludes all variables related to exchange; True includes all 
!> variables related to exchange; list includes variables for the specified exchange species.
!> Default True.
!> @n SurfaceMolalities: False excludes all variables related to surfaces; True includes all 
!> variables related to surfaces; list includes variables for the specified surface species.
!> Default True.
!> @n EquilibriumPhases: False excludes all variables related to equilibrium phases; True includes all 
!> variables related to equilibrium phases; list includes variables for the specified
!> equilibiurm phases. Default True.
!> @n Gases: False excludes all variables related to gases; True includes all 
!> variables related to gases; list includes variables for the specified gas components. Default True.
!> @n KineticReactants: False excludes all variables related to kinetic reactants; True includes all 
!> variables related to kinetic reactants; list includes variables for the specified kinetic 
!> reactants. Default True.
!> @n SolidSolutions: False excludes all variables related to solid solutions; True includes all 
!> variables related to solid solutions; list includes variables for the specified solid solutions
!> components. Default True.
!> @n CalculateValues: False excludes all calculate values; True includes all 
!> calculate values; list includes the specified calculate values. CALCLUATE_VALUES can be
!> used to calculate geochemical quantities not available in the other sets of variables. 
!> Default True.
!> @n SolutionActivities: False excludes all aqueous species; True includes all 
!> aqueous species; list includes only the specified aqueous species. Default False.
!> @n SolutionMolalities: False excludes all aqueous species; True includes all 
!> aqueous species; list includes only the specified aqueous species. Default False.
!> @n SaturationIndices: False excludes all saturation indices; True includes all 
!> saturation indices; list includes only the specified saturation indices. Default False.
    INTEGER FUNCTION bmif_add_output_vars(self, option, def)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_AddOutputVars(id, var, src) &
            BIND(C, NAME='RMF_BMI_AddOutputVars')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        CHARACTER(KIND=C_CHAR), INTENT(in) :: src
        END FUNCTION RMF_BMI_AddOutputVars
        END INTERFACE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: option
    CHARACTER(len=*), INTENT(in) :: def
    bmif_add_output_vars = RMF_BMI_AddOutputVars(self%bmiphreeqcrm_id, trim(option)//C_NULL_CHAR, trim(def)//C_NULL_CHAR)
    return
    END FUNCTION bmif_add_output_vars

    ! ====================================================
    ! Functions not implemented
    ! ====================================================

    INCLUDE "BMI_not_implemented.inc"

    ! ====================================================
    ! Utility functions
    ! ====================================================
        
    FUNCTION Lower(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i

    DO i = 1,LEN(s1)
        ch = s1(i:i)
        IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
        s2(i:i) = ch
    END DO
    END FUNCTION Lower

    INTEGER FUNCTION SUCCESS(i)
    implicit none
    integer, intent(in) :: i
    success = BMI_FAILURE
    if (i .ge. 0) success = BMI_SUCCESS
    END FUNCTION SUCCESS
    
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
#ifdef EXTEND_BMIPHREEQCRM
   
    !> abort the program.
    !> @a iresult will be interpreted as an IRESULT value and decoded; 
    !> @a err_str will be printed; and the reaction module will be destroyed.
    !> If using MPI, an MPI_abort message will be sent before the reaction
    !> module is destroyed. If the @a id is an invalid instance, abort will 
    !> return a value ofIBADINSTANCE, otherwise the program will exit with a 
    !> return code of 4.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param iresult    Integer treated as an IRESULT return code.
    !> @param err_str       String to be printed as an error message.
    !> @retval IRESULT   Program will exit before returning unless @a id is an 
    !> invalid reaction module id.
    !> @see
    !> @ref destroy,
    !> @ref errormessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> string = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1"
    !> status = brm%runstring(string)
    !> if (status .lt. 0) status = brm%abort(status, "runstring failed")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root or workers.

    INTEGER FUNCTION abort(self, iresult, err_str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: iresult
    CHARACTER(len=*), INTENT(in) :: err_str
    abort = RM_abort(self%bmiphreeqcrm_id, iresult, err_str)
    RETURN
    END FUNCTION abort

    !> Close the output and log files.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref openfiles,
    !> @ref setfileprefix.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%closefiles()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called only by root.

    INTEGER FUNCTION closefiles(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    closefiles = RM_closefiles(self%bmiphreeqcrm_id)
    RETURN
    END FUNCTION closefiles

    !> @a N sets of component concentrations are converted to SOLUTIONs numbered 
    !> 1-@a n in the Utility IPhreeqc. The solutions can be reacted and manipulated 
    !> with the methods of IPhreeqc. If solution concentration units (@ref setunitssolution) 
    !> are per liter, one liter of solution is created in the Utility instance; if solution
    !> concentration units are mass fraction, one kilogram of solution is created in the Utility 
    !> instance. The motivation for this method is the mixing of solutions in wells, where 
    !> it may be necessary to calculate solution properties (pH for example)
    !> or react the mixture to form scale minerals. The code fragments below make a mixture of
    !> concentrations and then calculate the pH of the mixture.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param c             Array of concentrations to be made SOLUTIONs in Utility IPhreeqc, 
    !> array size is (@a n, @a ncomps) where @a ncomps is the number of components 
    !> (@ref getcomponentcount).
    !> @param n             The number of sets of concentrations.
    !> @param tc            Array of temperatures to apply to the SOLUTIONs, in degree C. 
    !> Array of size @a n.
    !> @param p_atm         Array of pressures to apply to the SOLUTIONs, in atm. 
    !> Array of size n.
    !> @retval   id number of the utility IPhreeqc instance. 
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate (c_well(1,ncomps))
    !> do i = 1, ncomps
    !>   c_well(1,i) = 0.5 * c(1,i) + 0.5 * c(10,i)
    !> enddo
    !> allocate(tc(1), p_atm(1))
    !> tc(1) = 15.0
    !> p_atm(1) = 3.0
    !> iphreeqc_id = brm%concentrations2utility(c_well, 1, tc, p_atm)
    !> string = "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1"
    !> status = brm%runstring(iphreeqc_id, string)
    !> status = brm%setcurrentselectedoutputusernumber(iphreeqc_id, 5)
    !> status = brm%getselectedoutputValue(iphreeqc_id, 1, 1, vtype, pH, svalue)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called only by root.

    INTEGER FUNCTION concentrations2utility(self, c, n, tc, p_atm)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    real(kind=8), INTENT(in), DIMENSION(:,:) :: c
    INTEGER, INTENT(in) :: n
    real(kind=8), INTENT(in), DIMENSION(:) :: tc, p_atm
    concentrations2utility = RM_concentrations2utility(self%bmiphreeqcrm_id, c, n, tc, p_atm)
    return
    END FUNCTION concentrations2utility

    !> Provides a mapping from grid cells in the user's model to reaction cells in PhreeqcRM.
    !> The mapping is used to eliminate inactive cells and to use symmetry to decrease the 
    !> number of cells for which chemistry must be run. The array @a grid2chem of size 
    !> @a nxyz (the number of grid cells, @ref getgridcellcount) must contain the set of 
    !> all integers 0 <= @a i < @a count_chemistry, where @a count_chemistry is a number less 
    !> than or equal to @a nxyz. Inactive cells are assigned a negative integer.
    !> The mapping may be many-to-one to account for symmetry. Default is a one-to-one 
    !> mapping--all user grid cells are reaction cells (equivalent to @a grid2chem values 
    !> of 0,1,2,3,...,@a nxyz-1).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param grid2chem        An array of integers: Nonnegative is a reaction cell number (0 based), negative is an inactive cell. Array of size @a nxyz (number of grid cells).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> ! For demonstation, two equivalent rows by symmetry
    !> allocate(grid2chem(nxyz))
    !> do i = 1, nxyz/2
    !>   grid2chem(i) = i - 1
    !>   grid2chem(i+nxyz/2) = i - 1
    !> enddo
    !> status = brm%createmapping(grid2chem)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.

    INTEGER FUNCTION createmapping(self, grid2chem)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: grid2chem
    createmapping = RM_createmapping(self%bmiphreeqcrm_id, grid2chem)
    return
    END FUNCTION createmapping

    !> If @a e is negative, this method prints an error message corresponding to IRESULT @a e. 
    !> If @a e is non-negative, no action is taken.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param e                    An IRESULT value returned by one of the reaction-module methods.
    !> @retval IRESULT          0 is success, negative is failure (See @ref decodeerror).
    !> @par IRESULT definition:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> typedef enum {
    !>   IOK            =  0,  ! Success
    !>   IOUTOFMEMORY   = -1,  ! Failure, Out of memory
    !>   IBADVARTYPE    = -2,  ! Failure, Invalid VAR type
    !>   IINVALIDARG    = -3,  ! Failure, Invalid argument
    !>   IINVALIDROW    = -4,  ! Failure, Invalid row
    !>   IINVALIDCOL    = -5,  ! Failure, Invalid column
    !>   IBADINSTANCE   = -6,  ! Failure, Invalid rm instance id
    !>   IFAIL          = -7,  ! Failure, Unspecified
    !> } IRESULT;
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%createmapping(grid2chem)
    !> if (status < 0) status = brm%decodeerror(status)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Can be called by root and (or) workers.

    INTEGER FUNCTION decodeerror(self, e)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: e
    decodeerror = RM_decodeerror(self%bmiphreeqcrm_id, e)
    return
    END FUNCTION decodeerror

    !> destroys a reaction module.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref create.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%destroy()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.

    INTEGER FUNCTION destroy(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    destroy = RM_destroy(self%bmiphreeqcrm_id)
    return
    END FUNCTION destroy

    !> Writes the contents of all workers to file in _RAW formats, including SOLUTIONs and 
    !> all reactants.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param dump_on          Signal for writing the dump file: 1 true, 0 false.
    !> @param append           Signal to append to the contents of the dump file: 1 true, 0 false.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setdumpfilename.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>      
    !> dump_on = 1
    !> append = 0
    !> status = brm%setdumpfilename("advection_f90.dmp")
    !> status = brm%dumpmodule(dump_on, append)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root; workers must be in the loop of @ref mpiworker.

    INTEGER FUNCTION dumpmodule(self, dump_on, append)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: dump_on
    INTEGER, INTENT(in) :: append
    dumpmodule = RM_dumpmodule(self%bmiphreeqcrm_id, dump_on, append)
    return
    END FUNCTION dumpmodule

    !> Send an error message to the screen, the output file, and the log file.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param errstr           String to be printed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref logmessage,
    !> @ref openfiles,
    !> @ref outputmessage,
    !> @ref screenmessage,
    !> @ref warningmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%errormessage("Goodbye world")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers; root writes to output and log files.

    INTEGER FUNCTION errormessage(self, errstr)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: errstr
    errormessage = RM_errormessage(self%bmiphreeqcrm_id, errstr)
    return
    END FUNCTION errormessage

    !> Returns the number of items in the list of all elements in the InitialPhreeqc 
    !> instance. Elements are those that have been defined in a solution or any other 
    !> reactant (EQUILIBRIUM_PHASE, KINETICS, and others). The method can be called 
    !> multiple times and the list that is created is cummulative. The list is the set 
    !> of components that needs to be transported. By default the list includes water, 
    !> excess H and excess O (the H and O not contained in water); alternatively, the 
    !> list may be set to contain total H and total O (@ref setcomponenth2o),
    !> which requires transport results to be accurate to eight or nine significant digits.
    !> If multicomponent diffusion (MCD) is to be modeled, there is a capability to retrieve 
    !> aqueous species concentrations (@ref getspeciesconcentrations) and to set new 
    !> solution concentrations after MCD by using individual species concentrations
    !> (@ref speciesconcentrations2module). To use these methods the save-species 
    !> property needs to be turned on (@ref setspeciessaveon). If the save-species 
    !> property is on, findcomponents will generate a list of aqueous species 
    !> (@ref getspeciescount, @ref getspeciesnames), their diffusion coefficients 
    !> at 25 C (@ref getspeciesd25), their charge (@ref getspeciesz).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval              Number of components currently in the list, or IRESULT error 
    !> code (see @ref decodeerror).
    !> @see
    !> @ref getcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciesz,
    !> @ref setcomponenth2o.
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par The findcomponents method also generates lists of reactants--equilibrium phases,
    !> exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
    !> The lists are cumulative, including all reactants that were
    !> defined in the initial phreeqc instance at any time findcomponents was called.
    !> In addition, a list of phases is generated for which saturation indices may be calculated from the
    !> cumulative list of components.
    !> @see also
    !> @ref getequilibriumphasesnames,
    !> @ref getequilibriumphasescount,
    !> @ref getexchangenames,
    !> @ref getexchangespeciesnames,
    !> @ref getexchangespeciescount,
    !> @ref getgascomponentsnames,
    !> @ref getgascomponentscount,
    !> @ref getkineticreactionsnames,
    !> @ref getkineticreactionscount,
    !> @ref getsicount,
    !> @ref getsinames,
    !> @ref getsolidsolutioncomponentsnames,
    !> @ref getsolidsolutioncomponentscount,
    !> @ref getsolidsolutionnames,
    !> @ref getsurfacenames,
    !> @ref getsurfacespeciesnames,
    !> @ref getsurfacespeciescount,
    !> @ref getsurfacetypes.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: components(:)
    !> ncomps = brm%findcomponents()
    !> status = brm%getcomponents(components)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.

    INTEGER FUNCTION findcomponents(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
    class(bmi), intent(inout) :: self
    findcomponents = RM_findcomponents(self%bmiphreeqcrm_id)
    return
    END FUNCTION findcomponents

    !> Fills an array with the cell numbers in the user's numbering sytstem that map 
    !> to a cell in the PhreeqcRM numbering system. The mapping is defined by
    !> @ref createmapping.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param n             A cell number in the PhreeqcRM numbering system (0 <= n < 
    !> @ref getchemistrycellcount).
    !> @param list          Allocatable array to store the user cell numbers mapped to 
    !> PhreeqcRM cell @a n.
    !> @retval              IRESULT error code (see @ref decodeerror).
    !> @see
    !> @ref createmapping,
    !> @ref getchemistrycellcount,
    !> @ref getgridcellcount.
    !> @par C Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> if (brm%getbackwardmapping(cell_number, list) .eq. 0) then
    !>   if (fstr(1:l) .eq. "HYDRAULIC_K") then
    !>     my_basic_fortran_callback = K_ptr(list(1)+1)
    !>   endif
    !> endif
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.

    INTEGER FUNCTION getbackwardmapping(self, n, list)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in)    :: n
    INTEGER, INTENT(inout), allocatable :: list(:)
    getbackwardmapping = RM_getbackwardmapping(self%bmiphreeqcrm_id, n, list)
    return
    END FUNCTION getbackwardmapping

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Returns the number of chemistry cells in the reaction module. The number of chemistry 
    !> cells is defined by the set of non-negative integers in the mapping from user grid 
    !> cells (@ref createmapping). The number of chemistry cells is less than or equal to 
    !> the number of cells in the user's model.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval              Number of chemistry cells, or IRESULT error code 
    !> (see @ref decodeerror).
    !> @see
    !> @ref createmapping,
    !> @ref getgridcellcount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%createmapping(grid2chem)
    !> nchem = brm%getchemistrycellcount()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.

    INTEGER FUNCTION getchemistrycellcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getchemistrycellcount = RM_getchemistrycellcount(self%bmiphreeqcrm_id)
    return
    END FUNCTION getchemistrycellcount

    !> Returns a list of the names of the components identified by PhreeqcRM.
    !> The list contains all components (elements) found in solutions and reactants in the
    !> InitialPhreeqc instance by call(s) to @ref findcomponents.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param components       Allocatable, 1D character variable to receive the component names. 
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: components(:)
    !> status = brm%getcomponents(components)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    
    INTEGER FUNCTION getcomponents(self, components)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: components
    getcomponents = RM_getcomponents(self%bmiphreeqcrm_id, components)
    return 
    END FUNCTION getcomponents

    !> PRIVATE Returns the number of components in the reaction-module component list.
    !> The component list is generated by calls to @ref findcomponents.
    !> The return value from the last call to @ref findcomponents is equal to the return value from getcomponentcount.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 The number of components in the reaction-module component list, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getcomponents.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> ncomps1 = brm%getcomponentcount()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

    INTEGER FUNCTION getcomponentcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getcomponentcount = RM_getcomponentcount(self%bmiphreeqcrm_id)
    END FUNCTION getcomponentcount

    !> Transfer solution concentrations from each reaction cell to the concentration array 
    !> given in the argument list (@a c). Units of concentration for @a c are defined by 
    !> @ref setunitssolution. For concentration units of per liter, the solution volume 
    !> is used to calculate the concentrations for @a c. For mass fraction concentration units,
    !> the solution mass is used to calculate concentrations for @a c. Two options are 
    !> available for the volume and mass of solution that are used in converting to transport 
    !> concentrations: (1) the volume and mass of solution are calculated by PHREEQC, or
    !> (2) the volume of solution is the product of saturation (@ref setsaturationuser),
    !> porosity (@ref setporosity), and representative volume (@ref setrepresentativevolume), 
    !> and the mass of solution is volume times density as defined by @ref setdensityuser.
    !> @ref usesolutiondensityvolume determines which option is used. For option 1, the 
    !> databases that have partial molar volume definitions needed to calculate 
    !> solution volume accurately are phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param c         Array to receive the concentrations. Dimension of the array will be  
    !> set to (@a nxyz, @a ncomps), where @a nxyz is the number of user grid cells and 
    !> @a ncomps is the result of @ref findcomponents or @ref getcomponentcount.
    !> Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !>
    !> @see
    !> @ref findcomponents,
    !> @ref getcomponentcount,
    !> @ref getdensitycalculated,
    !> @ref getsaturationcalculated,
    !> @ref setconcentrations,
    !> @ref setdensityuser,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser,
    !> @ref setunitssolution,
    !> @ref usesolutiondensityvolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: c(:,:)
    !> status = brm%runcells()
    !> status = brm%getconcentrations(c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.

    INTEGER FUNCTION getconcentrations(self, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), allocatable :: c(:,:)
    getconcentrations = RM_getconcentrations(self%bmiphreeqcrm_id, c)
    return
    END FUNCTION getconcentrations
    
    !> Returns the user number of the current selected-output definition.
    !> @ref setcurrentselectedoutputusernumber or @ref setnthselectedoutput 
    !> specifies which of the selected-output definitions is used.
    !> @retval          User number of the the current selected-output definition,
    !> negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do isel = 1, brm%getselectedoutputcount()
    !>   status = brm%setnthselectedoutput(isel)
    !>   n_user = brm%getcurrentselectedoutputusernumber()
    !>   col = brm%getselectedoutputcolumncount()
    !>   allocate(selected_out(nxyz,col))
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !>   deallocate(selected_out)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getcurrentselectedoutputusernumber(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getcurrentselectedoutputusernumber = RM_getcurrentselectedoutputusernumber(self%bmiphreeqcrm_id)
    END FUNCTION getcurrentselectedoutputusernumber

    !> Transfer solution densities from the reaction cells to the array given in the 
    !> argument list (@a density). Densities are those calculated by the reaction module. 
    !> This method always  returns the calculated densities; @ref setdensityuser does 
    !> not affect the result. Only the following databases distributed with PhreeqcRM 
    !> have molar volume information needed to accurately calculate density: phreeqc.dat, 
    !> Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param density         Allocatable array to receive the densities. Dimension 
    !> of the array is set to @a nxyz, where @a nxyz is the number of user grid cells 
    !> (@ref getgridcellcount). Values for inactive cells are set to 1e30.
    !> @retval IRESULT     0 is success, negative is failure (See @ref decodeerror).
    !>
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: density(:)
    !> status = brm%runcells()
    !> status = brm%getdensitycalculated(density)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getdensitycalculated(self, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), dimension(:), allocatable  :: density
    getdensitycalculated = RM_getdensitycalculated(self%bmiphreeqcrm_id, density)
    return
    END FUNCTION getdensitycalculated
    
    INTEGER FUNCTION GetDensity(self, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), dimension(:), allocatable  :: density
    GetDensity = RM_getdensitycalculated(self%bmiphreeqcrm_id, density)
    return
    END FUNCTION GetDensity

    !> Returns an array with the ending cell numbers from the range of cell numbers 
    !> assigned to each worker.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param ec           Array to receive the ending cell numbers. Dimension of the array is
    !> the number of threads (OpenMP) or the number of processes (MPI).
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref create,
    !> @ref getmpitasks,
    !> @ref getstartcell,
    !> @ref getthreadcount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer, allocatable :: ec(:)
    !> status = brm%getendcell(ec)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getendcell(self, ec)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(inout), DIMENSION(:), allocatable :: ec
    getendcell = RM_getendcell(self%bmiphreeqcrm_id, ec)
    RETURN
    END FUNCTION getendcell

    !> Returns the number of equilibrium phases in the initial-phreeqc module.
    !> @ref findcomponents must be called before @ref getequilibriumphasescount.
    !> This method may be useful when generating selected output definitions related to
    !> equilibrium phases.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The number of equilibrium phases in the initial-phreeqc module.
    !> @see
    !> @ref findcomponents,
    !> @ref getequilibriumphasesnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getequilibriumphasesnames(names)
    !> do i = 1, brm%getequilibriumphasescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getequilibriumphasescount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getequilibriumphasescount = RM_getequilibriumphasescount(self%bmiphreeqcrm_id)
    END FUNCTION getequilibriumphasescount

    !> Retrieve a list of equilibrium phase names.
    !> The list includes all phases included in any EQUILIBRIUM_PHASES definitions in
    !> the initial-phreeqc module. @ref findcomponents must be called before 
    !> @ref getequilibriumphasesnames. This method may be useful when generating 
    !> selected output definitions related to equilibrium phases.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Array of equilibrium phase names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getequilibriumphasescount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%getequilibriumphasesnames(names)
    !> do i = 1, brm%getequilibriumphasescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getequilibriumphasesnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getequilibriumphasesnames = RM_getequilibriumphasesnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getequilibriumphasesnames


    ! INTEGER FUNCTION getequilibriumphasesname(self, num, name)
    ! USE ISO_C_BINDING
    ! IMPLICIT NONE
    ! INTERFACE
    ! INTEGER(KIND=C_INT) FUNCTION RMF_getequilibriumphasesname(id, num, name, l) &
    !     BIND(C, NAME='RMF_getequilibriumphasesname')
    ! USE ISO_C_BINDING
    ! IMPLICIT NONE
    ! INTEGER(KIND=C_INT), INTENT(inout) :: id, num, l
    ! CHARACTER(KIND=C_CHAR), INTENT(inout) :: name(*)
    ! END FUNCTION RMF_getequilibriumphasesname
    ! END INTERFACE
    ! INTEGER, INTENT(inout) :: id, num
    ! CHARACTER(len=*), INTENT(inout) :: name
    ! getequilibriumphasesname = RMF_getequilibriumphasesname(id, num, name, len(name))
    ! return
    ! END FUNCTION getequilibriumphasesname

    !> Returns a string containing error messages related to the last call to a PhreeqcRM 
    !> method to the character argument (@a errstr).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param errstr           The error string related to the last call to a PhreeqcRM method.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: errstr
    !> if (status .lt. 0) then
    !>   status = brm%geterrorstring(errstr)
    !>   write(*,"(A)") errstr
    !>   status = brm%destroy()
    !>   stop
    !> endif 
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION geterrorstring(self, errstr)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, INTENT(inout) :: errstr
    geterrorstring = RM_geterrorstring(self%bmiphreeqcrm_id, errstr)
    END FUNCTION geterrorstring

    !> Retrieves a list of exchange names.
    !> @ref findcomponents must be called before @ref getexchangenames.
    !> The exchange names array is the same length as the exchange species names array
    !> and provides the corresponding exchange site (for example, X corresponing to NaX).
    !> This method may be useful when generating selected output definitions related to exchangers.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array of  exchange names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getexchangespeciescount, @ref getexchangespeciesnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getexchangespeciesnames(names)
    !> do i = 1, brm%getexchangespeciescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getexchangenames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getexchangenames = RM_getexchangenames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getexchangenames

    !> Returns the number of exchange species in the initial-phreeqc module.
    !> @ref findcomponents must be called before @ref getexchangespeciescount.
    !> This method may be useful when generating selected output definitions related 
    !> to exchangers.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The number of exchange species found by call(s) to 
    !> @ref findcomponents.
    !> @see
    !> @ref findcomponents, @ref getexchangespeciesnames, @ref getexchangenames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names
    !> status = brm%getexchangespeciesnames(names)
    !> do i = 1, brm%getexchangespeciescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getexchangespeciescount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getexchangespeciescount = RM_getexchangespeciescount(self%bmiphreeqcrm_id)
    END FUNCTION getexchangespeciescount

    !> Retrieves a list of exchange species names.
    !> The list of exchange species (such as "NaX") is derived from the list of components
    !> (@ref findcomponents) and the list of all exchange names (such as "X")
    !> that are found by call(s) to @ref findcomponents. @ref findcomponents must 
    !> be called before @ref getexchangespeciesnames. This method may be useful 
    !> when generating selected output definitions related to exchangers.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array of exchange species names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getexchangespeciescount, @ref getexchangenames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names
    !> status = brm%getexchangespeciesnames(names)
    !> do i = 1, brm%getexchangespeciescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getexchangespeciesnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getexchangespeciesnames = RM_getexchangespeciesnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getexchangespeciesnames

    !> Returns the reaction-module file prefix to the character argument (@a prefix).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param prefix           Allocatable character string where the prefix is written.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setfileprefix.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: prefix
    !> status = brm%getfileprefix(prefix)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getfileprefix(self, prefix)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, INTENT(inout) :: prefix
    getfileprefix = RM_getfileprefix(self%bmiphreeqcrm_id, prefix)
    END FUNCTION getfileprefix

    !> Returns the number of gas phase components found by call(s) to @ref findcomponents.
    !> @ref findcomponents must be called before @ref getgascomponentscount.
    !> This method may be useful when generating selected output definitions related to
    !> gas phases.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The number of gas phase components in the initial-phreeqc module.
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentsnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getgascomponentsnames(names)
    !> do i = 1, brm%getgascomponentscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getgascomponentscount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getgascomponentscount = RM_getgascomponentscount(self%bmiphreeqcrm_id)
    END FUNCTION getgascomponentscount

    !> Retrieves a list of the gas component names.
    !> The list includes all gas components found by calls to @ref findcomponents.
    !> @ref findcomponents must be called before @ref getgascomponentsnames.
    !> This method may be useful when generating selected output definitions related 
    !> to gas phases.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array of gas component names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getgascomponentsnames(names)
    !> do i = 1, brm%getgascomponentscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getgascomponentsnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getgascomponentsnames = RM_getgascomponentsnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getgascomponentsnames

    !> Transfer moles of gas components from each reaction cell
    !> to the array given in the argument list (@a gas_moles).
    !>
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_moles        Allocatable array to receive the moles of gas components 
    !> for each cell. Dimension of the array is set to (@a nxyz, @a ngas_comps),
    !> where @a nxyz is the number of user grid cells and @a ngas_comps is the result
    !> of @ref getgascomponentscount. If a gas component is not defined for a cell,
    !> the number of moles is set to -1. Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !>
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascomppressures,
    !> @ref getgascompphi,
    !> @ref getgasphasevolume,
    !> @ref setgascompmoles,
    !> @ref setgasphasevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: gas_moles(:,:)
    !> status = brm%runcells()
    !> status = brm%getgascompmoles(gas_moles)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getgascompmoles(self, gas_moles)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:,:), allocatable, TARGET :: gas_moles    
    getgascompmoles = RM_getgascompmoles(self%bmiphreeqcrm_id, gas_moles)
    return
    END FUNCTION getgascompmoles

    !> Transfer pressures of gas components from each reaction cell
    !> to the array given in the argument list (@a gas_p).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_p          Allocatable array to receive the moles of gas components for 
    !> each cell. Dimension of the array is set to (@a nxyz, @a ngas_comps),
    !> where @a nxyz is the number of user grid cells and @a ngas_comps is the result
    !> of @ref getgascomponentscount. If a gas component is not defined for a cell,
    !> the pressure is set to -1. Values for inactive cells are set to 1e30.
    !> @retval IRESULT     0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascompmoles,
    !> @ref getgascompphi,
    !> @ref getgasphasevolume,
    !> @ref setgascompmoles,
    !> @ref setgasphasevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: gas_p(:,:)
    !> status = brm%runcells()
    !> status = brm%getgascomppressures(gas_p)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getgascomppressures(self, gas_p)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:,:), allocatable, TARGET :: gas_p
    getgascomppressures = RM_getgascomppressures(self%bmiphreeqcrm_id, gas_p )
    return
    END FUNCTION getgascomppressures

    !> Transfer fugacity coefficients (phi) of gas components from each reaction cell
    !> to the array given in the argument list (@a gas_phi). Fugacity of a gas component
    !> is equal to the pressure of the component times the fugacity coefficient.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_phi        Allocatable array to receive the fugacity coefficients
    !> of gas components for each cell. Dimension of the array is set to (@a nxyz, @a ngas_comps),
    !> where @a nxyz is the number of user grid cells and @a ngas_comps is the result
    !> of @ref getgascomponentscount. If a gas component is not defined for a cell,
    !> the fugacity coefficient is set to -1. Values for inactive cells are set to 1e30.
    !> @retval IRESULT     0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascompmoles,
    !> @ref getgascomppressures,
    !> @ref getgasphasevolume,
    !> @ref setgascompmoles,
    !> @ref setgasphasevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: gas_phi(:,:)
    !> status = brm%runcells()
    !> status = brm%getgascompphi(gas_phi)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getgascompphi(self, gas_phi)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:,:), allocatable, TARGET :: gas_phi
    getgascompphi = RM_getgascompphi(self%bmiphreeqcrm_id, gas_phi)
    return
    END FUNCTION getgascompphi

    !> Transfer volume of gas from each reaction cell
    !> to the array given in the argument list (@a gas_volume).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_volume       Array to receive the gas phase volumes.
    !> Dimension of the array is set to @a nxyz, where @a nxyz is the number of user 
    !> grid cells (@ref getgridcellcount). If a gas phase is not defined for a cell, 
    !> the volume is set to -1. Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascompmoles,
    !> @ref getgascompphi,
    !> @ref getgascomppressures,
    !> @ref setgascompmoles
    !> @ref setgasphasevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: gas_volume(:)
    !> status = brm%runcells()
    !> status = brm%getgasphasevolume(gas_volume)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getgasphasevolume(self, gas_volume)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:), allocatable, TARGET :: gas_volume
    getgasphasevolume = RM_getgasphasevolume(self%bmiphreeqcrm_id, gas_volume)
    return
    END FUNCTION getgasphasevolume

    !> Returns the gram formula weights (g/mol) for the components in the reaction-module 
    !> component list.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gfw              Array to receive the gram formula weights. Dimension of the 
    !> array is set to @a ncomps, where @a ncomps is the number of components in the 
    !> component list.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getcomponents,
    !> @ref getcomponentcount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: components(:)
    !> real(kind=8), allocatable   :: gfw(:)
    !> ncomps = brm%findcomponents()
    !> status = brm%getgfw(gfw)
    !> status = brm%getcomponents(components)
    !> do i = 1, ncomps
    !>   write(string,"(A10, F15.4)") components(i), gfw(i)
    !>   status = brm%outputmessage(string)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getgfw(self, gfw)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(inout), allocatable  :: gfw
    getgfw = RM_getgfw(self%bmiphreeqcrm_id, gfw)
    END FUNCTION getgfw

    !> Returns the number of grid cells in the user's model, which is defined in 
    !> the creation or initialization of the PhreeqcRM instance.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            Number of grid cells in the user's model, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref create,
    !> @ref createmapping,
    !> @ref initializeyaml.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> nxyz = brm%getgridcellcount()
    !> write(string1, "(A,I)") "Number of grid cells in the user's model: ", nxyz
    !> status = brm%outputmessage(trim(string1))
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getgridcellcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getgridcellcount = RM_getgridcellcount(self%bmiphreeqcrm_id)
    END FUNCTION getgridcellcount

    !> Returns an IPhreeqc id for the @a ith IPhreeqc instance in the reaction module.
    !> For the threaded version, there are @a nthreads + 2 IPhreeqc instances, where
    !> @a nthreads is defined in the constructor (@ref create) or initialization.
    !> The number of threads can be determined by @ref getthreadcount.
    !> The first @a nthreads (0 based) instances will be the workers, the
    !> next (@a nthreads) is the InitialPhreeqc instance, and the next (@a nthreads + 1) 
    !> is the Utility instance. Getting the IPhreeqc pointer for one of these instances 
    !> allows the user to use any of the IPhreeqc methods on that instance.
    !> For MPI, each process has exactly three IPhreeqc instances, one worker (number 0),
    !> one InitialPhreeqc instance (number 1), and one Utility instance (number 2).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param i         The number of the IPhreeqc instance to be retrieved (0 based).
    !> @retval          IPhreeqc id for the @a ith IPhreeqc instance, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref create,
    !> @ref getthreadcount.
    !> See IPhreeqc documentation for descriptions of IPhreeqc methods.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> ! Utility pointer is worker number nthreads + 1
    !> iphreeqc_id1 = brm%getiphreeqcid(brm%getthreadcount() + 1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getiphreeqcid(self, i)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: i
    getiphreeqcid = RM_getiphreeqcid(self%bmiphreeqcrm_id, i)
    END FUNCTION getiphreeqcid

    !> Transfer the concentration from each cell for one component to the array given in the 
    !> argument list (@a c). The concentrations are those resulting from the last call
    !> to @ref runcells. Units of concentration for @a c are defined by @ref setunitssolution.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param i                One-based index for the component to retrieve. Indices refer
    !> to the order produced by @ref getcomponents. The total number of components is given by
    !> @ref getcomponentcount.
    !> @param c                Allocatable array to receive the component concentrations.
    !> Dimension of the array is set to @a nxyz, where @a nxyz is the number of 
    !> user grid cells (@ref getgridcellcount). Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see    @ref findcomponents, @ref getcomponents, @ref getcomponentcount, 
    !> @ref getconcentrations.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable, dimension(:) :: c
    !> status = brm%runcells()
    !> status = brm%phreeqc_rm.getithconcentration(1, c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getithconcentration(self, i, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: i
    real(kind=8), INTENT(inout), DIMENSION(:), allocatable :: c
    getithconcentration = RM_getithconcentration(self%bmiphreeqcrm_id, i, c)
    return
    END FUNCTION getithconcentration

    !> Transfer the concentrations for one species from each cell to the array given in the
    !> argument list (@a c). The concentrations are those resulting from the last call
    !> to @ref runcells. Units of concentration for @a c are mol/L.
    !> To retrieve species concentrations, @ref setspeciessaveon must be set to @a true.
    !> This method is for use with multicomponent diffusion calculations.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param i                One-based index for the species to retrieve. Indices refer
    !> to the order given by @ref getspeciesnames. The total number of species is given
    !> by @ref getspeciescount.
    !> @param c                Allocatable array to receive the species concentrations.
    !> Dimension of the array is set to @a nxyz, where @a nxyz is the number of
    !> user grid cells (@ref getgridcellcount). Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see         @ref findcomponents, @ref getspeciescount, @ref getspeciesnames,
    !> @ref getspeciesconcentrations, @ref setspeciessaveon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable, dimension(:) :: c
    !> status = brm%runcells()
    !> status = brm%getithspeciesconcentration(1, c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getithspeciesconcentration(self, i, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: i
    real(kind=8), INTENT(inout), DIMENSION(:), allocatable :: c
    getithspeciesconcentration = RM_getithspeciesconcentration(self%bmiphreeqcrm_id, i, c)
    return
    END FUNCTION getithspeciesconcentration

    !> Returns the number of kinetic reactions found by call(s) to @ref findcomponents.
    !> @ref findcomponents must be called before @ref getkineticreactionscount.
    !> This method may be useful when generating selected output definitions related to
    !> kinetic reactions.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The number of kinetic reactions in the initial-phreeqc module.
    !> @see
    !> @ref findcomponents,
    !> @ref getkineticreactionsnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getkineticreactionsnames(names)
    !> do i = 1, brm%getkineticreactionscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getkineticreactionscount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getkineticreactionscount = RM_getkineticreactionscount(self%bmiphreeqcrm_id)
    END FUNCTION getkineticreactionscount

    !> Retrieves a list of kinetic reaction names.
    !> The list includes all kinetic reactions found by call(s) to @ref findcomponents.
    !> @ref findcomponents must be called before @ref getkineticreactionsnames.
    !> This method may be useful when generating selected output definitions related to 
    !> kinetic reactions.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array for kinetic reaction names
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getkineticreactionscount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getkineticreactionsnames(names)
    !> do i = 1, brm%getkineticreactionscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getkineticreactionsnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getkineticreactionsnames = RM_getkineticreactionsnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getkineticreactionsnames

    !> Returns the MPI task number. For the OPENMP version, the task number is always
    !> zero and the result of @ref getmpitasks is one. For the MPI version,
    !> the root task number is zero, and all workers have a task number greater than zero.
    !> The number of tasks can be obtained with @ref getmpitasks. The number of
    !> tasks and computer hosts are determined at run time by the mpiexec command, and the
    !> number of reaction-module processes is defined by the communicator used in
    !> constructing the reaction modules (@ref create).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 The MPI task number for a process, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref getmpitasks.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string1, "(A,I)") "MPI task number: ", brm%getmpimyself()
    !> status = brm%outputmessage(string1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.

    INTEGER FUNCTION getmpimyself(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getmpimyself = RM_getmpimyself(self%bmiphreeqcrm_id)
    END FUNCTION getmpimyself

    !> Returns the number of MPI processes (tasks) assigned to the reaction module.
    !> For the OPENMP version, the number of tasks is always
    !> one (although there may be multiple threads, @ref getthreadcount),
    !> and the task number returned by @ref getmpimyself is zero. For the MPI version, 
    !> the number of tasks and computer hosts are determined at run time by the mpiexec 
    !> command. An MPI communicator is used in constructing reaction modules for MPI. 
    !> The communicator may define a subset of the total number of MPI processes. 
    !> The root task number is zero, and all workers have a task number greater than zero.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 The number of MPI  processes assigned to the reaction module,
    !> negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getmpimyself,
    !> @ref create.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> mpi_tasks = brm%getmpitasks()
    !> write(string1, "(A,I)") "Number of MPI processes: ", mpi_tasks
    !> status = brm%outputmessage(string1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getmpitasks(self)
    USE ISO_C_BINDING
	class(bmi), intent(inout) :: self
    getmpitasks = RM_getmpitasks(self%bmiphreeqcrm_id)
    END FUNCTION getmpitasks

    !> Returns the user number for the @a nth selected-output definition.
    !> Definitions are sorted by user number. Phreeqc allows multiple selected-output
    !> definitions, each of which is assigned a nonnegative integer identifier by the
    !> user. The number of definitions can be obtained by @ref getselectedoutputcount.
    !> To cycle through all of the definitions, getnthselectedoutputusernumber
    !> can be used to identify the user number for each selected-output definition
    !> in sequence. @ref setcurrentselectedoutputusernumber is then used to select
    !> that user number for selected-output processing.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param n                The sequence number of the selected-output definition for 
    !> which the user number will be returned. Fortran, 1 based.
    !> @retval                 The user number of the @a nth selected-output definition, 
    !> negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   write(*,*) "Selected output sequence number: ", isel)
    !>   write(*,*) "Selected output user number:     ", n_user)
    !>   col = brm%getselectedoutputcolumncount()
    !>   allocate(selected_out(nxyz,col))
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !>   deallocate(selected_out)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getnthselectedoutputusernumber(self, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: n
    getnthselectedoutputusernumber = RM_getnthselectedoutputusernumber(self%bmiphreeqcrm_id, n)
    END FUNCTION getnthselectedoutputusernumber

    !> Transfer current porosities to the array given in the argument list (@a porosity).
    !> Porosity is not changed by PhreeqcRM; the values are either the default values
    !> or the values set by the last call to @ref setporosity.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param porosity           Array to receive the porosities. Dimension of the array 
    !> is set to @a nxyz, where @a nxyz is the number of user grid cells  
    !> (@ref getgridcellcount). Values for inactive cells are set to 1e30.
    !> @retval IRESULT          0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: porosity(:)
    !> status = brm%getporosity(porosity)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getporosity(self, porosity)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), dimension(:), allocatable :: porosity
    getporosity = RM_getporosity(self%bmiphreeqcrm_id, porosity)
    return
    END FUNCTION getporosity

    !> Transfer current pressures to the array given in the argument list (@a pressure).
    !> Pressure is not usually calculated by PhreeqcRM; the values are either the default values
    !> or the values set by the last call to @ref setpressure. Pressures can be calculated
    !> by PhreeqcRM if a fixed-volume GAS_PHASE is used.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param pressure        Array to receive the porosities. Dimension of the array is set 
    !> to @a nxyz, where @a nxyz is the number of user grid cells (@ref getgridcellcount). 
    !> Values for inactive cells are set to 1e30.
    !> @retval IRESULT          0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: pressure(:)
    !> status = brm%getpressure(pressure)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getpressure(self, pressure)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), dimension(:), allocatable :: pressure
    getpressure = RM_getpressure(self%bmiphreeqcrm_id, pressure)
    return
    END FUNCTION getpressure

    !> Returns an array of saturations (@a sat_calc) as calculated by the reaction module.
    !> Reactions will change the volume of solution in a cell. This method always returns 
    !> solution_volume/(rv * porosity); the method  @ref setsaturationuser has no effect 
    !> on the values returned. The transport code must decide whether to ignore or account 
    !> for this change in solution volume due to reactions. Following reactions, the cell 
    !> saturation is calculated as solution volume (@ref getsolutionvolume) divided by 
    !> the product of representative volume (@ref setrepresentativevolume) and the porosity 
    !> (@ref setporosity). The cell saturation returned by @a getsaturationcalculated 
    !> may be less than or greater than the saturation set by the transport code 
    !> (@ref setsaturationuser), and may be greater than or less than 1.0, even in 
    !> fully saturated simulations. Only the following databases distributed with PhreeqcRM
    !> have molar volume information needed to accurately calculate solution volume and 
    !> saturation: phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param sat_calc      Array to receive the saturations. Dimension of the array is set to 
    !> @a nxyz, where @a nxyz is the number of user grid cells (@ref getgridcellcount).
    !> Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getsolutionvolume,
    !> @ref setporosity,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: sat_calc(:)
    !> status = brm%runcells()
    !> status = brm%getsaturationcalculated(sat_calc)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getsaturationcalculated(self, sat_calc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:), allocatable :: sat_calc
    getsaturationcalculated = RM_getsaturationcalculated(self%bmiphreeqcrm_id, sat_calc)
    END FUNCTION getsaturationcalculated
    
    INTEGER FUNCTION GetSaturation(self, sat_calc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(inout), DIMENSION(:), allocatable :: sat_calc
    GetSaturation = RM_getsaturationcalculated(self%bmiphreeqcrm_id, sat_calc)
    END FUNCTION GetSaturation

    !> Populates an array with values from the current selected-output definition. 
    !> @ref setcurrentselectedoutputusernumber or @ref setnthselectedoutput determines 
    !> which of the selected-output definitions is used to populate the array.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param so           An array to contain the selected-output values. Size of the array 
    !> is set to (@a nxyz, @a col), where @a nxyz is the number of grid cells in the user's 
    !> model (@ref getgridcellcount), and @a col is the number of columns in the 
    !> selected-output definition (@ref getselectedoutputcolumncount).
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: headings
    !> real(kind=8), allocatable :: selected_out
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   status = brm%getselectedoutput(selected_out)
    !>   status = brm%getselectedoutputheadings(headings)
    !>   ! Print results
    !>   do i = 1, brm%getselectedoutputrowcount()
    !>     write(*,*) "Cell number ", i
    !>     write(*,*) "     Selected output: "
    !>     do j = 1, brm%getselectedoutputcolumncount()
    !>       write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ", trim(headings(j)),": ", selected_out(i,j)
    !>     enddo
    !>   enddo
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getselectedoutput(self, so)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:,:), INTENT(inout), allocatable :: so
    getselectedoutput = RM_getselectedoutput(self%bmiphreeqcrm_id, so)
    END FUNCTION getselectedoutput

    !> Returns the number of columns in the current selected-output definition. 
    !> @ref setcurrentselectedoutputusernumber or @ref setnthselectedoutput
    !> determines which of the selected-output definitions is used.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 Number of columns in the current selected-output 
    !> definition, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: selected_out(:,:)
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getselectedoutputcolumncount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getselectedoutputcolumncount = RM_getselectedoutputcolumncount(self%bmiphreeqcrm_id)
    END FUNCTION getselectedoutputcolumncount

    !> Returns the number of selected-output definitions. 
    !> @ref setcurrentselectedoutputusernumber or @ref setnthselectedoutput
    !> determines which of the selected-output definitions is used.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 Number of selected-output definitions, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: selected_out(:,:)
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getselectedoutputcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getselectedoutputcount = RM_getselectedoutputcount(self%bmiphreeqcrm_id)
    END FUNCTION getselectedoutputcount

   
    !> Returns the selected-output headings for the current selected-output file.
    !> @ref setcurrentselectedoutputusernumber or @ref setnthselectedoutput
    !> determines which of the selected-output definitions is used.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param headings         Allocatable, 1D character variable.
    !> to receive the headings. Character length and dimension will be allocated as needed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable, :: headings(:)
    !> do isel = 1, brm%getselectedoutputcount()
    !>   status = brm%setnthselectedoutput(isel)
    !>   status = brm%getselectedoutputheadings(headings)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getselectedoutputheadings(self, headings)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: headings
    getselectedoutputheadings = RM_getselectedoutputheadings(self%bmiphreeqcrm_id, headings)
    return
    END FUNCTION getselectedoutputheadings

    !> Returns the number of rows in the current selected-output definition. However, the method
    !> is included only for convenience; the number of rows is always equal to the number of
    !> grid cells in the user's model, and is equal to @ref getgridcellcount.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 Number of rows in the current selected-output definition, 
    !> negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: headings
    !> real(kind=8), allocatable :: selected_out
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   status = brm%getselectedoutput(selected_out)
    !>   status = brm%getselectedoutputheadings(headings)
    !>   ! Print results
    !>   do i = 1, brm%getselectedoutputrowcount()
    !>     write(*,*) "Cell number ", i
    !>     write(*,*) "     Selected output: "
    !>     do j = 1, brm%getselectedoutputcolumncount()
    !>       write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ", trim(headings(j)),": ", selected_out(i,j)
    !>     enddo
    !>   enddo
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

    INTEGER FUNCTION getselectedoutputrowcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getselectedoutputrowcount = RM_getselectedoutputrowcount(self%bmiphreeqcrm_id)
    END FUNCTION getselectedoutputrowcount

    !> Returns the number of phases in the found by @ref findcomponents for 
    !> which saturation indices can be calculated.
    !> @ref findcomponents must be called before @ref getsicount.
    !> This method may be useful when generating selected output definitions related to
    !> saturation indices.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval              The number of phases in the initial-phreeqc module for which 
    !> saturation indices can be calculated.
    !> @see
    !> @ref findcomponents,
    !> @ref getsinames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getsinames(names)
    !> do i = 1, brm%getsicount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsicount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getsicount = RM_getsicount(self%bmiphreeqcrm_id)
    END FUNCTION getsicount

    !> Retrieves a list of all phases for which saturation indices can be calculated.
    !> The list includes all phases that contain only elements included in the components in
    !> the initial-phreeqc module. The list assumes that all components are present to be 
    !> able to calculate the entire list of SIs; it may be that one or more components are 
    !> missing in any specific cell. @ref findcomponents must be called before 
    !> @ref getsinames. This method may be useful when generating selected output 
    !> definitions related to saturation indices.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array for saturation-index-phase names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getsicount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getsinames(names)
    !> do i = 1, brm%getsicount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsinames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsinames = RM_getsinames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsinames

    !> Returns the number of solid solution components in the initial-phreeqc module.
    !> @ref findcomponents must be called before @ref getsolidsolutioncomponentscount.
    !> This method may be useful when generating selected output definitions related to 
    !> solid solutions.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval         The number of solid solution components in the initial-phreeqc module.
    !> @see
    !> @ref findcomponents, @ref getsolidsolutioncomponentsnames,
    !> @ref getsolidsolutionnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getsolidsolutioncomponentsnames(names)
    !> do i = 1, brm%getsolidsolutioncomponentscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsolidsolutioncomponentscount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getsolidsolutioncomponentscount = RM_getsolidsolutioncomponentscount(self%bmiphreeqcrm_id)
    END FUNCTION getsolidsolutioncomponentscount

    !> Retrieves a list of the solid solution component names.
    !> The list includes all solid solution components found by call(s) to 
    !> @ref findcomponents. @ref findcomponents must be called before 
    !> @ref getsolidsolutioncomponentsnames. This method may be useful when 
    !> generating selected output definitions related to solid solutions.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names       Allocatable array for the solid solution compnent names
    !> @retval IRESULT 0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents, @ref getsolidsolutioncomponentscount, 
    !> @ref getsolidsolutionnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%getsolidsolutioncomponentsnames(names)
    !> do i = 1, brm%getsolidsolutioncomponentscount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsolidsolutioncomponentsnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsolidsolutioncomponentsnames = RM_getsolidsolutioncomponentsnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsolidsolutioncomponentsnames

    !> Retrieves a list of solid solution names.
    !> The list includes solid solution names found by call(s) to @ref findcomponents.
    !> The solid solution names array is the same length as the solid solution components 
    !> array and provides the corresponding name of solid solution containing the component.
    !> @ref findcomponents must be called before @ref getsolidsolutionnames.
    !> This method may be useful when generating selected output definitions related to 
    !> solid solutions.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array for the solid solution names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents, @ref getsolidsolutioncomponentscount, 
    !> @ref getsolidsolutioncomponentsnames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:), ss_names(:)
    !> status = brm%getsolidsolutioncomponentsnames(names)
    !> status = brm%getsolidsolutionnames(names)
    !> do i = 1, brm%getsolidsolutioncomponentscount()
    !>   write(*,*) names(i), ss_names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsolidsolutionnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsolidsolutionnames = RM_getsolidsolutionnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsolidsolutionnames


    !> Transfer solution volumes from the reaction cells to the array given in the 
    !> argument list (@a vol).
    !> Solution volumes are those calculated by the reaction module. Only the following 
    !> databases distributed with PhreeqcRM have molar volume information needed to 
    !< accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param vol            Allocatable array to receive the solution volumes. 
    !> Dimension of the array is set to (@a nxyz), where @a nxyz is the number of user 
    !> grid cells. Values for inactive cells are set to 1e30.
    !> @retval IRESULT    0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getsaturationcalculated.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: volume(:)allocate(volume)
    !> status = brm%runcells()
    !> status = brm%getsolutionvolume(volume)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getsolutionvolume(self, vol)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:), allocatable :: vol
    getsolutionvolume = RM_getsolutionvolume(self%bmiphreeqcrm_id, vol)
    END FUNCTION getsolutionvolume

    !> Transfer concentrations of aqueous species for each cell to the argument 
    !> (@a species_conc).
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> The list of aqueous species is determined by call(s) to  @ref findcomponents 
    !> and includes all aqueous species that can be made from the set of components.
    !> Solution volumes used to calculate mol/L are calculated by the reaction module.
    !> Only the following databases distributed with PhreeqcRM have molar volume information
    !> needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param species_conc     Allocatable array to receive the aqueous species 
    ! concentrations. Dimension of the array is set to (@a nxyz, @a nspecies),
    !> where @a nxyz is the number of user grid cells (@ref getgridcellcount),
    !> and @a nspecies is the number of aqueous species (@ref getspeciescount).
    !> Concentrations are moles per liter. Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: species_c(:,:)
    !> status = brm%setspeciessaveon(1)
    !> ncomps = brm%findcomponents()
    !> status = brm%runcells()
    !> status = brm%getspeciesconcentrations(species_c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getspeciesconcentrations(self, species_conc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:,:), allocatable :: species_conc
    getspeciesconcentrations = RM_getspeciesconcentrations(self%bmiphreeqcrm_id, species_conc)
    END FUNCTION getspeciesconcentrations

    !> Returns the number of aqueous species  in the reaction module.
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> The list of aqueous species is determined by @ref findcomponents and 
    !> includes all aqueous species that can be made from the set of components.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      The number of aqueous species, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setspeciessaveon(1)
    !> ncomps = findcomponents()
    !> status = brm%getspeciescount()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getspeciescount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getspeciescount = RM_getspeciescount(self%bmiphreeqcrm_id)
    END FUNCTION getspeciescount

    !> Transfers diffusion coefficients at 25C to the array argument (@a diffc).
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally 
    !> in the database file. Databases distributed with the reaction module that have 
    !> diffusion coefficients defined are phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param diffc          Allocatable array to receive the diffusion coefficients 
    !> at 25 C, m^2/s. Dimension of the array is set to  @a nspecies, where @a nspecies 
    !> is the number of aqueous species (@ref getspeciescount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: diffc(:)
    !> status = brm%setspeciessaveon(1)
    !> ncomps = brm%findcomponents()
    !> nspecies = brm%getspeciesd25(diffc)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getspeciesd25(self, diffc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:), allocatable :: diffc
    getspeciesd25 = RM_getspeciesd25(self%bmiphreeqcrm_id, diffc)
    END FUNCTION getspeciesd25

    !> Transfer log10 aqueous-species activity coefficients to the array argument 
    !> (@a species_log10gammas)
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> The list of aqueous species is determined by call(s) to  
    !> @ref findcomponents and includes all aqueous species that can be made 
    !> from the set of components.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param species_log10gammas  Allocatable array to receive the log10 aqueous 
    !> species activity coefficients. Dimension of the array is (@a nxyz, @a nspecies),
    !> where @a nxyz is the number of user grid cells (@ref getgridcellcount),
    !> and @a nspecies is the number of aqueous species (@ref getspeciescount).
    !> Values for inactive cells are set to 1e30.
    !> @retval IRESULT     0 is success, negative is failure 
    !> (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: species_log10gammas(:)
    !> status = brm%getspecieslog10gammas(species_log10gammas)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getspecieslog10gammas(self, species_log10gammas)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:,:), allocatable :: species_log10gammas
    getspecieslog10gammas = RM_getspecieslog10gammas(self%bmiphreeqcrm_id, species_log10gammas)
    END FUNCTION getspecieslog10gammas

    !> Transfer log10 aqueous-species molalities to the array argument 
    !> (@a species_log10molalities)
    !> To use this method @ref setspeciessaveon must be set to @a true.
    !> The list of aqueous species is determined by @ref findcomponents 
    !> and includes all aqueous species that can be made from the set of components.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param species_log10molalities  Array to receive the aqueous species log10 
    !> molalities. Dimension of the array is set to (@a nxyz, @a nspecies),
    !> where @a nxyz is the number of user grid cells (@ref getgridcellcount),
    !> and @a nspecies is the number of aqueous species (@ref getspeciescount).
    !> Values for inactive cells are set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: species_log10gammas(:)
    !> status = brm%getspecieslog10molalities(species_log10gammas)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION getspecieslog10molalities(self, species_log10molalities)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:,:), allocatable :: species_log10molalities
    getspecieslog10molalities = RM_getspecieslog10molalities(self%bmiphreeqcrm_id, species_log10molalities)
    END FUNCTION getspecieslog10molalities

    !> Provides a list of aqueous species names to the argument (@a names).
    !> This method is intended for use with multicomponent-diffusion transport calculations,
    !> and @ref setspeciessaveon must be set to @a true. The list of aqueous
    !> species is determined by @ref findcomponents and includes all
    !> aqueous species that can be made from the set of components.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable character array to receive the species names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: names(:)
    !> status = brm%setspeciessaveon(1)
    !> ncomps = brm%findcomponents()
    !> status = brm%getspeciesnames(names)
    !> do i = 1, brm%getspeciescount()
    !>   write(*,*) names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getspeciesnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getspeciesnames = RM_getspeciesnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getspeciesnames    


    !> Returns the value of the species-save property.
    !> By default, concentrations of aqueous species are not saved. Setting the species-
    !> save property to true allows aqueous species concentrations to be retrieved
    !> with @ref getspeciesconcentrations, and solution compositions to be set with
    !> @ref speciesconcentrations2module.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      0, species are not saved; 1, species are saved; negative is 
    !> failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module. 
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> save_on = brm%getspeciessaveon()
    !> if (save_on .ne. 0) then
    !>   write(*,*) "Reaction module is saving species concentrations"
    !> else
    !>   write(*,*) "Reaction module is not saving species concentrations"
    !> end
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getspeciessaveon(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getspeciessaveon = RM_getspeciessaveon(self%bmiphreeqcrm_id)
    END FUNCTION getspeciessaveon

    !> Transfers the charge of each aqueous species to the array argument (@a  z).
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param z                Allocatable array that receives the charge for each 
    !> aqueous species. Dimension of the array is @a nspecies, where @a nspecies is is 
    !> the number of aqueous species (@ref getspeciescount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref setspeciessaveon,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: z(:)
    !> status = brm%setspeciessaveon(1)
    !> ncomps = brm%findcomponents()
    !> nspecies = brm%getspeciescount()
    !> status = brm%getspeciesz(z)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getspeciesz(self, z)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), DIMENSION(:), allocatable :: z
    getspeciesz = RM_getspeciesz(self%bmiphreeqcrm_id, z)
    END FUNCTION getspeciesz

    !> Returns an array with the starting cell numbers from the range of cell numbers 
    !> assigned to each worker.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param sc          Allocatable array to receive the starting cell numbers. 
    !> Dimension of the array is the number of threads (OpenMP) or the number of 
    !> processes (MPI).
    !> @retval IRESULT 0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref create,
    !> @ref getendcell,
    !> @ref getmpitasks,
    !> @ref getthreadcount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer, allocatable :: sc(:)
    !> status = brm%getstartcell(sc)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION getstartcell(self, sc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, intent(inout), DIMENSION(:), allocatable :: sc
    getstartcell = RM_getstartcell(self%bmiphreeqcrm_id, sc)
    RETURN
    END FUNCTION getstartcell

    !> Retrieves the surface names (such as "Hfo") that corresponds with
    !> the surface species names.
    !> The lists of surface species names and surface names are the same length.
    !> @ref findcomponents must be called before @ref getsurfacenames.
    !> This method may be useful when generating selected output definitions related 
    !> to surfaces.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array for the surface names.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getsurfacespeciescount, @ref getsurfacespeciesnames, @ref getsurfacetypes.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: surface_species_names(:), types(:), surface_names(:)
    !> status = brm%getsurfacespeciesnames(surface_species_names)
    !> status = brm%getsurfacetypes(types)
    !> status = brm%getsurfacenames(surface_names)
    !> do i = 1, brm%getsurfacespeciescount()
    !>   write(*,*) surface_species_names(i), types(i), surface_names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsurfacenames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsurfacenames = RM_getsurfacenames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsurfacenames

    !> Returns the number of surface species (such as "Hfo_wOH") found by call(s) to 
    !> @ref findcomponents.
    !> @ref findcomponents must be called before @ref getsurfacespeciescount.
    !> This method may be useful when generating selected output definitions related 
    !> to surfaces.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The number of surface species in the reaction module.
    !> @see
    !> @ref findcomponents,
    !> @ref getsurfacespeciesnames, @ref getsurfacetypes, @ref getsurfacenames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: surface_species_names(:), types(:), surface_names(:)
    !> status = brm%getsurfacespeciesnames(surface_species_names)
    !> status = brm%getsurfacetypes(types)
    !> status = brm%getsurfacenames(surface_names)
    !> do i = 1, brm%getsurfacespeciescount()
    !>   write(*,*) surface_species_names(i), types(i), surface_names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsurfacespeciescount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getsurfacespeciescount = RM_getsurfacespeciescount(self%bmiphreeqcrm_id)
    END FUNCTION getsurfacespeciescount

    !> Retrieves a list of surface species names.
    !> The list of surface species (for example, "Hfo_wOH") is derived from 
    !> the list of components (@ref findcomponents) and the list of all surface 
    !> types (such as "Hfo_w") that are found by call(s) to @ref findcomponents.
    !> @ref findcomponents must be called before @ref getsurfacespeciesnames.
    !> This method may be useful when generating selected output definitions related 
    !> to surfaces.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array to receive the surface species names. 
    !> Dimension of the array is set to @a nspecies, where @a nspecies is returned by
    !> @ref getspeciescount.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getsurfacespeciescount, @ref getsurfacetypes, @ref getsurfacenames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: surface_species_names(:), types(:), surface_names(:)
    !> status = brm%getsurfacespeciesnames(surface_species_names)
    !> status = brm%getsurfacetypes(types)
    !> status = brm%getsurfacenames(surface_names)
    !> do i = 1, brm%getsurfacespeciescount()
    !>   write(*,*) surface_species_names(i), types(i), surface_names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsurfacespeciesnames(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsurfacespeciesnames = RM_getsurfacespeciesnames(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsurfacespeciesnames

    !> Retrieves the surface site types (such as "Hfo_w") that correspond with
    !> the surface species names.
    !> The lists of surface species names and surface species types are the same length.
    !> @ref findcomponents must be called before @ref getsurfacetypes.
    !> This method may be useful when generating selected output definitions related to surfaces.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param names            Allocatable array to receive surface types.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getsurfacespeciescount, @ref getsurfacespeciesnames, @ref getsurfacenames.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable :: surface_species_names(:), types(:), surface_names(:)
    !> status = brm%getsurfacespeciesnames(surface_species_names)
    !> status = brm%getsurfacetypes(types)
    !> status = brm%getsurfacenames(surface_names)
    !> do i = 1, brm%getsurfacespeciescount()
    !>   write(*,*) surface_species_names(i), types(i), surface_names(i)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getsurfacetypes(self, names)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: names
    getsurfacetypes = RM_getsurfacetypes(self%bmiphreeqcrm_id, names)
    return 
    END FUNCTION getsurfacetypes

    !> Returns an array of temperatures (@a temperature) from the reaction module.
    !> Reactions do not change the temperature, so the temperatures are either the
    !> temperatures at initialization, or the values set with the last call to
    !> @ref settemperature.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param temperature      Allocatable array to receive the temperatures.
    !> Dimension of the array is set to @a nxyz, where @a nxyz is the number of 
    !> user grid cells (@ref getgridcellcount). Values for inactive cells are 
    !> set to 1e30.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref settemperature.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: temperature(:)
    !> status = brm%gettemperature(temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION gettemperature(self, temperature)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), dimension(:), allocatable :: temperature
    gettemperature = RM_gettemperature(self%bmiphreeqcrm_id, temperature)
    return
    END FUNCTION gettemperature

    !> Returns the number of threads, which is equal to the number of workers used 
    !> to run in parallel with OPENMP.
    !> For the OPENMP version, the number of threads is set implicitly or explicitly 
    !> with @ref create. For the MPI version, the number of threads is always one 
    !> for each process.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval           Number of threads, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getmpitasks.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string1, "(A,I)") "Number of threads: ", brm%getthreadcount()
    !> status = brm%outputmessage(string1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers; result is always 1.
    INTEGER FUNCTION getthreadcount(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    getthreadcount = RM_getthreadcount(self%bmiphreeqcrm_id)
    END FUNCTION getthreadcount

    !> Returns the current simulation time in seconds. The reaction module does not 
    !> change the time value, so the returned value is equal to the default (0.0) or 
    !> the last time set by @ref settime.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The current simulation time in seconds.
    !> @see
    !> @ref gettimeconversion,
    !> @ref gettimestep,
    !> @ref settime,
    !> @ref settimeconversion,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
    !>       brm%gettime() * brm%gettimeconversion(), " days"
    !> status = brm%logmessage(string)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    real(kind=8) FUNCTION gettime(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    gettime = RM_gettime(self%bmiphreeqcrm_id)
    END FUNCTION gettime

    !> Returns a multiplier to convert time from seconds to another unit, as 
    !> specified by the user. The reaction module uses seconds as the time unit. 
    !> The user can set a conversion factor (@ref settimeconversion) and retrieve 
    !> it with gettimeconversion. The reaction module only uses the conversion 
    !> factor when printing the long version of cell chemistry 
    !> (@ref setprintchemistryon), which is rare. Default conversion factor is 1.0.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval                 Multiplier to convert seconds to another time unit.
    !> @see
    !> @ref gettime,
    !> @ref gettimestep,
    !> @ref settime,
    !> @ref settimeconversion,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
    !>       brm%gettime() * brm%gettimeconversion(), " days"
    !> status = brm%logmessage(string)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    real(kind=8) FUNCTION gettimeconversion(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    gettimeconversion = RM_gettimeconversion(self%bmiphreeqcrm_id)
    END FUNCTION gettimeconversion

    !> Returns the current simulation time step in seconds.
    !> This is the time over which kinetic reactions are integrated in a 
    !> call to @ref runcells. The reaction module does not change the time 
    !> step value, so the returned value is equal to the default (0.0) or the 
    !> last time step set by @ref settimestep.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval            The current simulation time step in seconds.
    !> @see
    !> @ref gettime,
    !> @ref gettimeconversion,
    !> @ref settime,
    !> @ref settimeconversion,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string, "(A32,F15.1,A)") "          Time step             ", &
    !>       brm%gettimestep() * brm%gettimeconversion(), " days"
    !> status = brm%logmessage(string)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    real(kind=8) FUNCTION gettimestep(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    gettimestep = RM_gettimestep(self%bmiphreeqcrm_id)
    END FUNCTION gettimestep

    ! INTEGER FUNCTION GetVarItemsize(self, var)
    ! USE ISO_C_BINDING
    ! IMPLICIT NONE
    ! INTERFACE
    !     INTEGER(KIND=C_INT) FUNCTION RMF_GetVarItemsize(id, var) &
    !         BIND(C, NAME='RMF_GetVarItemsize')
    !     USE ISO_C_BINDING
    !     IMPLICIT NONE
    !     INTEGER(KIND=C_INT), INTENT(inout) :: id
    !     CHARACTER(KIND=C_CHAR), INTENT(inout) :: var(*)
    !     END FUNCTION RMF_GetVarItemsize
    ! END INTERFACE
    ! INTEGER, INTENT(inout) :: id
    ! CHARACTER(len=*), INTENT(inout) :: var
    ! GetVarItemsize = RMF_GetVarItemsize(self, trim(var)//C_NULL_CHAR)
    ! return
    ! END FUNCTION GetVarItemsize

    ! INTEGER FUNCTION GetVarNbytes(self, var)
    ! USE ISO_C_BINDING
    ! IMPLICIT NONE
    ! INTERFACE
    !     INTEGER(KIND=C_INT) FUNCTION RMF_GetVarNbytes(id, var) &
    !         BIND(C, NAME='RMF_GetVarNbytes')
    !     USE ISO_C_BINDING
    !     IMPLICIT NONE
    !     INTEGER(KIND=C_INT), INTENT(inout) :: id
    !     CHARACTER(KIND=C_CHAR), INTENT(inout) :: var(*)
    !     END FUNCTION RMF_GetVarNbytes
    ! END INTERFACE
    ! INTEGER, INTENT(inout) :: id
    ! CHARACTER(len=*), INTENT(inout) :: var
    ! GetVarNbytes = RMF_GetVarNbytes(id, trim(var)//C_NULL_CHAR)
    ! return
    ! END FUNCTION GetVarNbytes
    
    !> Transfer current viscosities to the array given in the argument list (@a viscosity).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param viscosity            Allocatable array to receive the viscosities. Dimension of 
    !> the array is @a nxyz, where @a nxyz is the number of user grid cells 
    !> (@ref getgridcellcount). Values for inactive cells are set to 1e30.
    !> @retval IRESULT          0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: viscosity(:)
    !> status = brm%getviscosity(viscosity)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION getviscosity(self, viscosity)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), intent(inout), dimension(:), allocatable :: viscosity
    getviscosity = RM_getviscosity(self%bmiphreeqcrm_id, viscosity)
    return
    END FUNCTION getviscosity

#ifdef USE_YAML
    !> A YAML file can be used to initialize an instance of PhreeqcRM.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param yaml_name         String containing the YAML file name.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @par
    !> The file contains a YAML map of PhreeqcRM methods
    !> and the arguments corresponding to the methods.
    !> Note that the PhreeqcRM methods do not have the "" prefix
    !> and the id argument is not included.
    !> For example,
    !> @par
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> loaddatabase: phreeqc.dat
    !> runfile:
    !> workers: true
    !> initial_phreeqc: true
    !> utility: true
    !> chemistry_name: advect.pqi
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par
    !> @ref initializeyaml will read the YAML file and execute the specified methods with
    !> the specified arguments. Using YAML
    !> terminology, the argument(s) for a method may be a scalar, a sequence, or a map,
    !> depending if the argument is
    !> a single item, a single vector, or there are multiple arguments.
    !> In the case of a map, the name associated
    !> with each argument (for example "chemistry_name" above) is arbitrary.
    !> The names of the map keys for map
    !> arguments are not used in parsing the YAML file; only the order of
    !> the arguments is important.
    !> @par
    !> The following list gives the PhreeqcRM methods that can be specified in a YAML file
    !> and the arguments that are required. The arguments are described with C++ formats, which
    !> are sufficient to identify which arguments are YAML scalars (single bool/logical,
    !> int, double, string/character argument),
    !> sequences (single vector argument), or maps (multiple arguments).
    !> @par
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> closefiles(void);
    !> createmapping(std::vector< int >& grid2chem);
    !> dumpmodule();
    !> findcomponents();
    !> initialequilibriumphases2module(std::vector< int > equilibrium_phases);
    !> initialexchanges2module(std::vector< int > exchanges);
    !> initialgasphases2module(std::vector< int > gas_phases);
    !> InitialKineticss2Module(std::vector< int > kinetics);
    !> initialsolidsolutions2module(std::vector< int > solid_solutions);
    !> initialsolutions2module(std::vector< int > solutions);
    !> initialsurfaces2module(std::vector< int > surfaces);
    !> initialphreeqc2module(std::vector< int > initial_conditions1);
    !> initialphreeqc2module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
    !> initialphreeqccell2module(int n, std::vector< int > cell_numbers);
    !> loaddatabase(std::string database);
    !> openfiles(void);
    !> outputmessage(std::string str);
    !> runcells(void);
    !> runfile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
    !> runstring(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    !> screenmessage(std::string str);
    !> setcomponenth2o(bool tf);
    !> setconcentrations(std::vector< double > c);
    !> setcurrentselectedoutputusernumber(int n_user);
    !> setdensityuser(std::vector< double > density);
    !> setdumpfilename(std::string dump_name);
    !> seterrorhandlermode(int mode);
    !> seterroron(bool tf);
    !> setfileprefix(std::string prefix);
    !> setgascompmoles(std::vector< double > gas_moles);
    !> setgasphasevolume(std::vector< double > gas_volume);
    !> setpartitionuzsolids(bool tf);
    !> setporosity(std::vector< double > por);
    !> setpressure(std::vector< double > p);
    !> setprintchemistrymask(std::vector< int > cell_mask);
    !> setprintchemistryon(bool workers, bool initial_phreeqc, bool utility);
    !> setrebalancebycell(bool tf);
    !> setrebalancefraction(double f);
    !> setrepresentativevolume(std::vector< double > rv);
    !> setsaturationuser(std::vector< double > sat);
    !> setscreenon(bool tf);
    !> setselectedoutputon(bool tf);
    !> setspeciessaveon(bool save_on);
    !> settemperature(std::vector< double > t);
    !> settime(double time);
    !> settimeconversion(double conv_factor);
    !> settimestep(double time_step);
    !> setunitsexchange(int option);
    !> setunitsgasphase(int option);
    !> setunitskinetics(int option);
    !> setunitsppassemblage(int option);
    !> setunitssolution(int option);
    !> setunitsssassemblage(int option);
    !> setunitssurface(int option);
    !> speciesconcentrations2module(std::vector< double > species_conc);
    !> statesave(int istate);
    !> stateapply(int istate);
    !> statedelete(int istate);
    !> usesolutiondensityvolume(bool tf);
    !> warningmessage(std::string warnstr);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !>
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !>    id = create(nxyz, MPI_COMM_WORLD)
    !>    status = brm%initializeyaml("myfile.yaml")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initializeyaml(self, yaml_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: yaml_name
    initializeyaml = RM_initializeyaml(self%bmiphreeqcrm_id, yaml_name)
    END FUNCTION initializeyaml
#endif

    !> Fills an array (@a bc_conc) with concentrations from solutions in the 
    !> InitialPhreeqc instance.
    !> The method is used to obtain concentrations for boundary conditions. 
    !> If a negative value is used for a cell in @a bc1, then the highest numbered 
    !> solution in the InitialPhreeqc instance will be used for that cell. Concentrations 
    !> may be a mixture of two solutions, @a bc1 and @a bc2, with a mixing fraction 
    !> for @a bc1 1 of @a f1 and mixing fraction for @a bc2 of (1 - @a f1).
    !> A negative value for @a bc2 implies no mixing, and the associated value for 
    !> @a f1 is ignored. If @a bc2 and @a f1 are omitted, no mixing is used; 
    !> concentrations are derived from @a bc1 only.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param bc_conc       Allocatable array of concentrations extracted from the 
    !> InitialPhreeqc instance. The dimension of @a bc_conc is set to (@a n_boundary, @a ncomp),
    !> where @a ncomp is the number of components returned from @ref findcomponents 
    !> or @ref getcomponentcount.
    !> @param n_boundary    The number of boundary condition solutions that need to be filled.
    !> @param bc1           Array of solution index numbers that refer to solutions in the 
    !> InitialPhreeqc instance. Size is @a n_boundary.
    !> @param bc2           Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance
    !> and are defined to mix with @a bc1.
    !> Size is @a n_boundary. This argument and @a f1 are optional.
    !> @param f1            Array with fraction of @a bc1 that mixes with (1-@a f1) of @a bc2.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getcomponentcount.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: bc_conc(:,:)
    !> integer, allocatable :: bc1(:), bc2(:)
    !> real(kind=8), allocatable :: f1(:)
    !> nbound = 1
    !> allocate(bc1(nbound), bc2(nbound), f1(nbound))
    !> bc1 = 0           ! solution 0 from InitialPhreeqc instance
    !> bc2 = -1          ! no bc2 solution for mixing
    !> f1 = 1.0          ! mixing fraction for bc1
    !> status = brm%initialphreeqc2concentrations(bc_conc, nbound, bc1, bc2, f1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION initialphreeqc2concentrations(self, bc_conc, n_boundary, bc1, bc2, f1)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(INOUT), DIMENSION(:,:), allocatable :: bc_conc
    INTEGER, INTENT(in) :: n_boundary
    INTEGER, INTENT(in), DIMENSION(:) :: bc1
    INTEGER, INTENT(in), DIMENSION(:) , OPTIONAL :: bc2
    real(kind=8), INTENT(in), DIMENSION(:) , OPTIONAL :: f1
    initialphreeqc2concentrations = RM_initialphreeqc2concentrations(self%bmiphreeqcrm_id, bc_conc, n_boundary, bc1, bc2, f1)
    END FUNCTION initialphreeqc2concentrations

    !> Transfer solutions and reactants from the InitialPhreeqc instance to the 
    !> reaction-module workers, possibly with mixing. In its simplest form, @a ic1 is 
    !> used to select initial conditions, including solutions and reactants,
    !> for each cell of the model, without mixing. @a ic1 is dimensioned (@a nxyz, 7), 
    !> where @a nxyz is the number of grid cells in the user's model
    !> (@ref getgridcellcount). The dimension of 7 refers to solutions and reactants 
    !> in the following order: (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, 
    !> (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS. In Fortran, 
    !> ic1(100, 4) = 2, indicates that cell 99 (0 based) contains the SURFACE definition 
    !> with user number 2 that has been defined in the InitialPhreeqc instance (either 
    !> by @ref runfile or @ref runstring).
    !> @n@n
    !> It is also possible to mix solutions and reactants to obtain the initial conditions 
    !> for cells. For mixing, @a ic2 contains numbers for a second entity that mixes with 
    !> the entity defined in @a ic1. @a F1 contains the mixing fraction for @a ic1, whereas 
    !> (1 - @a f1) is the mixing fraction for @a ic2. In Fortran, ic1(100, 4) = 2, 
    !> initial_conditions2(100, 4) = 3, f1(100, 4) = 0.25 indicates that cell 99 (0 based) 
    !> contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface
    !> compositions have been defined in the InitialPhreeqc instance. If the user number 
    !> in @a ic2 is negative, no mixing occurs. If @a ic2 and @a f1 are omitted,
    !> no mixing is used, and initial conditions are derived solely from @a ic1.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param ic1            Array of solution and reactant index numbers that refer to 
    !> definitions in the InitialPhreeqc instance. Size is (@a nxyz,7). The order of 
    !> definitions is given above. Negative values are ignored, resulting in no definition 
    !> of that entity for that cell.
    !> @param ic2            Array of solution and reactant index numbers that refer to 
    !> definitions in the InitialPhreeqc instance. Nonnegative values of @a ic2 result in 
    !> mixing with the entities defined in @a ic1. Negative values result in no mixing.
    !> Size is (@a nxyz,7). The order of definitions is given above.
    !> Optional in Fortran; omitting results in no mixing.
    !> @param f1           Fraction of ic1 that mixes with (1-@a f1) of ic2.
    !> Size is (nxyz,7). The order of definitions is given above.
    !> Optional in Fortran; omitting results in no mixing.
    !> @retval IRESULT          0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
    !> ic1 = -1
    !> ic2 = -1
    !> f1 = 1.0
    !> do i = 1, nxyz
    !>   ic1(i,1) = 1       ! Solution 1
    !>   ic1(i,2) = -1      ! Equilibrium phases none
    !>   ic1(i,3) = 1       ! Exchange 1
    !>   ic1(i,4) = -1      ! Surface none
    !>   ic1(i,5) = -1      ! Gas phase none
    !>   ic1(i,6) = -1      ! Solid solutions none
    !>   ic1(i,7) = -1      ! Kinetics none
    !> enddo
    !> status = brm%initialphreeqc2module(ic1, ic2, f1)1))
    !> ! No mixing is defined, so the following is equivalent
    !> status = brm%initialphreeqc2module(ic1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialphreeqc2module(self, ic1, ic2, f1)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(inout), DIMENSION(:,:) :: ic1
    INTEGER, INTENT(in), DIMENSION(:,:), OPTIONAL :: ic2
    real(kind=8), INTENT(inout), DIMENSION(:,:), OPTIONAL :: f1
    initialphreeqc2module = RM_initialphreeqc2module(self%bmiphreeqcrm_id, ic1, ic2, f1)
    END FUNCTION initialphreeqc2module

    !> Transfer SOLUTION definitions from the InitialPhreeqc instance to the reaction-
    !> module workers. 
    !> @a solutions is used to select SOLUTION definitions for each 
    !> cell of the model. @a solutions is dimensioned @a nxyz, where @a nxyz is the 
    !> number of grid cells in the  user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param solutions    Array of SOLUTION index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values 
    !> are ignored, resulting in no transfer of a SOLUTION definition for that cell.
    !> (Note that all cells must have a SOLUTION definition, which could be defined 
    !> by other calls to @a initialsolutions2module, @ref initialphreeqc2module, 
    !> or @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see     
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(solutions(nxyz))
    !> solutions = 1
    !> status = brm%initialsolutions2module(solutions);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialsolutions2module(self, solutions)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: solutions
    initialsolutions2module = RM_initialsolutions2module(self%bmiphreeqcrm_id, solutions)
    END FUNCTION initialsolutions2module  

    !> Transfer EQUILIBRIUM_PHASES definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a equilibrium_phases is used to select EQUILIBRIUM_PHASES definitions for each 
    !> cell of the model. @a equilibrium_phases is dimensioned @a nxyz, where @a nxyz is 
    !> the number of grid cells in the  user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param equilibrium_phases Array of EQUILIBRIUM_PHASES index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values are 
    !> ignored, resulting in no transfer of an EQUILIBRIUM_PHASES definition for that cell.
    !> (Note that an EQUILIBRIUM_PHASES definition for a cell could be defined by other 
    !> calls to @a initialequilibriumphases2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT    0 is success, negative is failure (See @ref decodeerror).
    !> @see 
    !> @ref initialsolutions2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(equilibrium_phases(nxyz))
    !> equilibrium_phases = 1
    !> status = brm%initialequilibriumphases2module(equilibrium_phases);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialequilibriumphases2module(self, equilibrium_phases)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: equilibrium_phases
    initialequilibriumphases2module = RM_initialequilibriumphases2module(self%bmiphreeqcrm_id, equilibrium_phases)
    END FUNCTION initialequilibriumphases2module

    !> Transfer EXCHANGE definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a exchanges is used to select EXCHANGE definitions for each cell of the model.
    !> @a exchanges is dimensioned @a nxyz, where @a nxyz is the number of grid cells 
    !> in the user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param exchanges    Vector of EXCHANGE index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values 
    !> are ignored, resulting in no transfer of an EXCHANGE definition for that cell.
    !> (Note that an EXCHANGE definition for a cell could be defined by other 
    !> calls to @a initialexchanges2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see 
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(exchanges(nxyz))
    !> exchanges = 1
    !> status = brm%initialexchanges2module(exchanges);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialexchanges2module(self, exchanges)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: exchanges
    initialexchanges2module = RM_initialexchanges2module(self%bmiphreeqcrm_id, exchanges)
    END FUNCTION initialexchanges2module

    !> Transfer GAS_PHASE definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a gas_phases is used to select GAS_PHASE definitions for each cell of the model.
    !> @a gas_phases is dimensioned @a nxyz, where @a nxyz is the number of grid cells 
    !> in the user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_phases   Vector of GAS_PHASE index numbers that refer to
    !> definitions in the InitialPhreeqc instance.Size is @a nxyz. Negative values are 
    !> ignored, resulting in no transfer of a GAS_PHASE definition for that cell.
    !> (Note that an GAS_PHASE definition for a cell could be defined by other 
    !> calls to @a initialgasphases2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see 
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(gas_phases(nxyz))
    !> gas_phases = 1
    !> status = brm%initialgasphases2module(gas_phases);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialgasphases2module(self, gas_phases)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: gas_phases
    initialgasphases2module = RM_initialgasphases2module(self%bmiphreeqcrm_id, gas_phases)
    END FUNCTION initialgasphases2module

    !> Transfer SOLID_SOLUTIONS definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a solid_solutions is used to select SOLID_SOLUTIONS definitions for each cell 
    !> of the model. @a solid_solutions is dimensioned @a nxyz, where @a nxyz is the 
    !> number of grid cells in the user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param solid_solutions Array of SOLID_SOLUTIONS index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values 
    !> are ignored, resulting in no transfer of a SOLID_SOLUTIONS definition for that cell.
    !> (Note that an SOLID_SOLUTIONS definition for a cell could be defined by other 
    !> calls to @a initialsolidsolutions2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see  
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(solid_solutions(nxyz))
    !> solid_solutions = 1
    !> status = brm%initialsolidsolutions2module(solid_solutions);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialsolidsolutions2module(self, solid_solutions)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: solid_solutions
    initialsolidsolutions2module = RM_initialsolidsolutions2module(self%bmiphreeqcrm_id, solid_solutions)
    END FUNCTION initialsolidsolutions2module

    !> Transfer SURFACE definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a surfaces is used to select SURFACE definitions for each cell of the model.
    !> @a surfaces is dimensioned @a nxyz, where @a nxyz is the number of grid cells 
    !> in the user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param surfaces    Array of SURFACE index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values 
    !> are ignored, resulting in no transfer of a SURFACE definition for that cell.
    !> (Note that an SURFACE definition for a cell could be defined by other 
    !> calls to @a initialsurfaces2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see 
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialkinetics2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(surfaces(nxyz))
    !> surfaces = 1
    !> status = brm%initialsurfaces2module(surfaces);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialsurfaces2module(self, surfaces)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: surfaces
    initialsurfaces2module = RM_initialsurfaces2module(self%bmiphreeqcrm_id, surfaces)
    END FUNCTION initialsurfaces2module
    
    !> Transfer KINETICS definitions from the InitialPhreeqc instance to the 
    !> reaction-module workers.
    !> @a kinetics is used to select KINETICS definitions for each cell of the model.
    !> @a kinetics is dimensioned @a nxyz, where @a nxyz is the number of grid cells in the 
    !> user's model (@ref getgridcellcount).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param kinetics    Array of KINETICS index numbers that refer to
    !> definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values are 
    !> ignored, resulting in no transfer of a KINETICS definition for that cell.
    !> (Note that an KINETICS definition for a cell could be defined by other 
    !> calls to @a initialkinetics2module, @ref initialphreeqc2module, or 
    !> @ref initialphreeqccell2module.)
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see 
    !> @ref initialsolutions2module,
    !> @ref initialequilibriumphases2module,
    !> @ref initialexchanges2module,
    !> @ref initialgasphases2module,
    !> @ref initialsolidsolutions2module,
    !> @ref initialsurfaces2module,
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> dimension(kinetics(nxyz))
    !> kinetics = 1
    !> status = brm%initialkinetics2module(kinetics);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialkinetics2module(self, kinetics)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in), DIMENSION(:) :: kinetics
    initialkinetics2module = RM_initialkinetics2module(self%bmiphreeqcrm_id, kinetics)
    END FUNCTION initialkinetics2module  
    
    !> Fills an array (@a bc_conc) with aqueous species concentrations from solutions 
    !> in the InitialPhreeqc instance.
    !> This method is intended for use with multicomponent-diffusion transport calculations,
    !> and @ref setspeciessaveon must be set to @a true.
    !> The method is used to obtain aqueous species concentrations for boundary conditions. 
    !> If a negative value is used for a cell in @a bc1, then the highest numbered solution 
    !> in the InitialPhreeqc instance will be used for that cell. Concentrations may be a 
    !> mixture of two solutions, @a bc1 and @a bc2, with a mixing fraction for @a bc1 of
    !> @a f1 and mixing fraction for @a bc2 of (1 - @a f1). A negative value for @a bc2 
    !> implies no mixing, and the associated value for @a f1 is ignored. If @a bc2 and 
    !> @a f1 are omitted, no mixing is used; concentrations are derived from @a bc1 only.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param bc_conc        Array of aqueous concentrations extracted from the 
    !> InitialPhreeqc instance.
    !> The dimension of @a species_c is set to (@a n_boundary, @a nspecies), where 
    !> @a nspecies is the number of aqueous species returned from @ref getspeciescount.
    !> @param n_boundary     The number of boundary condition solutions that need to be filled.
    !> @param bc1            Array of solution index numbers that refer to solutions in 
    !> the InitialPhreeqc instance. Size is @a n_boundary.
    !> @param bc2            Array of solution index numbers that that refer to solutions 
    !> in the InitialPhreeqc instance and are defined to mix with @a bc1. Size is 
    !> @a n_boundary. Optional in Fortran.
    !> @param f1             Fraction of @a bc1 that mixes with (1-@a f1) of @a bc2.
    !> Size is @a n_boundary. Optional in Fortran.
    !> @retval IRESULT    0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciescount,
    !> @ref setspeciessaveon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: bc_conc(:,:)
    !> nbound = 1
    !> allocate(bc1(nbound), bc2(nbound), f1(nbound))
    !> bc1 = 0           ! solution 0 from InitialPhreeqc instance
    !> bc2 = -1          ! no bc2 solution for mixing
    !> f1 = 1.0          ! mixing fraction for bc1
    !> status = brm%initialphreeqc2speciesconcentrations(bc_conc, nbound, bc1, bc2, f1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION initialphreeqc2speciesconcentrations(self, bc_conc, n_boundary, bc1, bc2, f1)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:,:), intent(inout), allocatable :: bc_conc
    INTEGER, INTENT(in) :: n_boundary
    INTEGER, INTENT(in), DIMENSION(:) :: bc1
    INTEGER, INTENT(in), DIMENSION(:), OPTIONAL :: bc2
    real(kind=8), INTENT(in), DIMENSION(:), OPTIONAL :: f1
    initialphreeqc2speciesconcentrations = &
        RM_initialphreeqc2speciesconcentrations(self%bmiphreeqcrm_id, bc_conc, n_boundary, bc1, bc2, f1)
    END FUNCTION initialphreeqc2speciesconcentrations
    !> A cell numbered @a n_user in the InitialPhreeqc instance is selected to populate 
    !> a series of cells.
    !> All reactants with the number @a n_user are transferred along with the solution.
    !> If MIX @a n_user exists, it is used for the definition of the solution.
    !> If @a n_user is negative, @a n_user is redefined to be the largest solution or 
    !> MIX number in the InitialPhreeqc instance. All reactants for each cell in the 
    !> list @a cell_numbers are removed before the cell definition is copied from the 
    !> InitialPhreeqc instance to the workers.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param n_user           Cell number refers to a solution or MIX and associated 
    !> reactants in the InitialPhreeqc instance. A negative number indicates the largest 
    !> solution or MIX number in the InitialPhreeqc instance will be used.
    !> @param cell_numbers     A list of cell numbers in the user's grid-cell numbering 
    !> system that will be populated with cell @a n_user from the InitialPhreeqc instance.
    !> @param n_cell           The number of cell numbers in the @a cell_numbers list.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate (cell_numbers(2))
    !> cell_numbers(1) = 18
    !> cell_numbers(2) = 19
    !> ! n will be the largest SOLUTION number in InitialPhreeqc instance
    !> ! copies solution and reactants to cells 18 and 19
    !> status = brm%initialphreeqccell2module(-1, cell_numbers, 2)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION initialphreeqccell2module(self, n_user, cell_numbers, n_cell)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: n_user
    INTEGER, INTENT(in), DIMENSION(:) :: cell_numbers
    INTEGER, INTENT(in) :: n_cell
    initialphreeqccell2module = RM_initialphreeqccell2module(self%bmiphreeqcrm_id, n_user, cell_numbers, n_cell)
    END FUNCTION initialphreeqccell2module

    !> Load a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. 
    !> All definitions of the reaction module are cleared (SOLUTION_SPECIES, PHASES, 
    !> SOLUTIONs, etc.), and the database is read.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param db_name          String containing the database name.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref create.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%loaddatabase("phreeqc.dat")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION loaddatabase(self, db_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: db_name
    loaddatabase = RM_loaddatabase(self%bmiphreeqcrm_id, db_name)
    END FUNCTION loaddatabase

    !> Print a message to the log file.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param str              String to be printed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref errormessage,
    !> @ref openfiles,
    !> @ref outputmessage,
    !> @ref screenmessage,
    !> @ref warningmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
    !>       brm%gettime() * brm%gettimeconversion(), " days"
    !> status = brm%logmessage(string)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION logmessage(self, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: str
    logmessage = RM_logmessage(self%bmiphreeqcrm_id, str)
    END FUNCTION logmessage

    !> MPI only. Workers (processes with @ref getmpimyself > 0) must call 
    !> mpiworker to be able to respond to messages from the root to accept data, 
    !> perform calculations, and (or) return data. mpiworker contains a loop that 
    !> reads a message from root, performs a task, and waits for another message from 
    !> root. @ref setconcentrations, @ref runcells, and @ref getconcentrations
    !> are examples of methods that send a message from root to get the workers to 
    !> perform a task. The workers will respond to all methods that are designated 
    !> "workers must be in the loop of mpiworker" in the MPI section of the method 
    !> documentation. The workers will continue to respond to messages from root until 
    !> root calls @ref mpiworkerbreak.
    !> @n@n
    !> (Advanced) The list of tasks that the workers perform can be extended by using 
    !> @ref SetmpiworkerCallback. It is then possible to use the MPI processes to 
    !> perform other developer-defined tasks, such as transport calculations, without
    !> exiting from the mpiworker loop. Alternatively, root calls @ref mpiworkerbreak 
    !> to allow the workers to continue past a call to mpiworker. The workers perform 
    !> developer-defined calculations, and then mpiworker is called again to respond to
    !> requests from root to perform reaction-module tasks.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror). 
    !> mpiworker returns a value only when @ref mpiworkerbreak is called by root.
    !> @see
    !> @ref mpiworkerbreak,
    !> @ref SetmpiworkerCallback.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%mpiworker()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by all workers.
    INTEGER FUNCTION mpiworker(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    mpiworker = RM_mpiworker(self%bmiphreeqcrm_id)
    END FUNCTION mpiworker

    !> MPI only. This method is called by root to force workers (processes with 
    !>@ref getmpimyself > 0) to return from a call to @ref mpiworker.
    !> @ref mpiworker contains a loop that reads a message from root, performs a
    !> task, and waits for another message from root. The workers respond to all 
    !> methods that are designated "workers must be in the loop of mpiworker" in 
    !> the MPI section of the method documentation. The workers will continue to 
    !> respond to messages from root until root calls mpiworkerbreak.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref mpiworker,
    !> @ref SetmpiworkerCallback.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%mpiworkerbreak()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION mpiworkerbreak(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    mpiworkerbreak = RM_mpiworkerbreak(self%bmiphreeqcrm_id)
    END FUNCTION mpiworkerbreak

    !> Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
    !> based on the prefix defined by @ref setfileprefix.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref closefiles,
    !> @ref errormessage,
    !> @ref getfileprefix,
    !> @ref logmessage,
    !> @ref outputmessage,
    !> @ref setfileprefix,
    !> @ref warningmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setfileprefix("Advect_f90")
    !> status = brm%openfiles()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION openfiles(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    openfiles = RM_openfiles(self%bmiphreeqcrm_id)
    END FUNCTION openfiles

    !> Print a message to the output file.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param str              String to be printed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref errormessage,
    !> @ref logmessage,
    !> @ref screenmessage,
    !> @ref warningmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string1, "(A,I10)") "Number of threads:                                ", brm%getthreadcount()
    !> status = brm%outputmessage(string1)
    !> write(string1, "(A,I10)") "Number of MPI processes:                          ", brm%getmpitasks()
    !> status = brm%outputmessage(string1)
    !> write(string1, "(A,I10)") "MPI task number:                                  ", brm%getmpimyself()
    !> status = brm%outputmessage(string1)
    !> status = brm%getfileprefix(string)
    !> write(string1, "(A,A)") "File prefix:                                        ", string
    !> status = brm%outputmessage(trim(string1))
    !> write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", brm%getgridcellcount()
    !> status = brm%outputmessage(trim(string1))
    !> write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", brm%getchemistrycellcount()
    !> status = brm%outputmessage(trim(string1))
    !> write(string1, "(A,I10)") "Number of components for transport:               ", brm%getcomponentcount()
    !> status = brm%outputmessage(trim(string1))
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION outputmessage(self, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: str
    outputmessage = RM_outputmessage(self%bmiphreeqcrm_id, str)
    END FUNCTION outputmessage

    !> Runs a reaction step for all of the cells in the reaction module.
    !> Normally, tranport concentrations are transferred to the reaction cells 
    !> (@ref setconcentrations) before reaction calculations are run. The 
    !> length of time over which kinetic reactions are integrated is set
    !> by @ref settimestep. Other properties that may need to be updated 
    !> as a result of the transport calculations include porosity (@ref setporosity), 
    !> saturation (@ref setsaturationuser), temperature (@ref settemperature), 
    !> and pressure (@ref setpressure).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setconcentrations,
    !> @ref setporosity,
    !> @ref setpressure,
    !> @ref setsaturationuser,
    !> @ref settemperature,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setporosity(por)                ! If pore volume changes
    !> status = brm%setsaturationuser(sat)          ! If saturation changes
    !> status = brm%settemperature(temperature)     ! If temperature changes
    !> status = brm%setpressure(pressure)           ! If pressure changes
    !> status = brm%setconcentrations(c)            ! Transported concentrations
    !> status = brm%settimestep(time_step)          ! Time step for kinetic reactions
    !> status = brm%runcells()
    !> status = brm%getconcentrations(c)            ! Concentrations after reaction
    !> status = brm%getdensitycalculated(density)   ! Density after reaction
    !> status = brm%getsolutionvolume(volume)       ! Solution volume after reaction
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION runcells(self)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    runcells = RM_runcells(self%bmiphreeqcrm_id)
    END FUNCTION runcells

    !> Run a PHREEQC input file. The first three arguments determine which 
    !> IPhreeqc instances will run the file--the workers, the InitialPhreeqc instance, 
    !> and (or) the Utility instance. Input files that modify the thermodynamic database 
    !> should be run by all three sets of instances. Files with SELECTED_OUTPUT 
    !> definitions that will be used during the time-stepping loop need to be run by 
    !> the workers. Files that contain initial conditions or boundary conditions should
    !> be run by the InitialPhreeqc instance.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param workers          1, the workers will run the file; 0, the workers will 
    !> not run the file.
    !> @param initial_phreeqc  1, the InitialPhreeqc instance will run the file; 0, the 
    !> InitialPhreeqc will not run the file.
    !> @param utility          1, the Utility instance will run the file; 0, the Utility 
    !> instance will not run the file.
    !> @param chem_name        Name of the file to run.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref runstring.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%runfile(1, 1, 1, "advect.pqi")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION runfile(self, workers, initial_phreeqc, utility, chem_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: workers, initial_phreeqc, utility
    CHARACTER(len=*), INTENT(in) :: chem_name
    runfile = RM_runfile(self%bmiphreeqcrm_id, workers, initial_phreeqc, utility, chem_name)
    END FUNCTION runfile

    !> Run a PHREEQC input string. The first three arguments determine which
    !> IPhreeqc instances will run
    !> the string--the workers, the InitialPhreeqc instance, and (or) the Utility 
    !> instance. Input strings that modify the thermodynamic database should be run 
    !> by all three sets of instances. Strings with SELECTED_OUTPUT definitions that 
    !> will be used during the time-stepping loop need to be run by the workers. 
    !> Strings that contain initial conditions or boundary conditions should
    !> be run by the InitialPhreeqc instance.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param workers          1, the workers will run the string; 0, the workers will 
    !> not run the string.
    !> @param initial_phreeqc  1, the InitialPhreeqc instance will run the string; 0, 
    !> the InitialPhreeqc will not run the string.
    !> @param utility          1, the Utility instance will run the string; 0, the 
    !> Utility instance will not run the string.
    !> @param input_string     String containing PHREEQC input.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref runfile.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> string = "DELETE; -all"
    !> status = brm%runstring(1, 0, 1, string)  ! workers, initial_phreeqc, utility
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION runstring(self, workers, initial_phreeqc, utility, input_string)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: initial_phreeqc, workers, utility
    CHARACTER(len=*), INTENT(in) :: input_string
    runstring = RM_runstring(self%bmiphreeqcrm_id, workers, initial_phreeqc, utility, input_string)
    END FUNCTION runstring

    !> Print message to the screen.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param str              String to be printed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref errormessage,
    !> @ref logmessage,
    !> @ref outputmessage,
    !> @ref warningmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> write(string, "(A32,F15.1,A)") "Beginning reaction calculation  ", &
    !>       time * brm%gettimeconversion(), " days"
    !> status = brm%screenmessage(string)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION screenmessage(self, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: str
    screenmessage = RM_screenmessage(self%bmiphreeqcrm_id, str)
    END FUNCTION screenmessage

    !> Select whether to include H2O in the component list.
    !> The concentrations of H and O must be known accurately (8 to 10 
    !> significant digits) for the numerical method of PHREEQC to produce 
    !> accurate pH and pe values. Because most of the H and O are in the water 
    !> species, it may be more robust (require less accuracy in transport) to
    !> transport the excess H and O (the H and O not in water) and water.
    !> The default setting (@a true) is to include water, excess H, and excess O 
    !> as components. A setting of @a false will include total H and total O 
    !> as components. @a setcomponenth2o must be called before 
    !> @ref findcomponents.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf               0, total H and O are included in the component list; 1, 
    !> excess H, excess O, and water are included in the component list.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setcomponenth2o(0)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setcomponenth2o(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: tf
    setcomponenth2o = RM_setcomponenth2o(self%bmiphreeqcrm_id, tf)
    END FUNCTION setcomponenth2o

    !> Use the array of concentrations (@a c) to set the moles of components in each 
    !> reaction cell.
    !> The volume of water in a cell is the product of porosity (@ref setporosity), 
    !> saturation (@ref setsaturationuser), and reference volume 
    !> (@ref setrepresentativevolume). The moles of each component are determined 
    !> by the volume of water and per liter concentrations. If concentration units 
    !> (@ref setunitssolution) are mass fraction, the density (as specified by
    !>  @ref setdensityuser) is used to convert from mass fraction to per mass per liter.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param c                Array of component concentrations. Size of array is 
    !> (@a nxyz, @a ncomps), where @a nxyz is the number of grid cells in the user's 
    !> model (@ref getgridcellcount), and @a ncomps is the number of components as 
    !> determined by @ref findcomponents or @ref getcomponentcount.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setdensityuser,
    !> @ref setporosity,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser,
    !> @ref setunitssolution.
    !>
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(c(nxyz, ncomps))
    !> ...
    !> call advect_f90(c, bc_conc, ncomps, nxyz)
    !> status = brm%setporosity(por)               ! If porosity changes
    !> status = brm%setsaturationuser(sat)         ! If saturation changes
    !> status = brm%settemperature(temperature))   ! If temperature changes
    !> status = brm%setpressure(pressure)          ! If pressure changes
    !> status = brm%setconcentrations(c)           ! Transported concentrations
    !> status = brm%settimestep(time_step)         ! Time step for kinetic reactions
    !> status = brm%settime(time)                  ! Current time
    !> status = brm%runcells()
    !> status = brm%getconcentrations(c)           ! Concentrations after reaction
    !> status = brm%getdensitycalculated(density)  ! Density after reaction
    !> status = brm%getsolutionvolume(volume)      ! Solution volume after reaction
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setconcentrations(self, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:,:), INTENT(in) :: c
    setconcentrations = RM_setconcentrations(self%bmiphreeqcrm_id, c)
    END FUNCTION setconcentrations

    !> Select the current selected output by user number. 
    !> The user may define  multiple SELECTED_OUTPUT data blocks for the workers. 
    !> A user number is specified for each data block. The value of the argument 
    !> @a n_user selects which of the SELECTED_OUTPUT definitions will be used
    !> for selected-output operations.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setnthselectedoutput,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), allocatable :: selected_out(:,:)
    !> do isel = 1, brm%getselectedoutputcount()
    !>   n_user = brm%getnthselectedoutputusernumber(isel)
    !>   status = brm%setcurrentselectedoutputusernumber(n_user)
    !>   col = brm%getselectedoutputcolumncount()
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION setcurrentselectedoutputusernumber(self, n_user)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: n_user
    setcurrentselectedoutputusernumber = RM_setcurrentselectedoutputusernumber(self%bmiphreeqcrm_id, n_user)
    END FUNCTION setcurrentselectedoutputusernumber

    !> Set the density used for units conversion. 
    !> These density values are used when converting from transported mass fraction 
    !> concentrations (@ref setunitssolution) to produce per liter concentrations 
    !> during a call to @ref setconcentrations. They are also used when converting 
    !> from module concentrations to transport concentrations of mass fraction 
    !> (@ref getconcentrations), if @ref usesolutiondensityvolume is set to @a false.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param density          Array of densities. Size of array is @a nxyz, where @a nxyz 
    !> is the number of grid cells in the user's model (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getconcentrations,
    !> @ref setconcentrations,
    !> @ref setunitssolution,
    !> @ref usesolutiondensityvolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(density(nxyz))
    !> density = 1.0
    !> status = brm%setdensityuser(density)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setdensityuser(self, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: density
    setdensityuser = RM_setdensityuser(self%bmiphreeqcrm_id, density)
    END FUNCTION setdensityuser
    
    INTEGER FUNCTION SetDensity(self, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: density
    SetDensity = RM_setdensityuser(self%bmiphreeqcrm_id, density)
    END FUNCTION SetDensity

    !> Set the name of the dump file. It is the name used by @ref dumpmodule.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param dump_name        Name of dump file.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref dumpmodule.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setdumpfilename("advection_f90.dmp")
    !> dump_on = 1
    !> append = 0
    !> status = brm%dumpmodule(dump_on, append)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION setdumpfilename(self, dump_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: dump_name
    setdumpfilename = RM_setdumpfilename(self%bmiphreeqcrm_id, dump_name)
    END FUNCTION setdumpfilename

    !> Set the action to be taken when the reaction module encounters an error.
    !> Options are 0, return to calling program with an error return code (default);
    !> 1, throw an exception, in C++, the exception can be caught, for C and Fortran, 
    !> the program will exit; or 2, attempt to exit gracefully.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param mode             Error handling mode: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> id = create(nxyz, nthreads)
    !> status = brm%seterrorhandlermode(2)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION seterrorhandlermode(self, mode)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: mode
    seterrorhandlermode = RM_seterrorhandlermode(self%bmiphreeqcrm_id, mode)
    END FUNCTION seterrorhandlermode

    !> Set the property that controls whether error messages are generated and displayed.
    !> Messages include PHREEQC "ERROR" messages, and any messages written with 
    !> @ref errormessage.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf           @a 1, enable error messages; @a 0, disable error messages. 
    !> Default is 1.
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref errormessage,
    !> @ref screenmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%seterroron(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION seterroron(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: tf
    seterroron = RM_seterroron(self%bmiphreeqcrm_id, tf)
    END FUNCTION seterroron

    !> Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
    !> These files are opened by @ref openfiles.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param prefix           Prefix used when opening the output and log files.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref openfiles,
    !> @ref closefiles.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setfileprefix("Advect_f90")
    !> status = brm%openfiles()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION setfileprefix(self, prefix)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: prefix
    setfileprefix = RM_setfileprefix(self%bmiphreeqcrm_id, prefix)
    END FUNCTION setfileprefix

    !> Use the array of concentrations (@a gas_moles) to set the moles of
    !> gas components in each reaction cell.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_moles        Array of moles of gas components. Dimensions 
    !> of the array are (nxyz, ngas_comps), where ngas_comps is the result of 
    !> @ref getgascomponentscount, and @a nxyz is the number of user grid cells 
    !> (@ref getgridcellcount). If the number of moles is set to a negative number, 
    !> the gas component will not be defined for the GAS_PHASE of the reaction cell.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascompmoles,
    !> @ref getgascomppressures,
    !> @ref getgascompphi,
    !> @ref getgasphasevolume,
    !> @ref setgasphasevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> ngas_comps = brm%SetGasComponentsCount()
    !> allocate(gas_moles(nxyz, ngas_comps))
    !> ...
    !> status = brm%setgascompmoles(gas_moles)
    !> status = brm%runcells()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setgascompmoles(self, gas_moles)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:,:), INTENT(in) :: gas_moles
    setgascompmoles = RM_setgascompmoles(self%bmiphreeqcrm_id, gas_moles)
    END FUNCTION setgascompmoles

    !> Transfer volumes of gas phases from the array given in the argument list 
    !> (@a gas_volume) to each reaction cell.
    !> The gas-phase volume affects the pressures calculated for fixed-volume
    !> gas phases. If a gas-phase volume is defined with this method for a 
    !> GAS_PHASE in a cell, the gas phase is forced to be a fixed-volume gas phase.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param gas_volume       Array of gas-phase volumes. Dimension of the array 
    !> is (nxyz), where @a nxyz is the number of user grid cells (@ref getgridcellcount).
    !> If an element of the array is set to a negative number, the gas component will
    !> not be defined for the GAS_PHASE of the reaction cell.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getgascomponentscount,
    !> @ref getgascompmoles,
    !> @ref getgascomppressures,
    !> @ref getgascompphi,
    !> @ref getgasphasevolume,
    !> @ref setgascompmoles.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(gas_volume(nxyz))
    !> ...
    !> status = brm%setgasphasevolume(gas_volume)
    !> status = brm%runcells()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setgasphasevolume(self, gas_volume)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: gas_volume
    setgasphasevolume = RM_setgasphasevolume(self%bmiphreeqcrm_id, gas_volume)
    END FUNCTION setgasphasevolume

    !> Transfer the concentrations for one component given by the vector @a c 
    !> to each reaction cell. 
    !> Units of concentration for @a c are defined by @ref setunitssolution. 
    !> It is required that  @a setithconcentration be called for each component 
    !> in the system before @ref runcells is called.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param i                One-based index for the component to transfer. 
    !> Indices refer to the order produced by @ref getcomponents. The total number 
    !> of components is given by @ref getcomponentcount.
    !> @param c                Array of concentrations to transfer to the reaction cells.
    !> Dimension of the vector is @a nxyz, where @a nxyz is the number of
    !> user grid cells (@ref getgridcellcount). Values for inactive cells are ignored.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see                    @ref findcomponents, @ref getcomponentcount, 
    !> @ref getcomponents, @ref setconcentrations.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%phreeqc_rm.setithconcentration(i, c) ! repeat for all components
    !> ...
    !> status = brm%phreeqc_rm.runcells()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setithconcentration(self, i, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: i
    real(kind=8), INTENT(in), DIMENSION(:) :: c
    setithconcentration = RM_setithconcentration(self%bmiphreeqcrm_id, i, c)
    return
    END FUNCTION setithconcentration

    !> Transfer the concentrations for one aqueous species given by the vector
    !> @a c to each reaction cell.
    !> Units of concentration for @a c are mol/L. To set species concentrations, 
    !> @ref setspeciessaveon must be set to @a true. It is required that
    !> @a setithspeciesconcentration be called for each aqueous species in the 
    !> system before @ref runcells is called. This method is for use with 
    !> multicomponent diffusion calculations. 
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param i                One-based index for the species to transfer. Indices 
    !> refer to the order produced by @ref getspeciesnames. The total number of 
    !> species is given by @ref getspeciescount.
    !> @param c                Array of concentrations to transfer to the reaction cells.
    !> Dimension of the array is @a nxyz, where @a nxyz is the number of user grid 
    !> cells (@ref getgridcellcount). Values for inactive cells are ignored.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see                    @ref findcomponents, @ref getspeciescount, @ref getspeciesnames,
    !> @ref speciesconcentrations2module, @ref setspeciessaveon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setithspeciesconcentration(i, c) ! repeat for all species
    !> ...
    !> status = brm%runcells()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setithspeciesconcentration(self, i, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: i
    real(kind=8), INTENT(in), DIMENSION(:) :: c
    setithspeciesconcentration = RM_setithspeciesconcentration(self%bmiphreeqcrm_id, i, c)
    return
    END FUNCTION setithspeciesconcentration

    !> MPI only. Defines a callback function that allows additional tasks to 
    !> be done by the workers. The method @ref mpiworker contains a loop,
    !> where the workers receive a message (an integer), run a function 
    !> corresponding to that integer, and then wait for another message.
    !> SetmpiworkerCallback allows the developer to add another function
    !> that responds to additional integer messages by calling developer-defined 
    !> functions corresponding to those integers. @ref mpiworker calls the 
    !> callback function when the message number is not one of the PhreeqcRM 
    !> message numbers. Messages are unique integer numbers. PhreeqcRM uses integers 
    !> in a range beginning at 0. It is suggested that developers use message numbers 
    !> starting at 1000 or higher for their tasks. The callback function calls a 
    !> developer-defined function specified by the message number and then returns 
    !> to @ref mpiworker to wait for another message.
    !> @n@n
    !> For Fortran, the functions that are called from the callback function
    !> can use USE statements to find the data necessary to perform the tasks, and
    !> the only argument to the callback function is an integer message argument.
    !> @a SetmpiworkerCallback must be called by each worker before @ref mpiworker 
    !> is called.
    !> @n@n
    !> The motivation for this method is to allow the workers to perform other
    !> tasks, for instance, parallel transport calculations, within the structure
    !> of @ref mpiworker. The callback function can be used to allow the workers 
    !> to receive data, perform transport calculations, and (or) send results, without 
    !> leaving the loop of @ref mpiworker. Alternatively, it is possible for the 
    !> workers to return from @ref mpiworker by a call to @ref mpiworkerbreak 
    !> by root. The workers could then call subroutines to receive data, calculate 
    !> transport, and send data, and then resume processing PhreeqcRM messages from 
    !> root with another call to @ref mpiworker.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param fcn              A function that returns an integer and has an integer argument.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref mpiworker,
    !> @ref mpiworkerbreak.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> Code executed by root:
    !> status = do_something()
    !>
    !> Code executed by workers:
    !> status = brm%SetmpiworkerCallback(worker_tasks_f)
    !> status = brm%mpiworker()
    !>
    !> Code executed by root and workers:
    !> integer function do_something
    !>   implicit none
    !>   INCLUDE 'mpif.h'
    !>   integer status
    !>   integer i, method_number, mpi_myself, mpi_task, mpi_tasks, worker_number;
    !>   method_number = 1000
    !>   call MPI_Comm_size(MPI_COMM_WORLD, mpi_tasks, status)
    !>   call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
    !>   if (mpi_myself .eq. 0) then
    !>     CALL MPI_Bcast(method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status)
    !>     write(*,*) "I am root."
    !>     do i = 1, mpi_tasks-1
    !>       CALL MPI_Recv(worker_number, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, status)
    !>       write(*,*) "Recieved data from worker number ", worker_number, "."
    !>     enddo
    !>   else
    !>     CALL MPI_Send(mpi_myself, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status)
    !>   endif
    !>   do_something = 0
    !> end function do_something
    !>
    !> Code called by workers from method mpiworker:
    !> integer(kind=C_INT) function worker_tasks_f(method_number) BIND(C, NAME='worker_tasks_f')
    !>   USE ISO_C_BINDING
    !>   implicit none
    !>   interface
    !>     integer function do_something
    !>     end function do_something
    !>   end interface
    !>   integer(kind=c_int), INTENT(inout) :: method_number
    !>   integer :: status
    !>   if (method_number .eq. 1000) then
    !>     status = do_something()
    !>   endif
    !>   worker_tasks_f = 0
    !> end function worker_tasks_f
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by workers, before call to @ref mpiworker.
    INTEGER FUNCTION SetmpiworkerCallback(self, fcn)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    
    INTERFACE
    INTEGER(kind=c_int) FUNCTION fcn(method_number) BIND(C)
    USE ISO_C_BINDING
    !INTEGER, INTENT(in) :: method_number
    INTEGER(kind=c_int), INTENT(in) :: method_number
    END FUNCTION fcn
    END INTERFACE
    SetmpiworkerCallback = RM_SetmpiworkerCallback(self%bmiphreeqcrm_id, fcn)
    END FUNCTION SetmpiworkerCallback

    !> Specify the current selected output by sequence number. The user may define 
    !> multiple SELECTED_OUTPUT data blocks for the workers. A user number is specified 
    !> for each data block, and the blocks are stored in user-number order. The value of
    !> the argument @a n selects the sequence number of the SELECTED_OUTPUT definition 
    !> that will be used for selected-output operations.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param n             Sequence number of the SELECTED_OUTPUT data block that is to be used.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setselectedoutputon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do isel = 1, brm%getselectedoutputcount()
    !>   status = brm%setnthselectedoutput(isel)
    !>   n_user = brm%getcurrentselectedoutputusernumber()
    !>   col = brm%getselectedoutputcolumncount()
    !>   allocate(selected_out(nxyz,col))
    !>   status = brm%getselectedoutput(selected_out)
    !>   ! Process results here
    !>   deallocate(selected_out)
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION setnthselectedoutput(self, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: n
    setnthselectedoutput = RM_setnthselectedoutput(self%bmiphreeqcrm_id, n)
    END FUNCTION setnthselectedoutput

    !> Sets the property for partitioning solids between the saturated and 
    !> unsaturated parts of a partially saturated cell.
    !> The option is intended to be used by saturated-only flow codes 
    !> that allow a variable water table. The value has meaning only when 
    !> saturations less than 1.0 are encountered. The partially saturated cells
    !> may have a small water-to-rock ratio that causes reactions to proceed 
    !> differently relative to fully saturated cells. By setting  
    !> @a setpartitionuzsolids to true, the amounts of solids and gases are 
    !> partioned according to the saturation. If a cell has a saturation of 0.5, then
    !> the water interacts with only half of the solids and gases; the other half 
    !> is unreactive until the water table rises. As the saturation in a cell varies,
    !> solids and gases are transferred between the saturated and unsaturated 
    !> (unreactive) reservoirs of the cell. Unsaturated-zone flow and transport codes 
    !> will probably use the default (false), which assumes all gases and solids are 
    !> reactive regardless of saturation.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf             @a True, the fraction of solids and gases available for
    !> reaction is equal to the saturation;
    !> @a False (default), all solids and gases are reactive regardless of saturation.
    !> @retval IRESULT    0 is success, negative is failure (See @ref decodeerror).
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setpartitionuzsolids(0)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setpartitionuzsolids(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in)  :: tf
    setpartitionuzsolids = RM_setpartitionuzsolids(self%bmiphreeqcrm_id, tf)
    END FUNCTION setpartitionuzsolids

    !> Set the porosity for each reaction cell.
    !> The volume of water in a reaction cell is the product of the porosity, 
    !> the saturation (@ref setsaturationuser), and the representative volume 
    !> (@ref setrepresentativevolume).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param por              Array of porosities, unitless. Default is 0.1. Size 
    !> of array is @a nxyz, where @a nxyz is the number of grid cells in the user's 
    !> model (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getsaturationcalculated,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(por(nxyz))
    !> por = 0.2
    !> status = brm%setporosity(por)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setporosity(self, por)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: por
    setporosity = RM_setporosity(self%bmiphreeqcrm_id, por)
    END FUNCTION setporosity

    !> Set the pressure for each reaction cell. 
    !> Pressure effects are considered explicitly only in three of the databases 
    !> distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param p                Array of pressures, in atm. Size of array is @a nxyz, 
    !> where @a nxyz is the number of grid cells in the user's model (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref settemperature.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(pressure(nxyz))
    !> pressure = 2.0
    !> status = brm%setpressure(pressure)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setpressure(self, p)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: p
    setpressure = RM_setpressure(self%bmiphreeqcrm_id, p)
    END FUNCTION setpressure

    !> Enable or disable detailed output for each reaction cell.
    !> Printing for a cell will occur only when the printing is enabled with 
    !> @ref setprintchemistryon and the @a cell_mask value is 1.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param cell_mask        Array of integers. Size of array is @a nxyz, 
    !> where @a nxyz is the number of grid cells in the user's model 
    !> (@ref getgridcellcount). A value of 0 will disable printing detailed 
    !> output for the cell; a value of 1 will enable printing detailed output for a cell.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setprintchemistryon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(print_chemistry_mask(nxyz))
    !>   do i = 1, nxyz/2
    !>   print_chemistry_mask(i) = 1
    !>   print_chemistry_mask(i+nxyz/2) = 0
    !> enddo
    !> status = brm%setprintchemistrymask(print_chemistry_mask)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setprintchemistrymask(self, cell_mask)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, DIMENSION(:), INTENT(in) :: cell_mask
    setprintchemistrymask = RM_setprintchemistrymask(self%bmiphreeqcrm_id, cell_mask)
    END FUNCTION setprintchemistrymask

    !> Setting to enable or disable printing detailed output from reaction calculations 
    !> to the output file for a set of cells defined by @ref setprintchemistrymask. 
    !> The detailed output prints all of the output typical of a PHREEQC reaction calculation, 
    !> which includes solution descriptions and the compositions of all other reactants. 
    !> The output can be several hundred lines per cell, which can lead to a very large 
    !> output file (prefix.chem.txt, @ref openfiles). For the worker instances, the 
    !> output can be limited to a set of cells (@ref setprintchemistrymask) and, in 
    !> general, the amount of information printed can be limited by use of options in 
    !> the PRINT data block of PHREEQC (applied by using @ref runfile or @ref runstring).
    !> Printing the detailed output for the workers is generally used only for debugging, 
    !> and PhreeqcRM will run significantly faster when printing detailed output for the 
    !> workers is disabled.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param workers          0, disable detailed printing in the worker instances; 
    !> 1, enable detailed printing in the worker instances.
    !> @param initial_phreeqc  0, disable detailed printing in the InitialPhreeqc instance;
    !> 1, enable detailed printing in the InitialPhreeqc instances.
    !> @param utility          0, disable detailed printing in the Utility instance; 
    !> 1, enable detailed printing in the Utility instance.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setprintchemistrymask.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setprintchemistryon(0, 1, 0)  ! workers, initial_phreeqc, utility
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setprintchemistryon(self, workers, initial_phreeqc, utility)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: workers, initial_phreeqc, utility
    setprintchemistryon = RM_setprintchemistryon(self%bmiphreeqcrm_id, workers, initial_phreeqc, utility)
    END FUNCTION setprintchemistryon

    !> Set the load-balancing algorithm.
    !> PhreeqcRM attempts to rebalance the load of each thread or process such that each
    !> thread or process takes the same amount of time to run its part of a @ref runcells
    !> calculation. Two algorithms are available; one uses individual times for each cell and
    !> accounts for cells that were not run because
    !> saturation was zero (default), and
    !> the other assigns an average time to all cells.
    !> The methods are similar, but limited testing indicates the default method performs better.
    !>
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param method           0, indicates average times are used in rebalancing; 
    !> 1 indicates individual cell times are used in rebalancing (default).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setrebalancefraction.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setrebalancebycell(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setrebalancebycell(self, method)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in)  :: method
    setrebalancebycell = RM_setrebalancebycell(self%bmiphreeqcrm_id, method)
    END FUNCTION setrebalancebycell

    !> Sets the fraction of cells that are transferred among threads or processes 
    !> when rebalancing.
    !> PhreeqcRM attempts to rebalance the load of each thread or process such that 
    !> each thread or process takes the same amount of time to run its part of a 
    !> @ref runcells calculation. The rebalancing transfers cell calculations among 
    !> threads or processes to try to achieve an optimum balance. @a setrebalancefraction
    !> adjusts the calculated optimum number of cell transfers by a fraction from 
    !> 0 to 1.0 to determine the actual number of cell transfers. A value of zero 
    !> eliminates load rebalancing. A value less than 1.0 is suggested to slow the approach 
    !> to the optimum cell distribution and avoid possible oscillations when too many cells 
    !> are transferred at one iteration, requiring reverse transfers at the next iteration.
    !> Default is 0.5.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param f                Fraction from 0.0 to 1.0.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setrebalancebycell.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setrebalancefraction(0.5d0)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setrebalancefraction(self, f)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(in)  :: f
    setrebalancefraction = RM_setrebalancefraction(self%bmiphreeqcrm_id, f)
    END FUNCTION setrebalancefraction

    !> Set the representative volume of each reaction cell.
    !> By default the representative volume of each reaction cell is 1 liter. The volume 
    !> of water in a reaction cell is determined by the procuct of the representative 
    !> volume, the porosity (@ref setporosity), and the saturation 
    !> (@ref setsaturationuser). The numerical method of PHREEQC is more robust if 
    !> the water volume for a reaction cell is within a couple orders of magnitude of 1.0.
    !> Small water volumes caused by small porosities and (or) small saturations (and (or) 
    !> small representative volumes) may cause non-convergence of the numerical method.
    !> In these cases, a larger representative volume may help. Note that increasing the 
    !> representative volume also increases the number of moles of the reactants in the 
    !> reaction cell (minerals, surfaces, exchangers, and others).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param rv               Vector of representative volumes, in liters. Default 
    !> is 1.0 liter. Size of array is @a nxyz, where @a nxyz is the number of grid cells 
    !> in the user's model (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setporosity,
    !> @ref setsaturationuser.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> real(kind=8), dimension(:), allocatable   :: rv
    !> allocate(rv(nxyz))
    !> rv = 1.0
    !> status = brm%setrepresentativevolume(rv)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setrepresentativevolume(self, rv)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: rv
    setrepresentativevolume = RM_setrepresentativevolume(self%bmiphreeqcrm_id, rv)
    END FUNCTION setrepresentativevolume

    !> Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
    !> The volume of water in a cell is the product of porosity (@ref setporosity), 
    !> saturation (@a setsaturationuser), and representative volume 
    !> (@ref setrepresentativevolume). As a result of a reaction calculation,
    !> solution properties (density and volume) will change; the databases phreeqc.dat, 
    !> Amm.dat, and pitzer.dat have the molar volume data to calculate these changes.
    !> The methods @ref getdensitycalculated, @ref getsolutionvolume, and 
    !> @ref getsaturationcalculated can be used to account for these changes in the 
    !> succeeding transport calculation. @a setrepresentativevolume should be called 
    !> before initial conditions are defined for the reaction cells.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param sat              Array of saturations, unitless. Size of array is @a nxyz, 
    !> where @a nxyz is the number of grid cells in the user's model 
    !> (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getdensitycalculated,
    !> @ref getsaturationcalculated,
    !> @ref getsolutionvolume,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(sat(nxyz))
    !> sat = 1.0
    !> status = brm%setsaturationuser(sat)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setsaturationuser(self, sat)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: sat
    setsaturationuser = RM_setsaturationuser(self%bmiphreeqcrm_id, sat)
    END FUNCTION setsaturationuser
    
    INTEGER FUNCTION SetSaturation(self, sat)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: sat
    SetSaturation = RM_setsaturationuser(self%bmiphreeqcrm_id, sat)
    END FUNCTION SetSaturation

    !> Set the property that controls whether messages are written to the screen.
    !> Messages include information about rebalancing during @ref runcells, and
    !> any messages written with @ref screenmessage.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf            @a 1, enable screen messages; @a 0, disable screen messages. 
    !> Default is 1.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref runcells,
    !> @ref screenmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setscreenon(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION setscreenon(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: tf
    setscreenon = RM_setscreenon(self%bmiphreeqcrm_id, tf)
    END FUNCTION setscreenon

    !> Setting determines whether selected-output results are available to be retrieved
    !> with @ref getselectedoutput. @a 1 indicates that selected-output results
    !> will be accumulated during @ref runcells and can be retrieved with
    !> @ref getselectedoutput; @a 0 indicates that selected-output results will not
    !> be accumulated during @ref runcells.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf               0, disable selected output; 1, enable selected output.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref getcurrentselectedoutputusernumber,
    !> @ref getnthselectedoutputusernumber,
    !> @ref getselectedoutput,
    !> @ref getselectedoutputcolumncount,
    !> @ref getselectedoutputcount,
    !> @ref getselectedoutputheadings,
    !> @ref getselectedoutputrowcount,
    !> @ref setcurrentselectedoutputusernumber,
    !> @ref setnthselectedoutput.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setselectedoutputon(1)        ! enable selected output
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setselectedoutputon(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: tf
    setselectedoutputon = RM_setselectedoutputon(self%bmiphreeqcrm_id, tf)
    END FUNCTION setselectedoutputon

    !> Sets the value of the species-save property.
    !> This method enables use of PhreeqcRM with multicomponent-diffusion 
    !> transport calculations. By default, concentrations of aqueous species are 
    !> not saved. Setting the species-save property to 1 allows aqueous species 
    !> concentrations to be retrieved with @ref getspeciesconcentrations, and 
    !> solution compositions to be set with @ref speciesconcentrations2module.
    !> setspeciessaveon must be called before calls to @ref findcomponents.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param save_on          0, indicates species concentrations are not saved; 
    !> 1, indicates species concentrations are saved.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref speciesconcentrations2module.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> save_on = brm%setspeciessaveon(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers.
    INTEGER FUNCTION setspeciessaveon(self, save_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: save_on
    setspeciessaveon = RM_setspeciessaveon(self%bmiphreeqcrm_id, save_on)
    END FUNCTION setspeciessaveon

    !> Set the temperature for each reaction cell. If @a settemperature is not called,
    !> worker solutions will have temperatures as defined by initial conditions
    !> (@ref initialphreeqc2module and @ref initialphreeqccell2module).
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param t                Array of temperatures, in degrees C. Size of array is 
    !> @a nxyz, where @a nxyz is the number of grid cells in the user's model 
    !> (@ref getgridcellcount).
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setpressure.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> allocate(temperature(nxyz))
    !> temperature = 20.0
    !> status = brm%settemperature(temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION settemperature(self, t)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:), INTENT(in) :: t
    settemperature = RM_settemperature(self%bmiphreeqcrm_id, t)
    END FUNCTION settemperature

    !> Set current simulation time for the reaction module.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param time             Current simulation time, in seconds.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref settimeconversion,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%settime(time)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION settime(self, time)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(in) :: time
    settime = RM_settime(self%bmiphreeqcrm_id, time)
    END FUNCTION settime

    !> Set a factor to convert to user time units. Factor times seconds produces user time units.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param conv_factor      Factor to convert seconds to user time units.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref settime,
    !> @ref settimestep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%settimeconversion(dble(1.0 / 86400.0)) ! days
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION settimeconversion(self, conv_factor)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(in) :: conv_factor
    settimeconversion = RM_settimeconversion(self%bmiphreeqcrm_id, conv_factor)
    END FUNCTION settimeconversion

    !> Set current time step for the reaction module. This is the length
    !> of time over which kinetic reactions are integrated.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param time_step        Current time step, in seconds.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref settime,
    !> @ref settimeconversion.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%settimestep(time_step)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION settimestep(self, time_step)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), INTENT(in) :: time_step
    settimestep = RM_settimestep(self%bmiphreeqcrm_id, time_step)
    END FUNCTION settimestep

    !> Sets input units for exchangers.
    !> In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
    !> @a setunitsexchange specifies how the number of moles of exchange sites in 
    !> a reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is 
    !> porosity (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
    !>
    !> If a single EXCHANGE definition is used for cells with different initial porosity,
    !> the three options scale quite differently.
    !> For option 0, the number of moles of exchangers will be the same regardless of porosity.
    !> For option 1, the number of moles of exchangers will be vary directly with porosity 
    !> and inversely with rock volume.
    !> For option 2, the number of moles of exchangers will vary directly with rock volume 
    !>and inversely with porosity.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for exchangers: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitsexchange(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitsexchange(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitsexchange = RM_setunitsexchange(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitsexchange

    !> Set input units for gas phases.
    !> In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
    !> @a setunitsgasphase specifies how the number of moles of component gases in a 
    !> reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is 
    !> porosity (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
    !>
    !> If a single GAS_PHASE definition is used for cells with different initial porosity,
    !> the three options scale quite differently.
    !> For option 0, the number of moles of a gas component will be the same regardless of porosity.
    !> For option 1, the number of moles of a gas component will be vary directly with porosity 
    !> and inversely with rock volume.
    !> For option 2, the number of moles of a gas component will vary directly with rock 
    !> volume and inversely with porosity.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for gas phases: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitsgasphase(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitsgasphase(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitsgasphase = RM_setunitsgasphase(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitsgasphase

    !> Set input units for kinetic reactants.
    !> In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
    !> @a setunitskinetics specifies how the number of moles of kinetic reactants in a 
    !> reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity 
    !> (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
    !>
    !> If a single KINETICS definition is used for cells with different initial porosity,
    !> the three options scale quite differently.
    !> For option 0, the number of moles of kinetic reactants will be the same regardless of porosity.
    !> For option 1, the number of moles of kinetic reactants will be vary directly with 
    !> porosity and inversely with rock volume.
    !> For option 2, the number of moles of kinetic reactants will vary directly with 
    !> rock volume and inversely with porosity.
    !>
    !> Note that the volume of water in a cell in the reaction module is equal to the product of
    !> porosity (@ref setporosity), the saturation (@ref setsaturationuser), and 
    !> representative volume (@ref setrepresentativevolume), which is usually less than 
    !> 1 liter. It is important to write the RATES definitions for homogeneous (aqueous) 
    !> kinetic reactions to account for the current volume of water, often by calculating 
    !> the rate of reaction per liter of water and multiplying by the volume
    !> of water (Basic function SOLN_VOL).
    !>
    !> Rates that depend on surface area of solids, are not dependent
    !> on the volume of water. However, it is important to get the correct surface area 
    !> for the kinetic reaction. To scale the surface area with the number of moles, the 
    !> specific area (m^2 per mole of reactant) can be defined as a parameter 
    !> (KINETICS; -parm), which is multiplied by the number of moles of
    !> reactant (Basic function M) in RATES to obtain the surface area.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for kinetic reactants: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitskinetics(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitskinetics(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitskinetics = RM_setunitskinetics(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitskinetics

    !> Set input units for pure phase assemblages (equilibrium phases).
    !> In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
    !> @a setunitsppassemblage specifies how the number of moles of phases in a 
    !> reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is 
    !> porosity (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
    !>
    !> If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity,
    !>    the three options scale quite differently.
    !> For option 0, the number of moles of a mineral will be the same regardless of porosity.
    !> For option 1, the number of moles of a mineral will be vary directly with porosity 
    !> and inversely with rock volume.
    !> For option 2, the number of moles of a mineral will vary directly with rock volume and 
    !> inversely with porosity.
    !>
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for equilibrium phases: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitsppassemblage(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitsppassemblage(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitsppassemblage = RM_setunitsppassemblage(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitsppassemblage

    !> Solution concentration units used by the transport model.
    !> Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
    !> PHREEQC defines solutions by the number of moles of each
    !> element in the solution.
    !> @n@n
    !> To convert from mg/L to moles
    !> of element in the representative volume of a reaction cell, mg/L is converted 
    !> to mol/L and multiplied by the solution volume, which is the product of porosity 
    !> (@ref setporosity), saturation (@ref setsaturationuser), and representative 
    !> volume (@ref setrepresentativevolume). To convert from mol/L to moles
    !> of element in the representative volume of a reaction cell, mol/L is multiplied 
    !> by the solution volume. To convert from mass fraction to moles of element in the 
    !> representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied 
    !> by density (@ref setdensityuser) and multiplied by the solution volume.
    !> @n@n
    !> To convert from moles  of element in the representative volume of a reaction cell 
    !> to mg/L, the number of moles of an element is divided by the solution volume 
    !> resulting in mol/L, and then converted to mg/L. To convert from moles of element 
    !> in a cell to mol/L,  the number of moles of an element is divided by the
    !> solution volume resulting in mol/L.
    !> @n@n
    !> To convert from moles of element in a cell to mass fraction, the number of moles 
    !> of an element is converted to kg and divided by the total mass of the solution.
    !> Two options are available for the volume and mass of solution that are used in 
    !> converting to transport concentrations: (1) the volume and mass of solution are
    !> calculated by PHREEQC, or (2) the volume of solution is the product of porosity 
    !> (@ref setporosity), saturation (@ref setsaturationuser), and representative 
    !> volume (@ref setrepresentativevolume), and the mass of solution is volume times 
    !> density as defined by @ref setdensityuser. Which option is used is determined 
    !> by @ref usesolutiondensityvolume.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref setdensityuser,
    !> @ref setporosity,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser,
    !> @ref usesolutiondensityvolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitssolution(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitssolution(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitssolution = RM_setunitssolution(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitssolution

    !> Set input units for solid-solution assemblages.
    !> In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
    !> @a setunitsssassemblage specifies how the number of moles of solid-solution 
    !> components in a reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is 
    !> porosity (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for solid solutions: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> If a single SOLID_SOLUTION definition is used for cells with different initial porosity,
    !> the three options scale quite differently.
    !> For option 0, the number of moles of a solid-solution component will be the same 
    !> regardless of porosity.
    !> For option 1, the number of moles of a solid-solution component will be vary directly 
    !> with porosity and inversely with rock volume.
    !> For option 2, the number of moles of a solid-solution component will vary directly 
    !> with rock volume and inversely with porosity.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitsssassemblage(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitsssassemblage(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: option
    setunitsssassemblage = RM_setunitsssassemblage(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitsssassemblage

    !> Set input units for surfaces.
    !> In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
    !> @a setunitssurface specifies how the number of moles of surface sites in a 
    !> reaction cell (@a Mc) is calculated from the input value (@a Mp).
    !>
    !> Options are
    !> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the 
    !> representative volume (@ref setrepresentativevolume);
    !> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is 
    !> porosity (@ref setporosity); or
    !> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
    !>
    !> If a single SURFACE definition is used for cells with different initial porosity,
    !> the three options scale quite differently.
    !> For option 0, the number of moles of surface sites will be the same regardless of porosity.
    !> For option 1, the number of moles of surface sites will be vary directly with 
    !> porosity and inversely with rock volume.
    !> For option 2, the number of moles of surface sites will vary directly with 
    !> rock volume and inversely with porosity.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param option           Units option for surfaces: 0, 1, or 2.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref initialphreeqc2module,
    !> @ref initialphreeqccell2module,
    !> @ref setporosity,
    !> @ref setrepresentativevolume.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setunitssurface(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION setunitssurface(self, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    integer, intent(in) :: option
    setunitssurface = RM_setunitssurface(self%bmiphreeqcrm_id, option)
    END FUNCTION setunitssurface

    !> Set solution concentrations in the reaction cells based on the array 
    !> of aqueous species concentrations (@a species_conc).
    !> This method is intended for use with multicomponent-diffusion transport 
    !> calculations, and @ref setspeciessaveon must be set to @a true.
    !> The list of aqueous species is determined by @ref findcomponents and 
    !> includes all aqueous species that can be made from the set of components.
    !> The method determines the total concentration of a component by summing 
    !> the molarities of the individual species times the stoichiometric
    !> coefficient of the element in each species. Solution compositions in the 
    !> reaction cells are updated with these component concentrations.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param species_conc     Array of aqueous species concentrations. Dimension of 
    !> the array is (@a nxyz, @a nspecies), where @a nxyz is the number of user grid 
    !> cells (@ref getgridcellcount), and @a nspecies is the number of aqueous 
    !> species (@ref getspeciescount). Concentrations are moles per liter.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref findcomponents,
    !> @ref getspeciesconcentrations,
    !> @ref getspeciescount,
    !> @ref getspeciesd25,
    !> @ref getspecieslog10gammas,
    !> @ref getspecieslog10molalities,
    !> @ref getspeciesnames,
    !> @ref getspeciessaveon,
    !> @ref getspeciesz,
    !> @ref setspeciessaveon.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%setspeciessaveon(1)
    !> ncomps = brm%findcomponents()
    !> nspecies = brm%getspeciescount()
    !> nxyz = brm%getgridcellcount()
    !> allocate(species_c(nxyz, nspecies))
    !> ...
    !> status = brm%speciesconcentrations2module(species_c(1,1))
    !> status = brm%runcells()
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION speciesconcentrations2module(self, species_conc)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    real(kind=8), DIMENSION(:,:), INTENT(in) :: species_conc
    speciesconcentrations2module = RM_speciesconcentrations2module(self%bmiphreeqcrm_id, species_conc)
    END FUNCTION speciesconcentrations2module

    !> Save the state of the chemistry in all model cells, including SOLUTIONs,
    !> EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
    !> Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
    !> will be saved for each cell, if they have been defined in the worker IPhreeqc instances.
    !> The distribution of cells among the workers and the chemistry of fully or partially
    !> unsaturated cells are also saved. The state is saved in memory; use @ref dumpmodule 
    !> to save the state to file. PhreeqcRM can be reset to this state by using @ref stateapply.
    !> A state is identified by an integer, and multiple states can be saved.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param istate        Integer identifying the state that is saved.
    !> @retval IRESULT   0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref dumpmodule,
    !> @ref stateapply, and
    !> @ref statedelete.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%statesave(1)
    !> ...
    !> status = brm%stateapply(1)
    !> status = brm%statedelete(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION statesave(self, istate)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: istate
    statesave = RM_statesave(self%bmiphreeqcrm_id, istate)
    END FUNCTION statesave

    !> Reset the state of the module to a state previously saved with @ref statesave.
    !> The chemistry of all model cells are reset, including SOLUTIONs,
    !> EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
    !> MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
    !> will be reset for each cell, if they were defined in the worker IPhreeqc instances
    !> at the time the state was saved.
    !> The distribution of cells among the workers and the chemistry of fully or partially
    !> unsaturated cells are also reset to the saved state.
    !> The state to be applied is identified by an integer.
    !> @param self Fortran-supplied BMIPhreeqcRM instance..
    !> @param istate       Integer identifying the state that is to be applied.
    !> @retval IRESULT  0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref statesave and
    !> @ref statedelete.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%statesave(1)
    !> ...
    !> status = brm%stateapply(1)
    !> status = brm%statedelete(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION stateapply(self, istate)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: istate
    stateapply = RM_stateapply(self%bmiphreeqcrm_id, istate)
    END FUNCTION stateapply

    !> Delete a state previously saved with @ref statesave.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param istate      Integer identifying the state that is to be deleted.
    !> @retval IRESULT 0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref statesave and
    !> ref stateapply.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%statesave(1)
    !> ...
    !> status = brm%stateapply(1)
    !> status = brm%statedelete(1)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.
    INTEGER FUNCTION statedelete(self, istate)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: istate
    statedelete = RM_statedelete(self%bmiphreeqcrm_id, istate)
    END FUNCTION statedelete
    !> Determines the volume and density to use when converting from the reaction-module concentrations
    !> to transport concentrations (@ref getconcentrations).
    !> Two options are available to convert concentration units:
    !> (1) the density and solution volume calculated by PHREEQC are used, or
    !> (2) the specified density (@ref setdensityuser)
    !> and solution volume are defined by the product of
    !> saturation (@ref setsaturationuser), porosity (@ref setporosity),
    !> and representative volume (@ref setrepresentativevolume).
    !> Transport models that consider density-dependent flow will probably use the
    !> PHREEQC-calculated density and solution volume (default),
    !> whereas transport models that assume constant-density flow will probably use
    !> specified values of density and solution volume.
    !> Only the following databases distributed with PhreeqcRM have molar volume information
    !> needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
    !> Density is only used when converting to transport units of mass fraction.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param tf               @a True indicates that the solution density and volume as
    !> calculated by PHREEQC will be used to calculate concentrations.
    !> @a False indicates that the solution density set by @ref setdensityuser and 
    !> the volume determined by the product of  @ref setsaturationuser, @ref setporosity, 
    !> and @ref setrepresentativevolume, will be used to calculate concentrations retrieved 
    !> by @ref getconcentrations.
    !> @see
    !> @ref getconcentrations,
    !> @ref setdensityuser,
    !> @ref setporosity,
    !> @ref setrepresentativevolume,
    !> @ref setsaturationuser.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%usesolutiondensityvolume(0)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref mpiworker.

    INTEGER FUNCTION usesolutiondensityvolume(self, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    INTEGER, INTENT(in) :: tf
    usesolutiondensityvolume = RM_usesolutiondensityvolume(self%bmiphreeqcrm_id, tf)
    END FUNCTION usesolutiondensityvolume

    !> Print a warning message to the screen and the log file.
    !> @param self Fortran-supplied BMIPhreeqcRM instance.
    !> @param warn_str         String to be printed.
    !> @retval IRESULT      0 is success, negative is failure (See @ref decodeerror).
    !> @see
    !> @ref errormessage,
    !> @ref logmessage,
    !> @ref openfiles,
    !> @ref outputmessage,
    !> @ref screenmessage.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = brm%warningmessage("Parameter is out of range, using default")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and (or) workers; only root writes to the log file.
    INTEGER FUNCTION warningmessage(self, warn_str)
    USE ISO_C_BINDING
    IMPLICIT NONE
	class(bmi), intent(inout) :: self
    CHARACTER(len=*), INTENT(in) :: warn_str
    warningmessage = RM_warningmessage(self%bmiphreeqcrm_id, warn_str)
    END FUNCTION warningmessage

#endif 
!End EXTEND_BMIPHREEQCRM

    
    END MODULE BMIPhreeqcRM
