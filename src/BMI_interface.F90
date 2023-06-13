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
    !> You will need installed include directory (-I) added to the project) to 
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
#ifndef USE_MPI
    INTERFACE bmif_create
		module procedure bmif_create_default 
		module procedure bmif_create  
    END INTERFACE bmif_create
#endif

    !> @a bmif_get_value retrieves model variables. 
    !>
    !> Only variables in the list
    !> provided by @ref bmif_get_output_var_names can be retrieved. 
    INTERFACE bmif_get_value
        module procedure bmif_get_value_logical
        module procedure bmif_get_value_char
        module procedure bmif_get_value_char1
        module procedure bmif_get_value_double
        module procedure bmif_get_value_double1
        module procedure bmif_get_value_double2
        module procedure bmif_get_value_float ! not implemented
        module procedure bmif_get_value_int
        module procedure bmif_get_value_int1
        module procedure bmif_get_value_int2
    END INTERFACE bmif_get_value
    !> @a get_value_at_indices is not implemented. 
    !INTERFACE get_value_at_indices
    ! module procedure get_value_at_indices_double ! not implemented
    ! module procedure get_value_at_indices_float  ! not implemented
    ! module procedure get_value_at_indices_int    ! not implemented
    ! END INTERFACE get_value_at_indices
    !> @a bmif_set_value sets model variables. 
    !>
    !> Only variables in the list
    !> provided by @ref bmif_get_input_var_names can be set. 
   INTERFACE bmif_set_value
        module procedure bmif_set_value_b
        module procedure bmif_set_value_c
        module procedure bmif_set_value_double
        module procedure bmif_set_value_double1
        module procedure bmif_set_value_double2
        module procedure bmif_set_value_float ! not implemented
        module procedure bmif_set_value_int
        module procedure bmif_set_value_int1
        module procedure bmif_set_value_int2
    END INTERFACE bmif_set_value
    
    !> @a set_value_at_indices is not implemented. 
    INTERFACE set_value_at_indices
        module procedure set_value_at_indices_double ! not implemented
        module procedure set_value_at_indices_float  ! not implemented
        module procedure set_value_at_indices_int    ! not implemented
    END INTERFACE set_value_at_indices
    !> @a bmif_get_value_ptr retrieves pointers to model variables. 
    !>
    !> Only variables in the list
    !> provided by @ref bmif_get_pointable_var_names can be pointed to. 
    INTERFACE bmif_get_value_ptr
        module procedure bmif_get_value_ptr_logical
        module procedure bmif_get_value_ptr_double
        module procedure bmif_get_value_ptr_double1
        module procedure bmif_get_value_ptr_float ! not implemented
        module procedure bmif_get_value_ptr_integer
    END INTERFACE bmif_get_value_ptr
    
    CONTAINS
#ifndef USE_MPI
    ! ====================================================
    ! Initialize, run, finalize (IRF) 
    ! ====================================================
    !> @a bmif_create creates a BMIPhreeqcRM instance. 
    !> The default method, with no arguments, can be
    !> used if the instance is to be initialized with a YAML file.
    !> The YAML file must provide the number of cells in the user's 
    !> model through the YAMLSetGridCellCount method, and, optionally,
    !> the number of threads with YAMLThreadCount. 
    INTEGER FUNCTION bmif_create_default()
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RM_BMI_Create_default() &
        BIND(C, NAME='RM_BMI_Create_default')
    USE ISO_C_BINDING
    IMPLICIT NONE
    END FUNCTION RM_BMI_Create_default
    END INTERFACE
    bmif_create_default = RM_BMI_Create_default()
    return
    END FUNCTION bmif_create_default 
#endif    
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
    !> id = bmif_create(nxyz, nthreads)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.
    INTEGER FUNCTION bmif_create(nxyz, nthreads)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RM_BMI_Create(nxyz, nthreads) &
        BIND(C, NAME='RM_BMI_Create')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: nxyz
    INTEGER(KIND=C_INT), INTENT(in) :: nthreads
    END FUNCTION RM_BMI_Create
    END INTERFACE
    INTEGER, INTENT(in) :: nxyz
    INTEGER, INTENT(in) :: nthreads
    bmif_create = RM_BMI_Create(nxyz, nthreads)
    return
    END FUNCTION bmif_create    

    !> @a bmif_initialize uses a YAML file to initialize an instance of BMIPhreeqcRM. Same as
    !> @ref RM_InitializeYAML.
    !> @param id            The instance @a id returned from @ref bmif_create.
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
    !>  key: LoadDatabase
    !>   database: phreeqc.dat
    !> - key: RunFile
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
    !> @n CloseFiles();
    !> @n CreateMapping(std::vector< int >& grid2chem);
    !> @n DumpModule();
    !> @n FindComponents();
    !> @n InitialEquilibriumPhases2Module(std::vector< int > equilibrium_phases);
    !> @n InitialExchanges2Module(std::vector< int > exchanges);
    !> @n InitialGasPhases2Module(std::vector< int > gas_phases);
    !> @n InitialKineticss2Module(std::vector< int > kinetics);
    !> @n InitialSolidSolutions2Module(std::vector< int > solid_solutions);
    !> @n InitialSolutions2Module(std::vector< int > solutions);
    !> @n InitialSurfaces2Module(std::vector< int > surfaces);
    !> @n InitialPhreeqc2Module(std::vector< int > initial_conditions1);
    !> @n InitialPhreeqc2Module(std::vector< int > initial_conditions1, 
    !>    std::vector< int > initial_conditions2, std::vector< double > fraction1);
    !> @n InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
    !> @n LoadDatabase(std::string database);
    !> @n OpenFiles();
    !> @n OutputMessage(std::string str);
    !> @n RunCells();
    !> @n RunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
    !> @n RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    !> @n ScreenMessage(std::string str);
    !> @n SetComponentH2O(bool tf);
    !> @n SetConcentrations(std::vector< double > c);
    !> @n SetCurrentSelectedOutputUserNumber(int n_user);
    !> @n SetDensityUser(std::vector< double > density);
    !> @n SetDumpFileName(std::string dump_name);
    !> @n SetErrorHandlerMode(int mode);
    !> @n SetErrorOn(bool tf);
    !> @n SetFilePrefix(std::string prefix);
    !> @n SetGasCompMoles(std::vector< double > gas_moles);
    !> @n SetGasPhaseVolume(std::vector< double > gas_volume);
    !> @n SetPartitionUZSolids(bool tf);
    !> @n SetPorosity(std::vector< double > por);
    !> @n SetPressure(std::vector< double > p);
    !> @n SetPrintChemistryMask(std::vector< int > cell_mask);
    !> @n SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
    !> @n SetRebalanceByCell(bool tf);
    !> @n SetRebalanceFraction(double f);
    !> @n SetRepresentativeVolume(std::vector< double > rv);
    !> @n SetSaturationUser(std::vector< double > sat);
    !> @n SetScreenOn(bool tf);
    !> @n SetSelectedOutputOn(bool tf);
    !> @n SetSpeciesSaveOn(bool save_on);
    !> @n SetTemperature(std::vector< double > t);
    !> @n SetTime(double time);
    !> @n SetTimeConversion(double conv_factor);
    !> @n SetTimeStep(double time_step);
    !> @n SetUnitsExchange(int option);
    !> @n SetUnitsGasPhase(int option);
    !> @n SetUnitsKinetics(int option);
    !> @n SetUnitsPPassemblage(int option);
    !> @n SetUnitsSolution(int option);
    !> @n SetUnitsSSassemblage(int option);
    !> @n SetUnitsSurface(int option);
    !> @n SpeciesConcentrations2Module(std::vector< double > species_conc);
    !> @n StateSave(int istate);
    !> @n StateApply(int istate);
    !> @n StateDelete(int istate);
    !> @n UseSolutionDensityVolume(bool tf);
    !> @n WarningMessage(std::string warnstr);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !>
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> id = bmif_create(nxyz, nthreads)
    !> status = bmif_initializeYAML(id, "myfile.yaml")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
#ifdef USE_YAML
    INTEGER FUNCTION bmif_initialize(id, config_file)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: config_file
    bmif_initialize = success(RMF_BMI_Initialize(id, trim(config_file)//C_NULL_CHAR))
    return
    END FUNCTION bmif_initialize
#endif

    !> @a bmif_update runs a reaction step for all of the cells in the reaction module.
    !> Same as @ref RM_RunCells.
    !> @param id   The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_set_value(id, "Porosity", por)                ! If pore volume changes
    !> status = bmif_set_value(id, "SaturationUser", sat)              ! If saturation changes
    !> status = bmif_set_value(id, "Temperature", temperature)     ! If temperature changes
    !> status = bmif_set_value(id, "Pressure", pressure)           ! If pressure changes
    !> status = bmif_set_value(id, "Concentrations", c)            ! Transported concentrations
    !> status = bmif_set_value(id, "TimeStep", time_step)          ! Time step for kinetic reactions
    !> status = bmif_update(id)
    !> status = bmif_get_value(id, "Concentrations", c)            ! Concentrations after reaction
    !> status = bmif_get_value(id, "DensityCalculated", density)   ! Density after reaction
    !> status = bmif_get_value(id, "SolutionVolume", volume)       ! Solution volume after reaction
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
    INTEGER FUNCTION bmif_update(id)
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
    INTEGER, INTENT(in) :: id
    bmif_update = success(RMF_BMI_Update(id))
    return
    END FUNCTION bmif_update
     
    !> @a bmif_update_until is the same as @ref bmif_update, except the time step is calculated
    !> from the argument @a end_time. The time step is calculated to be @a end_time minus 
    !> the current time (@ref bmif_get_current_time).
    !> @param id       The instance @a id returned from @ref bmif_create.
    !> @param end_time Time at the end of the time step. 
    !> @see
    !> @ref bmif_initialize,
    !> @ref bmif_update.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_set_value(id, "Time", time)
    !> status = bmif_set_value(id, "Concentrations", c)
    !> status = bmif_update_until(time + 86400.0)
    !> status = bmif_get_value("Concentrations", c)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.    

    INTEGER FUNCTION bmif_update_until(id, end_time)
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
    INTEGER, INTENT(in) :: id
    real(kind=8), INTENT(in) :: end_time
    bmif_update_until = success(RMF_BMI_UpdateUntil(id, end_time))
    return
    END FUNCTION bmif_update_until

    !> @a bmif_finalize destroys a reaction module, same as @ref RM_Destroy.
    !> @param id  The instance @a id returned from @ref bmif_create.
    !> @retval    0 is success, 1 is failure.
    !> @see
    !> @ref bmif_create.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_finalize(id)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.
    INTEGER FUNCTION bmif_finalize(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_Destroy(id) &
        BIND(C, NAME='RMF_BMI_Destroy')
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER(KIND=C_INT), INTENT(in) :: id
    END FUNCTION RMF_BMI_Destroy
    END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER :: status
#ifdef USE_MPI
    !status = RM_MpiWorkerBreak(id)
#endif    
    bmif_finalize = success(RMF_BMI_Destroy(id))
    return
    END FUNCTION bmif_finalize
    
    ! ====================================================
    ! Exchange items
    ! ====================================================

    !> @a bmif_get_component_name returns the component name--BMIPhreeqcRM. 
    !>
    !> @param id             The instance @a id returned from @ref bmif_create.
    !> @param component_name Filled with "BMIPhreeqcRM", the name of the component.
    !> @retval               0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_get_component_name(id, component_name)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_component_name(id, component_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: component_name
    component_name = "BMI_PhreeqcRM"
    bmif_get_component_name = BMI_SUCCESS
    return
    END FUNCTION bmif_get_component_name

    !> @a bmif_get_input_item_count returns count of variables that
    !> can be set with @ref bmif_set_value.
    !> @param id     The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_input_item_count(id, count)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    !>
    INTEGER FUNCTION bmif_get_input_item_count(id, count)
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
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetInputItemCount(id)
    bmif_get_input_item_count = success(count)
    END FUNCTION bmif_get_input_item_count

    !> @a bmif_get_output_item_count returns count of output variables that can be
    !> retrieved with @ref bmif_get_value.
    !> @param id     The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_output_item_count(id, outputvarcount)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_output_item_count(id, count)
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
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetOutputItemCount(id)
    bmif_get_output_item_count = success(count)
    END FUNCTION bmif_get_output_item_count
    
    !> @a bmif_get_pointable_item_count returns count of variables for which
    !> pointers can be retrieved with @ref bmif_get_value_ptr.
    !> @param id     The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_pointable_item_count(id, pointablevarcount)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    !>
    INTEGER FUNCTION bmif_get_pointable_item_count(id, count)
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
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(inout) :: count
    INTEGER :: status
    count = RMF_BMI_GetPointableItemCount(id)
    bmif_get_pointable_item_count = success(count)
    END FUNCTION bmif_get_pointable_item_count

    !> Basic Model Interface method that returns a list of the variable names that can be set
    !> with @ref bmif_set_value.
    !> @param id          The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_input_var_names(id, inputvars)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_input_var_names(id, var_names)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "inputvarnames"
    status = RMF_BMI_GetNamesSize(id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_input_item_count(id, dim)
    nbytes = dim * itemsize
    dim = nbytes / itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_input_var_names = RMF_BMI_GetNames(id, trim(vartype)//C_NULL_CHAR, var_names(1))
    return
    END FUNCTION bmif_get_input_var_names

    !> @a bmif_get_output_var_names returns a list of the variable names that can be
    !> retrieved with @ref bmif_get_value.
    !> @param id          The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_output_var_names(id, var_names)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_output_var_names(id, var_names)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "outputvarnames"
    status = RMF_BMI_GetNamesSize(id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_output_item_count(id, dim)
    nbytes = dim * itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_output_var_names = RMF_BMI_GetNames(id, trim(vartype)//C_NULL_CHAR, var_names(1))
    return
    END FUNCTION bmif_get_output_var_names

    !> @a bmif_get_pointable_var_names returns a list of the variable names for which a pointer can be
    !> retrieved with @ref bmif_get_value_ptr.
    !> @param id          The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_pointable_var_names(id, var_names)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_pointable_var_names(id, var_names)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: var_names
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    vartype = "pointablevarnames"
    status = RMF_BMI_GetNamesSize(id, trim(vartype)//C_NULL_CHAR, itemsize)
    status = bmif_get_pointable_item_count(id, dim)
    nbytes = dim * itemsize
    dim = nbytes / itemsize
    if(allocated(var_names)) deallocate(var_names)
    allocate(character(len=itemsize) :: var_names(dim))
    bmif_get_pointable_var_names = RMF_BMI_GetNames(id, trim(vartype)//C_NULL_CHAR, var_names(1))
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
    !> @param id    The instance @a id returned from @ref bmif_create.
    !> @param var   Varaiable name. (Return value is the same regardless of @a name.)
    !> @param grid  1 (points). BMIPhreeqcRM cells derive meaning from the user's model. 
    !> @retval      0 is success, 1 is failure.
    integer FUNCTION bmif_get_var_grid(id, var, grid)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
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
    !> @param id    The instance @a id returned from @ref bmif_create.
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
    !>     status = bmif_get_var_units(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = bmif_get_var_type(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = bmif_get_var_itemsize(id, inputvars(i), itemsize)
    !>     write(*, "(5x, I15)") itemsize
    !>     status = bmif_get_var_nbytes(id, inputvars(i), nbytes)
    !>     write(*, "(5x, I15)") nbytes
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_type(id, var, vtype)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: vtype
    bmif_get_var_type = success(RMF_BMI_GetVarType(id, trim(var)//C_NULL_CHAR, vtype, len(vtype)))
    return
    END FUNCTION bmif_get_var_type
    
    !> @a bmif_get_var_units retrieves the units of a
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
!> @param id      The instance @a id returned from @ref bmif_create.
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
    !>     status = bmif_get_var_units(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = bmif_get_var_type(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = bmif_get_var_itemsize(id, inputvars(i), itemsize)
    !>     write(*, "(5x, I15)") itemsize
    !>     status = bmif_get_var_nbytes(id, inputvars(i), nbytes)
    !>     write(*, "(5x, I15)") nbytes
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_units(id, var, units)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: units
    bmif_get_var_units = success(RMF_BMI_GetVarUnits(id, trim(var)//C_NULL_CHAR, units, len(units)))
    return
    END FUNCTION bmif_get_var_units

    !> @a bmif_get_var_itemsize retrieves the size, in bytes, of a
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> Sizes may be the size of an integer, real(kind=8), 
    !> or a character length for string variables.
    !> @param id         The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_var_nbytes(id, "Temperature", nbytes)
    !> status = bmif_get_var_itemsize(id, "Temperature", item_size)
    !> dim    = nbytes/item_size
    !> allocate(bmi_temperature(dim))
    !> bmi_temperature = 25.0
    !> status = bmif_set_value("Temperature", bmi_temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_itemsize(id, var, itemsize)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    INTEGER, INTENT(out) :: itemsize
    itemsize = RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR)
    bmif_get_var_itemsize = success(bmif_get_var_itemsize)
    END FUNCTION bmif_get_var_itemsize

    !> @a bmif_get_var_nbytes retrieves the total number of bytes needed for a 
    !> variable that can be set with
    !> @ref bmif_set_value, retrieved with @ref bmif_get_value, or pointed to with
    !> @ref bmif_get_value_ptr.
    !> @param id      The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_var_nbytes(id, "Temperature", nbytes)
    !> status = bmif_get_var_itemsize(id, "Temperature", item_size)
    !> dim    = nbytes/item_size
    !> allocate(bmi_temperature(dim))
    !> bmi_temperature = 25.0
    !> status = bmif_set_value("Temperature", bmi_temperature)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_var_nbytes(id, var, nbytes)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    INTEGER, INTENT(out) :: nbytes
    nbytes  = RMF_BMI_GetVarNbytes(id, trim(var)//C_NULL_CHAR)
    bmif_get_var_nbytes = success(nbytes)
    END FUNCTION bmif_get_var_nbytes

    ! ====================================================
    ! Time information
    ! ====================================================

    !> @a bmif_get_current_time returns the current simulation time, in seconds. (Same as @ref RM_GetTime.)
    !> @param id      The instance @a id returned from @ref bmif_create.
    !> @param time    The current simulation time, in seconds.
    !> @retval        0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_end_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_get_current_time(id, time)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_current_time(id, time)
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
    INTEGER, INTENT(in) :: id
    real(kind=8), intent(inout) :: time
    time = RMF_BMI_GetCurrentTime(id)
    bmif_get_current_time = BMI_SUCCESS
    END FUNCTION bmif_get_current_time

    !> @a bmif_get_start_time returns the current simulation time, in seconds. 
    !> (Same as @ref bmif_get_current_time and @ref RM_GetTime.)
    !> @param id          The instance @a id returned from @ref bmif_create.
    !> @param start_time  The current simulation time, in seconds.
    !> @retval        0 is success, 1 is failure.
    INTEGER FUNCTION bmif_get_start_time(id, start_time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    real(kind=8), INTENT(inout) :: start_time
    bmif_get_start_time = bmif_get_current_time(id, start_time)
    END FUNCTION bmif_get_start_time

    !> @a bmif_get_end_time returns @ref bmif_get_current_time plus
    !> @ref bmif_get_time_step, in seconds.
    !> @param id        The instance @a id returned from @ref bmif_create.
    !> @param end_time  The end of the time step, in seconds.
    !> @retval          0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_get_end_time(id, end_time)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_end_time(id, end_time)
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
    INTEGER, INTENT(in) :: id
    real(kind=8), intent(inout) :: end_time
    end_time = RMF_BMI_GetEndTime(id)
    bmif_get_end_time = BMI_SUCCESS
    END FUNCTION bmif_get_end_time

    !> @a bmif_get_time_units returns the time units of PhreeqcRM.
    !> All time units are seconds for PhreeqcRM.
    !> @param id            The instance @a id returned from @ref bmif_create.
    !> @param time_units    Returns the string "seconds".
    !> @retval              0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_end_time,
    !> @ref bmif_get_time_step,
    !> @ref bmif_set_value,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep,
    !> @ref bmif_set_value.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(20) time_units
    !> status = bmif_get_time_units(id, time_units)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_time_units(id, time_units)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetTimeUnits(id, time_units, l) &
            BIND(C, NAME='RMF_BMI_GetTimeUnits')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id, l
        CHARACTER(KIND=C_CHAR), INTENT(inout) :: time_units(*)
        END FUNCTION RMF_BMI_GetTimeUnits
        END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: time_units
    bmif_get_time_units = success(RMF_BMI_GetTimeUnits(id, time_units, len(time_units)))
    return
    END FUNCTION bmif_get_time_units


    !> @a bmif_get_time_step returns the current simulation time step,
    !> in seconds. (Same as @ref RM_GetTimeStep.)
    !> @param id            The instance @a id returned from @ref bmif_create.
    !> @param time_step     The current simulation time step, in seconds.
    !> @retval              0 is success, 1 is failure.
    !> @see
    !> @ref bmif_get_current_time,
    !> @ref bmif_get_end_time,
    !> @ref bmif_set_value,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_get_time_step(id, time_step)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_get_time_step(id, time_step)
    USE ISO_C_BINDING
    IMPLICIT NONE
        INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetTimeStep(id) &
            BIND(C, NAME='RMF_BMI_GetTimeStep')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetTimeStep
        END INTERFACE
    INTEGER, INTENT(in) :: id
    real(kind=8), intent(inout) :: time_step
    time_step = RMF_BMI_GetTimeStep(id)
    bmif_get_time_step = BMI_SUCCESS
    END FUNCTION bmif_get_time_step
        
    ! ====================================================
    ! Getters, by type
    ! ====================================================

    !> @a bmif_get_value retrieves model variables. Only variables in the list
    !> provided by @ref bmif_get_output_var_names can be retrieved. 
    !> @param id     The instance @a id returned from @ref bmif_create.
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
    !> status = bmif_get_value(id, "DensityCalculated", bmi_density)
    !> status = bmif_get_value(id, "Components", bmi_comps)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
!Add NEW_VARIABLE to bmif_get_value Documentation
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_logical(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL, INTENT(inout) :: dest
    character(100) :: vartype
    integer :: status, dest_int
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "logical") then
        write(*,*) vartype, " logical"
        stop "Variable type error."
    endif
    bmif_get_value_logical = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest_int)
    dest = .true.
    if (dest_int .eq. 0) dest = .false.
    return
    END FUNCTION bmif_get_value_logical

    !> \overload
    INTEGER FUNCTION bmif_get_value_char(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: dest
    character(BMI_MAX_TYPE_NAME) :: vartype
    integer :: itemsize, status
    CHARACTER(len=:), allocatable :: temp
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "character") then
        write(*,*) vartype, " character"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(id, var, itemsize)
    allocate(character(len=itemsize) :: temp)
    status = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, temp)
    if (len(dest) .gt. 0) then
        dest = temp
    else
        status = RM_ErrorMessage(id, "Variable length is zero")   
        status = -1
    endif
    bmif_get_value_char = success(status)
    return
    END FUNCTION bmif_get_value_char
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_char_alloc(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=:), allocatable, INTENT(inout) :: dest
    character(BMI_MAX_TYPE_NAME) :: vartype
    integer :: itemsize, status
    CHARACTER(len=:), allocatable :: temp
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "character") then
        write(*,*) vartype, " character"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(id, var, itemsize)
    allocate(character(len=itemsize) :: temp)
    status = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, temp)
    dest = temp
    bmif_get_value_char_alloc = success(status)
    return
    END FUNCTION bmif_get_value_char_alloc

    !> \overload
    INTEGER FUNCTION bmif_get_value_char1(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "character(len=:),allocatable,dimension(:)") then
        write(*,*) vartype, " character(len=:),allocatable,dimension(:)"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(id, var, itemsize)
    status = bmif_get_var_nbytes(id, var, nbytes)
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
    bmif_get_value_char1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_char1

    !> \overload
    INTEGER FUNCTION bmif_get_value_double(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8"
        stop "Variable type error."
    endif
    bmif_get_value_double = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
    END FUNCTION bmif_get_value_double

    !> \overload
    INTEGER FUNCTION bmif_get_value_double1(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), allocatable, dimension(:), INTENT(inout) :: dest
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8 1d"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(id, var, itemsize)
    status = bmif_get_var_nbytes(id, var, nbytes)
    dim = nbytes / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
    endif
    bmif_get_value_double1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_double1

    !> \overload
    INTEGER FUNCTION bmif_get_value_double2(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), allocatable, INTENT(inout) :: dest(:,:)
    character(100) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        write(*,*) vartype, " real(kind=8 2d"
        stop "Variable type error."
    endif
    varname = Lower(var)
    need_alloc = .true.
    if (varname .eq. "concentrations") then
        dim1 = RM_GetGridCellCount(id)
        dim2 = RM_GetComponentCount(id)
    else if (varname .eq. "selectedoutput") then
        dim1 = RM_GetSelectedOutputRowCount(id)
        dim2 = RM_GetSelectedOutputColumnCount(id)
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
    bmif_get_value_double2 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION bmif_get_value_double2

    !> \overload
    INTEGER FUNCTION bmif_get_value_int(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer") then
        write(*,*) vartype, " integer"
        stop "Variable type error."
    endif
    bmif_get_value_int = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
    END FUNCTION bmif_get_value_int

    !> \overload
    INTEGER FUNCTION bmif_get_value_int1(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, allocatable, INTENT(inout) :: dest(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:)") then
        write(*,*) vartype, " integer 1d"
        stop "Variable type error."
    endif
    status = bmif_get_var_itemsize(id, var, itemsize)
    status = bmif_get_var_nbytes(id, var, nbytes)  
    dim = nbytes / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
    endif
    bmif_get_value_int1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
    END FUNCTION bmif_get_value_int1

    !> \overload
    INTEGER FUNCTION bmif_get_value_int2(id, var, dest)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, allocatable, INTENT(inout) :: dest(:,:)
    character(100) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer") then
        write(*,*) vartype, " integer 2d"
        stop "Variable type error."
    endif
    varname = Lower(varname)
    need_alloc = .true.
    stop "Unknown 2d variable"
    bmif_get_value_int2 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION bmif_get_value_int2
   
    !> @a bmif_get_value_ptr retrieves pointers to model variables. Only variables in the list
    !> provided by @ref bmif_get_pointable_var_names can be pointed to. 
    !> @param id     The instance @a id returned from @ref bmif_create.
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
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.     
!Add NEW_VARIABLE to bmif_get_value_ptr Documentation  
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_double(id, var, ptr)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=c_double), pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer :: status
    status = RMF_BMI_GetValuePtr(id, trim(var)//C_NULL_CHAR, src)
    call C_F_POINTER(src, ptr)
    bmif_get_value_ptr_double = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_double
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_double1(id, var, ptr)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=c_double), pointer, INTENT(inout) :: ptr(:)
    type (c_ptr) :: src
    integer nbytes, itemsize, dim, status
    status = bmif_get_var_nbytes(id, var, nbytes)
    status = bmif_get_var_itemsize(id, var, itemsize)
    dim = nbytes/itemsize
    status = RMF_BMI_GetValuePtr(id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr, [dim])
    bmif_get_value_ptr_double1 = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_double1
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_integer(id, var, ptr)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer status
    status = RMF_BMI_GetValuePtr(id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr)
    bmif_get_value_ptr_integer = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_integer

    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_logical(id, var, ptr)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    logical(kind=1), pointer, INTENT(inout) :: ptr
    type (c_ptr) :: src
    integer status
    status = RMF_BMI_GetValuePtr(id, trim(var)//C_NULL_CHAR, src)
    call c_f_pointer(src, ptr)
    bmif_get_value_ptr_logical = success(status)
    return 
    END FUNCTION bmif_get_value_ptr_logical    

    ! ====================================================
    ! Setters, by type
    ! ====================================================

    !> @a bmif_set_value sets model variables. Only variables in the list
    !> provided by @ref bmif_get_input_var_names can be set. 
    !> @param id     The instance id returned from @ref bmif_create.
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
    !> status = bmif_set_value(id, "Temperature", tc)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
    !> \overload
    INTEGER FUNCTION bmif_set_value_b(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL, INTENT(in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim, src_int
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "logical") then
        stop "Variable type error."
    endif
    src_int = 1
    if(.not. src) src_int = 0
    bmif_set_value_b = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src_int)
    return
    END FUNCTION bmif_set_value_b

    !> \overload
    INTEGER FUNCTION bmif_set_value_c(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "character") then
        stop "Variable type error."
    endif
    bmif_set_value_c = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, trim(src)//C_NULL_CHAR)
    return
    END FUNCTION bmif_set_value_c

    !> \overload
    INTEGER FUNCTION bmif_set_value_int(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src
    integer :: src_copy
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer") then
        stop "Variable type error."
    endif
    src_copy = src
    if (var .eq. "NthSelectedOutput") src_copy = src - 1
    bmif_set_value_int = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src_copy)
    return
    END FUNCTION bmif_set_value_int

    !> \overload
    INTEGER FUNCTION bmif_set_value_int1(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(id, var, nbytes)
    status =  bmif_get_var_itemsize(id, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_int1 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src(1))
    return
    END FUNCTION bmif_set_value_int1

    !> \overload
    INTEGER FUNCTION bmif_set_value_int2(id, var, src)
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
    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT( in) :: src(:,:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "integer,allocatable,dimension(:,:)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(id, var, nbytes)
    status = bmif_get_var_itemsize(id, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src,1)*size(src,2)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_int2 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src(1,1))
    return
    END FUNCTION bmif_set_value_int2

    !> \overload
    INTEGER FUNCTION bmif_set_value_double(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    bmif_set_value_double = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src)
    return
    END FUNCTION bmif_set_value_double

    !> \overload
    INTEGER FUNCTION bmif_set_value_double1(id, var, src)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src(:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(id, var, nbytes)
    status = bmif_get_var_itemsize(id, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_double1 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src(1))
    return
    END FUNCTION bmif_set_value_double1

    !> \overload
    INTEGER FUNCTION bmif_set_value_double2(id, var, src)
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
    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT( in) :: src(:,:)
    character(100) :: vartype
    integer :: nbytes, status, dim, itemsize
    status = bmif_get_var_type(id, var, vartype)
    if (vartype .ne. "real(kind=8)") then
        stop "Variable type error."
    endif
    status = bmif_get_var_nbytes(id, var, nbytes)
    status = bmif_get_var_itemsize(id, var, itemsize)
    dim = nbytes / itemsize
    if (dim .ne. size(src,1)*size(src,2)) then
        stop "Variable dimension error"
    endif
    bmif_set_value_double2 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, src(1,1))
    return
    END FUNCTION bmif_set_value_double2

    ! ====================================================
    ! Grid information
    ! ====================================================

    !> @a bmif_get_grid_rank returns a rank of 1 for grid 0. 
    !> BMIPhreeqcRM has a 1D series of
    !> cells; any grid or spatial information must
    !> be found in the user's model.
    !> @param id     The instance id returned from @ref bmif_create.
    !> @param grid   Grid number, only grid 0 is considered.
    !> @param rank   Rank of 1 is returned for grid 0; 0 for
    !> all other values of @a grid.
    !> @retval       0 is success, 1 is failure.    
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_grid_rank(id, grid, rank)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_rank(id, grid, rank)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id, grid
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
    !> @param id    The instance id returned from @ref bmif_create.
    !> @param grid  Grid number, only grid 0 is considered.
    !> @param ngrid Same value as @ref RM_GetGridCellCount 
    !> or @ref bmif_get_value "GridCellCount" is returned for grid 0; 
    !> 0 for all other values of @a grid.
    !> @retval       0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_grid_size(id, grid, ngrid)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_size(id, grid, ngrid)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id, grid
    INTEGER, INTENT(inout) :: ngrid
    if (grid .eq. 0) then
        bmif_grid_size = success(bmif_get_value(id, "GridCellCount", ngrid))
    else
        ngrid = 0
        bmif_grid_size = BMI_FAILURE
    endif
    END FUNCTION bmif_grid_size

    !> @a bmif_grid_type defines the grid to be points. No grid
    !> information is available in BMIPhreeqcRM; all grid 
    !> information must be found in the user's model.
    !> @param id    The instance id returned from @ref bmif_create.
    !> @param grid  Grid number, only grid 0 is considered.
    !> @param str   "Points" is returned for grid 0;
    !> "Undefined grid identifier" is returned for all other 
    !> values of @a grid.
    !> @retval      0 is success, 1 is failure.
    !> @par Fortran example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = bmif_grid_type(id, grid, str)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION bmif_grid_type(id, grid, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id, grid
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
!> the PhreeqcRM method FindComponents. FindComponents generates SELECTED_OUTPUT 333 and
!> USER_PUNCH 333 data blocks that make the variables accessible. Variables will
!> only be accessible if the system includes the given reactant; for example, no
!> gas variables will be created if there are no GAS_PHASEs in the model. 
!>
!> @param id    The instance id returned from @ref bmif_create.
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
    INTEGER FUNCTION bmif_add_output_vars(id, option, def)
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
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: option
    CHARACTER(len=*), INTENT(in) :: def
    bmif_add_output_vars = RMF_BMI_AddOutputVars(id, trim(option)//C_NULL_CHAR, trim(def)//C_NULL_CHAR)
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
    
    END MODULE BMIPhreeqcRM
