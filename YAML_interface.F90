#if defined(USE_YAML) 
    !> @file YAML_interface.F90
    !> @brief YAMLPhreeqcRM module definition 
    !>
    !
    !*MODULE YAMLPhreeqcRM Helper module for building YAML initialization files. 
    !> @brief Fortran documentation for using YAML to initialize instances 
	!> of BMIPhreeqcRM and PhreeqcRM.
    !> @par ""
    !> "USE YAMLPhreeqcRM" defines a module that can be used in preprocessors or
	!> Graphical User Interfaces to store initialization data for BMIPhreeqcRM or
	!> PhreeqcRM instances. PhreeqcRM methods and data can be stored in a YAML file.
	!> After an instance of BMIPhreeqcRM or PhreeqcRM has been created, the method
	!> bmif_initialize or RM_InitializeYAML can be used to run the specified methods
	!> with the specified data to define properties and initial conditions for 
	!> the instance.
    !>    
MODULE YAMLPhreeqcRM
    contains
!> Creates a YAMLPhreeqcRM instance with a YAML document that is ready to 
!> for writing data for initiation of a PhreeqcRM instance.
!> @retval id   Id of the new YAMLPhreeqcRM instance.
!> @see
!> @ref DestroyYAMLPhreeqcRM.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! Create YAMLPhreeqcRM document
!> id = CreateYAMLPhreeqcRM() 
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION CreateYAMLPhreeqcRM()
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION CreateYAMLPhreeqcRM_F() &
			BIND(C, NAME='CreateYAMLPhreeqcRM_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
		END FUNCTION CreateYAMLPhreeqcRM_F
    END INTERFACE
	CreateYAMLPhreeqcRM = CreateYAMLPhreeqcRM_F()
    END FUNCTION CreateYAMLPhreeqcRM
!> Deletes the YAMLPhreeqcRM instance and all data.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @see
!> @ref YAMLClear.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> YAML_filename = "AdvectBMI_f90.yaml"
!> status = WriteYAMLDoc(id, YAML_filename)
!> status = YAMLClear(id)  
!> status = DestroyYAMLPhreeqcRM(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly 
    INTEGER FUNCTION DestroyYAMLPhreeqcRM(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION DestroyYAMLPhreeqcRM_F(id) &
			BIND(C, NAME='DestroyYAMLPhreeqcRM_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION DestroyYAMLPhreeqcRM_F
    END INTERFACE
    integer, intent(in) :: id
	DestroyYAMLPhreeqcRM = DestroyYAMLPhreeqcRM_F(id)
    END FUNCTION DestroyYAMLPhreeqcRM
!> Writes YAML document to file.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param file_name     Name of file to write YAML document.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @see
!> @ref DestroyYAMLPhreeqcRM, @ref YAMLClear.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> YAML_filename = "AdvectBMI_f90.yaml"
!> status = WriteYAMLDoc(id, YAML_filename)
!> status = YAMLClear(id)  
!> status = DestroyYAMLPhreeqcRM(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION WriteYAMLDoc(id, file_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION WriteYAMLDoc_F(id, file_name) &
			BIND(C, NAME='WriteYAMLDoc_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: file_name(*)
		END FUNCTION WriteYAMLDoc_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: file_name
	WriteYAMLDoc = WriteYAMLDoc_F(id, trim(file_name)//C_NULL_CHAR)
    END FUNCTION WriteYAMLDoc
    
!> Clears all definitions from the YAML document.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @see
!> @ref DestroyYAMLPhreeqcRM.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> YAML_filename = "AdvectBMI_f90.yaml"
!> status = WriteYAMLDoc(id, YAML_filename)
!> status = YAMLClear(id)  
!> status = DestroyYAMLPhreeqcRM(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLClear(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLClear_F(id) &
			BIND(C, NAME='YAMLClear_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION YAMLClear_F
    END INTERFACE
    integer, intent(in) :: id
	YAMLClear = YAMLClear_F(id)
    END FUNCTION YAMLClear
!> Inserts data into the YAML document to select sets of output variables.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance. Sets of variables can be included or excluded with
!> multiple calls to this method. All calls must precede the final call to
!> @ref YAMLFindComponents. FindComponents generates SELECTED_OUTPUT 333 and
!> USER_PUNCH 333 data blocks that make the variables accessible. Variables will
!> only be accessible if the system includes the given reactant; for example, no
!> gas variables will be created if there are no GAS_PHASEs in the model. 
!> @param id        The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option    A string value, among those listed below, that selects sets of variables 
!> that can be retieved by the bmif_get_value method. 
!> @param def A string value that can be "false", "true", or a list of items to be included as
!> accessible variables. A value of "false", excludes all variables of the given type; a 
!> value of "true" includes all variables of the given type for the current system; a list
!> specifies a subset of items of the given type. 
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
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
    INTEGER FUNCTION YAMLAddOutputVars(id, option, def)
    USE ISO_C_BINDING
    IMPLICIT NONE
		INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLAddOutputVars_F(id, var, src) &
			BIND(C, NAME='YAMLAddOutputVars_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
		INTEGER(KIND=C_INT), INTENT(in) :: id
		CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
		CHARACTER(KIND=C_CHAR), INTENT(in) :: src
		END FUNCTION YAMLAddOutputVars_F
		END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: option
    CHARACTER(len=*), INTENT(in) :: def
    character(100) :: vartype
    integer :: bytes, nbytes, status, dim
    YAMLAddOutputVars = YAMLAddOutputVars_F(id, trim(option)//C_NULL_CHAR, trim(def)//C_NULL_CHAR)
    return
    END FUNCTION YAMLAddOutputVars  
!> Inserts data into the YAML document for the PhreeqcRM method CloseFiles.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> CloseFiles closes the output and log files.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLCloseFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLCloseFiles(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLCloseFiles_F(id) &
			BIND(C, NAME='YAMLCloseFiles_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION YAMLCloseFiles_F
    END INTERFACE
    integer, intent(in) :: id
	YAMLCloseFiles = YAMLCloseFiles_F(id)
    END FUNCTION YAMLCloseFiles
!> Inserts data into the YAML document for the PhreeqcRM method CreateMapping.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param grid2chem     Integer array of mapping from user's model grid to cells
!> for which chemistry will be run. 
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> CreateMapping
!> provides a mapping from grid cells in the user's model to reaction cells for which chemistry needs to be run.
!> The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells
!> for which chemistry must be run.
!> The array @a grid2chem of size @a nxyz (the number of grid cells)
!> must contain the set of all integers 0 <= @a i < @a count_chemistry,
!> where @a count_chemistry is a number less than or equal to @a nxyz.
!> Inactive cells are assigned a negative integer.
!> The mapping may be many-to-one to account for symmetry.
!> Default is a one-to-one mapping--all user grid cells are reaction cells
!> (equivalent to @a grid2chem values of 0,1,2,3,...,nxyz-1).
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! Demonstation of mapping, two equivalent rows by symmetry
!> ! zero-based indexing
!> integer, allocatable, dimension(:) :: grid2chem
!> allocate(grid2chem(nxyz))
!> grid2chem = -1
!> do i = 1, nxyz / 2 
!> 	 grid2chem(i) = i - 1
!> 	 grid2chem(i + nxyz / 2) = i - 1
!> enddo
!> status = YAMLCreateMapping(id, grid2chem)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLCreateMapping(id, grid2chem)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLCreateMapping_F(id, grid2chem, l) &
			BIND(C, NAME='YAMLCreateMapping_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: grid2chem
        integer(kind=C_INT), intent(in) :: l
		END FUNCTION YAMLCreateMapping_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: grid2chem
	YAMLCreateMapping = YAMLCreateMapping_F(id, grid2chem(1), size(grid2chem))
    END FUNCTION YAMLCreateMapping
!> Inserts data into the YAML document for the PhreeqcRM method DumpModule.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param dump_on          Signal for writing the dump file, true or false.
!> @param append           Signal to append to the contents of the dump file, true or false.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> DumpModule writes the contents of all workers to file in _RAW formats (see appendix of PHREEQC version 3 manual),
!> including SOLUTIONs and all reactants.
!> @see                    @ref YAMLSetDumpFileName.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> logical dump_on, append
!> dump_on = .true.
!> append = .false.
!> status = YAMLSetDumpFileName(id, "Advect_cpp.dmp")
!> status = YAMLDumpModule(id, dump_on, append)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLDumpModule(id, dump_on, append)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLDumpModule_F(id, idump_on, iappend) &
			BIND(C, NAME='YAMLDumpModule_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer, intent(in) :: id
        integer(kind=C_INT), intent(in) :: idump_on, iappend
		END FUNCTION YAMLDumpModule_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: dump_on, append  
    integer :: idump_on, iappend
	idump_on = 0
	iappend = 0
    if (dump_on) idump_on = 1
    if (append) iappend = 1
	YAMLDumpModule = YAMLDumpModule_F(id, idump_on, iappend)
    END FUNCTION YAMLDumpModule
!> Inserts data into the YAML document for the PhreeqcRM method FindComponents.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id    The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval      Zero indicates success, negative indicates failure.
!> @par
!> FindComponents accumulates a list of elements. Elements are those that have been
!> defined in a solution or any other reactant
!> (EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
!> This method can be called multiple times and the list that is created is cummulative.
!> The list is the set of components that needs to be transported. By default the list
!> includes water, excess H and excess O (the H and O not contained in water);
!> alternatively, the list may be set to contain total H and total O (@ref YAMLSetComponentH2O),
!> which requires transport results to be accurate to eight or nine significant digits.
!> If multicomponent diffusion (MCD) is to be modeled,
!> there is a capability to retrieve aqueous species concentrations 
!> and to set new solution concentrations after
!> MCD by using individual species concentrations
!> (@ref YAMLSpeciesConcentrations2Module).
!> To use these methods, the save-species property needs to be turned on (@ref YAMLSetSpeciesSaveOn).
!> If the save-species property is on, FindComponents will generate
!> a list of aqueous species,
!> their diffusion coefficients at 25 C,
!> and their charge.
!> @see
!> @ref YAMLSetComponentH2O,
!> @ref YAMLSetSpeciesSaveOn,
!> @ref YAMLSpeciesConcentrations2Module.
!> @par 
!> The FindComponents method also generates lists of reactants--equilibrium phases,
!> exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
!> The lists are cumulative, including all reactants that were
!> defined in the initial phreeqc instance at any time FindComponents was called.
!> In addition, a list of phases is generated for which saturation indices may be calculated from the
!> cumulative list of components.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLRunFile(id, yaml_file)
!> status = YAMLFindComponents(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly    
    INTEGER FUNCTION YAMLFindComponents(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLFindComponents_F(id) &
			BIND(C, NAME='YAMLFindComponents_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION YAMLFindComponents_F
    END INTERFACE
    integer, intent(in) :: id
	YAMLFindComponents = YAMLFindComponents_F(id)
    END FUNCTION YAMLFindComponents    
    
!> Inserts data into the YAML document for the PhreeqcRM method InitialSolutions2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param solutions   Vector of SOLUTION index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval            Zero indicates success, negative indicates failure.
!> @par
!> InitialSolutions2Module transfers SOLUTION definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a solutions is a vector of SOLUTION index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialSolutions2Module(id, solutions)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialSolutions2Module_F(id, solutions, dim) &
			BIND(C, NAME='YAMLInitialSolutions2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: solutions(*)
		END FUNCTION YAMLInitialSolutions2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: solutions
	integer :: dim
	dim = size(solutions)
	YAMLInitialSolutions2Module = YAMLInitialSolutions2Module_F(id, solutions, dim)
    END FUNCTION YAMLInitialSolutions2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialEquilibriumPhases2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                   The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param equilibrium_phases   Vector of EQUILIBRIUM_PHASES index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                     Zero indicates success, negative indicates failure.
!> @par
!> InitialEquilibriumPhases2Module transfers EQUILIBRIUM_PHASES definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a equilibrium_phases is a vector of EQUILIBRIUM_PHASES index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialEquilibriumPhases2Module(id, equilibrium_phases)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialEquilibriumPhases2Module_F(id, equilibrium_phases, dim) &
			BIND(C, NAME='YAMLInitialEquilibriumPhases2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: equilibrium_phases(*)
		END FUNCTION YAMLInitialEquilibriumPhases2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: equilibrium_phases
	integer :: dim
	dim = size(equilibrium_phases)
	YAMLInitialEquilibriumPhases2Module = YAMLInitialEquilibriumPhases2Module_F(id, equilibrium_phases, dim)
    END FUNCTION YAMLInitialEquilibriumPhases2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialExchanges2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                   The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param exchanges            Vector of EXCHANGE index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                     Zero indicates success, negative indicates failure.
!> @par
!> InitialExchanges2Module transfers EXCHANGE definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a exchanges is a vector of EXCHANGE index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialExchanges2Module(id, exchanges)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialExchanges2Module_F(id, exchanges, dim) &
			BIND(C, NAME='YAMLInitialExchanges2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: exchanges(*)
		END FUNCTION YAMLInitialExchanges2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: exchanges
	integer :: dim
	dim = size(exchanges)
	YAMLInitialExchanges2Module = YAMLInitialExchanges2Module_F(id, exchanges, dim)
    END FUNCTION YAMLInitialExchanges2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialSurfaces2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                  The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param surfaces            Vector of SURFACE index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                     Zero indicates success, negative indicates failure.
!> @par
!> InitialSurfaces2Module transfers SURFACE definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a surfaces is a vector of SURFACE index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialSurfaces2Module(id, surfaces)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialSurfaces2Module_F(id, surfaces, dim) &
			BIND(C, NAME='YAMLInitialSurfaces2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: surfaces(*)
		END FUNCTION YAMLInitialSurfaces2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: surfaces
	integer :: dim
	dim = size(surfaces)
	YAMLInitialSurfaces2Module = YAMLInitialSurfaces2Module_F(id, surfaces, dim)
    END FUNCTION YAMLInitialSurfaces2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialGasPhases2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                  The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param gas_phases          Vector of GAS_PHASE index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                    Zero indicates success, negative indicates failure.
!> @par
!> InitialGasPhases2Module transfers GAS_PHASE definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a gas_phases is a vector of GAS_PHASE index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialGasPhases2Module(id, gas_phases)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialGasPhases2Module_F(id, gas_phases, dim) &
			BIND(C, NAME='YAMLInitialGasPhases2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: gas_phases(*)
		END FUNCTION YAMLInitialGasPhases2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: gas_phases
	integer :: dim
	dim = size(gas_phases)
	YAMLInitialGasPhases2Module = YAMLInitialGasPhases2Module_F(id, gas_phases, dim)
    END FUNCTION YAMLInitialGasPhases2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialSolidSolutions2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                  The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param solid_solutions     Vector of SOLID_SOLUTIONS index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                    Zero indicates success, negative indicates failure.
!> @par
!> InitialSolidSolutions2Module transfers SOLID_SOLUTIONS definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a solid_solutions is a vector of SOLID_SOLUTIONS index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialSolidSolutions2Module(id, solid_solutions)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialSolidSolutions2Module_F(id, solid_solutions, dim) &
			BIND(C, NAME='YAMLInitialSolidSolutions2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id, dim
        integer(kind=C_INT), intent(in) :: solid_solutions(*)
		END FUNCTION YAMLInitialSolidSolutions2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: solid_solutions
	integer :: dim
	dim = size(solid_solutions)
	YAMLInitialSolidSolutions2Module = YAMLInitialSolidSolutions2Module_F(id, solid_solutions, dim)
    END FUNCTION YAMLInitialSolidSolutions2Module
	
!> Inserts data into the YAML document for the PhreeqcRM method InitialKinetics2Module.
!> When the YAML document is written to file it can be processed by the method 
!> InitializeYAML to initialize a PhreeqcRM instance.
!> @param id                  The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param kinetics            Vector of KINETICS index numbers that is dimensioned @a nxyz, 
!> where @a nxyz is the number of grid cells in the user's model. 
!> @retval                    Zero indicates success, negative indicates failure.
!> @par
!> InitialKinetics2Module transfers KINETICS definitions from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a kinetics is a vector of KINETICS index numbers that refer to
!> definitions in the InitialPhreeqc instance.
    INTEGER FUNCTION YAMLInitialKinetics2Module(id, kinetics)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialKinetics2Module_F(id, kinetics, dim) &
			BIND(C, NAME='YAMLInitialKinetics2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: kinetics(*)
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLInitialKinetics2Module_F
    END INTERFACE
    integer, intent(in) :: id
	integer :: dim
    integer, allocatable, dimension(:), intent(in) :: kinetics
	dim = size(kinetics)
	YAMLInitialKinetics2Module = YAMLInitialKinetics2Module_F(id, kinetics, dim)
    END FUNCTION YAMLInitialKinetics2Module
    
!> Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id   The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param ic1  Vector of solution and reactant index numbers that refer to
!> definitions in the InitialPhreeqc instance.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc 
!> instance to the reaction-module workers.
!> @a ic1 is used to select initial conditions, including solutions and reactants,
!> for each cell of the model, without mixing.
!> @a ic1 is dimensioned 7 times @a nxyz, where @a nxyz is the 
!> number of grid cells in the user's model. 
!> The dimension of 7 refers to solutions and reactants in the following order:
!> (0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
!> (5) SOLID_SOLUTIONS, and (6) KINETICS.
!> The definition initial_solution1[3*nxyz + 99] = 2, indicates that
!> cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 
!> in the InitialPhreeqc instance.
!> 
!> Size is 7 times @a nxyz. The order of definitions is given above.
!> Negative values are ignored, resulting in no definition of that entity for that cell.
!> @see                        @ref YAMLInitialPhreeqcCell2Module, 
!> @ref YAMLInitialPhreeqc2Module_mix.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> integer, allocatable, dimension(:,:) :: ic1
!> allocate(ic1(nxyz,7))
!> do i = 1, nxyz
!> 	ic1(i,1) = 1     ! Solution 1
!> 	ic1(i,2) = -1    ! Equilibrium phases none
!> 	ic1(i,3) = 1     ! Exchange 1
!> 	ic1(i,4) = -1    ! Surface none
!> 	ic1(i,5) = -1    ! Gas phase none
!> 	ic1(i,6) = -1    ! Solid solutions none
!> 	ic1(i,7) = -1    ! Kinetics none
!> enddo
!> status = YAMLInitialPhreeqc2Module_mix(id, ic1) 
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLInitialPhreeqc2Module(id, ic1)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialPhreeqc2Module_F(id, ic1, dim) &
			BIND(C, NAME='YAMLInitialPhreeqc2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: ic1
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLInitialPhreeqc2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:,:), intent(in) :: ic1
    integer :: l
    l = size(ic1,1)*size(ic1,2)
	YAMLInitialPhreeqc2Module = YAMLInitialPhreeqc2Module_F(id, ic1(1,1), l)
    END FUNCTION YAMLInitialPhreeqc2Module  
!> Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param ic1    Vector of solution and reactant index numbers that refer to
!> definitions in the InitialPhreeqc instance.
!> Size is 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model.
!> The order of reactants is given below and in the example.
!> Negative values are ignored, resulting in no definition of that entity for that cell.
!> @param ic2    Vector of solution and reactant index numbers that refer to
!> definitions in the InitialPhreeqc instance.
!> Nonnegative values of @a ic2 result in mixing with the entities defined in @a ic1.
!> Negative values result in no mixing.
!> Size is 7 times @a nxyz. 
!> @param f1           Fraction of @a ic1 that mixes with (1 - @a f1)
!> of @a ic2.
!> Size is 7 times @a nxyz. 
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc instance to
!> the reaction-module workers, possibly with mixing.
!> In its simplest form, @a  ic1 is used to select initial conditions, including solutions and reactants,
!> for each cell of the model, without mixing.
!> The dimension of 7 refers to solutions and reactants in the following order:
!> (0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
!> (5) SOLID_SOLUTIONS, and (6) KINETICS.
!> The definition ic1[3*nxyz + 99] = 2, indicates that
!> cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 
!> in the InitialPhreeqc instance (either by RunFile or RunString).
!> @n@n
!> It is also possible to mix solutions and reactants to obtain the initial conditions 
!> for cells. For mixing,
!> @a initials_conditions2 contains numbers for a second entity that mixes with 
!> the entity defined in @a ic1.
!> @a f1 contains the mixing fraction for @a ic1,
!> whereas (1 - @a f1) is the mixing fraction for @a ic2.
!> The definitions ic1[3*nxyz + 99] = 2, initial_solution2[3*nxyz + 99] = 3,
!> f1[3*nxyz + 99] = 0.25 indicates that
!> cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3,
!> where the surface compositions have been defined in the InitialPhreeqc instance.
!> If the user number in @a ic2 is negative, no mixing occurs.
!> 
!> @see                        @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLInitialPhreeqc2Module.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> integer, allocatable, dimension(:,:) :: ic1
!> integer, allocatable, dimension(:,:) :: ic2
!> real(kind=8), allocatable, dimension(:,:) :: f1
!> allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
!> ic1 = -1
!> ic2 = -1
!> f1 = 1.0d0
!> do i = 1, nxyz
!> 	ic1(i,1) = 1     ! Solution 1
!> 	ic1(i,2) = -1    ! Equilibrium phases none
!> 	ic1(i,3) = 1     ! Exchange 1
!> 	ic1(i,4) = -1    ! Surface none
!> 	ic1(i,5) = -1    ! Gas phase none
!> 	ic1(i,6) = -1    ! Solid solutions none
!> 	ic1(i,7) = -1    ! Kinetics none
!> enddo
!> status = YAMLInitialPhreeqc2Module_mix(id, ic1, ic2, f1) 
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLInitialPhreeqc2Module_mix(id, ic1, ic2, f1)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialPhreeqc2Module_mix_F(id, ic1, ic2, f1, dim) &
			BIND(C, NAME='YAMLInitialPhreeqc2Module_mix_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: ic1
        integer(kind=C_INT), intent(in) :: ic2
        real(kind=C_DOUBLE), intent(in) :: f1
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLInitialPhreeqc2Module_mix_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:,:), intent(in) :: ic1
    integer, allocatable, dimension(:,:), intent(in) :: ic2
    real(kind=8), allocatable, dimension(:,:), intent(in) :: f1
    integer :: l1, l2, l3
    l1 = size(ic1,1)*size(ic1,2)
    l2 = size(ic2,1)*size(ic2,2)
    l3 = size(f1,1)*size(f1,2)
    if ((l1 .ne. l2) .or. (l1 .ne. l3)) then
        stop "Dimension error in YAMLInitialPhreeqc2Module"
    endif
	YAMLInitialPhreeqc2Module_mix = YAMLInitialPhreeqc2Module_mix_F(id, ic1(1,1), ic2(1,1), f1(1,1), l1)
    END FUNCTION YAMLInitialPhreeqc2Module_mix  
!> Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqcCell2Module.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n                  Number that refers to a solution or MIX and associated 
!> reactants in the InitialPhreeqc instance.
!> @param cell_numbers       A vector of grid-cell numbers (user's grid-cell numbering system) that
!> will be populated with cell @a n from the InitialPhreeqc instance.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> InitialPhreeqcCell2Module uses a cell numbered @a n in the InitialPhreeqc instance to 
!> populate a series of transport cells.
!> All reactants with the number @a n are transferred along with the solution.
!> If MIX @a n exists, it is used for the definition of the solution.
!> If @a n is negative, @a n is redefined to be the largest solution or MIX number 
!> in the InitialPhreeqc instance.
!> All reactants for each cell in the list @a cell_numbers are removed before the cell
!> definition is copied from the InitialPhreeqc instance to the workers.
!> 
!> @see                      @ref YAMLInitialPhreeqc2Module,
!> @ref YAMLInitialPhreeqc2Module_mix.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> integer, allocatable, dimension(:) :: module_cells
!> allocate(module_cells(2))
!> module_cells(1) = 18
!> module_cells(2) = 19
!> status = YAMLInitialPhreeqcCell2Module(id, -1, module_cells)
!> </PRE>
!> </CODE>
!> @endhtmlonly  
    INTEGER FUNCTION YAMLInitialPhreeqcCell2Module(id, n, cell_numbers)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLInitialPhreeqcCell2Module_F(id, n, cell_numbers, dim) &
			BIND(C, NAME='YAMLInitialPhreeqcCell2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
        integer(kind=C_INT), intent(in) :: cell_numbers
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLInitialPhreeqcCell2Module_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
    integer, allocatable, dimension(:), intent(in) :: cell_numbers
    integer :: l
    l = size(cell_numbers)
	YAMLInitialPhreeqcCell2Module = YAMLInitialPhreeqcCell2Module_F(id, n, cell_numbers(1), l)
    END FUNCTION YAMLInitialPhreeqcCell2Module    
!> Inserts data into the YAML document for the PhreeqcRM method LoadDatabase.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param file_name         String containing the database name.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> LoadDatabase loads a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
!> of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLLoadDatabase(id, "phreeqc.dat") 
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLLoadDatabase(id, file_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLLoadDatabase_F(id, file_name) &
			BIND(C, NAME='YAMLLoadDatabase_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: file_name(*)
		END FUNCTION YAMLLoadDatabase_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: file_name
	YAMLLoadDatabase = YAMLLoadDatabase_F(id, trim(file_name)//C_NULL_CHAR)
    END FUNCTION YAMLLoadDatabase
!> Inserts data into the YAML document for the PhreeqcRM method LogMessage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param str              String to be printed.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> LogMessage prints a message to the log file.
!> @see                    @ref YAMLOutputMessage, @ref YAMLScreenMessage, @ref YAMLWarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLLogMessage(id, "Finished section 1 of initialization")
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLLogMessage(id, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLLogMessage_F(id, str) &
			BIND(C, NAME='YAMLLogMessage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: str(*)
		END FUNCTION YAMLLogMessage_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: str
	YAMLLogMessage = YAMLLogMessage_F(id, trim(str)//C_NULL_CHAR)
    END FUNCTION YAMLLogMessage
!> Inserts data into the YAML document for the PhreeqcRM method OpenFiles.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> OpenFiles opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
!> based on the prefix defined by @ref YAMLSetFilePrefix.
!> @see                    @ref YAMLSetFilePrefix, @ref YAMLCloseFiles,
!> @ref YAMLLogMessage, @ref YAMLOutputMessage, and @ref YAMLWarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetFilePrefix(id, "Advect_cpp")
!> status = YAMLOpenFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly    
    INTEGER FUNCTION YAMLOpenFiles(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLOpenFiles_F(id) &
			BIND(C, NAME='YAMLOpenFiles_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION YAMLOpenFiles_F
    END INTERFACE
    integer, intent(in) :: id
	YAMLOpenFiles = YAMLOpenFiles_F(id)
    END FUNCTION YAMLOpenFiles   
!> Inserts data into the YAML document for the PhreeqcRM method OutputMessage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param str              String to be printed.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> OutputMessage prints a message to the output file.
!> @see                    @ref YAMLLogMessage, @ref YAMLScreenMessage, @ref YAMLWarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLOutputMessage(id, "Finished section 1 of initialization")
!> </PRE>
!> </CODE>
!> @endhtmlonly    
    INTEGER FUNCTION YAMLOutputMessage(id, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLOutputMessage_F(id, str) &
			BIND(C, NAME='YAMLOutputMessage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: str(*)
		END FUNCTION YAMLOutputMessage_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: str
	YAMLOutputMessage = YAMLOutputMessage_F(id, trim(str)//C_NULL_CHAR)
    END FUNCTION YAMLOutputMessage
!> Inserts data into the YAML document for the PhreeqcRM method RunCells.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> RunCells runs reactions for all cells in the reaction module.
!> During initialization, RunCells can be used to equilibrate each solution with all 
!> reactants in a cell while
!> using a time step of zero (@ref YAMLSetTimeStep) to avoid kinetic reactions.
!> Other properties that may need to be initialized before RunCells is invoked
!> include porosity (@ref YAMLSetPorosity),
!> saturation (@ref YAMLSetSaturationUser),
!> temperature (@ref YAMLSetTemperature), and pressure (@ref YAMLSetPressure).
!> 
!> @see                    @ref YAMLSetPorosity,
!> @ref YAMLSetPressure, @ref YAMLSetSaturationUser, @ref YAMLSetTemperature, @ref YAMLSetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetTimeStep(id, 0.0)
!> status = YAMLRunCells(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLRunCells(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLRunCells_F(id) &
			BIND(C, NAME='YAMLRunCells_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
		END FUNCTION YAMLRunCells_F
    END INTERFACE
    integer, intent(in) :: id
	YAMLRunCells = YAMLRunCells_F(id)
    END FUNCTION YAMLRunCells       
!> Inserts data into the YAML document for the PhreeqcRM method RunFile.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param workers          @a True, the workers will run the file; @a False, the workers will not run the file.
!> @param initial_phreeqc  @a True, the InitialPhreeqc instance will run the file; @a False, the InitialPhreeqc will not run the file.
!> @param utility          @a True, the Utility instance will run the file; @a False, the Utility instance will not run the file.
!> @param file_name        Name of the file to run.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> RunFile runs a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
!> the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
!> files that modify the thermodynamic database should be run by all three sets of instances.
!> Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
!> be run by the workers. Files that contain initial conditions or boundary conditions should
!> be run by the InitialPhreeqc instance.
!> @see                    @ref YAMLRunString.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLRunFile(id, .true., .true., .true., "advect.pqi")
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLRunFile(id, workers, initial_phreeqc, utility, file_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLRunFile_F(id, iworkers, iinitial_phreeqc, iutility, file_name) &
			BIND(C, NAME='YAMLRunFile_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: iworkers, iinitial_phreeqc, iutility
        character(KIND=C_CHAR), intent(in) :: file_name(*)
		END FUNCTION YAMLRunFile_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: workers, initial_phreeqc, utility
    integer :: iworkers, iinitial_phreeqc, iutility
    character(len=*), intent(in) :: file_name
	iworkers = 0
	iinitial_phreeqc = 0
	iutility = 0
	if (workers) iworkers = 1
	if (initial_phreeqc) iinitial_phreeqc = 1
	if (utility) iutility = 1

	YAMLRunFile = YAMLRunFile_F(id, iworkers, iinitial_phreeqc, iutility, trim(file_name)//C_NULL_CHAR)
    END FUNCTION YAMLRunFile   
!> Inserts data into the YAML document for the PhreeqcRM method RunString.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance. 
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param workers          @a True, the workers will run the string; @a False, the workers will not run the string.
!> @param initial_phreeqc  @a True, the InitialPhreeqc instance will run the string; @a False, the InitialPhreeqc will not run the string.
!> @param utility          @a True, the Utility instance will run the string; @a False, the Utility instance will not run the string.
!> @param input_string     String containing PHREEQC input.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> RunString runs a PHREEQC input string. The first three arguments determine which
!> IPhreeqc instances will run
!> the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
!> strings that modify the thermodynamic database should be run by all three sets of instances.
!> Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
!> be run by the workers. Strings that contain initial conditions or boundary conditions should
!> be run by the InitialPhreeqc instance.
!> @see                    @ref YAMLRunFile.
!> @par C++ Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> input = "DELETE; -all"
!> status = YAMLRunString(id, .true., .true., .true., "advect.pqi")
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLRunString(id, workers, initial_phreeqc, utility, input_string)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLRunString_F(id, iworkers, iinitial_phreeqc, iutility, input_string) &
			BIND(C, NAME='YAMLRunString_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: iworkers, iinitial_phreeqc, iutility
        character(KIND=C_CHAR), intent(in) :: input_string(*)
		END FUNCTION YAMLRunString_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: workers, initial_phreeqc, utility
    integer :: iworkers, iinitial_phreeqc, iutility
    character(len=*), intent(in) :: input_string
	iworkers = 0
	iinitial_phreeqc = 0
	iutility = 0
	if (workers) iworkers = 1
	if (initial_phreeqc) iinitial_phreeqc = 1
	if (utility) iutility = 1
	YAMLRunString = YAMLRunString_F(id, iworkers, iinitial_phreeqc, iutility, trim(input_string)//C_NULL_CHAR)
    END FUNCTION YAMLRunString 
!> Inserts data into the YAML document for the PhreeqcRM method ScreenMessage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param str              String to be printed.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> ScreenMessage prints a message to the screen.
!> @see                    @ref YAMLLogMessage, @ref YAMLOutputMessage, @ref YAMLWarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> string = "Beginning to process YAML for initial conditions"
!> status = YAMLScreenMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLScreenMessage(id, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLScreenMessage_F(id, str) &
			BIND(C, NAME='YAMLScreenMessage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: str(*)
		END FUNCTION YAMLScreenMessage_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: str
	YAMLScreenMessage = YAMLScreenMessage_F(id, trim(str)//C_NULL_CHAR)
    END FUNCTION YAMLScreenMessage
!> Inserts data into the YAML document for the PhreeqcRM method SetComponentH2O.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf               @a True (default), excess H, excess O, and water are included in the component list;
!> @a False, total H and O are included in the component list.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetComponentH2O selects whether to include H2O in the component list.
!> The concentrations of H and O must be known
!> accurately (8 to 10 significant digits) for the numerical method of
!> PHREEQC to produce accurate pH and pe values.
!> Because most of the H and O are in the water species,
!> it may be more robust (require less accuracy in transport) to
!> transport the excess H and O (the H and O not in water) and water.
!> The default setting (@a true) is to include water, excess H, and excess O as components.
!> A setting of @a false will include total H and total O as components.
!> YAMLSetComponentH2O must be called before @ref YAMLFindComponents.
!> @see                    @ref YAMLFindComponents.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetComponentH2O(id, .false.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetComponentH2O(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetComponentH2O_F(id, itf) &
			BIND(C, NAME='YAMLSetComponentH2O_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetComponentH2O_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
    integer             :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetComponentH2O = YAMLSetComponentH2O_F(id, itf)
    END FUNCTION YAMLSetComponentH2O
!> Inserts data into the YAML document for the PhreeqcRM method SetConcentrations.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param c               Vector of component concentrations. Size of vector is @a ncomps times @a nxyz,
!> where @a ncomps is the number of components as determined
!> by FindComponents or GetComponentCount and
!> @a nxyz is the number of grid cells in the user's model.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> The only way to use this method is to have pre-calculated PHREEQC solution concentrations,
!> which is not common. Concentrations are normally initialized
!> with @ref YAMLInitialPhreeqc2Module or @ref YAMLInitialPhreeqcCell2Module.
!> @see                    @ref YAMLSetDensityUser, @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume,
!> @ref YAMLSetSaturationUser, @ref YAMLSetUnitsSolution.
!> @par C++ Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetConcentrations(id, c)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetConcentrations(id, c)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetConcentrations_F(id, c, dim) &
			BIND(C, NAME='YAMLSetConcentrations_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: c
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetConcentrations_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:,:), intent(in) :: c
    integer :: dim
    dim = size(c,1)*size(c,2)
	YAMLSetConcentrations = YAMLSetConcentrations_F(id, c(1,1), dim)
    END FUNCTION YAMLSetConcentrations   
!> Inserts data into the YAML document for the PhreeqcRM method SetCurrentSelectedOutputUserNumber.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetCurrentSelectedOutputUserNumber selects the current selected output by user number.
!> The user may define multiple SELECTED_OUTPUT
!> data blocks for the workers. A user number is specified for each data block. The value of
!> the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
!> for selected-output operations.
!> @see
!> @ref YAMLSetNthSelectedOutput,
!> @ref YAMLSetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetCurrentSelectedOutputUserNumber(id, n_user)
!> </PRE>
!> </CODE>
!> @endhtmlonly    
    INTEGER FUNCTION YAMLSetCurrentSelectedOutputUserNumber(id, n_user)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetCurrentSelectedOutputUserNumber_F(id, n) &
			BIND(C, NAME='YAMLSetCurrentSelectedOutputUserNumber_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetCurrentSelectedOutputUserNumber_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n_user
	YAMLSetCurrentSelectedOutputUserNumber = YAMLSetCurrentSelectedOutputUserNumber_F(id, n_user)
    END FUNCTION YAMLSetCurrentSelectedOutputUserNumber    
!> Inserts data into the YAML document for the PhreeqcRM method SetDensityUser.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> 
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param density          Vector of densities. Size of vector is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetDensityUser sets the density for each reaction cell. These density values are used
!> when converting from transported mass-fraction concentrations (@ref YAMLSetUnitsSolution) to
!> produce per liter concentrations during a call to SetConcentrations.
!> They are also used when converting from reaction-cell concentrations to transport concentrations,
!> if UseSolutionDensityVolume is set to @a false.
!> 
!> @see
!> @ref YAMLSetUnitsSolution, @ref YAMLUseSolutionDensityVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:)   :: density
!> allocate(density(nxyz))
!> por = 1.0d0
!> status = YAMLSetDensityUser(id, por)
!> </PRE>
!> </CODE>
!> @endhtmlonly    
    INTEGER FUNCTION YAMLSetDensityUser(id, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetDensityUser_F(id, density, dim) &
			BIND(C, NAME='YAMLSetDensityUser_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: density
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetDensityUser_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: density
	YAMLSetDensityUser = YAMLSetDensityUser_F(id, density(1), size(density))
    END FUNCTION YAMLSetDensityUser
!> Inserts data into the YAML document for the PhreeqcRM method SetDumpFileName.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param file_name        Name of dump file.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetDumpFileName	sets the name of the dump file. It is the name used by the method DumpModule.
!> @see                    @ref YAMLDumpModule.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> logical dump_on, append
!> status = YAMLSetDumpFileName(id, "Advect_cpp.dmp")
!> dump_on = .true.
!> append = .false.
!> status = YAMLDumpModule(id, dump_on, append)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetDumpFileName(id, file_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetDumpFileName_F(id, file_name) &
			BIND(C, NAME='YAMLSetDumpFileName_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: file_name(*)
		END FUNCTION YAMLSetDumpFileName_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: file_name
	YAMLSetDumpFileName = YAMLSetDumpFileName_F(id, trim(file_name)//C_NULL_CHAR)
    END FUNCTION YAMLSetDumpFileName
!> Inserts data into the YAML document for the PhreeqcRM method SetErrorHandlerMode.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param mode             Error handling mode: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetErrorHandlerMode sets the action to be taken when the reaction module encounters an error.
!> Options are 0, return to calling program with an error return code (default);
!> 1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
!> 2, attempt to exit gracefully.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetErrorHandlerMode(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetErrorHandlerMode(id, mode)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetErrorHandlerMode_F(id, n) &
			BIND(C, NAME='YAMLSetErrorHandlerMode_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetErrorHandlerMode_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: mode
	YAMLSetErrorHandlerMode = YAMLSetErrorHandlerMode_F(id, mode)
    END FUNCTION YAMLSetErrorHandlerMode  
!> Inserts data into the YAML document for the PhreeqcRM method SetErrorOn.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf  @a True, enable error messages; @a False, disable error messages. Default is true.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetErrorOn sets the property that controls whether error messages are generated and displayed.
!> Messages include PHREEQC "ERROR" messages, and
!> any messages written with the method ErrorMessage.
!> 
!> @see                  @ref YAMLLogMessage, @ref YAMLOutputMessage, @ref YAMLScreenMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetErrorOn(id, .true.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetErrorOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetErrorOn_F(id, itf) &
			BIND(C, NAME='YAMLSetErrorOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetErrorOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
    integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetErrorOn = YAMLSetErrorOn_F(id, itf)
    END FUNCTION YAMLSetErrorOn   
!> Inserts data into the YAML document for the PhreeqcRM method SetFilePrefix.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param prefix           Prefix used when opening the output and log files.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetFilePrefix sets the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
!> These files are opened by the method OpenFiles.
!> @see                    @ref YAMLOpenFiles, @ref YAMLCloseFiles.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetFilePrefix(id, "Advect_cpp")
!> status = YAMLOpenFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetFilePrefix(id, prefix)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetFilePrefix_F(id, prefix) &
			BIND(C, NAME='YAMLSetFilePrefix_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: prefix(*)
		END FUNCTION YAMLSetFilePrefix_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: prefix
	YAMLSetFilePrefix = YAMLSetFilePrefix_F(id, trim(prefix)//C_NULL_CHAR)
    END FUNCTION YAMLSetFilePrefix
!> Inserts data into the YAML document for the PhreeqcRM method SetGasCompMoles.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param  gas_moles               Vector of moles of gas components.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetGasCompMoles transfers moles of gas components from
!> the vector given in the argument list (@a gas_moles) to each reaction cell.
!> 
!> Dimension of the vector is set to @a ngas_comps times @a nxyz,
!> where, @a ngas_comps is the result of GetGasComponentsCount,
!> and @a nxyz is the number of user grid cells.
!> If the number of moles is set to a negative number, the gas component will
!> not be defined for the GAS_PHASE of the reaction cell.
!> 
!> @see
!> @ref YAMLFindComponents,
!> @ref YAMLSetGasPhaseVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:,:) :: gas_moles
!> allocate(gas_moles(nxyz*ngas))
!> status = YAMLSetGasCompMoles(id, gas_moles)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetGasCompMoles(id, gas_moles)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetGasCompMoles_F(id, gas_moles, dim) &
			BIND(C, NAME='YAMLSetGasCompMoles_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: gas_moles
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetGasCompMoles_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:,:), intent(in) :: gas_moles
    integer :: dim
    dim = size(gas_moles,1)*size(gas_moles,2)
	YAMLSetGasCompMoles = YAMLSetGasCompMoles_F(id, gas_moles(1,1), dim)
    END FUNCTION YAMLSetGasCompMoles   
!> Inserts data into the YAML document for the PhreeqcRM method SetGasPhaseVolume.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param  gas_volume               Vector of volumes for each gas phase.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetGasPhaseVolume transfers volumes of gas phases from
!> the vector given in the argument list (@a gas_volume) to each reaction cell.
!> The gas-phase volume affects the gas-component pressures calculated for fixed-volume
!> gas phases. If a gas-phase volume is defined with this methood
!> for a GAS_PHASE in a cell,
!> the gas phase is forced to be a fixed-volume gas phase.
!> 
!> Dimension of the vector is @a nxyz,
!> where @a nxyz is the number of user grid cells.
!> If the volume is set to a negative number for a cell, the gas-phase volume for that cell is
!> not changed.
!> 
!> @see
!> @ref YAMLFindComponents,
!> @ref YAMLSetGasCompMoles.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:) :: gas_volume
!> allocate(gas_volume(nxyz))
!> status = YAMLSetGasPhaseVolume(id, gas_volume)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetGasPhaseVolume(id, gas_volume)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetGasPhaseVolume_F(id, gas_volume, dim) &
			BIND(C, NAME='YAMLSetGasPhaseVolume_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: gas_volume
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetGasPhaseVolume_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: gas_volume
	YAMLSetGasPhaseVolume = YAMLSetGasPhaseVolume_F(id, gas_volume(1), size(gas_volume))
    END FUNCTION YAMLSetGasPhaseVolume
!> Inserts data into the YAML document to define the number of cells in the user's model.
!> Once the YAML document is written, the number of model cells can be extracted
!> with the method GetGridCellCountYAML. GetGridCellCountYAML is NOT a PhreeqcRM 
!> method; it is a global method and must be used BEFORE the PhreeqcRM instance
!> is created. SetGridCellCount will be ignored once the PhreeqcRM instance exists.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n           Number of cells for the PhreeqcRM instance. The number of cells
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> can be used in the creation of the PhreeqcRM instance. The PhreeqcRM constructor
!> takes two arguments. GetGridCellCountYAML
!> provides the value for the first argument. If the YAML file does not contain the
!> node "SetGridCellCount:", GetGridCellCountYAML will return zero.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> integer nxyz
!> nxyz = 40
!> status = YAMLSetGridCellCount(nxyz)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetGridCellCount(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetGridCellCount_F(id, n) &
			BIND(C, NAME='YAMLSetGridCellCount_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetGridCellCount_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetGridCellCount = YAMLSetGridCellCount_F(id, n)
    END FUNCTION YAMLSetGridCellCount  
!> Inserts data into the YAML document for the PhreeqcRM method SetNthSelectedOutput.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n           Sequence number of the SELECTED_OUTPUT data block that is to be used.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetNthSelectedOutput specifies the current selected output by sequence number (one-based). 
!> The user may define multiple SELECTED_OUTPUT
!> data blocks for the workers. A user number is specified for each data block, and the blocks are
!> stored in user-number order. The value of
!> the argument @a n selects the sequence number of the SELECTED_OUTPUT definition that will be used
!> for selected-output operations. 
!> @see
!> @ref YAMLSetCurrentSelectedOutputUserNumber,
!> @ref YAMLSetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetCurrentSelectedOutput(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetNthSelectedOutput(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetNthSelectedOutput_F(id, n) &
			BIND(C, NAME='YAMLSetNthSelectedOutput_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetNthSelectedOutput_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetNthSelectedOutput = YAMLSetNthSelectedOutput_F(id, n)
    END FUNCTION YAMLSetNthSelectedOutput  
!> Inserts data into the YAML document for the PhreeqcRM method SetPartitionUZSolids.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf       @a True, the fraction of solids and gases available for
!> reaction is equal to the saturation;
!> @a False (default), all solids and gases are reactive regardless of saturation.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetPartitionUZSolids sets the property for partitioning solids between the 
!> saturated and unsaturated parts of a partially saturated cell.
!> @par
!> The option is intended to be used by saturated-only
!> flow codes that allow a variable water table.
!> The value has meaning only when saturations
!> less than 1.0 are encountered. The partially saturated cells
!> may have a small water-to-rock ratio that causes
!> reactions to proceed differently relative to fully saturated cells.
!> By setting  @a SetPartitionUZSolids to true, the
!> amounts of solids and gases are partioned according to the saturation.
!> If a cell has a saturation of 0.5, then
!> the water interacts with only half of the solids and gases; the other half is unreactive
!> until the water table rises. As the saturation in a cell varies,
!> solids and gases are transferred between the
!> saturated and unsaturated (unreactive) reservoirs of the cell.
!> Unsaturated-zone flow and transport codes will probably use the default (false),
!> which assumes all gases and solids are reactive regardless of saturation.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetPartitionUZSolids(id, .false.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetPartitionUZSolids(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPartitionUZSolids_F(id, itf) &
			BIND(C, NAME='YAMLSetPartitionUZSolids_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetPartitionUZSolids_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
	integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetPartitionUZSolids = YAMLSetPartitionUZSolids_F(id, itf)
    END FUNCTION YAMLSetPartitionUZSolids       
!> Inserts data into the YAML document for the PhreeqcRM method SetPorosity.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param por              Vector of porosities, unitless. Default is 0.1.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetPorosity sets the porosity for each reaction cell.
!> The volume of water in a reaction cell is the product of porosity, saturation
!> (SetSaturationUser), and representative volume (SetRepresentativeVolume).
!> Size of vector is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model.
!> @see                    @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:) :: por
!> allocate(por(nxyz))
!> por = 0.2d0
!> status = YAMLSetPorosity(id, por)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetPorosity(id, por)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPorosity_F(id, por, dim) &
			BIND(C, NAME='YAMLSetPorosity_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: por
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetPorosity_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: por
	YAMLSetPorosity = YAMLSetPorosity_F(id, por(1), size(por))
    END FUNCTION YAMLSetPorosity
!> Inserts data into the YAML document for the PhreeqcRM method SetPressure.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param p                Vector of pressures, in atm. Size of vector is @a nxyz,
!> where @a nxyz is the number of grid cells in the user's model.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetPressure sets the pressure for each reaction cell. Pressure effects are 
!> considered only in three of the
!> databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
!> @see                    @ref YAMLSetTemperature.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:) :: pressure
!> allocate(por(nxyz))
!> pressure = 2.0d0
!> status = YAMLSetPressure(id, pressure)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetPressure(id, p)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPressure_F(id, p, dim) &
			BIND(C, NAME='YAMLSetPressure_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: p
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetPressure_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: p
	YAMLSetPressure = YAMLSetPressure_F(id, p(1), size(p))
    END FUNCTION YAMLSetPressure
!> Inserts data into the YAML document for the PhreeqcRM method SetPrintChemistryMask.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param mask        Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model. A value of 0 will
!> disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetPrintChemistryMask enables or disables detailed output for each reaction cell.
!> Printing for a reaction cell will occur only when the
!> printing is enabled with SetPrintChemistryOn and the @a mask value is 1.
!> @see                    @ref YAMLSetPrintChemistryOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> integer, allocatable, dimension(:) :: print_chemistry_mask;
!> allocate(print_chemistry_mask(nxyz))
!> print_chemistry_mask = 0
!> do i = 1, nxyz / 2 
!> 	print_chemistry_mask(i) = 1
!> enddo
!> status = YAMLSetPrintChemistryMask(id, print_chemistry_mask)  
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetPrintChemistryMask(id, mask)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPrintChemistryMask_F(id, mask, dim) &
			BIND(C, NAME='YAMLSetPrintChemistryMask_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: mask
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetPrintChemistryMask_F
    END INTERFACE
    integer, intent(in) :: id
    integer, allocatable, dimension(:), intent(in) :: mask
	YAMLSetPrintChemistryMask = YAMLSetPrintChemistryMask_F(id, mask(1), size(mask))
    END FUNCTION YAMLSetPrintChemistryMask  
!> Inserts data into the YAML document for the PhreeqcRM method SetPrintChemistryOn.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param workers          @a True, enable detailed printing in the worker instances;
!> @a False, disable detailed printing in the worker instances.
!> @param initial_phreeqc  @a True, enable detailed printing in the InitialPhreeqc instance;
!> @a False, disable detailed printing in the InitialPhreeqc instance.
!> @param utility          @a True, enable detailed printing in the Utility instance;
!> @a False, disable detailed printing in the Utility instance.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetPrintChemistryOn
!> sets the property that enables or disables printing detailed output from reaction calculations
!> to the output file for a set of cells defined by SetPrintChemistryMask.
!> The detailed output prints all of the output typical of a PHREEQC reaction calculation,
!> which includes solution descriptions and the compositions of all other reactants.
!> The output can be several hundred lines per cell, which can lead to a very
!> large output file (prefix.chem.txt opened by the method OpenFiles).
!> For the worker instances, the output can be limited to a set of cells
!> (method SetPrintChemistryMask) and, in general, the
!> amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
!> (applied by using methods RunFile or RunString).
!> Printing the detailed output for the workers is generally used only for debugging,
!> and PhreeqcRM will run significantly faster
!> when printing detailed output for the workers is disabled.
!> @see                    @ref YAMLOpenFiles, @ref YAMLRunFile, @ref YAMLRunString, 
!> @ref YAMLSetPrintChemistryMask.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetPrintChemistryOn(id, .false., .true., .false.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetPrintChemistryOn(id, workers, initial_phreeqc, utility)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPrintChemistryOn_F(id, iworkers, iinitial_phreeqc, iutility) &
			BIND(C, NAME='YAMLSetPrintChemistryOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: iworkers, iinitial_phreeqc, iutility
		END FUNCTION YAMLSetPrintChemistryOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: workers, initial_phreeqc, utility
	integer :: iworkers, iinitial_phreeqc, iutility
	iworkers = 0
	iinitial_phreeqc = 0
	iutility = 0
	if (workers) iworkers = 1
	if (initial_phreeqc) iinitial_phreeqc = 1
	if (utility) iutility = 1
	YAMLSetPrintChemistryOn = YAMLSetPrintChemistryOn_F(id, iworkers, iinitial_phreeqc, iutility)
    END FUNCTION YAMLSetPrintChemistryOn
!> Inserts data into the YAML document for the PhreeqcRM method SetRebalanceByCell.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf           @a True, indicates individual cell times are used in rebalancing (default);
!> @a False, indicates average times are used in rebalancing.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetRebalanceByCell
!> sets the load-balancing algorithm.
!> PhreeqcRM attempts to rebalance the load of each thread or process such that each
!> thread or process takes the same amount of time to run its part of a RunCells
!> calculation. Two algorithms are available; one uses individual times for each cell and
!> accounts for cells that were not run because
!> saturation was zero (default), and
!> the other assigns an average time to all cells.
!> The methods are similar, but limited testing indicates the default method performs better.
!> @see                    @ref YAMLSetRebalanceFraction.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetRebalanceByCell(id, .true.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetRebalanceByCell(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetRebalanceByCell_F(id, itf) &
			BIND(C, NAME='YAMLSetRebalanceByCell_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetRebalanceByCell_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
	integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetRebalanceByCell = YAMLSetRebalanceByCell_F(id, itf)
    END FUNCTION YAMLSetRebalanceByCell  
!> Inserts data into the YAML document for the PhreeqcRM method SetRebalanceFraction.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param f                Fraction from 0.0 to 1.0.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetRebalanceFraction sets the fraction of cells that are transferred 
!> among threads or processes when rebalancing.
!> PhreeqcRM attempts to rebalance the load of each thread or process such that each
!> thread or process takes the same amount of time to run its part of a RunCells
!> calculation. The rebalancing transfers cell calculations among threads or processes to
!> try to achieve an optimum balance. @a SetRebalanceFraction
!> adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
!> determine the actual number of cell transfers. A value of zero eliminates
!> load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
!> distribution and avoid possible oscillations
!> when too many cells are transferred at one iteration, requiring reverse transfers 
!> at the next iteration. Default is 0.5.
!> 
!> @see                    @ref YAMLSetRebalanceByCell.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetRebalanceFraction(id, 0.5)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetRebalanceFraction(id, f)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetRebalanceFraction_F(id, f) &
			BIND(C, NAME='YAMLSetRebalanceFraction_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: f
		END FUNCTION YAMLSetRebalanceFraction_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), intent(in) :: f
	YAMLSetRebalanceFraction = YAMLSetRebalanceFraction_F(id, f)
    END FUNCTION YAMLSetRebalanceFraction  
!> Inserts data into the YAML document for the PhreeqcRM method SetRepresentativeVolume.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param rv              Vector of representative volumes, in liters. Default is 1.0 liter.
!> Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetRepresentativeVolume
!> sets the representative volume of each reaction cell.
!> By default the representative volume of each reaction cell is 1 liter.
!> The volume of water in a reaction cell is determined by the product of the representative volume,
!> the porosity (SetPorosity), and the saturation (SetSaturationUser).
!> The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
!> within a couple orders of magnitude of 1.0.
!> Small water volumes caused by small porosities and (or) small saturations 
!> (and (or) small representative volumes)
!> may cause non-convergence of the numerical method.
!> In these cases, a larger representative volume may help. Note
!> that increasing the representative volume also increases
!> the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
!> and others), which are defined as moles per representative volume.
!> @a SetRepresentativeVolume should be called before initial conditions 
!> are defined for the reaction cells.
!> 
!> @see                    @ref YAMLSetPorosity, @ref YAMLSetSaturationUser.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> real(kind=8), allocatable, dimension(:)   :: rv
!> allocate(rv(nxyz))
!> rv = 1.0d0
!> status = YAMLSetRepresentativeVolume(id, rv)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetRepresentativeVolume(id, rv)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetRepresentativeVolume_F(id, rv, dim) &
			BIND(C, NAME='YAMLSetRepresentativeVolume_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: rv
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetRepresentativeVolume_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: rv
	YAMLSetRepresentativeVolume = YAMLSetRepresentativeVolume_F(id, rv(1), size(rv))
    END FUNCTION YAMLSetRepresentativeVolume
!> Inserts data into the YAML document for the PhreeqcRM method SetSaturationUser.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param sat              Vector of saturations, unitless. Default 1.0. Size of vector is @a nxyz,
!> where @a nxyz is the number of grid cells in the user's model.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetSaturationUser
!> sets the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
!> The volume of water in a cell is the product of porosity (SetPorosity), saturation (SetSaturationUser),
!> and representative volume (SetRepresentativeVolume). As a result of a reaction calculation,
!> solution properties (density and volume) will change;
!> the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar 
!> volume data to calculate these changes. The methods GetDensityCalculated,
!> GetSolutionVolume, and GetSaturationCalculated can be used to account
!> for these changes in the succeeding transport calculation.
!> 
!> @see
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>   
!> real(kind=8), allocatable, dimension(:) :: sat 
!> allocate(sat(nxyz))
!> sat = 1.0d0
!> status = YAMLSetSaturationUser(id, sat)  
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetSaturationUser(id, sat)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSaturationUser_F(id, sat, dim) &
			BIND(C, NAME='YAMLSetSaturationUser_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: sat
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetSaturationUser_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: sat
	YAMLSetSaturationUser = YAMLSetSaturationUser_F(id, sat(1), size(sat))
    END FUNCTION YAMLSetSaturationUser
!> Inserts data into the YAML document for the PhreeqcRM method SetScreenOn.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf  @a True, enable screen messages; @a False, disable screen messages. Default is true.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetScreenOn
!> sets the property that controls whether messages are written to the screen.
!> Messages include information about rebalancing during RunCells, and
!> any messages written with ScreenMessage.
!> 
!> @see                    @ref YAMLRunCells, @ref YAMLScreenMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetScreenOn(id, .true.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetScreenOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetScreenOn_F(id, itf) &
			BIND(C, NAME='YAMLSetScreenOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetScreenOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
	integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetScreenOn = YAMLSetScreenOn_F(id, itf)
    END FUNCTION YAMLSetScreenOn  
!> Inserts data into the YAML document for the PhreeqcRM method SetSelectedOutputOn.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf  @a True, enable selected output; @a False, disable selected output.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetSelectedOutputOn
!> sets the property that controls whether selected-output results are available to be retrieved
!> with GetSelectedOutput. @a True indicates that selected-output results
!> will be accumulated during RunCells and can be retrieved with GetSelectedOutput;
!> @a False indicates that selected-output results will not
!> be accumulated during RunCells.
!> 
!> @see
!> @ref YAMLSetCurrentSelectedOutputUserNumber,
!> @ref YAMLSetNthSelectedOutput,
!> @ref YAMLSetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetSelectedOutputOn(id, .true.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetSelectedOutputOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSelectedOutputOn_F(id, itf) &
			BIND(C, NAME='YAMLSetSelectedOutputOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLSetSelectedOutputOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
	integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLSetSelectedOutputOn = YAMLSetSelectedOutputOn_F(id, itf)
    END FUNCTION YAMLSetSelectedOutputOn  
!> Inserts data into the YAML document for the PhreeqcRM method SetSpeciesSaveOn.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param save_on          @a True indicates species concentrations are saved;
!> @a False indicates species concentrations are not saved.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetSpeciesSaveOn
!> sets the value of the species-save property.
!> This method enables or disables use of PhreeqcRM with multicomponent-diffusion transport calculations.
!> By default, concentrations of aqueous species are not saved.
!> Setting the species-save property to @a true allows
!> aqueous species concentrations to be retrieved
!> with GetSpeciesConcentrations, and solution compositions to be set with
!> SpeciesConcentrations2Module.
!> SetSpeciesSaveOn must be called before calls to FindComponents.
!> 
!> @see                    @ref YAMLFindComponents,
!> @ref YAMLSpeciesConcentrations2Module.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetSpeciesSaveOn(id, .true.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetSpeciesSaveOn(id, save_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSpeciesSaveOn_F(id, save_on) &
			BIND(C, NAME='YAMLSetSpeciesSaveOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: save_on
		END FUNCTION YAMLSetSpeciesSaveOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: save_on
	integer :: itf
	itf = 0
	if (save_on) itf = 1
	YAMLSetSpeciesSaveOn = YAMLSetSpeciesSaveOn_F(id, itf)
    END FUNCTION YAMLSetSpeciesSaveOn 
!> Inserts data into the YAML document for the PhreeqcRM method SetTemperature.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tc                Vector of temperatures, in degrees C.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetTemperature
!> sets the temperature for each reaction cell. If SetTemperature is not called,
!> worker solutions will have temperatures as defined by initial conditions
!> (InitialPhreeqc2Module and InitialPhreeqcCell2Module).
!> Size of vector is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model.
!> @see                    @ref YAMLInitialPhreeqc2Module,
!> @ref YAMLInitialPhreeqcCell2Module, @ref YAMLSetPressure.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>    
!> real(kind=8), allocatable, dimension(:) :: tc
!> allocate(tc(nxyz))
!> por = 0.2d0
!> status = YAMLSetTemperature(id, tc)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetTemperature(id, tc)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTemperature_F(id, tc, dim) &
			BIND(C, NAME='YAMLSetTemperature_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: tc
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetTemperature_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:), intent(in) :: tc
	YAMLSetTemperature = YAMLSetTemperature_F(id, tc(1), size(tc))
    END FUNCTION YAMLSetTemperature
!> Inserts data into the YAML document for the PhreeqcRM method SetTime.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param time             Current simulation time, in seconds.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetTime
!> sets current simulation time for the reaction module.
!> @see                    @ref YAMLSetTimeStep, @ref YAMLSetTimeConversion.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetTime(id, 0.0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetTime(id, time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTime_F(id, time) &
			BIND(C, NAME='YAMLSetTime_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: time
		END FUNCTION YAMLSetTime_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), intent(in) :: time
	YAMLSetTime = YAMLSetTime_F(id, time)
    END FUNCTION YAMLSetTime  
!> Inserts data into the YAML document for the PhreeqcRM method SetTimeConversion.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param conv_factor      Factor to convert seconds to user time units.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetTimeConversion
!> Set a factor to convert from seconds to user time units. Factor times seconds produces user time units
!> that is used in some PhreeqcRM printing.
!> 
!> @see                    @ref YAMLSetTime, @ref YAMLSetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> conv_factor = 1.0d0 / 86400.d0
!> status = YAMLSetTimeConversion(id, conv_factor)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetTimeConversion(id, conv_factor)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTimeConversion_F(id, conv_factor) &
			BIND(C, NAME='YAMLSetTimeConversion_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: conv_factor
		END FUNCTION YAMLSetTimeConversion_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), intent(in) :: conv_factor
	YAMLSetTimeConversion = YAMLSetTimeConversion_F(id, conv_factor)
    END FUNCTION YAMLSetTimeConversion 
!> Inserts data into the YAML document for the PhreeqcRM method SetTimeStep.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param time_step        Time step, in seconds.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetTimeStep
!> sets current time step for the reaction module. This is the length
!> of time over which kinetic reactions are integrated.
!> 
!> @see                    @ref YAMLSetTime, @ref YAMLSetTimeConversion.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> time_step = 86400.d0
!> status = YAMLSetTimeStep(time_step)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetTimeStep(id, time_step)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTimeStep_F(id, time_step) &
			BIND(C, NAME='YAMLSetTimeStep_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: time_step
		END FUNCTION YAMLSetTimeStep_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), intent(in) :: time_step
	YAMLSetTimeStep = YAMLSetTimeStep_F(id, time_step)
    END FUNCTION YAMLSetTimeStep  
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsExchange.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option           Units option for exchangers: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsExchange
!> sets input units for exchangers.
!> In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
!> SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
!> @par
!> If a single EXCHANGE definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of exchangers will be the same regardless of porosity.
!> For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsExchange(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsExchange(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsExchange_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsExchange_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsExchange_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsExchange = YAMLSetUnitsExchange_F(id, option)
    END FUNCTION YAMLSetUnitsExchange 
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsGasPhase.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option           Units option for gas phases: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsGasPhase
!> sets input units for gas phases.
!> In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
!> @a SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> @par
!> If a single GAS_PHASE definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of a gas component will be the same regardless of porosity.
!> For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsGasPhase(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsGasPhase(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsGasPhase_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsGasPhase_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsGasPhase_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsGasPhase = YAMLSetUnitsGasPhase_F(id, option)
    END FUNCTION YAMLSetUnitsGasPhase 
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsKinetics.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option           Units option for kinetic reactants: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsKinetics
!> sets input units for kinetic reactants.
!> @par
!> In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
!> @a SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> @par
!> If a single KINETICS definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of kinetic reactants will be the same regardless of porosity.
!> For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.
!> @par
!> Note that the volume of water in a cell in the reaction module is equal to the product of
!> porosity (SetPorosity), the saturation (SetSaturationUser), and representative volume (SetRepresentativeVolume),
!> which is usually less than 1 liter. It is important to write the RATES
!> definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
!> water, often by calculating the rate of reaction per liter of water and multiplying by the volume
!> of water (Basic function SOLN_VOL).
!> @par
!> Rates that depend on surface area of solids, are not dependent
!> on the volume of water. However, it is important to get the correct surface area for the kinetic
!> reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant)
!> can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of
!> reactant (Basic function M) in RATES to obtain the surface area.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsKinetics(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsKinetics(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsKinetics_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsKinetics_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsKinetics_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsKinetics = YAMLSetUnitsKinetics_F(id, option)
    END FUNCTION YAMLSetUnitsKinetics
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsPPassemblage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option           Units option for equilibrium phases: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsPPassemblage
!> sets input units for pure phase assemblages (equilibrium phases).
!> In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
!> @a SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
!> @par
!> If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of a mineral will be the same regardless of porosity.
!> For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsPPassemblage(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsPPassemblage(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsPPassemblage_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsPPassemblage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsPPassemblage_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsPPassemblage = YAMLSetUnitsPPassemblage_F(id, option)
    END FUNCTION YAMLSetUnitsPPassemblage
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsSolution.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsSolution
!> sets solution concentration units used by the transport model.
!> Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
!> PHREEQC defines solutions by the number of moles of each
!> element in the solution.
!> @par
!> To convert from mg/L to moles
!> of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
!> multiplied by the solution volume,
!> which is the product of porosity (SetPorosity), saturation (SetSaturationUser),
!> and representative volume (SetRepresentativeVolume).
!> To convert from mol/L to moles
!> of element in the representative volume of a reaction cell, mol/L is
!> multiplied by the solution volume.
!> To convert from mass fraction to moles
!> of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
!> (SetDensityUser) and
!> multiplied by the solution volume.
!> @par
!> To convert from moles
!> of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
!> solution volume resulting in mol/L, and then converted to mg/L.
!> To convert from moles
!> of element in a cell to mol/L,  the number of moles of an element is divided by the
!> solution volume resulting in mol/L.
!> To convert from moles
!> of element in a cell to mass fraction, the number of moles of an element is converted to kg and divided
!> by the total mass of the solution.
!> Two options are available for the volume and mass of solution
!> that are used in converting to transport concentrations: (1) the volume and mass of solution are
!> calculated by PHREEQC, or (2) the volume of solution is the product of porosity (SetPorosity),
!> saturation (SetSaturationUser), and representative volume (SetRepresentativeVolume),
!> and the mass of solution is volume times density as defined by SetDensityUser.
!> Which option is used is determined by UseSolutionDensityVolume.
!> 
!> @see                    @ref YAMLSetDensityUser, @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume,
!> @ref YAMLSetSaturationUser,
!> @ref YAMLUseSolutionDensityVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsSolution(id, 2)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsSolution(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSolution_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsSolution_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsSolution_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsSolution = YAMLSetUnitsSolution_F(id, option)
    END FUNCTION YAMLSetUnitsSolution
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsSSassemblage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option        Units option for solid solutions: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsSSassemblage
!> sets input units for solid-solution assemblages.
!> In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
!> @a SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@ P)*RV.
!> @par
!> If a single SOLID_SOLUTION definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of a solid-solution component will be the same regardless of porosity.
!> For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsSSassemblage(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsSSassemblage(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSSassemblage_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsSSassemblage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsSSassemblage_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsSSassemblage = YAMLSetUnitsSSassemblage_F(id, option)
    END FUNCTION YAMLSetUnitsSSassemblage
!> Inserts data into the YAML document for the PhreeqcRM method SetUnitsSurface.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param option        Units option for surfaces: 0, 1, or 2.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SetUnitsSurface
!> sets input units for surfaces.
!> In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
!> @a SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> @par
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> @par
!> If a single SURFACE definition is used for cells with different initial porosity,
!>    the three options scale quite differently.
!> For option 0, the number of moles of surface sites will be the same regardless of porosity.
!> For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume.
!> For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.
!> 
!> @see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSetUnitsSurface(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSetUnitsSurface(id, option)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSurface_F(id, option) &
			BIND(C, NAME='YAMLSetUnitsSurface_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: option
		END FUNCTION YAMLSetUnitsSurface_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: option
	YAMLSetUnitsSurface = YAMLSetUnitsSurface_F(id, option)
    END FUNCTION YAMLSetUnitsSurface
!> Inserts data into the YAML document for the PhreeqcRM method SpeciesConcentrations2Module.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param species_conc     Vector of aqueous species concentrations. Dimension of the 
!> array is @a nspecies times @a nxyz,
!> where  @a nspecies is the number of aqueous species,
!> and @a nxyz is the number of user grid cells.
!> Concentrations are moles per liter.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> SpeciesConcentrations2Module
!> sets solution concentrations in the reaction cells
!> based on the vector of aqueous species concentrations (@a species_conc).
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and SetSpeciesSaveOn must be set to @a true.
!> The list of aqueous species is determined by FindComponents and includes all
!> aqueous species that can be made from the set of components.
!> The method determines the total concentration of a component
!> by summing the molarities of the individual species times the stoichiometric
!> coefficient of the element in each species.
!> Solution compositions in the reaction cells are updated with these component concentrations.
!> Usually, accurate concentrations will not be known to use YAMLSpeciesConcentrations2Module during
!> initialization.
!> 
!> @see                    @ref YAMLFindComponents,
!> @ref YAMLSetSpeciesSaveOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLSpeciesConcentrations2Module(id, c)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLSpeciesConcentrations2Module(id, species_conc)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSpeciesConcentrations2Module_F(id, c, dim) &
			BIND(C, NAME='YAMLSpeciesConcentrations2Module_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: c
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSpeciesConcentrations2Module_F
    END INTERFACE
    integer, intent(in) :: id
    real(kind=8), allocatable, dimension(:,:), intent(in) :: species_conc
    integer :: dim
    dim = size(species_conc,1)*size(species_conc,2)
	YAMLSpeciesConcentrations2Module = YAMLSpeciesConcentrations2Module_F(id, species_conc(1,1), dim)
    END FUNCTION YAMLSpeciesConcentrations2Module   
!> Inserts data into the YAML document for the PhreeqcRM method StateSave.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n      Integer identifying the state that is saved.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> StateSave
!> saves the state of the chemistry in all model cells, including SOLUTIONs,
!> EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
!> Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
!> will be saved for each cell, if they have been defined in the worker IPhreeqc instances.
!> The distribution of cells among the workers and the chemistry of fully or partially
!> unsaturated cells are also saved. The state is saved in memory; use DumpModule to save the state
!> to file. PhreeqcRM can be reset to this state by using StateApply.
!> A state is identified by an integer, and multiple states can be saved.
!> 
!> @see                    @ref YAMLDumpModule,
!> @ref YAMLStateApply, and
!> @ref YAMLStateDelete.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLStateSave(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLStateSave(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLStateSave_F(id, n) &
			BIND(C, NAME='YAMLStateSave_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLStateSave_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLStateSave = YAMLStateSave_F(id, n)
    END FUNCTION YAMLStateSave
!> Inserts data into the YAML document for the PhreeqcRM method StateApply.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n     Integer identifying the state that is to be applied.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> StateApply
!> resets the state of the module to a state previously saved with StateSave.
!> The chemistry of all model cells are reset, including SOLUTIONs,
!> EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
!> MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
!> will be reset for each cell, if they were defined in the worker IPhreeqc instances
!> at the time the state was saved.
!> The distribution of cells among the workers and the chemistry of fully or partially
!> unsaturated cells are also reset to the saved state.
!> The state to be applied is identified by an integer.
!> 
!> @see                    @ref YAMLStateSave and
!> @ref YAMLStateDelete.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLStateApply(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLStateApply(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLStateApply_F(id, n) &
			BIND(C, NAME='YAMLStateApply_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLStateApply_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLStateApply = YAMLStateApply_F(id, n)
    END FUNCTION YAMLStateApply
!> Inserts data into the YAML document for the PhreeqcRM method StateDelete.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n     Integer identifying the state that is to be deleted.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> StateDelete
!> deletes a state previously saved with StateSave.
!> 
!> @see                    @ref YAMLStateSave and
!> @ref YAMLStateApply.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLStateDelete(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLStateDelete(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLStateDelete_F(id, n) &
			BIND(C, NAME='YAMLStateDelete_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLStateDelete_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLStateDelete = YAMLStateDelete_F(id, n)
    END FUNCTION YAMLStateDelete
    
!> Inserts data into the YAML document to define the number of threads to use
!> with PhreeqcRM calculations..
!> Once the YAML document is written, the number threads to use can be extracted
!> when bmif_initialize is called. The data for ThreadCount will be ignored 
!> if the PhreeqcRM instance has already been initialized.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param n           Number of threads to use for multiprocessing in PhreeqcRM instance. 
!> A value of zero will cause PhreeqcRM to use the number of logical processors available
!> on the computer.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLThreadCount(0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLThreadCount(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLThreadCount_F(id, n) &
			BIND(C, NAME='YAMLThreadCount_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLThreadCount_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLThreadCount = YAMLThreadCount_F(id, n)
    END FUNCTION YAMLThreadCount
    
!> Inserts data into the YAML document for the PhreeqcRM method UseSolutionDensityVolume.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param tf          @a True indicates that the solution density and volume as
!> calculated by PHREEQC will be used to calculate concentrations.
!> @a False indicates that the solution density set by SetDensityUser and the volume determined by the
!> product of  SetSaturationUser, SetPorosity, and SetRepresentativeVolume,
!> will be used to calculate concentrations retrieved by GetConcentrations.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> UseSolutionDensityVolume
!> determines the volume and density to use when converting from the reaction-cell concentrations
!> to transport concentrations (GetConcentrations).
!> Two options are available to convert concentration units:
!> (1) the density and solution volume calculated by PHREEQC are used, or
!> (2) the specified density (SetDensityUser)
!> and solution volume are determined by the product of
!> saturation (SetSaturationUser), porosity (SetPorosity),
!> and representative volume (SetRepresentativeVolume).
!> Transport models that consider density-dependent flow will probably use the
!> PHREEQC-calculated density and solution volume (default),
!> whereas transport models that assume constant-density flow will probably use
!> specified values of density and solution volume.
!> Only the following databases distributed with PhreeqcRM have molar-volume information
!> needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
!> Density is only used when converting to or from transport units of mass fraction.
!> 
!> @see                    @ref YAMLSetDensityUser,
!> @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = YAMLUseSolutionDensityVolume(id, .false.)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLUseSolutionDensityVolume(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLUseSolutionDensityVolume_F(id, itf) &
			BIND(C, NAME='YAMLUseSolutionDensityVolume_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: itf
		END FUNCTION YAMLUseSolutionDensityVolume_F
    END INTERFACE
    integer, intent(in) :: id
    logical, intent(in) :: tf
	integer :: itf
	itf = 0
	if (tf) itf = 1
	YAMLUseSolutionDensityVolume = YAMLUseSolutionDensityVolume_F(id, itf)
    END FUNCTION YAMLUseSolutionDensityVolume 
!> Inserts data into the YAML document for the PhreeqcRM method WarningMessage.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @param str          String to be printed.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
!> @par
!> WarningMessage
!> prints a warning message to the screen and the log file.
!> 
!> @see                    @ref YAMLOpenFiles, @ref YAMLLogMessage,
!> @ref YAMLOutputMessage, @ref YAMLScreenMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = WarningMessage(id, "Need to check these definitions.")
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLWarningMessage(id, str)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLWarningMessage_F(id, str) &
			BIND(C, NAME='YAMLWarningMessage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        character(KIND=C_CHAR), intent(in) :: str(*)
		END FUNCTION YAMLWarningMessage_F
    END INTERFACE
    integer, intent(in) :: id
    character(len=*), intent(in) :: str
	YAMLWarningMessage = YAMLWarningMessage_F(id, trim(str)//C_NULL_CHAR)
    END FUNCTION YAMLWarningMessage
    
END MODULE YAMLPhreeqcRM 
 
#endif
