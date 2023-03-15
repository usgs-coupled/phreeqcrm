#ifdef USE_YAML
MODULE YAML_interface
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
!> @param grid2chem        A vector of integers: Nonnegative is a reaction-cell number (0 based),
!> negative is an inactive cell. Vector is of size @a nxyz (number of grid cells).
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
!> status = YAMLSetDumpFileName("Advect_cpp.dmp")
!> status = YAMLDumpModule(dump_on, append)
!> </PRE>
!> </CODE>
!> @endhtmlonly
    INTEGER FUNCTION YAMLDumpModule(id, dump_on, append)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLDumpModule_F(id, dump_on, append) &
			BIND(C, NAME='YAMLDumpModule_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer, intent(in) :: id
        logical(kind=C_INT), intent(in) :: dump_on, append
		END FUNCTION YAMLDumpModule_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: dump_on, append   
	YAMLDumpModule = YAMLDumpModule_F(id, dump_on, append)
    END FUNCTION YAMLDumpModule
!> Inserts data into the YAML document for the PhreeqcRM method FindComponents.
!> When the YAML document is written to file it can be processed by the method InitializeYAML to
!> initialize a PhreeqcRM instance.
!> @param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
!> @retval IRM_RESULT   Zero indicates success, negative indicates failure.
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
!> double precision, allocatable, dimension(:,:) :: f1
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
    double precision, allocatable, dimension(:,:), intent(in) :: f1
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
!> @ref @ref YAMLInitialPhreeqc2Module_mix.
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
!> status = YAMLLogMessage("Finished section 1 of initialization");
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
!> status = YAMLSetFilePrefix("Advect_cpp");
!> status = YAMLOpenFiles();
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
!> status = YAMLOutputMessage("Finished section 1 of initialization");
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
    
    INTEGER FUNCTION YAMLRunFile(id, workers, initial_phreeqc, utility, file_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLRunFile_F(id, workers, initial_phreeqc, utility, file_name) &
			BIND(C, NAME='YAMLRunFile_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: workers, initial_phreeqc, utility
        character(KIND=C_CHAR), intent(in) :: file_name(*)
		END FUNCTION YAMLRunFile_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: workers, initial_phreeqc, utility
    character(len=*), intent(in) :: file_name
	YAMLRunFile = YAMLRunFile_F(id, workers, initial_phreeqc, utility, trim(file_name)//C_NULL_CHAR)
    END FUNCTION YAMLRunFile   
    
    INTEGER FUNCTION YAMLRunString(id, workers, initial_phreeqc, utility, input_string)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLRunString_F(id, workers, initial_phreeqc, utility, input_string) &
			BIND(C, NAME='YAMLRunString_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: workers, initial_phreeqc, utility
        character(KIND=C_CHAR), intent(in) :: input_string(*)
		END FUNCTION YAMLRunString_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: workers, initial_phreeqc, utility
    character(len=*), intent(in) :: input_string
	YAMLRunString = YAMLRunString_F(id, workers, initial_phreeqc, utility, trim(input_string)//C_NULL_CHAR)
    END FUNCTION YAMLRunString 
    
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
    
    
    INTEGER FUNCTION YAMLSetComponentH2O(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetComponentH2O_F(id, tf) &
			BIND(C, NAME='YAMLSetComponentH2O_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetComponentH2O_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetComponentH2O = YAMLSetComponentH2O_F(id, tf)
    END FUNCTION YAMLSetComponentH2O
    
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
    double precision, allocatable, dimension(:,:), intent(in) :: c
    integer :: dim
    dim = size(c,1)*size(c,2)
	YAMLSetConcentrations = YAMLSetConcentrations_F(id, c(1,1), dim)
    END FUNCTION YAMLSetConcentrations   
    
    INTEGER FUNCTION YAMLSetCurrentSelectedOutputUserNumber(id, n)
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
    integer, intent(in) :: n
	YAMLSetCurrentSelectedOutputUserNumber = YAMLSetCurrentSelectedOutputUserNumber_F(id, n)
    END FUNCTION YAMLSetCurrentSelectedOutputUserNumber    
    
    INTEGER FUNCTION YAMLSetDensity(id, density)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetDensity_F(id, density, dim) &
			BIND(C, NAME='YAMLSetDensity_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: density
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetDensity_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, allocatable, dimension(:), intent(in) :: density
	YAMLSetDensity = YAMLSetDensity_F(id, density(1), size(density))
    END FUNCTION YAMLSetDensity
    
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
    
    INTEGER FUNCTION YAMLSetErrorHandlerMode(id, n)
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
    integer, intent(in) :: n
	YAMLSetErrorHandlerMode = YAMLSetErrorHandlerMode_F(id, n)
    END FUNCTION YAMLSetErrorHandlerMode  
    
    INTEGER FUNCTION YAMLSetErrorOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetErrorOn_F(id, tf) &
			BIND(C, NAME='YAMLSetErrorOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetErrorOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetErrorOn = YAMLSetErrorOn_F(id, tf)
    END FUNCTION YAMLSetErrorOn   
    
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
    
    INTEGER FUNCTION YAMLSetGasCompMoles(id, gmoles)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetGasCompMoles_F(id, gmoles, dim) &
			BIND(C, NAME='YAMLSetGasCompMoles_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: gmoles
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetGasCompMoles_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, allocatable, dimension(:,:), intent(in) :: gmoles
    integer :: dim
    dim = size(gmoles,1)*size(gmoles,2)
	YAMLSetGasCompMoles = YAMLSetGasCompMoles_F(id, gmoles(1,1), dim)
    END FUNCTION YAMLSetGasCompMoles   
    
    INTEGER FUNCTION YAMLSetGasPhaseVolume(id, vol)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetGasPhaseVolume_F(id, vol, dim) &
			BIND(C, NAME='YAMLSetGasPhaseVolume_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: vol
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetGasPhaseVolume_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, allocatable, dimension(:), intent(in) :: vol
	YAMLSetGasPhaseVolume = YAMLSetGasPhaseVolume_F(id, vol(1), size(vol))
    END FUNCTION YAMLSetGasPhaseVolume
    
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
    
    INTEGER FUNCTION YAMLSetPartitionUZSolids(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPartitionUZSolids_F(id, tf) &
			BIND(C, NAME='YAMLSetPartitionUZSolids_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetPartitionUZSolids_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetPartitionUZSolids = YAMLSetPartitionUZSolids_F(id, tf)
    END FUNCTION YAMLSetPartitionUZSolids       
    
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
    double precision, allocatable, dimension(:), intent(in) :: por
	YAMLSetPorosity = YAMLSetPorosity_F(id, por(1), size(por))
    END FUNCTION YAMLSetPorosity
    
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
    double precision, allocatable, dimension(:), intent(in) :: p
	YAMLSetPressure = YAMLSetPressure_F(id, p(1), size(p))
    END FUNCTION YAMLSetPressure
    
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
    
    INTEGER FUNCTION YAMLSetPrintChemistryOn(id, workers, initial_phreeqc, utility)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetPrintChemistryOn_F(id, workers, initial_phreeqc, utility) &
			BIND(C, NAME='YAMLSetPrintChemistryOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: workers, initial_phreeqc, utility
		END FUNCTION YAMLSetPrintChemistryOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: workers, initial_phreeqc, utility
	YAMLSetPrintChemistryOn = YAMLSetPrintChemistryOn_F(id, workers, initial_phreeqc, utility)
    END FUNCTION YAMLSetPrintChemistryOn
    
    INTEGER FUNCTION YAMLSetRebalanceByCell(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetRebalanceByCell_F(id, tf) &
			BIND(C, NAME='YAMLSetRebalanceByCell_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetRebalanceByCell_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetRebalanceByCell = YAMLSetRebalanceByCell_F(id, tf)
    END FUNCTION YAMLSetRebalanceByCell  
    
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
    double precision, intent(in) :: f
	YAMLSetRebalanceFraction = YAMLSetRebalanceFraction_F(id, f)
    END FUNCTION YAMLSetRebalanceFraction  
    
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
    double precision, allocatable, dimension(:), intent(in) :: rv
	YAMLSetRepresentativeVolume = YAMLSetRepresentativeVolume_F(id, rv(1), size(rv))
    END FUNCTION YAMLSetRepresentativeVolume
    
    INTEGER FUNCTION YAMLSetSaturation(id, sat)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSaturation_F(id, sat, dim) &
			BIND(C, NAME='YAMLSetSaturation_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: sat
        integer(kind=C_INT), intent(in) :: dim
		END FUNCTION YAMLSetSaturation_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, allocatable, dimension(:), intent(in) :: sat
	YAMLSetSaturation = YAMLSetSaturation_F(id, sat(1), size(sat))
    END FUNCTION YAMLSetSaturation
    
    INTEGER FUNCTION YAMLSetScreenOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetScreenOn_F(id, tf) &
			BIND(C, NAME='YAMLSetScreenOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetScreenOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetScreenOn = YAMLSetScreenOn_F(id, tf)
    END FUNCTION YAMLSetScreenOn  
    
    INTEGER FUNCTION YAMLSetSelectedOutputOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSelectedOutputOn_F(id, tf) &
			BIND(C, NAME='YAMLSetSelectedOutputOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetSelectedOutputOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetSelectedOutputOn = YAMLSetSelectedOutputOn_F(id, tf)
    END FUNCTION YAMLSetSelectedOutputOn  
    
    INTEGER FUNCTION YAMLSetSpeciesSaveOn(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetSpeciesSaveOn_F(id, tf) &
			BIND(C, NAME='YAMLSetSpeciesSaveOn_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLSetSpeciesSaveOn_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLSetSpeciesSaveOn = YAMLSetSpeciesSaveOn_F(id, tf)
    END FUNCTION YAMLSetSpeciesSaveOn 
    
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
    double precision, allocatable, dimension(:), intent(in) :: tc
	YAMLSetTemperature = YAMLSetTemperature_F(id, tc(1), size(tc))
    END FUNCTION YAMLSetTemperature
    
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
    double precision, intent(in) :: time
	YAMLSetTime = YAMLSetTime_F(id, time)
    END FUNCTION YAMLSetTime  
    
    INTEGER FUNCTION YAMLSetTimeConversion(id, conv)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTimeConversion_F(id, conv) &
			BIND(C, NAME='YAMLSetTimeConversion_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: conv
		END FUNCTION YAMLSetTimeConversion_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, intent(in) :: conv
	YAMLSetTimeConversion = YAMLSetTimeConversion_F(id, conv)
    END FUNCTION YAMLSetTimeConversion 
    
    INTEGER FUNCTION YAMLSetTimeStep(id, time)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetTimeStep_F(id, time) &
			BIND(C, NAME='YAMLSetTimeStep_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        real(kind=C_DOUBLE), intent(in) :: time
		END FUNCTION YAMLSetTimeStep_F
    END INTERFACE
    integer, intent(in) :: id
    double precision, intent(in) :: time
	YAMLSetTimeStep = YAMLSetTimeStep_F(id, time)
    END FUNCTION YAMLSetTimeStep  
    
    INTEGER FUNCTION YAMLSetUnitsExchange(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsExchange_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsExchange_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsExchange_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsExchange = YAMLSetUnitsExchange_F(id, n)
    END FUNCTION YAMLSetUnitsExchange 
    
    INTEGER FUNCTION YAMLSetUnitsGasPhase(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsGasPhase_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsGasPhase_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsGasPhase_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsGasPhase = YAMLSetUnitsGasPhase_F(id, n)
    END FUNCTION YAMLSetUnitsGasPhase 
    
    INTEGER FUNCTION YAMLSetUnitsKinetics(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsKinetics_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsKinetics_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsKinetics_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsKinetics = YAMLSetUnitsKinetics_F(id, n)
    END FUNCTION YAMLSetUnitsKinetics
    
    INTEGER FUNCTION YAMLSetUnitsPPassemblage(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsPPassemblage_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsPPassemblage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsPPassemblage_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsPPassemblage = YAMLSetUnitsPPassemblage_F(id, n)
    END FUNCTION YAMLSetUnitsPPassemblage
    
    INTEGER FUNCTION YAMLSetUnitsSolution(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSolution_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsSolution_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsSolution_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsSolution = YAMLSetUnitsSolution_F(id, n)
    END FUNCTION YAMLSetUnitsSolution
    
    INTEGER FUNCTION YAMLSetUnitsSSassemblage(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSSassemblage_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsSSassemblage_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsSSassemblage_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsSSassemblage = YAMLSetUnitsSSassemblage_F(id, n)
    END FUNCTION YAMLSetUnitsSSassemblage
    
    INTEGER FUNCTION YAMLSetUnitsSurface(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLSetUnitsSurface_F(id, n) &
			BIND(C, NAME='YAMLSetUnitsSurface_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        integer(kind=C_INT), intent(in) :: n
		END FUNCTION YAMLSetUnitsSurface_F
    END INTERFACE
    integer, intent(in) :: id
    integer, intent(in) :: n
	YAMLSetUnitsSurface = YAMLSetUnitsSurface_F(id, n)
    END FUNCTION YAMLSetUnitsSurface
    
    INTEGER FUNCTION YAMLSpeciesConcentrations2Module(id, c)
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
    double precision, allocatable, dimension(:,:), intent(in) :: c
    integer :: dim
    dim = size(c,1)*size(c,2)
	YAMLSpeciesConcentrations2Module = YAMLSpeciesConcentrations2Module_F(id, c(1,1), dim)
    END FUNCTION YAMLSpeciesConcentrations2Module   
    
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
    
    INTEGER FUNCTION YAMLUseSolutionDensityVolume(id, tf)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION YAMLUseSolutionDensityVolume_F(id, tf) &
			BIND(C, NAME='YAMLUseSolutionDensityVolume_F')
		USE ISO_C_BINDING
		IMPLICIT NONE
        integer(kind=C_INT), intent(in) :: id
        logical(kind=C_INT), intent(in) :: tf
		END FUNCTION YAMLUseSolutionDensityVolume_F
    END INTERFACE
    integer, intent(in) :: id
    logical(kind=4), intent(in) :: tf
	YAMLUseSolutionDensityVolume = YAMLUseSolutionDensityVolume_F(id, tf)
    END FUNCTION YAMLUseSolutionDensityVolume 
    
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
    
END MODULE YAML_interface 
 
#endif
