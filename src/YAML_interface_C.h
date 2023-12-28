/*! @file YAML_interface_C.h
	@brief C header file for YAMLPhreeqcRM.
*/
#ifdef USE_YAML
#ifndef INC_YAML_interface_C_H
#define INC_YAML_interface_C_H

#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif
	/**
	Creates a YAMLPhreeqcRM instance with a YAML document that is ready to
	for writing data for initiation of a PhreeqcRM instance.
	@retval id   Id of the new YAMLPhreeqcRM instance.

	@see
	@ref DestroyYAMLPhreeqcRM.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	id = CreateYAMLPhreeqcRM();
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT int CreateYAMLPhreeqcRM(void);

	/**
	Deletes the YAMLPhreeqcRM instance and all data.
	@param id     The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval int   Zero indicates success, negative indicates failure.
	@see
	@ref YAMLClear.
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	n = DestroyYAMLPhreeqcRM(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT int DestroyYAMLPhreeqcRM(int id);

	/**
	Clears all definitions from the YAML document.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.

	@see
	@ref DestroyYAMLPhreeqcRM.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLClear(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLClear(int id);

	/**
	Writes YAML document to file.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param file_name   Name of file to write YAML document.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.

	@see
	@ref DestroyYAMLPhreeqcRM,
	@ref YAMLClear.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = WriteYAMLDoc(id, "AdvectBMI_f90.yaml");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT WriteYAMLDoc(int id, const char* file_name);

	/**

	Inserts data into the YAML document to select sets of output variables. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance. Sets of variables can be
	included or excluded with multiple calls to this method. All calls must
	precede the final call to @ref YAMLFindComponents. FindComponents generates
	SELECTED_OUTPUT 333 and USER_PUNCH 333 data blocks that make the variables
	accessible. Variables will only be accessible if the system includes the
	given reactant; for example, no gas variables will be created if there are
	no GAS_PHASEs in the model.

	@param id			The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option		A string value, among those listed below, that selects
	sets of variables that can be retieved by the bmif_get_value method.
	@param def A string value that can be "false", "true", or a list of items
	to be included as accessible variables. A value of "false", excludes all
	variables of the given type; a value of "true" includes all variables of
	the given type for the current system; a list specifies a subset of items
	of the given type.
	@retval IRM_RESULT 	Zero indicates success, negative indicates failure.
	<p>
	Values for the the parameter @a option:
	</p>
	<p>
	@n@a AddOutputVars: False excludes all variables; True causes the settings
	for each variable group to determine the variables that will be defined.
	Default True;

	@n@a SolutionProperties: False excludes all solution property variables;
	True includes variables pH, pe, alkalinity, ionic strength, water mass,
	charge balance, percent error, and specific conductance. Default True.

	@n@a SolutionTotalMolalities: False excludes all total element and element
	redox state variables;
	True includes all elements and element redox state variables for the system
	defined for the
	calculation; list restricts variables to the specified elements and redox
	states.
	Default True.

	@n@a ExchangeMolalities: False excludes all variables related to exchange;
	True includes all variables related to exchange; list includes variables
	for the specified exchange species. Default True.

	@n@a SurfaceMolalities: False excludes all variables related to surfaces;
	True includes all variables related to surfaces; list includes variables
	for the specified surface species. Default True.

	@n@a EquilibriumPhases: False excludes all variables related to equilibrium
	phases; True includes all variables related to equilibrium phases; list
	includes variables for the specified equilibiurm phases. Default True.

	@n@a Gases: False excludes all variables related to gases; True includes all
	variables related to gases; list includes variables for the specified gas
	components. Default True.

	@n@a KineticReactants: False excludes all variables related to kinetic
	reactants; True includes all variables related to kinetic reactants; list
	includes variables for the specified kinetic reactants. Default True.

	@n@a SolidSolutions: False excludes all variables related to solid
	solutions; True includes all variables related to solid solutions; list
	includes variables for the specified solid solutions components. Default
	True.

	@n@a CalculateValues: False excludes all calculate values; True includes all
	calculate values; list includes the specified calculate values.
	CALCLUATE_VALUES can be used to calculate geochemical quantities not
	available in the other sets of variables. Default True.

	@n@a SolutionActivities: False excludes all aqueous species; True includes
	all aqueous species; list includes only the specified aqueous species.
	Default False.

	@n@a SolutionMolalities: False excludes all aqueous species; True includes
	all aqueous species; list includes only the specified aqueous species.
	Default False.

	@n@a SaturationIndices: False excludes all saturation indices; True includes
	all saturation indices; list includes only the specified saturation
	indices. Default False.
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLAddOutputVars(id, "SolutionMolalities", "True");
	status = YAMLAddOutputVars(id, "SaturationIndices", "Calcite Dolomite");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLAddOutputVars(int id, char* option, char*
		def);

	/**
	Inserts data into the YAML document for the PhreeqcRM method CloseFiles.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	CloseFiles closes the output and log files.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLCloseFiles(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLCloseFiles(int id);

	/**
	Inserts data into the YAML document for the PhreeqcRM method CreateMapping.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param grid2chem     Integer array of mapping from user's model grid to
	cells for which chemistry will be run.
	@param dim     Dimension of the @a grid2chem vector.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	@a CreateMapping provides a mapping from grid cells in the user's model to
	reaction cells for which chemistry needs to be run. The mapping is used to
	eliminate inactive cells and to use symmetry to decrease the number of
	cells for which chemistry must be run. The array @a grid2chem of size @a
	nxyz (the number of grid cells) must contain the set of all integers 0 <=
	@a i < @a count_chemistry, where @a count_chemistry is a number less than
	or equal to @a nxyz. Inactive cells are assigned a negative integer. The
	mapping may be many-to-one to account for symmetry. Default is a
	one-to-one mapping--all user grid cells are reaction cells (equivalent to
	@a grid2chem values of 0,1,2,3,...,nxyz-1).
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	grid2chem = (int*)malloc(nxyz*sizeof(int));
	for(i=0; i < nxyz/2; i++) grid2chem[i] = i;
	status = YAMLCreateMapping(id, grid2chem, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLCreateMapping(int id, int* grid2chem, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method DumpModule.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param dump_on      Signal for writing the dump file, 1 or 0.
	@param append       Signal to append to the contents of the dump file,
	1 or 0.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a DumpModule writes the contents of all workers to file in _RAW formats (see
	appendix of PHREEQC version 3 manual), including SOLUTIONs and all
	reactants.
	</p>
	@see
	@ref YAMLSetDumpFileName.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetDumpFileName(id, "Advect_c.dmp");
	status = YAMLDumpModule(id, 1, 0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLDumpModule(int id, int dump_on, int append);

	/**
	Inserts data into the YAML document for the PhreeqcRM method FindComponents.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a FindComponents accumulates a list of elements. Elements are those that
	have been defined in a solution or any other reactant (EQUILIBRIUM_PHASE,
	KINETICS, and others), including charge imbalance. This method can be
	called multiple times and the list that is created is cummulative. The
	list is the set of components that needs to be transported. By default the
	list includes water, excess H and excess O (the H and O not contained in
	water); alternatively, the list may be set to contain total H and total O
	(@ref YAMLSetComponentH2O), which requires transport results to be
	accurate to eight or nine significant digits. If multicomponent diffusion
	(MCD) is to be modeled, there is a capability to retrieve aqueous species
	concentrations and to set new solution concentrations after MCD by using
	individual species concentrations (@ref YAMLSpeciesConcentrations2Module).
	To use these methods, the save-species property needs to be turned on
	(@ref YAMLSetSpeciesSaveOn). If the save-species property is on,
	FindComponents will generate a list of aqueous species, their diffusion
	coefficients at 25 C, and their charge.
	</p>
	@see
	@ref YAMLSetComponentH2O,
	@ref YAMLSetSpeciesSaveOn,
	@ref YAMLSpeciesConcentrations2Module.
	<p>
	The @a FindComponents method also generates lists of reactants--equilibrium
	phases, exchangers, gas components, kinetic reactants, solid solution
	components, and surfaces. The lists are cumulative, including all
	reactants that were defined in the initial phreeqc instance at any time
	@a FindComponents was called. In addition, a list of phases is generated for
	which saturation indices may be calculated from the cumulative list of
	components.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLFindComponents(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLFindComponents(int id);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialSolutions2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param solutions   Vector of SOLUTION index numbers that is dimensioned
	@a nxyz, where @a nxyz is the number of grid cells in the user's model.
	@param dim         Dimension of the @a solutions vector.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.
	<p>
	@a InitialSolutions2Module transfers SOLUTION definitions from the
	InitialPhreeqc instance to the reaction-module workers. @a solutions is a
	vector of SOLUTION index numbers that refer to definitions in the
	InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolutions2Module(int id, int*
		solutions, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialEquilibriumPhases2Module. When the YAML document is written to file
	it can be processed by the method InitializeYAML or the BMI method
	BMI_Initialize or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id                  The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param equilibrium_phases  Vector of EQUILIBRIUM_PHASES index numbers that
	is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	@param dim                 Dimension of the @a equilibrium_phases vector.
	@retval IRM_RESULT         Zero indicates success, negative indicates
	failure.
	<p>
	@a InitialEquilibriumPhases2Module transfers EQUILIBRIUM_PHASES definitions
	from the InitialPhreeqc instance to the reaction-module workers. @a
	equilibrium_phases is a vector of EQUILIBRIUM_PHASES index numbers that
	refer to definitions in the InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialEquilibriumPhases2Module(int id, int*
		equilibrium_phases, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialExchanges2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id             The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param exchanges      Vector of EXCHANGE index numbers that is
	dimensioned @a nxyz, where @a nxyz is the number of grid cells in the
	user's model.
	@param dim            Dimension of the @a exchanges vector.
	@retval IRM_RESULT    Zero indicates success, negative indicates failure.
	<p>
	@a InitialExchanges2Module transfers EXCHANGE definitions from the
	InitialPhreeqc instance to the reaction-module workers. @a exchanges is a
	vector of EXCHANGE index numbers that refer to definitions in the
	InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialExchanges2Module(int id, int*
		exchanges, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialSurfaces2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param surfaces     Vector of SURFACE index numbers that is dimensioned @a
	nxyz, where @a nxyz is the number of grid cells in the user's model.
	@param dim     Dimension of the @a surfaces vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a InitialSurfaces2Module transfers SURFACE definitions from the
	InitialPhreeqc instance to the reaction-module workers. @a surfaces is a
	vector of SURFACE index numbers that refer to definitions in the
	InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialSurfaces2Module(int id, int* surfaces,
		int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialGasPhases2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param gas_phases   Vector of GAS_PHASE index numbers that is dimensioned
	@a nxyz, where @a nxyz is the number of grid cells in the user's model.
	@param dim          Dimension of the @a gas_phases vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a InitialGasPhases2Module transfers GAS_PHASE definitions from the
	InitialPhreeqc instance to the reaction-module workers. @a gas_phases is a
	vector of GAS_PHASE index numbers that refer to definitions in the
	InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialGasPhases2Module(int id, int*
		gas_phases, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialSolidSolutions2Module. When the YAML document is written to file it
	can be processed by the method InitializeYAML or the BMI method
	BMI_Initialize or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id               The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param solid_solutions  Vector of SOLID_SOLUTIONS index numbers that is
	dimensioned @a nxyz, where @a nxyz is the number of grid cells in the
	user's model.
	@param dim              Dimension of the @a solid_solutions vector.
	@retval IRM_RESULT      Zero indicates success, negative indicates failure.
	<p>
	@a InitialSolidSolutions2Module transfers SOLID_SOLUTIONS definitions from
	the InitialPhreeqc instance to the reaction-module workers. @a
	solid_solutions is a vector of SOLID_SOLUTIONS index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolidSolutions2Module(int id, int*
		solid_solutions, int dim);
	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialKinetics2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param kinetics     Vector of KINETICS index numbers that is dimensioned
	@a nxyz, where @a nxyz is the number of grid cells in the user's model.
	@param dim          Dimension of the @a kinetics vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a InitialKinetics2Module transfers KINETICS definitions from the
	InitialPhreeqc instance to the reaction-module workers. @a kinetics is a
	vector of KINETICS index numbers that refer to definitions in the
	InitialPhreeqc instance.
	</p>
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialKinetics2Module(int id, int* kinetics,
		int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialPhreeqc2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param ic1         Vector of solution and reactant index numbers that
	refer to definitions in the InitialPhreeqc instance.
	@param dim         Dimension of the @a grid2chem vector.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.
	<p>
	@a InitialPhreeqc2Module transfers solutions and reactants from the
	InitialPhreeqc instance to the reaction-module workers. @a ic1 is used to
	select initial conditions, including solutions and reactants, for each
	cell of the model, without mixing. @a ic1 is dimensioned 7 times @a nxyz,
	where @a nxyz is the number of grid cells in the user's model. The
	dimension of 7 refers to solutions and reactants in the following order:
	(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4)
	GAS_PHASE, (5) SOLID_SOLUTIONS, and (6) KINETICS. The definition
	initial_solution1[3*nxyz + 99] = 2, indicates that cell 99 (0 based)
	contains the SURFACE definition (index 3) defined by SURFACE 2 in the
	InitialPhreeqc instance. Size is 7 times @a nxyz. The order of definitions
	is given above. Negative values are ignored, resulting in no definition of
	that entity for that cell.
	</p>
	@see
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLInitialPhreeqc2Module_mix.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	ic1=(int*)malloc(nxyz*7*sizeof(int));
	for(i=0; i< nxyz; i++)
	{
		ic1(i)  = 1;             // Solution 1
		ic1(1*nxyz + i) = -1;    // Equilibrium phases none
		ic1(2*nxyz + i) = 1;     // Exchange 1
		ic1(3*nxyz + i) = -1;    // Surface none
		ic1(4*nxyz + i) = -1;    // Gas phase none
		ic1(5*nxyz + i) = -1;    // Solid solutions none
		ic1(6*nxyz + i) = -1;    // Kinetics none
	}
	status = YAMLInitialPhreeqc2Module_mix(id, ic1, 7*nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module(int id, int* ic1, int
		dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialPhreeqc2Module. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param ic1         Vector of solution and reactant index numbers that refer
	to
	definitions in the InitialPhreeqc instance. Size is 7 times @a nxyz, where
	@a nxyz is the number of grid cells in the user's model. The order of
	reactants is given below and in the example. Negative values are ignored,
	resulting in no definition of that entity for that cell.
	@param ic2          Vector of solution and reactant index numbers that
	refer to
	definitions in the InitialPhreeqc instance. Nonnegative values of @a ic2
	result in mixing with the entities defined in @a ic1. Negative values
	result in no mixing. Size is 7 times @a nxyz.
	@param f1           Fraction of @a ic1 that mixes with (1 - @a f1) of @a
	ic2. Size is 7 times @a nxyz.
	@param dim          Dimension of the vectors.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.<
	<p>
	@a InitialPhreeqc2Module transfers solutions and reactants from the
	InitialPhreeqc instance to the reaction-module workers, possibly with
	mixing. In its simplest form, @a  ic1 is used to select initial
	conditions, including solutions and reactants, for each cell of the model,
	without mixing. The dimension of 7 refers to solutions and reactants in
	the following order: (0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE,
	(3) SURFACE, (4) GAS_PHASE, (5) SOLID_SOLUTIONS, and (6) KINETICS. The
	definition ic1[3*nxyz + 99] = 2, indicates that cell 99 (0 based) contains
	the SURFACE definition (index 3) defined by SURFACE 2 in the
	InitialPhreeqc instance (either by RunFile or RunString).
	</p>
	<p>
	It is also possible to mix solutions and reactants to obtain the initial
	conditions for cells. For mixing, @a initials_conditions2 contains numbers
	for a second entity that mixes with the entity defined in @a ic1. @a f1
	contains the mixing fraction for @a ic1, whereas (1 - @a f1) is the mixing
	fraction for @a ic2. The definitions ic1[3*nxyz + 99] = 2,
	initial_solution2[3*nxyz + 99] = 3, f1[3*nxyz + 99] = 0.25 indicates that
	cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3,
	where the surface compositions have been defined in the InitialPhreeqc
	instance. If the user number in @a ic2 is negative, no mixing occurs.
	</p>
	@see
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLInitialPhreeqc2Module.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	ic1=(int*)malloc(nxyz*7*sizeof(int));
	ic2=(int*)malloc(nxyz*7*sizeof(int));
	if1=(double*)malloc(nxyz*7*sizeof(double));
	for(i=0; i< nxyz; i++)
	{
		ic1(i)  = 1;             // Solution 1
		ic1(1*nxyz + i) = -1;    // Equilibrium phases none
		ic1(2*nxyz + i) = 1;     // Exchange 1
		ic1(3*nxyz + i) = -1;    // Surface none
		ic1(4*nxyz + i) = -1;    // Gas phase none
		ic1(5*nxyz + i) = -1;    // Solid solutions none
		ic1(6*nxyz + i) = -1;    // Kinetics none
	}
	for(i=0; i < 7*nxyz; i++) ic2[i] = -1;
	for(i=0; i < 7*nxyz; i++) f1[i] = 1;
	status = YAMLInitialPhreeqc2Module_mix(id, ic1, ic2, f1, 7*nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module_mix(int id, int* ic1,
		int* ic2, double* f1, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	InitialPhreeqcCell2Module. When the YAML document is written to file it can
	be processed by the method InitializeYAML or the BMI method BMI_Initialize
	or the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param n            Number that refers to a solution or MIX and associated
	reactants in the InitialPhreeqc instance.
	@param cell_numbers A vector of grid-cell numbers.
	@param dim          Dimension of the vector cell_numbers.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a InitialPhreeqcCell2Module uses a cell numbered @a n in the InitialPhreeqc
	instance to populate a series of transport cells. All reactants with the
	number @a n are transferred along with the solution. If MIX @a n exists,
	it is used for the definition of the solution. If @a n is negative, @a n
	is redefined to be the largest solution or MIX number in the
	InitialPhreeqc instance. All reactants for each cell in the list @a
	cell_numbers are removed before the cell definition is copied from the
	InitialPhreeqc instance to the workers.

	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqc2Module_mix.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	module_cells = (int*)malloc(2*sizeof(int));
	module_cells[0] = 18;
	module_cells[1] = 19;
	status = YAMLInitialPhreeqcCell2Module(id, -1, module_cells, 2);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqcCell2Module(int id,
		int n, int* cell_numbers, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method LoadDatabase.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param database   String containing the database name.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.
	<p>
	@a LoadDatabase loads a database for all IPhreeqc instances--workers,
	InitialPhreeqc, and Utility. All definitions of the reaction module are
	cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is
	read.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLLoadDatabase(id, "phreeqc.dat");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLLoadDatabase(int id, const char* database);

	/**
	Inserts data into the YAML document for the PhreeqcRM method LogMessage.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param str          String to be printed.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a LogMessage prints a message to the log file.
	</p>
	@see
	@ref YAMLOutputMessage,
	@ref YAMLScreenMessage,
	@ref YAMLWarningMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLLogMessage(id, "Finished section 1 of initialization");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLLogMessage(int id, const char* str);

	/**
	Inserts data into the YAML document for the PhreeqcRM method OpenFiles. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a OpenFiles opens the output and log files. Files are named prefix.chem.txt
	and prefix.log.txt based on the prefix defined by @ref YAMLSetFilePrefix.
	</p>
	@see
	@ref YAMLSetFilePrefix,
	@ref YAMLCloseFiles,
	@ref YAMLLogMessage,
	@ref YAMLOutputMessage, and
	@ref YAMLWarningMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetFilePrefix(id, "Advect_c");
	status = YAMLOpenFiles(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLOpenFiles(int id);

	/**
	Inserts data into the YAML document for the PhreeqcRM method OutputMessage.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param str          String to be printed.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a OutputMessage prints a message to the output file.
	</p>
	@see
	@ref YAMLLogMessage,
	@ref YAMLScreenMessage,
	@ref YAMLWarningMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLOutputMessage(id, "Finished section 1 of initialization");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLOutputMessage(int id, const char* str);

	/**
	Inserts data into the YAML document for the PhreeqcRM method RunCells. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a RunCells runs reactions for all cells in the reaction module. During
	initialization, RunCells can be used to equilibrate each solution with all
	reactants in a cell while using a time step of zero (@ref YAMLSetTimeStep)
	to avoid kinetic reactions. Other properties that may need to be
	initialized before RunCells is invoked include porosity (@ref
	YAMLSetPorosity), saturation (@ref YAMLSetSaturationUser), temperature
	(@ref YAMLSetTemperature), and pressure (@ref YAMLSetPressure).
	</p>
	@see
	@ref YAMLSetPorosity,
	@ref YAMLSetPressure,
	@ref YAMLSetSaturationUser,
	@ref YAMLSetTemperature,
	@ref YAMLSetTimeStep.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLRunCells(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLRunCells(int id);

	/**
	Inserts data into the YAML document for the PhreeqcRM method RunFile. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id               The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param workers          1, the workers will run the file; 0, the workers
	will not run the file.
	@param initial_phreeqc  1, the InitialPhreeqc instance will run the file;
	0, the InitialPhreeqc will not run the file.
	@param utility          1, the Utility instance will run the file; 0, the
	Utility instance will not run the file.
	@param file_name        Name of the file to run.
	@retval IRM_RESULT      Zero indicates success, negative indicates failure.
	<p>
	@a RunFile runs a PHREEQC input file. The first three arguments determine
	which IPhreeqc instances will run the file--the workers, the
	InitialPhreeqc instance, and (or) the Utility instance. Input files that
	modify the thermodynamic database should be run by all three sets of
	instances. Files with SELECTED_OUTPUT definitions that will be used during
	the time-stepping loop need to be run by the workers. Files that contain
	initial conditions or boundary conditions should be run by the
	InitialPhreeqc instance.
	</p>
	@see
	@ref YAMLRunString.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLRunFile(id, 1, 1, 1, "advect.pqi");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLRunFile(int id, int workers, int
		initial_phreeqc,
		int utility, const char* file_name);

	/**
	Inserts data into the YAML document for the PhreeqcRM method RunString. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id               The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param workers          1, the workers will run the string; 0, the workers
	will not run the string.
	@param initial_phreeqc  1, the InitialPhreeqc instance will run the
	string; 0, the InitialPhreeqc will not run the string.
	@param utility          1, the Utility instance will run the string; 0,
	the Utility instance will not run the string.
	@param input_string     String containing PHREEQC input.
	@retval IRM_RESULT      Zero indicates success, negative indicates failure.
	<p>
	@a RunString runs a PHREEQC input string. The first three arguments determine
	which IPhreeqc instances will run the string--the workers, the
	InitialPhreeqc instance, and (or) the Utility instance. Input strings that
	modify the thermodynamic database should be run by all three sets of
	instances. Strings with SELECTED_OUTPUT definitions that will be used
	during the time-stepping loop need to be run by the workers. Strings that
	contain initial conditions or boundary conditions should be run by the
	InitialPhreeqc instance.
	</p>
	@see
	@ref YAMLRunFile.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLRunString(id, 1, 1, 1, "DELETE; -all");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLRunString(int id, int workers, int
		initial_phreeqc,
		int utility, const char* input_string);

	/**
	Inserts data into the YAML document for the PhreeqcRM method ScreenMessage.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param str          String to be printed.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a ScreenMessage prints a message to the screen.
	</p>
	@see
	@ref YAMLLogMessage,
	@ref YAMLOutputMessage,
	@ref YAMLWarningMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLScreenMessage(id, "Beginning to process YAML for initial
	conditions");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLScreenMessage(int id, const char* str);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetComponentH2O. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf           1 (default), excess H, excess O, and water are
	included in the component list; 0, total H and O are included in the
	component list.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetComponentH2O selects whether to include H2O in the component list. The
	concentrations of H and O must be known accurately (8 to 10 significant
	digits) for the numerical method of PHREEQC to produce accurate pH and pe
	values. Because most of the H and O are in the water species, it may be
	more robust (require less accuracy in transport) to transport the excess H
	and O (the H and O not in water) and water. The default setting (1) is to
	include water, excess H, and excess O as components. A setting of 0 will
	include total H and total O as components. YAMLSetComponentH2O must be
	called before @ref YAMLFindComponents.
	</p>
	@see
	@ref YAMLFindComponents.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetComponentH2O(id, 0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetComponentH2O(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetConcentrations. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param c            Vector of component concentrations. Size of vector is
	@a ncomps times @a nxyz, where @a ncomps is the number of components as
	determined by FindComponents or GetComponentCount and @a nxyz is the
	number of grid cells in the user's model.
	@param dim          Dimension of the c vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	The only way to use this method is to have pre-calculated PHREEQC solution
	concentrations, which is not common. Concentrations are normally
	initialized with @ref YAMLInitialPhreeqc2Module or @ref
	YAMLInitialPhreeqcCell2Module.
	</p>
	@see
	@ref YAMLSetDensityUser,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser,
	@ref YAMLSetUnitsSolution.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetConcentrations(id, c, ncomps*nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetConcentrations(int id, double* c, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetCurrentSelectedOutputUserNumber. When the YAML document is written to
	file it can be processed by the method InitializeYAML or the BMI method
	BMI_Initialize or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param n_user       User number of the SELECTED_OUTPUT data block that is
	to be used.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetCurrentSelectedOutputUserNumber selects the current selected output by
	user number. The user may define multiple SELECTED_OUTPUT data blocks for
	the workers. A user number is specified for each data block. The value of
	the argument @a n_user selects which of the SELECTED_OUTPUT definitions
	will be used for selected-output operations.
	</p>
	@see
	@ref YAMLSetNthSelectedOutput,
	@ref YAMLSetSelectedOutputOn.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetCurrentSelectedOutputUserNumber(id, n_user);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetCurrentSelectedOutputUserNumber(int id,
		int n_user);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetDensityUser. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param density      Vector of densities. Size of vector is @a nxyz, where
	@a nxyz is the number
	of grid cells in the user's model.
	@param dim          Dimension of the @a density vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetDensityUser sets the density for each reaction cell. These density
	values are used when converting from transported mass-fraction
	concentrations (@ref YAMLSetUnitsSolution) to produce per liter
	concentrations during a call to SetConcentrations. They are also used when
	converting from reaction-cell concentrations to transport concentrations,
	if UseSolutionDensityVolume is set to 0 (false).
	</p>
	@see
	@ref YAMLSetUnitsSolution,
	@ref YAMLUseSolutionDensityVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	density = (double*)malloc(nxyz*sizeof(double));
	for(i = 0; i < nxyz; i++) density[i] = 1.0;
	status = YAMLSetDensityUser(id, density);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetDensityUser(int id, double* density, int
		dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetDumpFileName. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param file_name    Name of dump file.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetDumpFileName	sets the name of the dump file. It is the name used by
	the method DumpModule.
	</p>
	@see
	@ref YAMLDumpModule.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetDumpFileName(id, "Advect_c.dmp");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetDumpFileName(int id, const char*
		file_name);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetErrorHandlerMode. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param mode         Error handling mode: 0, 1, or 2.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetErrorHandlerMode sets the action to be taken when the reaction module
	encounters an error. Options are 0, return to calling program with an
	error return code (default); 1, throw an exception, in C++, the exception
	can be caught, for C and Fortran, the program will exit; or 2, attempt to
	exit gracefully.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetErrorHandlerMode(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorHandlerMode(int id, int mode);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetErrorOn.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf           1, enable error messages; 0, disable error messages.
	Default is 1.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetErrorOn sets the property that controls whether error messages are
	generated and displayed. Messages include PHREEQC "ERROR" messages, and
	any messages written with the method ErrorMessage.
	</p>
	@see
	@ref YAMLLogMessage,
	@ref YAMLOutputMessage,
	@ref YAMLScreenMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetErrorOn(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorOn(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetFilePrefix.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param prefix       Prefix used when opening the output and log files.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetFilePrefix sets the prefix for the output (prefix.chem.txt) and log
	(prefix.log.txt) files. These files are opened by the method OpenFiles.
	</p>
	@see
	@ref YAMLOpenFiles,
	@ref YAMLCloseFiles.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetFilePrefix(id, "Advect_c");
	status = YAMLOpenFiles(id);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetFilePrefix(int id, const char* prefix);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetGasCompMoles. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param  gas_moles   Vector of moles of gas components.
	@param dim          Dimension of vector @a gas_moles.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetGasCompMoles transfers moles of gas components from the vector given
	in the argument list (@a gas_moles) to each reaction cell. Dimension of
	the vector must be @a ngas_comps times @a nxyz, where, @a ngas_comps is
	the result of GetGasComponentsCount, and @a nxyz is the number of user
	grid cells. If the number of moles is set to a negative number, the gas
	component will not be defined for the GAS_PHASE of the reaction cell.
	</p>
	@see
	@ref YAMLFindComponents,
	@ref YAMLSetGasPhaseVolume.
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	gas_moles = (double*)malloc(nxyz*ngas*sizeof(double));
	status = YAMLSetGasCompMoles(id, gas_moles, nxyz*ngas);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetGasCompMoles(int id, double* gas_moles,
		int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetGasPhaseVolume. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param  gas_volume   Vector of volumes for each gas phase.
	@param dim           Dimension of the vector @a gas_volume.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	@a SetGasPhaseVolume transfers volumes of gas phases from the vector given
	in the argument list (@a gas_volume) to each reaction cell. The gas-phase
	volume affects the gas-component pressures calculated for fixed-volume
	gas phases. If a gas-phase volume is defined with this methood for a
	GAS_PHASE in a cell, the gas phase is forced to be a fixed-volume gas
	phase. Dimension of the vector is @a nxyz, where @a nxyz is the number of
	user grid cells. If the volume is set to a negative number for a cell,
	the gas-phase volume for that cell is not changed.
	</p>
	@see
	@ref YAMLFindComponents,
	@ref YAMLSetGasCompMoles.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	gas_volume = (double*)malloc(nxyz*sizeof(double));
	status = YAMLSetGasPhaseVolume(id, gas_volume, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetGasPhaseVolume(int id, double* gas_volume,
		int dim);

	/**
	Inserts data into the YAML document to define the number of cells in the
	user's model. Once the YAML document is written, the number of model
	cells can be extracted with the method GetGridCellCountYAML.
	GetGridCellCountYAML is NOT a PhreeqcRM method; it is a global method and
	must be used BEFORE the PhreeqcRM instance is created. SetGridCellCount
	will be ignored once the PhreeqcRM instance exists.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param count        Number of cells for the PhreeqcRM instance. The number
	of cells
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a YAMLSetGridCellCount can be used in the creation of the PhreeqcRM
	instance. The PhreeqcRM constructor takes two arguments.
	GetGridCellCountYAML provides the value for the first argument. If the
	YAML file does not contain the node "SetGridCellCount:",
	GetGridCellCountYAML will return zero.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetGridCellCount(nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetGridCellCount(int id, int count);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetNthSelectedOutput. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param n            Sequence number of the SELECTED_OUTPUT data block
	that is to be used (zero-based).
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetNthSelectedOutput specifies the current selected output by sequence
	number (one-based). The user may define multiple SELECTED_OUTPUT data
	blocks for the workers. A user number is specified for each data block,
	and the blocks are stored in user-number order. The value of the argument
	@a n selects the sequence number of the SELECTED_OUTPUT definition that
	will be used for selected-output operations.
	</p>
	@see
	@ref YAMLSetCurrentSelectedOutputUserNumber,
	@ref YAMLSetSelectedOutputOn.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetCurrentSelectedOutput(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetNthSelectedOutput(int id, int n);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetPartitionUZSolids. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize or
	the BMI method BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf           1, the fraction of solids and gases available for
	reaction is equal to the saturation; 0 (default), all solids and gases
	are reactive regardless of saturation.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetPartitionUZSolids sets the property for partitioning solids between
	the saturated and unsaturated parts of a partially saturated cell. The
	option is intended to be used by saturated-only flow codes that allow a
	variable water table. The value has meaning only when saturations less
	than 1.0 are encountered. The partially saturated cells may have a small
	water-to-rock ratio that causes reactions to proceed differently relative
	to fully saturated cells. By setting  @a SetPartitionUZSolids to true,
	the amounts of solids and gases are partioned according to the
	saturation. If a cell has a saturation of 0.5, then the water interacts
	with only half of the solids and gases; the other half is unreactive
	until the water table rises. As the saturation in a cell varies, solids
	and gases are transferred between the saturated and unsaturated
	(unreactive) reservoirs of the cell. Unsaturated-zone flow and transport
	codes will probably use the default (false), which assumes all gases and
	solids are reactive regardless of saturation.
	</p>
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetPartitionUZSolids(id, 0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetPartitionUZSolids(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetPorosity.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param por          Vector of porosities, unitless. Default is 0.1.
	@param dim          Dimension of the vector @ por.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetPorosity sets the porosity for each reaction cell. The volume of water
	in a reaction cell is the product of porosity, saturation
	(SetSaturationUser), and representative volume (SetRepresentativeVolume).
	Size of vector is @a nxyz, where @a nxyz is the number of grid cells in
	the user's model.
	</p>
	@see
	@ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	por = (double*)maloc(nxyz*sizeof(double));
	for(i = 0; i < nxyz; i++) por[i] = 0.2;
	status = YAMLSetPorosity(id, por, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetPorosity(int id, double* por, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetPressure.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param p            Vector of pressures, in atm. Size of vector is @a
	nxyz, where @a nxyz is the number of grid cells in the user's model.
	@param dim Dimension of the vector @a p.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetPressure sets the pressure for each reaction cell. Pressure effects
	are considered only in three of the databases distributed with PhreeqcRM:
	phreeqc.dat, Amm.dat, and pitzer.dat.
	</p>
	@see
	@ref YAMLSetTemperature.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	pressure = (double*)maloc(nxyz*sizeof(double));
	for(i = 0; i < nxyz; i++) pressure[i] = 2.0;
	status = YAMLSetPressure(id, pressure, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetPressure(int id, double* p, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetPrintChemistryMask. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param mask         Vector of integers. Size of vector is @a nxyz, where
	@a nxyz is the number of grid cells in the user's model. A value of 0
	will disable printing detailed output for the cell; a value of 1 will
	enable printing detailed output for a cell.
	@param dim          Dimension of the vector @ mask.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetPrintChemistryMask enables or disables detailed output for each
	reaction cell. Printing for a reaction cell will occur only when the
	printing is enabled with SetPrintChemistryOn and the @a mask value is 1.
	</p>
	@see
	@ref YAMLSetPrintChemistryOn.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	print_chemistry_mask = (int*)malloc(nxyz*sizeof(int));
	for(i=0; i < nxyz; i=1) print_chemistry_mask[i] = 1;
	status = YAMLSetPrintChemistryMask(id, print_chemistry_mask, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryMask(int id, int* mask,
		int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetPrintChemistryOn. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id               The instance id returned from @ref
	CreateYAMLPhreeqcRM.
	@param workers          1, enable detailed printing in the worker
	instances; 0, disable detailed printing in the worker instances.
	@param initial_phreeqc  1, enable detailed printing in the InitialPhreeqc
	instance; 0, disable detailed printing in the InitialPhreeqc instance.
	@param utility          1, enable detailed printing in the Utility
	instance; 0, disable detailed printing in the Utility instance.
	@retval IRM_RESULT      Zero indicates success, negative indicates failure.
	<p>
	@a SetPrintChemistryOn sets the property that enables or disables printing
	detailed output from reaction calculations to the output file for a set
	of cells defined by SetPrintChemistryMask. The detailed output prints all
	of the output typical of a PHREEQC reaction calculation, which includes
	solution descriptions and the compositions of all other reactants. The
	output can be several hundred lines per cell, which can lead to a very
	large output file (prefix.chem.txt opened by the method OpenFiles). For
	the worker instances, the output can be limited to a set of cells (method
	SetPrintChemistryMask) and, in general, the amount of information printed
	can be limited by use of options in the PRINT data block of PHREEQC
	(applied by using methods RunFile or RunString). Printing the detailed
	output for the workers is generally used only for debugging, and
	PhreeqcRM will run significantly faster when printing detailed output for
	the workers is disabled.
	</p>
	@see
	@ref YAMLOpenFiles,
	@ref YAMLRunFile,
	@ref YAMLRunString,
	@ref YAMLSetPrintChemistryMask.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetPrintChemistryOn(id, 0, 1, 0)
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryOn(int id, int workers, int
		initial_phreeqc,
		int utility);
	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetRebalanceByCell. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf          1, indicates individual cell times are used in
	rebalancing (default); 0, indicates average times are used in rebalancing.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetRebalanceByCell sets the load-balancing algorithm. PhreeqcRM attempts
	to rebalance the load of each thread or process such that each thread or
	process takes the same amount of time to run its part of a RunCells
	calculation. Two algorithms are available; one uses individual times for
	each cell and accounts for cells that were not run because saturation was
	zero (default), and the other assigns an average time to all cells. The
	methods are similar, but limited testing indicates the default method
	performs better.
	</p>
	@see
	@ref YAMLSetRebalanceFraction.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetRebalanceByCell(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceByCell(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetRebalanceFraction. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param f            Fraction from 0.0 to 1.0.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetRebalanceFraction sets the fraction of cells that are transferred
	among threads or processes when rebalancing. PhreeqcRM attempts to
	rebalance the load of each thread or process such that each thread or
	process takes the same amount of time to run its part of a RunCells
	calculation. The rebalancing transfers cell calculations among threads or
	processes to try to achieve an optimum balance. @a SetRebalanceFraction
	adjusts the calculated optimum number of cell transfers by a fraction
	from 0 to 1.0 to determine the actual number of cell transfers. A value
	of zero eliminates load rebalancing. A value less than 1.0 is suggested
	to slow the approach to the optimum cell distribution and avoid possible
	oscillations when too many cells are transferred at one iteration,
	requiring reverse transfers at the next iteration. Default is 0.5.
	</p>
	@see
	@ref YAMLSetRebalanceByCell.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetRebalanceFraction(id, 0.5);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceFraction(int id, double f);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetRepresentativeVolume. When the YAML document is written to file it can
	be processed by the method InitializeYAML or the BMI method BMI_Initialize
	to initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param rv           Vector of representative volumes, in liters. Default is
	1.0 liter.
	Size of array is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model.
	@param dim          Dimension of the vector @a rv.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetRepresentativeVolume sets the representative volume of each reaction
	cell. By default the representative volume of each reaction cell is 1
	liter. The volume of water in a reaction cell is determined by the
	product of the representative volume, the porosity (SetPorosity), and the
	saturation (SetSaturationUser). The numerical method of PHREEQC is more
	robust if the water volume for a reaction cell is within a couple orders
	of magnitude of 1.0. Small water volumes caused by small porosities and
	(or) small saturations (and (or) small representative volumes) may cause
	non-convergence of the numerical method. In these cases, a larger
	representative volume may help. Note that increasing the representative
	volume also increases the number of moles of the reactants in the
	reaction cell (minerals, surfaces, exchangers, and others), which are
	defined as moles per representative volume. @a SetRepresentativeVolume
	should be called before initial conditions are defined for the reaction
	cells.
	</p>
	@see
	@ref YAMLSetPorosity,
	@ref YAMLSetSaturationUser.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	rv = (double*)malloc(nxyz*sizeof(double));
	for(i=0; i < nxyz; i++) rv[i] = 1.0;
	status = YAMLSetRepresentativeVolume(id, rv, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetRepresentativeVolume(int id, double* rv,
		int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetSaturationUser. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param sat          Vector of saturations, unitless. Default 1.0. Size of
	vector is @a nxyz, where @a nxyz is the number of grid cells in the
	user's model.
	@param dim          Dimension of the @ sat vector.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetSaturationUser sets the saturation of each reaction cell. Saturation
	is a fraction ranging from 0 to 1. The volume of water in a cell is the
	product of porosity (SetPorosity), saturation (SetSaturationUser), and
	representative volume (SetRepresentativeVolume). As a result of a
	reaction calculation, solution properties (density and volume) will
	change; the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar
	volume data to calculate these changes. The methods GetDensityCalculated,
	GetSolutionVolume, and GetSaturationCalculated can be used to account for
	these changes in the succeeding transport calculation.
	</p>
	@see
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	sat = (double*)malloc(nxyz*sizeof(double));
	for(i=0; i < nxyz; i++) sat[i] = 1.0;
	status = YAMLSetSaturationUser(id, sat, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetSaturationUser(int id, double* sat, int
		dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetScreenOn.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf          1, enable screen messages; 0, disable screen messages.
	Default is true.
	@retval IRM_RESULT Zero indicates success, negative indicates failure.
	<p>
	@a SetScreenOn
	sets the property that controls whether messages are written to the
	screen. Messages include information about rebalancing during RunCells,
	and any messages written with ScreenMessage.
	</p>
	@see
	@ref YAMLRunCells,
	@ref YAMLScreenMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetScreenOn(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetScreenOn(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetSelectedOutputOn. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf           1, enable selected output; 0, disable selected output.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetSelectedOutputOn sets the property that controls whether
	selected-output results are available to be retrieved with
	GetSelectedOutput. 1 indicates that selected-output results will be
	accumulated during RunCells and can be retrieved with GetSelectedOutput;
	0 indicates that selected-output results will not be accumulated during
	RunCells.
	</p>
	@see
	@ref YAMLSetCurrentSelectedOutputUserNumber,
	@ref YAMLSetNthSelectedOutput,
	@ref YAMLSetSelectedOutputOn.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetSelectedOutputOn(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetSelectedOutputOn(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetSpeciesSaveOn. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param save_on      1 indicates species concentrations are saved; 0
	indicates species concentrations are not saved.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetSpeciesSaveOn sets the value of the species-save property. This method
	enables or disables use of PhreeqcRM with multicomponent-diffusion
	transport calculations. By default, concentrations of aqueous species are
	not saved. Setting the species-save property to @a true allows aqueous
	species concentrations to be retrieved with GetSpeciesConcentrations, and
	solution compositions to be set with SpeciesConcentrations2Module.
	SetSpeciesSaveOn must be called before calls to FindComponents.
	</p>
	@see
	@ref YAMLFindComponents,
	@ref YAMLSpeciesConcentrations2Module.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetSpeciesSaveOn(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetSpeciesSaveOn(int id, int save_on);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetTemperature. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tc           Vector of temperatures, in degrees C.
	@param dim          Dimension of vector @a tc.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetTemperature sets the temperature for each reaction cell. If
	@a SetTemperature is not called, worker solutions will have temperatures as
	defined by initial conditions (InitialPhreeqc2Module and
	InitialPhreeqcCell2Module). Size of vector is @a nxyz, where @a nxyz is
	the number of grid cells in the user's model.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPressure.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	tc = (double*)malloc(nxyz*sizeof(double));
	for(i=0; i < nxyz; i++) tc[i] = 1.0;
	status = YAMLSetTemperature(id, tc, nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetTemperature(int id, double* tc, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetTime. When
	the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param time         Current simulation time, in seconds.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetTime sets current simulation time for the reaction module.
	</p>
	@see
	@ref YAMLSetTimeStep,
	@ref YAMLSetTimeConversion.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetTime(id, 0.0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetTime(int id, double time);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetTimeConversion. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param conv_factor  Factor to convert seconds to user time units.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetTimeConversion sets a factor to convert from seconds to user time
	units. Factor times seconds produces user time units that is used in some
	PhreeqcRM printing.
	</p>
	@see
	@ref YAMLSetTime,
	@ref YAMLSetTimeStep.
	@par Fortran Example:

	@htmlonly
	<CODE>
	<PRE>
	conv_factor = 1.0 / 86400.;
	status = YAMLSetTimeConversion(id, conv_factor);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeConversion(int id, double conv_factor);

	/**
	Inserts data into the YAML document for the PhreeqcRM method SetTimeStep.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param time_step    Time step, in seconds.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetTimeStep sets current time step for the reaction module. This is the
	length of time over which kinetic reactions are integrated.
	</p>
	@see
	@ref YAMLSetTime,
	@ref YAMLSetTimeConversion.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	time_step = 86400.;
	status = YAMLSetTimeStep(id, time_step);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeStep(int id, double time_step);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsExchange. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option       Units option for exchangers: 0, 1, or 2.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsExchange sets input units for exchangers. In PHREEQC input,
	exchangers are defined by moles of exchange sites (@a Mp).
	SetUnitsExchange specifies how the number of moles of exchange sites in a
	reaction cell (@a Mc) is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a
	Mp*(1-P)*RV.
	</p>
	<p>
	If a single EXCHANGE definition is used for cells with different initial
	porosity, the three options scale quite differently. For option 0, the
	number of moles of exchangers will be the same regardless of porosity. For
	option 1, the number of moles of exchangers will be vary directly with
	porosity and inversely with rock volume. For option 2, the number of moles
	of exchangers will vary directly with rock volume and inversely with
	porosity.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsExchange(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsExchange(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsGasPhase. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option       Units option for gas phases: 0, 1, or 2.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsGasPhase
	sets input units for gas phases.
	In PHREEQC input, gas phases are defined by moles of component gases (@a
	Mp).
	@a SetUnitsGasPhase specifies how the number of moles of component gases in
	a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV, @a Mc = @a Mp*(1-@a
	P)*RV.
	</p>
	<p>
	If a single GAS_PHASE definition is used for cells with different initial
	porosity, the three options scale quite differently. For option 0, the
	number of moles of a gas component will be the same regardless of porosity.
	For option 1, the number of moles of a gas component will be vary directly
	with porosity and inversely with rock volume. For option 2, the number of
	moles of a gas component will vary directly with rock volume and inversely
	with porosity.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsGasPhase(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsGasPhase(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsKinetics. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option       Units option for kinetic reactants: 0, 1, or 2.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsKinetics
	sets input units for kinetic reactants.
	</p>
	<p>
	In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a
	Mp). @a SetUnitsKinetics specifies how the number of moles of kinetic
	reactants in a reaction cell (@a Mc) is calculated from the input value (@a
	Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV, @a Mc = @a Mp*(1-@a
	P)*RV.
	</p>
	<p>
	If a single KINETICS definition is used for cells with different initial
	porosity, the three options scale quite differently. For option 0, the
	number of moles of kinetic reactants will be the same regardless of
	porosity. For option 1, the number of moles of kinetic reactants will be
	vary directly with porosity and inversely with rock volume. For option 2,
	the number of moles of kinetic reactants will vary directly with rock
	volume and inversely with porosity.
	</p>
	<p>
	Note that the volume of water in a cell in the reaction module is equal to
	the product of porosity (SetPorosity), the saturation (SetSaturationUser),
	and representative volume (SetRepresentativeVolume), which is usually less
	than 1 liter. It is important to write the RATES definitions for
	homogeneous (aqueous) kinetic reactions to account for the current volume
	of water, often by calculating the rate of reaction per liter of water and
	multiplying by the volume of water (Basic function SOLN_VOL).
	</p>
	<p>
	Rates that depend on surface area of solids, are not dependent on the
	volume of water. However, it is important to get the correct surface area
	for the kinetic reaction. To scale the surface area with the number of
	moles, the specific area (m^2 per mole of reactant) can be defined as a
	parameter (KINETICS; -parm), which is multiplied by the number of moles of
	reactant (Basic function M) in RATES to obtain the surface area.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsKinetics(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsKinetics(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsPPassemblage. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option       Units option for equilibrium phases: 0, 1, or 2.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsPPassemblage sets input units for pure phase assemblages
	(equilibrium phases). In PHREEQC input, equilibrium phases are defined by
	moles of each phase (@a Mp). @a SetUnitsPPassemblage specifies how the
	number of moles of phases in a reaction cell (@a Mc) is calculated from the
	input value (@a Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV, @a Mc = @a
	Mp*(1-P)*RV.
	</p>
	<p>
	If a single EQUILIBRIUM_PHASES definition is used for cells with different
	initial porosity, the three options scale quite differently. For option 0,
	the number of moles of a mineral will be the same regardless of porosity.
	For option 1, the number of moles of a mineral will be vary directly with
	porosity and inversely with rock volume. For option 2, the number of moles
	of a mineral will vary directly with rock volume and inversely with
	porosity.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsPPassemblage(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsPPassemblage(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsSolution. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option       Units option for solutions: 1, 2, or 3, default is 1,
	mg/L.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsSolution sets solution concentration units used by the transport
	model. Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs. PHREEQC
	defines solutions by the number of moles of each element in the solution.
	</p>
	<p>
	To convert from mg/L to moles of element in the representative volume of a
	reaction cell, mg/L is converted to mol/L and multiplied by the solution
	volume, which is the product of porosity (SetPorosity), saturation
	(SetSaturationUser), and representative volume (SetRepresentativeVolume).
	To convert from mol/L to moles of element in the representative volume of a
	reaction cell, mol/L is multiplied by the solution volume. To convert from
	mass fraction to moles of element in the representative volume of a
	reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
	(SetDensityUser) and multiplied by the solution volume.
	</p>
	<p>
	To convert from moles of element in the representative volume of a reaction
	cell to mg/L, the number of moles of an element is divided by the solution
	volume resulting in mol/L, and then converted to mg/L. To convert from
	moles of element in a cell to mol/L,  the number of moles of an element is
	divided by the solution volume resulting in mol/L. To convert from moles of
	element in a cell to mass fraction, the number of moles of an element is
	converted to kg and divided by the total mass of the solution. Two options
	are available for the volume and mass of solution that are used in
	converting to transport concentrations: (1) the volume and mass of solution
	are calculated by PHREEQC, or (2) the volume of solution is the product of
	porosity (SetPorosity), saturation (SetSaturationUser), and representative
	volume (SetRepresentativeVolume), and the mass of solution is volume times
	density as defined by SetDensityUser. Which option is used is determined by
	UseSolutionDensityVolume.
	</p>
	@see
	@ref YAMLSetDensityUser,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser,
	@ref YAMLUseSolutionDensityVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsSolution(id, 2);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSolution(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsSSassemblage. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option        Units option for solid solutions: 0, 1, or 2.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsSSassemblage sets input units for solid-solution assemblages. In
	PHREEQC, solid solutions are defined by moles of each component (@a Mp). @a
	SetUnitsSSassemblage specifies how the number of moles of solid-solution
	components in a reaction cell (@a Mc) is calculated from the input value
	(@a Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV, @a Mc = @a Mp*(1-@
	P)*RV.
	</p>
	<p>
	If a single SOLID_SOLUTION definition is used for cells with different
	initial porosity, the three options scale quite differently. For option 0,
	the number of moles of a solid-solution component will be the same
	regardless of porosity. For option 1, the number of moles of a
	solid-solution component will be vary directly with porosity and inversely
	with rock volume. For option 2, the number of moles of a solid-solution
	component will vary directly with rock volume and inversely with porosity.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsSSassemblage(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSSassemblage(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SetUnitsSurface. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param option        Units option for surfaces: 0, 1, or 2.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	@a SetUnitsSurface sets input units for surfaces. In PHREEQC input, surfaces
	are defined by moles of surface sites (@a Mp). @a SetUnitsSurface specifies
	how the number of moles of surface sites in a reaction cell (@a Mc) is
	calculated from the input value (@a Mp).
	</p>
	<p>
	Options are 0, @a Mp is mol/L of RV (default), @a Mc = @a Mp*RV, where RV
	is the representative volume (SetRepresentativeVolume); 1, @a Mp is mol/L
	of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity
	(SetPorosity); or 2, @a Mp is mol/L of rock in the RV, @a Mc = @a Mp*(1-@a
	P)*RV.
	</p>
	<p>
	If a single SURFACE definition is used for cells with different initial
	porosity, the three options scale quite differently. For option 0, the
	number of moles of surface sites will be the same regardless of porosity.
	For option 1, the number of moles of surface sites will be vary directly
	with porosity and inversely with rock volume. For option 2, the number of
	moles of surface sites will vary directly with rock volume and inversely
	with porosity.
	</p>
	@see
	@ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSetUnitsSurface(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSurface(int id, int option);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	SpeciesConcentrations2Module. When the YAML document is written to file it
	can be processed by the method InitializeYAML or the BMI method
	BMI_Initialize to initialize a PhreeqcRM instance.

	@param id            The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param species_conc  Vector of aqueous species concentrations. Dimension of
	the array is @a nspecies times @a nxyz, where  @a nspecies is the number of
	aqueous species, and @a nxyz is the number of user grid cells.
	Concentrations are moles per liter.
	@param dim           Dimension of vector @a species_conc.
	@retval IRM_RESULT   Zero indicates success, negative indicates failure.
	<p>
	@a SpeciesConcentrations2Module sets solution concentrations in the reaction
	cells based on the vector of aqueous species concentrations (@a
	species_conc). This method is intended for use with
	multicomponent-diffusion transport calculations, and SetSpeciesSaveOn must
	be set to @a true. The list of aqueous species is determined by
	FindComponents and includes all aqueous species that can be made from the
	set of components. The method determines the total concentration of a
	component by summing the molarities of the individual species times the
	stoichiometric coefficient of the element in each species. Solution
	compositions in the reaction cells are updated with these component
	concentrations. Usually, accurate concentrations will not be known to use
	YAMLSpeciesConcentrations2Module during initialization.
	</p>
	@see
	@ref YAMLFindComponents,
	@ref YAMLSetSpeciesSaveOn.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLSpeciesConcentrations2Module(id, c, nspec*nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLSpeciesConcentrations2Module(int id, double*
		species_conc, int dim);

	/**
	Inserts data into the YAML document for the PhreeqcRM method StateSave.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param istate       Integer identifying the state that is saved.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a StateSave saves the state of the chemistry in all model cells, including
	SOLUTIONs, EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS,
	SOLID_SOLUTIONs, and SURFACEs. Although not generally used, MIXes,
	REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs will be saved for
	each cell, if they have been defined in the worker IPhreeqc instances. The
	distribution of cells among the workers and the chemistry of fully or
	partially unsaturated cells are also saved. The state is saved in memory;
	use DumpModule to save the state to file. PhreeqcRM can be reset to this
	state by using StateApply. A state is identified by an integer, and
	multiple states can be saved.
	</p>
	@see
	@ref YAMLDumpModule,
	@ref YAMLStateApply, and
	@ref YAMLStateDelete.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLStateSave(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLStateSave(int id, int istate);

	/**
	Inserts data into the YAML document for the PhreeqcRM method StateApply.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param istate       Integer identifying the state that is to be applied.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a StateApply resets the state of the module to a state previously saved with
	StateSave. The chemistry of all model cells are reset, including SOLUTIONs,
	EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and
	SURFACEs. MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
	will be reset for each cell, if they were defined in the worker IPhreeqc
	instances at the time the state was saved. The distribution of cells among
	the workers and the chemistry of fully or partially unsaturated cells are
	also reset to the saved state. The state to be applied is identified by an
	integer.
	</p>
	@see
	@ref YAMLStateSave and
	@ref YAMLStateDelete.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLStateApply(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLStateApply(int id, int istate);

	/**
	Inserts data into the YAML document for the PhreeqcRM method StateDelete.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML or the BMI method BMI_Initialize to initialize a PhreeqcRM
	instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param istate            Integer identifying the state that is to be deleted.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a StateDelete deletes a state previously saved with StateSave.
	</p>
	@see
	@ref YAMLStateSave and
	@ref YAMLStateApply.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLStateDelete(id, 1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLStateDelete(int id, int istate);

	/**
	Inserts data into the YAML document to define the number of threads to use
	with PhreeqcRM calculations. Once the YAML document is written, the number
	threads to use can be extracted when bmif_initialize is called. The data
	for ThreadCount will be ignored if the PhreeqcRM instance has already been
	initialized.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param nthreads            Number of threads to use for multiprocessing in
	PhreeqcRM instance. A value of zero will cause PhreeqcRM to use the number
	of logical processors available on the computer.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLThreadCount(id, 0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLThreadCount(int id, int nthreads);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	UseSolutionDensityVolume. When the YAML document is written to file it can
	be processed by the method InitializeYAML or the BMI method BMI_Initialize
	to initialize a PhreeqcRM instance.

	@param id          The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param tf          1 indicates that the solution density and volume as
	calculated by PHREEQC will be used to calculate concentrations. 0 indicates
	that the solution density set by SetDensityUser and the volume determined
	by the product of  SetSaturationUser, SetPorosity, and
	SetRepresentativeVolume, will be used to calculate concentrations retrieved
	by GetConcentrations.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a UseSolutionDensityVolume determines the volume and density to use when
	converting from the reaction-cell concentrations to transport
	concentrations (GetConcentrations). Two options are available to convert
	concentration units: (1) the density and solution volume calculated by
	PHREEQC are used, or (2) the specified density (SetDensityUser) and
	solution volume are determined by the product of saturation
	(SetSaturationUser), porosity (SetPorosity), and representative volume
	(SetRepresentativeVolume). Transport models that consider density-dependent
	flow will probably use the PHREEQC-calculated density and solution volume
	(default), whereas transport models that assume constant-density flow will
	probably use specified values of density and solution volume. Only the
	following databases distributed with PhreeqcRM have molar-volume
	information needed to accurately calculate density and solution volume:
	phreeqc.dat, Amm.dat, and pitzer.dat. Density is only used when converting
	to or from transport units of mass fraction.
	</p>
	@see
	@ref YAMLSetDensityUser,
	@ref YAMLSetPorosity,
	@ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = YAMLUseSolutionDensityVolume(id, 0);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLUseSolutionDensityVolume(int id, int tf);

	/**
	Inserts data into the YAML document for the PhreeqcRM method
	WarningMessage. When the YAML document is written to file it can be
	processed by the method InitializeYAML or the BMI method BMI_Initialize to
	initialize a PhreeqcRM instance.

	@param id           The instance id returned from @ref CreateYAMLPhreeqcRM.
	@param str          String to be printed.
	@retval IRM_RESULT  Zero indicates success, negative indicates failure.
	<p>
	@a WarningMessage
	prints a warning message to the screen and the log file.
	</p>
	@see
	@ref YAMLOpenFiles,
	@ref YAMLLogMessage,
	@ref YAMLOutputMessage,
	@ref YAMLScreenMessage.

	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	status = WarningMessage(id, "Need to check these definitions.");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT YAMLWarningMessage(int id, const char* str);

#if defined(__cplusplus)
}
#endif

#endif // INC_YAML_interface_C_H
#endif // USE_YAML
