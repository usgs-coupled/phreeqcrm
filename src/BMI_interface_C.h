#ifdef SKIP
/*! @file BMI_interface_C.h
	@brief C header file for BMIPhreeqcRM.
*/
#if !defined(BMI_INTERFACE_C_H_INCLUDED)
#define BMI_INTERFACE_C_H_INCLUDED
#include "IrmResult.h"
#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif	
	/**
    @a BMI_AddOutputVars allows selection of sets of variables that can be retieved
    by the @a BMI_GetValue methods. Sets of variables can be included or excluded with
    multiple calls to this method. All calls must precede the final call to
    the PhreeqcRM method FindComponents. FindComponents generates SELECTED_OUTPUT 333 and
    USER_PUNCH 333 data blocks that make the variables accessible. Variables will
    only be accessible if the system includes the given reactant; for example, no
    gas variables will be Created if there are no GAS_PHASEs in the model.
    @param id Id number returned by @ref BMI_Create.
    @param option A string value, among those listed below, that includes or
    excludes variables from @ref BMI_GetOutputVarName, @a BMI_GetValue methods,
    and other BMI methods.
    @param def A string value that can be "false", "true", or a list of items to be included as
    accessible variables. A value of "false", excludes all variables of the given type; a
    value of "true" includes all variables of the given type for the current system; a list
    specifies a subset of items of the given type.
    @retval    0 is success, 0 is failure.
    <p>
    Values for the the parameter @a option:
    </p>
    @n@a AddOutputVars: False excludes all variables; True causes the settings for each variable group
    to determine the variables that will be defined. Default True;
    @n@a SolutionProperties: False excludes all solution property variables; True includes variables pH, pe,
    alkalinity, ionic strength, water mass, charge balance, percent error, and specific conductance.
    Default True.
    @n@a SolutionTotalMolalities: False excludes all total element and element redox state variables;
    True includes all elements and element redox state variables for the system defined for the
    calculation; list restricts variables to the specified elements and redox states.
    Default True.
    @n@a ExchangeMolalities: False excludes all variables related to exchange; True includes all
    variables related to exchange; list includes variables for the specified exchange species.
    Default True.
    @n@a SurfaceMolalities: False excludes all variables related to surfaces; True includes all
    variables related to surfaces; list includes variables for the specified surface species.
    Default True.
    @n@a EquilibriumPhases: False excludes all variables related to equilibrium phases; True includes all
    variables related to equilibrium phases; list includes variables for the specified
    equilibiurm phases. Default True.
    @n@a Gases: False excludes all variables related to gases; True includes all
    variables related to gases; list includes variables for the specified gas components. Default True.
    @n@a KineticReactants: False excludes all variables related to kinetic reactants; True includes all
    variables related to kinetic reactants; list includes variables for the specified kinetic
    reactants. Default True.
    @n@a SolidSolutions: False excludes all variables related to solid solutions; True includes all
    variables related to solid solutions; list includes variables for the specified solid solutions
    components. Default True.
    @n@a CalculateValues: False excludes all calculate values; True includes all
    calculate values; list includes the specified calculate values. CALCLUATE_VALUES can be
    used to calculate geochemical quantities not available in the other sets of variables.
    Default True.
    @n@a SolutionActivities: False excludes all aqueous species; True includes all
    aqueous species; list includes only the specified aqueous species. Default False.
    @n@a SolutionMolalities: False excludes all aqueous species; True includes all
    aqueous species; list includes only the specified aqueous species. Default False.
    @n@a SaturationIndices: False excludes all saturation indices; True includes all
    saturation indices; list includes only the specified saturation indices. Default False.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_AddOutputVars(id, "SolutionMolalities", "True");
    status = BMI_AddOutputVars(id, "SaturationIndices", "Calcite Dolomite");
    </PRE>
    </CODE>
    @endhtmlonly
	*/
	IRM_DLL_EXPORT IRM_RESULT BMI_AddOutputVars(int id, char* option, char* def);
	/**
	@a BMI_Create Creates a reaction module. If the code is compiled with
	the preprocessor directive USE_OPENMP, the reaction module is multithreaded.
	If the code is compiled with the preprocessor directive USE_MPI, the reaction
	module will use MPI and multiple processes. If neither preprocessor directive is used,
	the reaction module will be serial (unparallelized).
	@param nxyz                         The number of grid cells in the user's model.
	@param nthreads (or @a comm, MPI)   When using OPENMP, the argument (@a nthreads)
	is the number of worker threads to be used.
	If @a nthreads <= 0, the number of threads is set equal to the number of
	processors of the computer.
	When using MPI, the argument (@a comm) is the MPI communicator to use within
	the reaction module.
	@retval Id of the BMIPhreeqcRM instance, negative is failure.
	@see
	@ref BMI_Finalize.
	@par C example:
	@htmlonly
	<CODE>
	<PRE>
	nxyz = 40;
	nthreads = 3;
	id = BMI_Create(nxyz, nthreads);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root and workers.
	*/
	IRM_DLL_EXPORT int        BMI_Create(int nxyz, int nthreads);
	//IRM_DLL_EXPORT int        BMI_Create_default();
    /**
    @a BMI_Destroy Destroys a reaction module; same as @ref BMI_Finalize.
    @param id Id number returned by @ref BMI_Create.
    @retval    0 is success, 0 is failure.
    @see
    @ref BMI_Create.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_Destroy(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root and workers.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_Destroy(int id);
    /**
    @a BMI_Finalize Destroys a reaction module.
    @param id Id number returned by @ref BMI_Create.
    @retval    0 is success, 0 is failure.
    @see
    @ref BMI_Create.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_Finalize(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root and workers.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_Finalize(int id);
    /**
    @a BMI_GetComponentName returns the component name--"BMIPhreeqcRM".
    @param id Id number returned by @ref BMI_Create.
    @param component_name Returns "BMIPhreeqcRM", the name of the component.
    @param l Length of string buffer @a component_name.
    @retval               0 is success, 1 is failure; negative indicates buffer is too small.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_GetComponentName(id, component_name);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetComponentName(int id, char* component_name, int l);
    /**
    @a BMI_GetCurrentTime returns the current simulation time, in seconds. 
    @param id Id number returned by @ref BMI_Create.
    @retval        The current simulation time, in seconds.
    @see
    @ref BMI_GetEndTime,
    @ref BMI_GetTimeStep,
    @ref BMI_GetTime.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    now = BMI_GetCurrentTime(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT double     BMI_GetCurrentTime(int id);
    /**
    @a BMI_GetEndTime returns @ref BMI_GetCurrentTime plus
    @ref BMI_GetTimeStep, in seconds.
    @param id Id number returned by @ref BMI_Create.
    @retval          The end of the time step, in seconds.
    @see
    @ref BMI_GetCurrentTime,
    @ref BMI_GetTimeStep.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_GetEndTime(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT double     BMI_GetEndTime(int id);
    /**
    @a BMI_GetGridRank returns a rank of 1 for grid 0.
    BMIPhreeqcRM has a 1D series of
    cells; any grid or spatial information must
    be found in the user's model.
    @param id Id number returned by @ref BMI_Create.
    @param grid   Grid number, only grid 0 is considered.
    @retval       Rank of 1 is returned for grid 0; 0 for
    all other values of @a grid.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    rank = BMI_GetGridRank(id, grid)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetGridRank(int id, int grid);
    /**
    @ref BMI_GetGridSize returns the number of cells specified
    at creation of the BMIPhreeqcRM instance.
    @param id Id number returned by @ref BMI_Create.
    @param grid  Grid number, only grid 0 is considered.
    @retval    Number of cells. Same value as @ref BMI_GetValueInt(id, "GridCellCount")
    is returned for grid 0;
    0 for all other values of @a grid.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    nxyz = BMI_GetGridSize(id, grid);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetGridSize(int id, int grid);
    /**
    @a BMI_GetGridType defines the grid to be points. No grid
    information is available in BMIPhreeqcRM; all grid
    information must be found in the user's model.
    @param id    Id number returned by @ref BMI_Create.
    @param grid  Grid number, only grid 0 is considered.
    @param str   "Points" is returned for grid 0;
    "Undefined grid identifier" is returned for all other
    values of @a grid.
    @param l Length of string buffer @a str.
    @retval      0 is success, 1 is failure, negative indicates the buffer is too small.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_GetGridType(id, grid, str, l)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetGridType(int id, int grid, char* str, int l);
    /**
    @a BMI_GetInputItemCount returns count of variables that
    can be set with @a BMI_SetValue methods.
    @param id Id number returned by @ref BMI_Create.
    @retval   Number of input variables that can be set with @a BMI_SetValue methods.
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = BMI_GetInputItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetInputItemCount(int id);
	//IRM_DLL_EXPORT int        BMI_GetInputVarNamesSize(int id);
    /**
    @a BMI_GetInputVarName returns the ith variable name that can be set
    with @a BMI_SetValue methods.
    @param id Id number returned by @ref BMI_Create.
    @param i 0-based index of variable name to retrieve.
    @param name   Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256];
    status = BMI_GetInputVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetInputVarName(int id, int i, char* name, int l);
	//IRM_DLL_EXPORT int        BMI_GetNames(int id, const char* type, char* dest);
	//IRM_DLL_EXPORT int        BMI_GetNamesSize(int id, const char* type, int* dest);
    /**
    @a BMI_GetOutputItemCount returns count of output variables that can be
    retrieved with @a BMI_GetValue methods.
    @param id Id number returned by @ref BMI_Create.
    @retval   Number of output variables that can be retrieved with @a BMI_GetValue methods.
    @see
    @ref BMI_GetOutputVarName,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = BMI_GetOutputItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetOutputItemCount(int id);
	//IRM_DLL_EXPORT int        BMI_GetOutputVarNamesSize(int id);
    /**
    @a BMI_GetOutputVarName returns ith variable name for which a pointer can be
    retrieved with @a BMI_GetValue methods.
    @param id Id number returned by @ref BMI_Create.
    @param i 0-based index of variable name to retrieve.
    @param name    Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256]
    status = BMI_GetOutputVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetOutputVarName(int id, int i, char* name, int l);
    /**
    @a BMI_GetPointableItemCount returns count of pointable variables that can be
    retrieved with @ref BMI_GetValuePtr.
    @param id Id number returned by @ref BMI_Create.
    @retval   Number of output variables that can be retrieved with @ref BMI_GetValuePtr.
    @see
    @ref BMI_GetPointableVarName,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = BMI_GetPointableItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int		  BMI_GetPointableItemCount(int id);
	//IRM_DLL_EXPORT int        BMI_GetPointableVarNamesSize(int id);
    /**
    @a BMI_GetPointableVarName returns ith variable name for which a pointer can be
    retrieved with @ref BMI_GetValuePtr.
    @param id Id number returned by @ref BMI_Create.
    @param i 0-based index of variable name to retrieve.
    @param name    Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetPointableItemCount,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256];
    status = BMI_GetPointableVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetPointableVarName(int id, int i, char* name, int l);
    /**
    @a BMI_GetStartTime returns the current simulation time, in seconds.
    (Same as @ref BMI_GetCurrentTime.)
    @param id Id number returned by @ref BMI_Create.
    @retval       The current simulation time, in seconds.
    */
	IRM_DLL_EXPORT double     BMI_GetStartTime(int id);
    /**
    @a BMI_GetTime returns the current simulation time, in seconds.
    (Same as @ref BMI_GetCurrentTime.)
    @param id Id number returned by @ref BMI_Create.
    @retval       The current simulation time, in seconds.
    */
    IRM_DLL_EXPORT double     BMI_GetTime(int id);
    /**
    @a BMI_GetTimeStep returns the current simulation time step,
    in seconds.
    @param id Id number returned by @ref BMI_Create.
    @retval              The current simulation time step, in seconds.
    @see
    @ref BMI_GetCurrentTime,
    @ref BMI_GetEndTime.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    time_step = BMI_GetTimeStep(id)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT double     BMI_GetTimeStep(int id);
    /**
    @a BMI_GetTimeUnits returns the time units of PhreeqcRM.
    All time units are seconds for PhreeqcRM.
    @param id Id number returned by @ref BMI_Create.
    @param units    Returns the string "seconds".
    @param l    Length of the string buffer @a units.
    @retval              0 is success, 1 failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetCurrentTime,
    @ref BMI_GetEndTime,
    @ref BMI_GetTimeStep.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char time_units[256]
    status = BMI_GetTimeUnits(id, time_units, 256)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetTimeUnits(int id, char* units, int l);
	// GetValue
    /**
    @a BMI_GetValueInt retrieves int model variables. Only variables in the list
    provided by @ref BMI_GetOutputVarName can be retrieved.
    @param id Id number returned by @ref BMI_Create.
    @param var    Name of the variable to retrieve.
    @param dest   Variable in which to place results.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var).
    </p>
    @n "ComponentCount"
    @n "CurrentSelectedOutputUserNumber"
    @n "GridCellCount"
    @n "SelectedOutputColumnCount"
    @n "SelectedOutputCount"
    @n "SelectedOutputOn"
    @n "SelectedOutputRowCount".
    @see
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_GetValueInt(id, "ComponentCount", &count);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueInt(int id, char* var, int* dest);
    /**
    @a BMI_GetValueDouble retrieves model variables. Only variables in the list
    provided by @ref BMI_GetOutputVarName can be retrieved. 
    @param id Id number returned by @ref BMI_Create.
    @param var    Name of the variable to retrieve.
    @param dest   Variable in which to place results.
    @retval       0 is success, 1 is failure.
    <p>
    Variables in addition to the ones listed below may be retrieved by this method,
    depending on variables selected by @ref BMI_AddOutputVars. All variables added
    by @ref BMI_AddOutputVars will be double arrays of size equal to the number of 
    model cells [@ref BMI_GetValueInt(id, "GridCellCount")].
    </p>
    <p>
    Variable names for the second argument (@a var).
    </p>
    @n "Concentrations"
    @n "DensityCalculated"
    @n "Gfw"
    @n "Porosity"
    @n "Pressure"
    @n "SaturationCalculated"
    @n "SelectedOutput"
    @n "SolutionVolume"
    @n "Temperature"
    @n "Time"
    @n "TimeStep"
    @n "Viscosity".
    @see
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    density = (double *)malloc(nxyz*sizeof(double));
    status = BMI_GetValueDouble(id, "DensityCalculated", density);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueDouble(int id, char* var, double* dest);
    /**
    @a BMI_GetValueChar retrieves char model variables. Only variables in the list
    provided by @ref BMI_GetOutputVarName can be retrieved.
    @param id Id number returned by @ref BMI_Create.
    @param var    Name of the variable to retrieve.
    @param dest   Variable in which to place results.
    @param l      Length of the string buffer @a dest.
    @retval       0 is success, 1 is failure; negative indicates buffer is too small.
    <p>
    The buffer length must be at least one character greater than the value 
    returned by @ref BMI_GetVarNbytes to allow for null termination. 
    "ErrorString" and "FilePrefix" return single strings.
    "Components" and "SelectedOutputHeadings" retrieve a string that is a 
    concatenated list of components or selected-output headings.
    The length of each item in a list is given by @ref BMI_GetVarItemsize. 
    The concatenated list must be processed to extract each component or heading 
    and a null termination must be appended.
    Alternatively, the components can be retrieved one at a time with 
    @a RM_GetComponent or @a RM_GetSelectedOutputHeading.
    </p>
    <p>
    Variable names for the second argument (@a var).
    </p>
    @n "Components"
    @n "ErrorString"
    @n "FilePrefix"
    @n "SelectedOutputHeadings".
    @see
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char string[256];
    status = BMI_GetValueChar(id, "FilePrefix", string);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueChar(int id, char* var, char* dest, int l);
	// GetValuePtr
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrInt(int id, char* var, int** dest);
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrDouble(int id, char* var, double** dest);
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrChar(int id, char* var, char** dest);
    /**
    @a BMI_GetValuePtr retrieves pointers to model variables. Only variables in the list
    provided by @ref BMI_GetPointableVarName can be pointed to.
    @param id Id number returned by @ref BMI_Create.
    @param var    Name of the variable to retrieve.
    @retval       Pointer to an up-to-date copy of the variable's data.
    <p>
    The following list gives the name in the second argument (@a var) and the
    data type the pointer:
    </p>
    @n "ComponentCount"
    @n "Concentrations"
    @n "DensityCalculated"
    @n "Gfw"
    @n "GridCellCount"
    @n "Porosity"
    @n "Pressure"
    @n "SaturationCalculated"
    @n "SelectedOutputOn"
    @n "SolutionVolume"
    @n "Temperature"
    @n "Time"
    @n "TimeStep"
    @n "Viscosity"
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT void* BMI_GetValuePtr(int id, char* var);
    /**
    @a BMI_GetVarGrid returns a value of 1, indicating points.
    BMIPhreeqcRM does not have a grid of its own. The cells
    of BMIPhreeqcRM are associated with the user's model grid,
    and all spatial characterists are assigned by the user's
    model.
    @param id Id number returned by @ref BMI_Create.
    @param var   Varaiable name. (Return value is the same regardless of value of @ var.)
    @retval      1 (points). BMIPhreeqcRM cells derive meaning from the user's model.
    */
    IRM_DLL_EXPORT int        BMI_GetVarGrid(int id, char* var);
    /**
    @a BMI_GetVarItemsize retrieves the size, in bytes, of a
    variable that can be set with
    @a BMI_SetValue methods, retrieved with @a BMI_GetValue methods, or pointed to with
    @ref BMI_GetValuePtr.
    Sizes may be the size of an integer, double,
    or a character length for string variables.
    @param id Id number returned by @ref BMI_Create.
    @param name        Name of the variable to retrieve the item size.
    @retval           Size, in bytes, of one element of the variable.
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetPointableVarName,
    @ref BMI_GetPointableItemCount,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    for(i = 0; i < GetInputVarCount(id); i++)
    {
        itemsize = GetVarItemsize(id, name);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetVarItemsize(int id, char* name);
    /**
    @a BMI_GetVarNbytes retrieves the total number of bytes needed for a
    variable that can be set with
    @a BMI_SetValue methods, retrieved with @a BMI_GetValue methods, or pointed to with
    @ref BMI_GetValuePtr.
    @param id Id number returned by @ref BMI_Create.
    @param name     Name of the variable to retrieve the number of bytes needed to
    retrieve or store the variable.
    @retval        Total number of bytes needed for the variable.
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetPointableVarName,
    @ref BMI_GetPointableItemCount,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    for(i = 0; i < GetInputVarCount(id); i++)
    {
        nbytes = GetVarNbytes(id, name);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT int        BMI_GetVarNbytes(int id, char* name);
    /**
    @a BMI_GetVarType retrieves the type of a variable that can be set with
    @a BMI_SetValue methods, retrieved with @a BMI_GetValue methods, or pointed to with
    @ref BMI_GetValuePtr.
    Types are "char", "double", or "int",
    or an array of these types.
    @param id Id number returned by @ref BMI_Create.
    @param name   Name of the variable to retrieve the type.
    @param vtype Type of the variable.
    @param l Length of string buffer @a vtype.
    @retval      0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetPointableVarName,
    @ref BMI_GetPointableItemCount,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char string[256];
    for(i = 0; i < GetInputVarCount(id); i++)
    {
        status = GetVarType(id, i, string, 256);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetVarType(int id, char* name, char* vtype, int l);
    /**
    @a BMI_GetVarType retrieves the units of a variable that can be set with
    @a BMI_SetValue methods, retrieved with @a BMI_GetValue methods, or pointed to with
    @ref BMI_GetValuePtr.
    @param id Id number returned by @ref BMI_Create.
    @param name   Name of the variable to retrieve the type.
    @param units Units of the variable.
    @param l Length of string buffer @a units.
    @retval      0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetOutputVarName,
    @ref BMI_GetOutputItemCount,
    @ref BMI_GetPointableVarName,
    @ref BMI_GetPointableItemCount,
    @ref BMI_GetValuePtr,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char string[256];
    for(i = 0; i < GetInputVarCount(id); i++)
    {
        status = GetVarUnits(id, i, string, 256);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_GetVarUnits(int id, char* name, char* units, int l);
	/**
	@a BMI_Initialize uses a YAML file to initialize an instance of BMIPhreeqcRM. 
    @param id Id number returned by @ref BMI_Create.
    @param config_file   String containing the YAML file name.
    @retval              0 is success, 1 is failure.
    <p>
    The file contains a YAML map of PhreeqcRM methods
    and the arguments corresponding to the methods.
    For example,
    </p>
    @htmlonly
    <CODE>
    <PRE>
    - key: LoadDatabase
      database: phreeqc.dat
    - key: RunFile
      workers: true
      initial_phreeqc: true
      utility: true
      chemistry_name: advect.pqi
    </PRE>
    </CODE>
    @endhtmlonly
    <p>
    @a BMI_Initialize will read the YAML file and execute the specified methods with
    the specified arguments. Using YAML
    terminology, the argument(s) for a method may be a scalar, a sequence, or a map,
    depending if the argument is
    a single item, a single vector, or there are multiple arguments.
    In the case of a map, the name associated
    with each argument (for example "chemistry_name" above) is arbitrary.
    The names of the map keys for map
    arguments are not used in parsing the YAML file; only the order of
    the arguments is important.
    </p>
    <p>
    The following list gives the PhreeqcRM methods that can be specified in a YAML file
    and the arguments that are required. The arguments are described with C++ formats, which
    are sufficient to identify which arguments are YAML scalars (single bool/logical,
    int, double, string/character argument),
    sequences (single vector argument), or maps (multiple arguments).
    </p>
    @htmlonly
    <CODE>
    <PRE>
    CloseFiles();
    CreateMapping(std::vector< int >& grid2chem);
    DumpModule();
    FindComponents();
    InitialEquilibriumPhases2Module(std::vector< int > equilibrium_phases);
    InitialExchanges2Module(std::vector< int > exchanges);
    InitialGasPhases2Module(std::vector< int > gas_phases);
    InitialKineticss2Module(std::vector< int > kinetics);
    InitialSolidSolutions2Module(std::vector< int > solid_solutions);
    InitialSolutions2Module(std::vector< int > solutions);
    InitialSurfaces2Module(std::vector< int > surfaces);
    InitialPhreeqc2Module(std::vector< int > initial_conditions1);
    InitialPhreeqc2Module(std::vector< int > initial_conditions1, 
    std::vector< int > initial_conditions2, std::vector< double > fraction1);
    InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
    LoadDatabase(std::string database);
    OpenFiles();
    OutputMessage(std::string str);
    RunCells();
    RunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
    RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    ScreenMessage(std::string str);
    SetComponentH2O(bool tf);
    SetConcentrations(std::vector< double > c);
    SetCurrentSelectedOutputUserNumber(int n_user);
    SetDensityUser(std::vector< double > density);
    SetDumpFileName(std::string dump_name);
    SetErrorHandlerMode(int mode);
    SetErrorOn(bool tf);
    SetFilePrefix(std::string prefix);
    SetGasCompMoles(std::vector< double > gas_moles);
    SetGasPhaseVolume(std::vector< double > gas_volume);
    SetPartitionUZSolids(bool tf);
    SetPorosity(std::vector< double > por);
    SetPressure(std::vector< double > p);
    SetPrintChemistryMask(std::vector< int > cell_mask);
    SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
    SetRebalanceByCell(bool tf);
    SetRebalanceFraction(double f);
    SetRepresentativeVolume(std::vector< double > rv);
    SetSaturationUser(std::vector< double > sat);
    SetScreenOn(bool tf);
    SetSelectedOutputOn(bool tf);
    SetSpeciesSaveOn(bool save_on);
    SetTemperature(std::vector< double > t);
    SetTime(double time);
    SetTimeConversion(double conv_factor);
    SetTimeStep(double time_step);
    SetUnitsExchange(int option);
    SetUnitsGasPhase(int option);
    SetUnitsKinetics(int option);
    SetUnitsPPassemblage(int option);
    SetUnitsSolution(int option);
    SetUnitsSSassemblage(int option);
    SetUnitsSurface(int option);
    SpeciesConcentrations2Module(std::vector< double > species_conc);
    StateSave(int istate);
    StateApply(int istate);
    StateDelete(int istate);
    UseSolutionDensityVolume(bool tf);
    WarningMessage(std::string warnstr);
    </PRE>
    </CODE>
    @endhtmlonly
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    id = BMI_Create(nxyz, nthreads);
    status = BMI_InitializeYAML(id, "myfile.yaml");
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_Initialize(int id, char* config_file);
    //
	// SetValue methods
    //
    /**
    @a BMI_SetValueChar sets model character variables. Only variables in the list
    provided by @ref BMI_GetInputVarName can be set.
    @param id Id number returned by @ref BMI_Create.
    @param name    Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "FilePrefix"
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_SetValueChar(id, "FilePrefix", "my_prefix");
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueChar(int id, char* name, const char* src);
    /**
    @a BMI_SetValueDouble sets model double variables. Only variables in the list
    provided by @ref BMI_GetInputVarName can be set.
    @param id Id number returned by @ref BMI_Create.
    @param name    Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "Time"
    @n "TimeStep".
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_SetValueDouble(id, "TimeStep", 86400.0);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueDouble(int id, char* name, double src);
    /**
    @a BMI_SetValueDoubleArray sets model double array variables. Only variables in the list
    provided by @ref BMI_GetInputVarName can be set.
    @param id Id number returned by @ref BMI_Create.
    @param name   Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "Concentrations"
    @n "DensityUser"
    @n "Porosity"
    @n "Pressure"
    @n "SaturationUser"
    @n "Temperature".
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    tc = (double *)malloc(nxyz*sizeof(double));
    for(i=0; i < nxyz; i++) tc[i] = 28.0e0;
    status = BMI_SetValueDoubleArray(id, "Temperature", tc);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueDoubleArray(int id, char* name, double* src);
    /**
    @a BMI_SetValueInt sets model int variables. Only variables in the list
    provided by @ref BMI_GetInputVarName can be set.
    @param id Id number returned by @ref BMI_Create.
    @param name   Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "NthSelectedOutput"
    @n "SelectedOutputOn".
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_SetValueInt(id, "SelectedOutputOn", 1);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueInt(int id, char* name, int src);
#ifdef SKIP
    /**
    @a BMI_SetValueIntArray sets model int array variables. Only variables in the list
    provided by @ref BMI_GetInputVarName can be set.
    @param id Id number returned by @ref BMI_Create.
    @param name   Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "Concentrations"
    @n "DensityUser"
    @n "FilePrefix"
    @n "NthSelectedOutput"
    @n "Porosity"
    @n "Pressure"
    @n "SaturationUser"
    @n "SelectedOutputOn"
    @n "Temperature"
    @n "Time"
    @n "TimeStep".
    @see
    @ref BMI_GetInputVarName,
    @ref BMI_GetInputItemCount,
    @ref BMI_GetVarItemsize,
    @ref BMI_GetVarNbytes,
    @ref BMI_GetVarType,
    @ref BMI_GetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    tc = (double *)malloc(nxyz*sizeof(double));
    for(i=0; i < nxyz; i++) tc[i] = 28.0e0;
    status = BMI_SetValue(id, "Temperature", tc);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueIntArray(int id, char* name, int* src);
#endif
    /**
    @a BMI_Update runs a reaction step for all of the cells in the reaction module.
    @param id Id number returned by @ref BMI_Create.
    @retval     0 is success, 1 is failure.
    <p>
    Tranported concentrations are transferred to the reaction cells
    (@ref BMI_SetValueDoubleArray "Concentrations") before
    reaction calculations are run. The length of time over which kinetic
    reactions are integrated is set
    by @ref BMI_SetValueDouble "TimeStep". Other properties that may need to be updated
    as a result of the transport
    calculations include 
    porosity,
    pressure,
    saturation,
    temperature.
    </p>
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_SetValue(id, "Porosity", por);                ! If pore volume changes
    status = BMI_SetValue(id, "SaturationUser", sat);          ! If saturation changes
    status = BMI_SetValue(id, "Temperature", temperature);     ! If temperature changes
    status = BMI_SetValue(id, "Pressure", pressure);           ! If pressure changes
    status = BMI_SetValue(id, "Concentrations", c);            ! Transported concentrations
    status = BMI_SetValue(id, "TimeStep", time_step);          ! Time step for kinetic reactions
    status = BMI_Update(id);
    status = BMI_GetValue(id, "Concentrations", c);            ! Concentrations after reaction
    status = BMI_GetValue(id, "DensityCalculated", density);   ! Density after reaction
    status = BMI_GetValue(id, "SolutionVolume", volume);       ! Solution volume after reaction
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_Update(int id);
    /**
    @a BMI_UpdateUntil is the same as @ref BMI_Update, except the time step is calculated
    from the argument @a end_time. The time step is calculated to be @a end_time minus
    the current time (@ref BMI_GetCurrentTime).
    @param id Id number returned by @ref BMI_Create..
    @param end_time Time at the end of the time step.
    @see
    @ref BMI_Initialize,
    @ref BMI_Update.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = BMI_SetValue(id, "Time", time);
    status = BMI_SetValue(id, "Concentrations", c);
    status = BMI_UpdateUntil(id, time + 86400.0);
    status = BMI_GetValue(id, "Concentrations", c);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
	IRM_DLL_EXPORT IRM_RESULT BMI_UpdateUntil(int id, double end_time);
    /**
    @a BMI_GetValueAtIndices is not implemented
    */
    IRM_DLL_EXPORT void BMI_GetValueAtIndices(int id, char* name, void* dest, int* inds, int count);
    /**
    @a BMI_SetValueAtIndices is not implemented.
    */
    IRM_DLL_EXPORT void BMI_SetValueAtIndices(int id, char* name, int* inds, int count, void* src);
    /**
    @a BMI_GetGridShape is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridShape(int id, const int grid, int* shape);
    /**
    @a BMI_GetGridSpacing is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridSpacing(int id, const int grid, double* spacing);
    /**
    @a BMI_GetGridOrigin is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridOrigin(int id, const int grid, double* origin);
    /**
    @a BMI_GetGridX is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridX(int id, const int grid, double* x);
    /**
    @a BMI_GetGridY is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridY(int id, const int grid, double* y);
    /**
    @a BMI_GetGridZ is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridZ(int id, const int grid, double* z);
    /**
    @a BMI_GetGridNodeCount is not implemented.
    */
    IRM_DLL_EXPORT int BMI_GetGridNodeCount(int id, const int grid);
    /**
    @a BMI_GetGridEdgeCount is not implemented.
    */
    IRM_DLL_EXPORT int BMI_GetGridEdgeCount(int id, const int grid);
    /**
    @a BMI_GetGridFaceCount is not implemented.
    */
    IRM_DLL_EXPORT int BMI_GetGridFaceCount(int id, const int grid);
    /**
    @a BMI_GetGridEdgeNodes is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridEdgeNodes(int id, const int grid, int* edge_nodes);
    /**
    @a BMI_GetGridFaceEdges is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridFaceEdges(int id, const int grid, int* face_edges);
    /**
    @a BMI_GetGridFaceNodes is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridFaceNodes(int id, const int grid, int* face_nodes);
    /**
    @a BMI_GetGridNodesPerFace is not implemented.
    */
    IRM_DLL_EXPORT void BMI_GetGridNodesPerFace(int id, const int grid, int* nodes_per_face);
#endif // SKIP
#if defined(__cplusplus)
}

#endif // BMI_INTERFACE_C_H_INCLUDED

#endif // SKIP