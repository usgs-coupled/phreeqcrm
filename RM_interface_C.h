/*! @file RM_interface_C.h
	@brief C header file for module BMIPhreeqcRM. RM_BmiCreate creates a module with 
    access to all BMI methods. For backward compatibility, the deprecated method RM_Create
    creates an old-style PhreeqcRM instance, which does not have access to BMI methods.
*/
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "IrmResult.h"
#ifndef RM_INTERFACE_C_H
#define RM_INTERFACE_C_H

#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif
    /**
    @a RM_BmiCreate Creates a BMI reaction module, which allows use of all of
    the RM_Bmi methods. 
    If the code is compiled with
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
    @ref RM_BmiFinalize.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    nxyz = 40;
    nthreads = 3;
    id = RM_BmiCreate(nxyz, nthreads);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root and workers.
    */
    IRM_DLL_EXPORT int        RM_BmiCreate(int nxyz, int nthreads);
    /**
    @a RM_BmiDestroy Destroys a BMI reaction module; same as @ref RM_BmiFinalize.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval    0 is success, 0 is failure.
    @see
    @ref RM_BmiCreate.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiDestroy(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root and workers.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiDestroy(int id);
    /**
    @a RM_BmiAddOutputVars allows selection of sets of variables that can be retieved
    by the @a RM_BmiGetValue methods. Sets of variables can be included or excluded with
    multiple calls to this method. All calls must precede the final call to
    the PhreeqcRM method FindComponents. FindComponents generates SELECTED_OUTPUT 333 and
    USER_PUNCH 333 data blocks that make the variables accessible. Variables will
    only be accessible if the system includes the given reactant; for example, no
    gas variables will be Created if there are no GAS_PHASEs in the model.
    @param id Id number returned by @ref RM_BmiCreate.
    @param option A string value, among those listed below, that includes or
    excludes variables from @ref RM_BmiGetOutputVarName, @a RM_BmiGetValue methods,
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
    status = RM_BmiAddOutputVars(id, "SolutionMolalities", "True");
    status = RM_BmiAddOutputVars(id, "SaturationIndices", "Calcite Dolomite");
    </PRE>
    </CODE>
    @endhtmlonly
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiAddOutputVars(int id, char* option, char* def);
    /**
    @a RM_BmiFinalize Destroys a reaction module.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval    0 is success, 0 is failure.
    @see
    @ref RM_BmiCreate.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiFinalize(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root and workers.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiFinalize(int id);
    /**
    @a RM_BmiGetComponentName returns the component name--"BMIPhreeqcRM".
    @param id Id number returned by @ref RM_BmiCreate.
    @param component_name Returns "BMIPhreeqcRM", the name of the component.
    @param l Length of string buffer @a component_name.
    @retval               0 is success, 1 is failure; negative indicates buffer is too small.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiGetComponentName(id, component_name);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetComponentName(int id, char* component_name, int l);
    /**
    @a RM_BmiGetCurrentTime returns the current simulation time, in seconds.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval        The current simulation time, in seconds.
    @see
    @ref RM_BmiGetEndTime,
    @ref RM_BmiGetTimeStep,
    @ref RM_BmiGetTime.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    now = RM_BmiGetCurrentTime(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT double     RM_BmiGetCurrentTime(int id);
    /**
    @a RM_BmiGetEndTime returns @ref RM_BmiGetCurrentTime plus
    @ref RM_BmiGetTimeStep, in seconds.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval          The end of the time step, in seconds.
    @see
    @ref RM_BmiGetCurrentTime,
    @ref RM_BmiGetTimeStep.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiGetEndTime(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT double     RM_BmiGetEndTime(int id);
    /**
    @a RM_BmiGetGridRank returns a rank of 1 for grid 0.
    BMIPhreeqcRM has a 1D series of
    cells; any grid or spatial information must
    be found in the user's model.
    @param id Id number returned by @ref RM_BmiCreate.
    @param grid   Grid number, only grid 0 is considered.
    @retval       Rank of 1 is returned for grid 0; 0 for
    all other values of @a grid.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    rank = RM_BmiGetGridRank(id, grid)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT int        RM_BmiGetGridRank(int id, int grid);
    /**
    @ref RM_BmiGetGridSize returns the number of cells specified
    at creation of the BMIPhreeqcRM instance.
    @param id Id number returned by @ref RM_BmiCreate.
    @param grid  Grid number, only grid 0 is considered.
    @retval    Number of cells. Same value as @ref RM_BmiGetValueInt(id, "GridCellCount")
    is returned for grid 0;
    0 for all other values of @a grid.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    nxyz = RM_BmiGetGridSize(id, grid);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT int        RM_BmiGetGridSize(int id, int grid);
    /**
    @a RM_BmiGetGridType defines the grid to be points. No grid
    information is available in BMIPhreeqcRM; all grid
    information must be found in the user's model.
    @param id    Id number returned by @ref RM_BmiCreate.
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
    status = RM_BmiGetGridType(id, grid, str, l)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetGridType(int id, int grid, char* str, int l);
    /**
    @a RM_BmiGetInputItemCount returns count of variables that
    can be set with @a RM_BmiSetValue methods.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval   Number of input variables that can be set with @a RM_BmiSetValue methods.
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = RM_BmiGetInputItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT int        RM_BmiGetInputItemCount(int id);
    //IRM_DLL_EXPORT int        RM_BmiGetInputVarNamesSize(int id);
    /**
    @a RM_BmiGetInputVarName returns the ith variable name that can be set
    with @a RM_BmiSetValue methods.
    @param id Id number returned by @ref RM_BmiCreate.
    @param i 0-based index of variable name to retrieve.
    @param name   Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256];
    status = RM_BmiGetInputVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetInputVarName(int id, int i, char* name, int l);
    //IRM_DLL_EXPORT int        RM_BmiGetNames(int id, const char* type, char* dest);
    //IRM_DLL_EXPORT int        RM_BmiGetNamesSize(int id, const char* type, int* dest);
    /**
    @a RM_BmiGetOutputItemCount returns count of output variables that can be
    retrieved with @a RM_BmiGetValue methods.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval   Number of output variables that can be retrieved with @a RM_BmiGetValue methods.
    @see
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = RM_BmiGetOutputItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT int        RM_BmiGetOutputItemCount(int id);
    //IRM_DLL_EXPORT int        RM_BmiGetOutputVarNamesSize(int id);
    /**
    @a RM_BmiGetOutputVarName returns ith variable name for which a pointer can be
    retrieved with @a RM_BmiGetValue methods.
    @param id Id number returned by @ref RM_BmiCreate.
    @param i 0-based index of variable name to retrieve.
    @param name    Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256]
    status = RM_BmiGetOutputVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetOutputVarName(int id, int i, char* name, int l);
    /**
    @a RM_BmiGetPointableItemCount returns count of pointable variables that can be
    retrieved with @ref RM_BmiGetValuePtr.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval   Number of output variables that can be retrieved with @ref RM_BmiGetValuePtr.
    @see
    @ref RM_BmiGetPointableVarName,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    count = RM_BmiGetPointableItemCount(id);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT int		  RM_BmiGetPointableItemCount(int id);
    //IRM_DLL_EXPORT int        RM_BmiGetPointableVarNamesSize(int id);
    /**
    @a RM_BmiGetPointableVarName returns ith variable name for which a pointer can be
    retrieved with @ref RM_BmiGetValuePtr.
    @param id Id number returned by @ref RM_BmiCreate.
    @param i 0-based index of variable name to retrieve.
    @param name    Retrieved variable name.
    @param l Length of buffer for @a name.
    @retval            0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetPointableItemCount,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char name[256];
    status = RM_BmiGetPointableVarName(id, 0, name, 256);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetPointableVarName(int id, int i, char* name, int l);
    /**
    @a RM_BmiGetStartTime returns the current simulation time, in seconds.
    (Same as @ref RM_BmiGetCurrentTime.)
    @param id Id number returned by @ref RM_BmiCreate.
    @retval       The current simulation time, in seconds.
    */
    IRM_DLL_EXPORT double     RM_BmiGetStartTime(int id);
    /**
    @a RM_BmiGetTime returns the current simulation time, in seconds.
    (Same as @ref RM_BmiGetCurrentTime.)
    @param id Id number returned by @ref RM_BmiCreate.
    @retval       The current simulation time, in seconds.
    */
    IRM_DLL_EXPORT double     RM_BmiGetTime(int id);
    /**
    @a RM_BmiGetTimeStep returns the current simulation time step,
    in seconds.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval              The current simulation time step, in seconds.
    @see
    @ref RM_BmiGetCurrentTime,
    @ref RM_BmiGetEndTime.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    time_step = RM_BmiGetTimeStep(id)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT double     RM_BmiGetTimeStep(int id);
    /**
    @a RM_BmiGetTimeUnits returns the time units of PhreeqcRM.
    All time units are seconds for PhreeqcRM.
    @param id Id number returned by @ref RM_BmiCreate.
    @param units    Returns the string "seconds".
    @param l    Length of the string buffer @a units.
    @retval              0 is success, 1 failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetCurrentTime,
    @ref RM_BmiGetEndTime,
    @ref RM_BmiGetTimeStep.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char time_units[256]
    status = RM_BmiGetTimeUnits(id, time_units, 256)
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetTimeUnits(int id, char* units, int l);
    // GetValue
    /**
    @a RM_BmiGetValueInt retrieves int model variables. Only variables in the list
    provided by @ref RM_BmiGetOutputVarName can be retrieved.
    @param id Id number returned by @ref RM_BmiCreate.
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
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiGetValueInt(id, "ComponentCount", &count);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValueInt(int id, char* var, int* dest);
    /**
    @a RM_BmiGetValueDouble retrieves model variables. Only variables in the list
    provided by @ref RM_BmiGetOutputVarName can be retrieved.
    @param id Id number returned by @ref RM_BmiCreate.
    @param var    Name of the variable to retrieve.
    @param dest   Variable in which to place results.
    @retval       0 is success, 1 is failure.
    <p>
    Variables in addition to the ones listed below may be retrieved by this method,
    depending on variables selected by @ref RM_BmiAddOutputVars. All variables added
    by @ref RM_BmiAddOutputVars will be double arrays of size equal to the number of
    model cells [@ref RM_BmiGetValueInt(id, "GridCellCount")].
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
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    density = (double *)malloc(nxyz*sizeof(double));
    status = RM_BmiGetValueDouble(id, "DensityCalculated", density);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValueDouble(int id, char* var, double* dest);
    /**
    @a RM_BmiGetValueChar retrieves char model variables. Only variables in the list
    provided by @ref RM_BmiGetOutputVarName can be retrieved.
    @param id Id number returned by @ref RM_BmiCreate.
    @param var    Name of the variable to retrieve.
    @param dest   Variable in which to place results.
    @param l      Length of the string buffer @a dest.
    @retval       0 is success, 1 is failure; negative indicates buffer is too small.
    <p>
    The buffer length must be at least one character greater than the value
    returned by @ref RM_BmiGetVarNbytes to allow for null termination.
    "ErrorString" and "FilePrefix" return single strings.
    "Components" and "SelectedOutputHeadings" retrieve a string that is a
    concatenated list of components or selected-output headings.
    The length of each item in a list is given by @ref RM_BmiGetVarItemsize.
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
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    char string[256];
    status = RM_BmiGetValueChar(id, "FilePrefix", string);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValueChar(int id, char* var, char* dest, int l);
    // GetValuePtr
    //IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValuePtrInt(int id, char* var, int** dest);
    //IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValuePtrDouble(int id, char* var, double** dest);
    //IRM_DLL_EXPORT IRM_RESULT RM_BmiGetValuePtrChar(int id, char* var, char** dest);
    /**
    @a RM_BmiGetValuePtr retrieves pointers to model variables. Only variables in the list
    provided by @ref RM_BmiGetPointableVarName can be pointed to.
    @param id Id number returned by @ref RM_BmiCreate.
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
    IRM_DLL_EXPORT void* RM_BmiGetValuePtr(int id, char* var);
    /**
    @a RM_BmiGetVarGrid returns a value of 1, indicating points.
    BMIPhreeqcRM does not have a grid of its own. The cells
    of BMIPhreeqcRM are associated with the user's model grid,
    and all spatial characterists are assigned by the user's
    model.
    @param id Id number returned by @ref RM_BmiCreate.
    @param var   Varaiable name. (Return value is the same regardless of value of @ var.)
    @retval      1 (points). BMIPhreeqcRM cells derive meaning from the user's model.
    */
    IRM_DLL_EXPORT int        RM_BmiGetVarGrid(int id, char* var);
    /**
    @a RM_BmiGetVarItemsize retrieves the size, in bytes, of a
    variable that can be set with
    @a RM_BmiSetValue methods, retrieved with @a RM_BmiGetValue methods, or pointed to with
    @ref RM_BmiGetValuePtr.
    Sizes may be the size of an integer, double,
    or a character length for string variables.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name        Name of the variable to retrieve the item size.
    @retval           Size, in bytes, of one element of the variable.
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetPointableVarName,
    @ref RM_BmiGetPointableItemCount,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
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
    IRM_DLL_EXPORT int        RM_BmiGetVarItemsize(int id, char* name);
    /**
    @a RM_BmiGetVarNbytes retrieves the total number of bytes needed for a
    variable that can be set with
    @a RM_BmiSetValue methods, retrieved with @a RM_BmiGetValue methods, or pointed to with
    @ref RM_BmiGetValuePtr.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name     Name of the variable to retrieve the number of bytes needed to
    retrieve or store the variable.
    @retval        Total number of bytes needed for the variable.
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetPointableVarName,
    @ref RM_BmiGetPointableItemCount,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
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
    IRM_DLL_EXPORT int        RM_BmiGetVarNbytes(int id, char* name);
    /**
    @a RM_BmiGetVarType retrieves the type of a variable that can be set with
    @a RM_BmiSetValue methods, retrieved with @a RM_BmiGetValue methods, or pointed to with
    @ref RM_BmiGetValuePtr.
    Types are "char", "double", or "int",
    or an array of these types.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name   Name of the variable to retrieve the type.
    @param vtype Type of the variable.
    @param l Length of string buffer @a vtype.
    @retval      0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetPointableVarName,
    @ref RM_BmiGetPointableItemCount,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarUnits.
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
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetVarType(int id, char* name, char* vtype, int l);
    /**
    @a RM_BmiGetVarType retrieves the units of a variable that can be set with
    @a RM_BmiSetValue methods, retrieved with @a RM_BmiGetValue methods, or pointed to with
    @ref RM_BmiGetValuePtr.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name   Name of the variable to retrieve the type.
    @param units Units of the variable.
    @param l Length of string buffer @a units.
    @retval      0 is success, 1 is failure; negative indicates buffer is too small.
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetOutputVarName,
    @ref RM_BmiGetOutputItemCount,
    @ref RM_BmiGetPointableVarName,
    @ref RM_BmiGetPointableItemCount,
    @ref RM_BmiGetValuePtr,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType.
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
    IRM_DLL_EXPORT IRM_RESULT RM_BmiGetVarUnits(int id, char* name, char* units, int l);
    /**
    @a RM_BmiInitialize uses a YAML file to initialize an instance of BMIPhreeqcRM.
    @param id Id number returned by @ref RM_BmiCreate.
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
    @a RM_BmiInitialize will read the YAML file and execute the specified methods with
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
    id = RM_BmiCreate(nxyz, nthreads);
    status = RM_BmiInitializeYAML(id, "myfile.yaml");
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiInitialize(int id, char* config_file);
    //
    // SetValue methods
    //
    /**
    @a RM_BmiSetValueChar sets model character variables. Only variables in the list
    provided by @ref RM_BmiGetInputVarName can be set.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name    Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "FilePrefix"
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiSetValueChar(id, "FilePrefix", "my_prefix");
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiSetValueChar(int id, char* name, const char* src);
    /**
    @a RM_BmiSetValueDouble sets model double variables. Only variables in the list
    provided by @ref RM_BmiGetInputVarName can be set.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name    Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "Time"
    @n "TimeStep".
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiSetValueDouble(id, "TimeStep", 86400.0);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiSetValueDouble(int id, char* name, double src);
    /**
    @a RM_BmiSetValueDoubleArray sets model double array variables. Only variables in the list
    provided by @ref RM_BmiGetInputVarName can be set.
    @param id Id number returned by @ref RM_BmiCreate.
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
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    tc = (double *)malloc(nxyz*sizeof(double));
    for(i=0; i < nxyz; i++) tc[i] = 28.0e0;
    status = RM_BmiSetValueDoubleArray(id, "Temperature", tc);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiSetValueDoubleArray(int id, char* name, double* src);
    /**
    @a RM_BmiSetValueInt sets model int variables. Only variables in the list
    provided by @ref RM_BmiGetInputVarName can be set.
    @param id Id number returned by @ref RM_BmiCreate.
    @param name   Name of variable to set.
    @param src    Data to use to set the variable.
    @retval       0 is success, 1 is failure.
    <p>
    Variable names for the second argument (@a var):
    </p>
    @n "NthSelectedOutput"
    @n "SelectedOutputOn".
    @see
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiSetValueInt(id, "SelectedOutputOn", 1);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiSetValueInt(int id, char* name, int src);
#ifdef SKIP
    /**
    @a RM_BmiSetValueIntArray sets model int array variables. Only variables in the list
    provided by @ref RM_BmiGetInputVarName can be set.
    @param id Id number returned by @ref RM_BmiCreate.
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
    @ref RM_BmiGetInputVarName,
    @ref RM_BmiGetInputItemCount,
    @ref RM_BmiGetVarItemsize,
    @ref RM_BmiGetVarNbytes,
    @ref RM_BmiGetVarType,
    @ref RM_BmiGetVarUnits.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    tc = (double *)malloc(nxyz*sizeof(double));
    for(i=0; i < nxyz; i++) tc[i] = 28.0e0;
    status = RM_BmiSetValue(id, "Temperature", tc);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiSetValueIntArray(int id, char* name, int* src);
#endif
    /**
    @a RM_BmiUpdate runs a reaction step for all of the cells in the reaction module.
    @param id Id number returned by @ref RM_BmiCreate.
    @retval     0 is success, 1 is failure.
    <p>
    Tranported concentrations are transferred to the reaction cells
    (@ref RM_BmiSetValueDoubleArray "Concentrations") before
    reaction calculations are run. The length of time over which kinetic
    reactions are integrated is set
    by @ref RM_BmiSetValueDouble "TimeStep". Other properties that may need to be updated
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
    status = RM_BmiSetValue(id, "Porosity", por);                ! If pore volume changes
    status = RM_BmiSetValue(id, "SaturationUser", sat);          ! If saturation changes
    status = RM_BmiSetValue(id, "Temperature", temperature);     ! If temperature changes
    status = RM_BmiSetValue(id, "Pressure", pressure);           ! If pressure changes
    status = RM_BmiSetValue(id, "Concentrations", c);            ! Transported concentrations
    status = RM_BmiSetValue(id, "TimeStep", time_step);          ! Time step for kinetic reactions
    status = RM_BmiUpdate(id);
    status = RM_BmiGetValue(id, "Concentrations", c);            ! Concentrations after reaction
    status = RM_BmiGetValue(id, "DensityCalculated", density);   ! Density after reaction
    status = RM_BmiGetValue(id, "SolutionVolume", volume);       ! Solution volume after reaction
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiUpdate(int id);
    /**
    @a RM_BmiUpdateUntil is the same as @ref RM_BmiUpdate, except the time step is calculated
    from the argument @a end_time. The time step is calculated to be @a end_time minus
    the current time (@ref RM_BmiGetCurrentTime).
    @param id Id number returned by @ref RM_BmiCreate..
    @param end_time Time at the end of the time step.
    @see
    @ref RM_BmiInitialize,
    @ref RM_BmiUpdate.
    @par C example:
    @htmlonly
    <CODE>
    <PRE>
    status = RM_BmiSetValue(id, "Time", time);
    status = RM_BmiSetValue(id, "Concentrations", c);
    status = RM_BmiUpdateUntil(id, time + 86400.0);
    status = RM_BmiGetValue(id, "Concentrations", c);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @a RM_MpiWorker.
    */
    IRM_DLL_EXPORT IRM_RESULT RM_BmiUpdateUntil(int id, double end_time);
    /**
    @a RM_BmiGetValueAtIndices is not implemented
    */
    IRM_DLL_EXPORT void RM_BmiGetValueAtIndices(int id, char* name, void* dest, int* inds, int count);
    /**
    @a RM_BmiSetValueAtIndices is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiSetValueAtIndices(int id, char* name, int* inds, int count, void* src);
    /**
    @a RM_BmiGetGridShape is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridShape(int id, const int grid, int* shape);
    /**
    @a RM_BmiGetGridSpacing is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridSpacing(int id, const int grid, double* spacing);
    /**
    @a RM_BmiGetGridOrigin is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridOrigin(int id, const int grid, double* origin);
    /**
    @a RM_BmiGetGridX is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridX(int id, const int grid, double* x);
    /**
    @a RM_BmiGetGridY is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridY(int id, const int grid, double* y);
    /**
    @a RM_BmiGetGridZ is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridZ(int id, const int grid, double* z);
    /**
    @a RM_BmiGetGridNodeCount is not implemented.
    */
    IRM_DLL_EXPORT int RM_BmiGetGridNodeCount(int id, const int grid);
    /**
    @a RM_BmiGetGridEdgeCount is not implemented.
    */
    IRM_DLL_EXPORT int RM_BmiGetGridEdgeCount(int id, const int grid);
    /**
    @a RM_BmiGetGridFaceCount is not implemented.
    */
    IRM_DLL_EXPORT int RM_BmiGetGridFaceCount(int id, const int grid);
    /**
    @a RM_BmiGetGridEdgeNodes is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridEdgeNodes(int id, const int grid, int* edge_nodes);
    /**
    @a RM_BmiGetGridFaceEdges is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridFaceEdges(int id, const int grid, int* face_edges);
    /**
    @a RM_BmiGetGridFaceNodes is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridFaceNodes(int id, const int grid, int* face_nodes);
    /**
    @a RM_BmiGetGridNodesPerFace is not implemented.
    */
    IRM_DLL_EXPORT void RM_BmiGetGridNodesPerFace(int id, const int grid, int* nodes_per_face);




/**
Abort the program. @a Result will be interpreted as
an IRM_RESULT value and decoded; @a err_str will be printed; and the reaction module
will be destroyed. If using MPI, an MPI_Abort message will be sent before the reaction
module is destroyed. If the @a id is an invalid instance, RM_Abort will return a value of
IRM_BADINSTANCE, otherwise the program will exit with a return code of 4.
@param id            The instance id returned from @ref RM_Create.
@param result        Integer treated as an IRM_RESULT return code.
@param err_str       String to be printed as an error message.
@retval IRM_RESULT   Program will exit before returning unless @a id is an invalid reaction module id.
@see                 
@ref RM_Destroy, 
@ref RM_ErrorMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm);
Utilities::strcpy_safe(str, MAX_LENGTH, "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1");
status = RunString(iphreeqc_id, str);
if (status != 0) status = RM_Abort(id, status, "IPhreeqc RunString failed");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root or workers.
*/
IRM_DLL_EXPORT IRM_RESULT RM_Abort(int id, int result, const char * err_str);
/**
Close the output and log files.
@param id            The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
@see                 
@ref RM_OpenFiles, 
@ref RM_SetFilePrefix.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_CloseFiles(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called only by root.
 */
IRM_DLL_EXPORT IRM_RESULT        RM_CloseFiles(int id);
/**
@a N sets of component concentrations are converted to SOLUTIONs numbered 1-@a n in the Utility IPhreeqc.
The solutions can be reacted and manipulated with the methods of IPhreeqc. If solution concentration units
(@ref RM_SetUnitsSolution) are per liter, one liter of solution is created in the Utility instance; if solution
concentration units are mass fraction, one kilogram of solution is created in the Utility instance.
The motivation for this
method is the mixing of solutions in wells, where it may be necessary to calculate solution properties
(pH for example) or react the mixture to form scale minerals. The code fragments below make a mixture of
concentrations and then calculate the pH of the mixture.
@param id            The instance @a id returned from @ref RM_Create.
@param c             Array of concentrations to be made SOLUTIONs in Utility IPhreeqc. Array storage is 
n * ncomps.
@param n             The number of sets of concentrations.
@param tc            Array of temperatures to apply to the SOLUTIONs, in degree C. Array of size @a n.
@param p_atm         Array of pressures to apply to the SOLUTIONs, in atm. Array of size n.
@retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
c_well = (double *) malloc((size_t) ((size_t) (1 * ncomps * sizeof(double))));
for (i = 0; i < ncomps; i++)
{
  c_well[i] = 0.5 * c[0 + nxyz*i] + 0.5 * c[9 + nxyz*i];
}
tc = (double *) malloc((size_t) (1 * sizeof(double)));
p_atm = (double *) malloc((size_t) (1 * sizeof(double)));
tc[0] = 15.0;
p_atm[0] = 3.0;
iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm);
Utilities::strcpy_safe(str, MAX_LENGTH, "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1");
status = RunString(iphreeqc_id, str);
status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5);
status = GetSelectedOutputValue2(iphreeqc_id, 1, 0, &vtype, &pH, svalue, 100);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called only by root.
 */
IRM_DLL_EXPORT int        RM_Concentrations2Utility(int id, double *c, int n, double *tc, double *p_atm);
/**
Creates a reaction module without BMI methods. This method is <b>deprecated</b> and included only for
backward compatibility. Use RM_BmiCreate to create a reaction module with access to all
RM_Bmi methods.
<p>
If the code is compiled with
the preprocessor directive USE_OPENMP, the reaction module is multithreaded.
If the code is compiled with the preprocessor directive USE_MPI, the reaction
module will use MPI and multiple processes. If neither preprocessor directive is used,
the reaction module will be serial (unparallelized).
</p>
@param nxyz                   The number of grid cells in the user's model.
@param nthreads (or @a comm, MPI)       When using OPENMP, the argument (@a nthreads) is the number of worker threads to be used.
If @a nthreads <= 0, the number of threads is set equal to the number of processors of the computer.
When using MPI, the argument (@a comm) is the MPI communicator to use within the reaction module.
@retval Id of the PhreeqcRM instance, negative is failure (See @ref RM_DecodeError).
@see                 
@ref RM_Destroy.
@par C Example:
@htmlonly
<CODE>
<PRE>
nxyz = 40;
#ifdef USE_MPI
  id = RM_Create(nxyz, MPI_COMM_WORLD);
  if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
  {
    exit(4);
  }
  if (mpi_myself > 0)
  {
    status = RM_MpiWorker(id);
    status = RM_Destroy(id);
    return;
  }
#else
  nthreads = 3;
  id = RM_Create(nxyz, nthreads);
#endif
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and workers.
 */
#ifdef USE_MPI
IRM_DLL_EXPORT int RM_Create(int nxyz, MPI_Comm comm);
#else
IRM_DLL_EXPORT int RM_Create(int nxyz, int nthreads);
#endif
/**
Provides a mapping from grid cells in the user's model to reaction cells in PhreeqcRM.
The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells 
for which chemistry must be run.
The array @a grid2chem of size @a nxyz (the number of grid cells, @ref RM_GetGridCellCount) 
must contain the set of all integers 0 <= @a i < @a count_chemistry, 
where @a count_chemistry is a number less than or equal to @a nxyz.
Inactive cells are assigned a negative integer. 
The mapping may be many-to-one to account for symmetry.
Default is a one-to-one mapping--all user grid cells are reaction cells (equivalent to @a grid2chem values of 0,1,2,3,...,@a nxyz-1).
@param id               The instance @a id returned from @ref RM_Create.
@param grid2chem        An array of integers: Nonnegative is a reaction cell number (0 based), negative is an inactive cell. Array of size @a nxyz (number of grid cells).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
// For demonstation, two equivalent rows by symmetry
grid2chem = (int *) malloc((size_t) (nxyz * sizeof(int)));
for (i = 0; i < nxyz/2; i++)
{
  grid2chem[i] = i;
  grid2chem[i+nxyz/2] = i;
}
status = RM_CreateMapping(id, grid2chem);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_CreateMapping (int id, int *grid2chem);
/**
If @a e is negative, this method prints an error message corresponding to IRM_RESULT @a e. If @a e is non-negative, no action is taken.
@param id                   The instance @a id returned from @ref RM_Create.
@param e                    An IRM_RESULT value returned by one of the reaction-module methods.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
@par IRM_RESULT definition:
@htmlonly
<CODE>
<PRE>
typedef enum {
  IRM_OK            =  0,  //Success
  IRM_OUTOFMEMORY   = -1,  //Failure, Out of memory
  IRM_BADVARTYPE    = -2,  //Failure, Invalid VAR type
  IRM_INVALIDARG    = -3,  //Failure, Invalid argument
  IRM_INVALIDROW    = -4,  //Failure, Invalid row
  IRM_INVALIDCOL    = -5,  //Failure, Invalid column
  IRM_BADINSTANCE   = -6,  //Failure, Invalid rm instance id
  IRM_FAIL          = -7,  //Failure, Unspecified
} IRM_RESULT;
</PRE>
</CODE>
@endhtmlonly
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_CreateMapping(id, grid2chem);
if (status < 0) status = RM_DecodeError(id, status);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Can be called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_DecodeError (int id, int e);
/**
Destroys a reaction module.
@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_Destroy(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_Destroy(int id);
/**
Writes the contents of all workers to file in _RAW formats, including SOLUTIONs and all reactants.
@param id               The instance @a id returned from @ref RM_Create.
@param dump_on          Signal for writing the dump file: 1 true, 0 false.
@param append           Signal to append to the contents of the dump file: 1 true, 0 false.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetDumpFileName.
@par C Example:
@htmlonly
<CODE>
<PRE>
dump_on = 1;
append = 0;
status = RM_SetDumpFileName(id, "advection_c.dmp");
status = RM_DumpModule(id, dump_on, append);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root; workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_DumpModule(int id, int dump_on, int append);
/**
Send an error message to the screen, the output file, and the log file.
@param id               The instance @a id returned from @ref RM_Create.
@param errstr           String to be printed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see               
@ref RM_LogMessage,     
@ref RM_OpenFiles,  
@ref RM_OutputMessage, 
@ref RM_ScreenMessage, 
@ref RM_WarningMessage.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_ErrorMessage(id, "Goodbye world");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; root writes to output and log files.
 */
IRM_DLL_EXPORT IRM_RESULT RM_ErrorMessage(int id, const char *errstr);
/**
Returns the number of items in the list of all elements in the InitialPhreeqc instance.
Elements are those that have been defined in a solution or any other reactant (EQUILIBRIUM_PHASE, KINETICS, and others).
The method can be called multiple times and the list that is created is cummulative.
The list is the set of components that needs to be transported.
By default the list
includes water, excess H and excess O (the H and O not contained in water);
alternatively, the list may be set to contain total H and total O (@ref RM_SetComponentH2O),
which requires transport results to be accurate to eight or nine significant digits.
If multicomponent diffusion (MCD) is to be modeled, there is a capability to retrieve aqueous species concentrations
(@ref RM_GetSpeciesConcentrations) and to set new solution concentrations after MCD by using individual species concentrations
(@ref RM_SpeciesConcentrations2Module). To use these methods the save-species property needs to be turned on (@ref RM_SetSpeciesSaveOn).
If the save-species property is on, RM_FindComponents will generate
a list of aqueous species (@ref RM_GetSpeciesCount, @ref RM_GetSpeciesName), their diffusion coefficients at 25 C (@ref RM_GetSpeciesD25),
their charge (@ref RM_GetSpeciesZ).
@param id            The instance @a id returned from @ref RM_Create.
@retval              Number of components currently in the list, or IRM_RESULT error code (see @ref RM_DecodeError).
@see                 
@ref RM_GetComponent, 
@ref RM_GetSpeciesConcentrations,
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas, 
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesZ, 
@ref RM_SetComponentH2O,
@ref RM_SetSpeciesSaveOn, 
@ref RM_SpeciesConcentrations2Module. 
@par The RM_FindComponents method also generates lists of reactants--equilibrium phases,
exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
The lists are cumulative, including all reactants that were
defined in the initial phreeqc instance at any time FindComponents was called.
In addition, a list of phases is generated for which saturation indices may be calculated from the
cumulative list of components.
@see 
@ref RM_GetEquilibriumPhasesName,
@ref RM_GetEquilibriumPhasesCount,
@ref RM_GetExchangeName,
@ref RM_GetExchangeSpeciesName,
@ref RM_GetExchangeSpeciesCount,
@ref RM_GetGasComponentsName,
@ref RM_GetGasComponentsCount,
@ref RM_GetKineticReactionsName,
@ref RM_GetKineticReactionsCount,
@ref RM_GetSICount,
@ref RM_GetSIName,
@ref RM_GetSolidSolutionComponentsName,
@ref RM_GetSolidSolutionComponentsCount,
@ref RM_GetSolidSolutionName,
@ref RM_GetSurfaceName,
@ref RM_GetSurfaceSpeciesName,
@ref RM_GetSurfaceSpeciesCount,
@ref RM_GetSurfaceType.
@par C Example:
@htmlonly
<CODE>
<PRE>
// Get list of components
ncomps = RM_FindComponents(id);
components = (char **) malloc((size_t) (ncomps * sizeof(char *)));
for (i = 0; i < ncomps; i++)
{
  components[i] = (char *) malloc((size_t) (100 * sizeof(char *)));
  status = RM_GetComponent(id, i, components[i], 100);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT int        RM_FindComponents(int id);
/**
Fills an array with the cell numbers in the user's numbering sytstem that map to a cell in the
PhreeqcRM numbering system. The mapping is defined by @ref RM_CreateMapping.

@param id            The instance @a id returned from @ref RM_Create.
@param n             A cell number in the PhreeqcRM numbering system (0 <= n < @ref RM_GetChemistryCellCount).
@param list          Array to store the user cell numbers mapped to PhreeqcRM cell @a n.
@param size          Input, the allocated size of @a list; it is an error if the array is too small. 
                     Output, the number of cells mapped to cell @a n.
@retval              IRM_RESULT error code (see @ref RM_DecodeError).

@see                 
@ref RM_CreateMapping, 
@ref RM_GetChemistryCellCount, 
@ref RM_GetGridCellCount.

@par C Example:
@htmlonly
<CODE>
<PRE>
if (RM_GetBackwardMapping(rm_id, rm_cell_number, list, &size) == 0)
{
  if (strcmp(str, "HYDRAULIC_K") == 0)
  {
    return data->K_ptr[list[0]];
  }
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetBackwardMapping(int id, int n, int *list, int *size);
/**
Returns the number of chemistry cells in the reaction module. The number of chemistry cells is defined by
the set of non-negative integers in the mapping from user grid cells (@ref RM_CreateMapping).
The number of chemistry cells is less than or equal to the number of cells in the user's model.
@param id            The instance @a id returned from @ref RM_Create.
@retval              Number of chemistry cells, or IRM_RESULT error code (see @ref RM_DecodeError).
@see                 
@ref RM_CreateMapping, 
@ref RM_GetGridCellCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_CreateMapping(id, grid2chem);
nchem = RM_GetChemistryCellCount(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int        RM_GetChemistryCellCount(int id);
/**
Retrieves an item from the reaction-module component list that was generated by calls to @ref RM_FindComponents.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the component to be retrieved. C, 0 based.
@param chem_name        The string value associated with component @a num.
@param l                The length of the maximum number of characters for @a chem_name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetComponentCount
@par C Example:
@htmlonly
<CODE>
<PRE>
// Get list of components
ncomps = RM_FindComponents(id);
components = (char **) malloc((size_t) (ncomps * sizeof(char *)));
for (i = 0; i < ncomps; i++)
{
  components[i] = (char *) malloc((size_t) (100 * sizeof(char *)));
  status = RM_GetComponent(id, i, components[i], 100);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetComponent(int id, int num, char *chem_name, int l);
/**
Returns the number of components in the reaction-module component list.
The component list is generated by calls to @ref RM_FindComponents.
The return value from the last call to @ref RM_FindComponents is equal to the return value from RM_GetComponentCount.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of components in the reaction-module component list, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetComponent.
@par C Example:
@htmlonly
<CODE>
<PRE>
ncomps1 = RM_GetComponentCount(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT int RM_GetComponentCount(int id);

/**
Transfer solution concentrations from each reaction cell
to the concentration array given in the argument list (@a c).
Units of concentration for @a c are defined by @ref RM_SetUnitsSolution.
For concentration units of per liter,
the solution volume is used to calculate the concentrations for @a c.
For mass fraction concentration units,
the solution mass is used to calculate concentrations for @a c.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or
(2) the volume of solution is the product of saturation (@ref RM_SetSaturationUser),
porosity (@ref RM_SetPorosity), and representative volume (@ref RM_SetRepresentativeVolume),
and the mass of solution is volume times density as defined by @ref RM_SetDensityUser.
@ref RM_UseSolutionDensityVolume determines which option is used.
For option 1, the databases that have partial molar volume definitions needed
to accurately calculate solution volume are
phreeqc.dat, Amm.dat, and pitzer.dat.

@param id               The instance @a id returned from @ref RM_Create.
@param c                Array to receive the concentrations. Dimension of the array is @a nxyz * @a ncomps,
where @a nxyz is the number of user grid cells and @a ncomps is the result of @ref RM_FindComponents or @ref RM_GetComponentCount.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see                    
@ref RM_FindComponents, 
@ref RM_GetComponentCount, 
@ref RM_GetSaturationCalculated,
@ref RM_SetConcentrations, 
@ref RM_SetDensityUser, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser,
@ref RM_SetUnitsSolution, 
@ref RM_UseSolutionDensityVolume.
@par C Example:
@htmlonly
<CODE>
<PRE>
c = (double *) malloc((size_t) (ncomps * nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetConcentrations(int id, double *c);
/**
Transfer the concentration from each cell for one component to the array given in the
argument list (@a c). The concentrations are those resulting from the last call
to @ref RM_RunCells. Units of concentration for @a c are defined by @ref RM_SetUnitsSolution.
@param id               The instance @a id returned from @ref RM_Create.
@param i                Zero-based index for the component to retrieve. Indices refer
to the order produced by @ref RM_GetComponent. The total number of components is given by
@ref RM_GetComponentCount.
@param c                Allocated array to receive the component concentrations.
Dimension of the array must be at least @a nxyz, where @a nxyz is the number of
user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see    @ref RM_FindComponents,
@ref RM_GetComponent,
@ref RM_GetComponentCount,
@ref RM_GetConcentrations.
@par C Example:
@htmlonly
<CODE>
<PRE>
c = (double*)malloc(nxyz * sizeof(double));
status = RM_RunCells(id);
status = RM_GetIthConcentration(id, 0, c)
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetIthConcentration(int id, int i, double* c);
/**
Transfer the concentrations for one species from each cell to the array given in the
argument list (@a c). The concentrations are those resulting from the last call
to @ref RM_RunCells. Units of concentration for @a c are mol/L.
To retrieve species concentrations, @ref RM_SetSpeciesSaveOn must be set to 1.
This method is for use with multicomponent diffusion calculations.
@param id               The instance @a id returned from @ref RM_Create.
@param i                Zero-based index for the species to retrieve. Indices refer
to the order given by @ref RM_GetSpeciesName. The total number of species is given
by @ref RM_GetSpeciesCount.
@param c                Allocated array to receive the species concentrations.
Dimension of the array must be at least @a nxyz, where @a nxyz is the number of
user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see         @ref RM_FindComponents,
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesConcentrations,
@ref RM_SetSpeciesSaveOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
c = (double*)malloc(nxyz*sizeof(double));
status = RM_RunCells(id);
status = RM_GetIthSpeciesConcentration(id, 0, c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetIthSpeciesConcentration(int id, int i, double* c);
/**
Transfer the concentrations for one component given by the vector @a c
to each reaction cell.
Units of concentration for @a c are defined by @ref RM_SetUnitsSolution.
It is required that  @a RM_SetIthConcentration be called for each component
in the system before @ref RM_RunCells is called.
@param id               The instance @a id returned from @ref RM_Create.
@param i                Zero-based index for the component to transfer.
Indices refer to the order produced by @ref RM_GetComponent. The total number
of components is given by @ref RM_GetComponentCount.
@param c                Array of concentrations to transfer to the reaction cells.
Dimension of the vector is @a nxyz, where @a nxyz is the number of
user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are ignored.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetComponentCount,
@ref RM_GetComponent,
@ref RM_SetConcentrations.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetIthConcentration(id, i, c); ! repeat for all components
...
status = phreeqc_rm.RunCells(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_SetIthConcentration(int id, int i, double* c);
/**
Transfer the concentrations for one aqueous species given by the vector
@a c to each reaction cell.
Units of concentration for @a c are mol/L. To set species concentrations,
@ref RM_SetSpeciesSaveOn must be set to @a true. It is required that
@a RM_SetIthSpeciesConcentration be called for each aqueous species in the
system before @ref RM_RunCells is called. This method is for use with
multicomponent diffusion calculations.
@param id               The instance @a id returned from @ref RM_Create.
@param i                Zero-based index for the species to transfer. Indices
refer to the order produced by @ref RM_GetSpeciesName. The total number of
species is given by @ref RM_GetSpeciesCount.
@param c                Array of concentrations to transfer to the reaction cells.
Dimension of the array is @a nxyz, where @a nxyz is the number of user grid
cells (@ref RM_GetGridCellCount). Values for inactive cells are ignored.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesName,
@ref RM_SpeciesConcentrations2Module,
@ref RM_SetSpeciesSaveOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetIthSpeciesConcentration(id, i, c); ! repeat for all species
...
status = RM_RunCells(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_SetIthSpeciesConcentration(int id, int i, double* c);

/**
Returns the user number of the current selected-output definition.
@ref RM_SetCurrentSelectedOutputUserNumber or @ref RM_SetNthSelectedOutput specifies which of the
selected-output definitions is used.
@retval                 User number of the the current selected-output definition,
negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
	{
	  status = RM_SetNthSelectedOutputUser(id, isel);
	  n_user = RM_GetCurrentSelectedOutputUserNumber(id);
	  col = RM_GetSelectedOutputColumnCount(id);
	  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
	  status = RM_GetSelectedOutput(id, selected_out);
	  // Process results here
	  free(selected_out);
	}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
*/
IRM_DLL_EXPORT int RM_GetCurrentSelectedOutputUserNumber(int id);

/**
Transfer solution densities from the reaction cells to the array given in the argument list (@a density).
Densities are those calculated by the reaction module. This method always 
returns the calculated densities; @ref RM_SetDensityUser does not affect the result.
Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate density:
phreeqc.dat, Amm.dat, and pitzer.dat.

@param id                   The instance @a id returned from @ref RM_Create.
@param density              Array to receive the densities. Dimension of the array is @a nxyz,
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).

@par C Example:
@htmlonly
<CODE>
<PRE>
density = (double *) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetDensityCalculated(id, density);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetDensityCalculated(int id, double *density);
/**
Deprecated equivalent of RM_GetDensityCalculated.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetDensity(int id, double* density);


/**
Returns an array with the ending cell numbers from the range of cell numbers assigned to each worker.
@param id               The instance @a id returned from @ref RM_Create.
@param ec               Array to receive the ending cell numbers. Dimension of the array is 
                        the number of threads (OpenMP) or the number of processes (MPI).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create, 
@ref RM_GetMpiTasks, 
@ref RM_GetStartCell, 
@ref RM_GetThreadCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
n = RM_GetThreadCount(id) * RM_GetMpiTasks(id);
ec = (int *) malloc((size_t) (n * sizeof(int)));
status = RM_GetEndCell(id, ec);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetEndCell(int id, int *ec);
/**
Returns the number of equilibrium phases in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetEquilibriumPhasesCount.
This method may be useful when generating selected output definitions related to
equilibrium phases.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of equilibrium phases in the initial-phreeqc module.
@see
@ref RM_FindComponents,
@ref RM_GetEquilibriumPhasesName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -equilibrium_phases\n");
for (i = 0; i < RM_GetEquilibriumPhasesCount(id); i++)
{
status = RM_GetEquilibriumPhasesName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int        RM_GetEquilibriumPhasesCount(int id);
/**
Retrieves an item from the equilibrium phase list.
The list includes all phases included in any EQUILIBRIUM_PHASES definitions in
the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetEquilibriumPhasesName.
This method may be useful when generating selected output definitions related to equilibrium phases.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the equilibrium phase name to be retrieved. (0 basedindex.)
@param name             The equilibrium phase name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetEquilibriumPhasesCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -equilibrium_phases\n");
for (i = 0; i < RM_GetEquilibriumPhasesCount(id); i++)
{
status = RM_GetEquilibriumPhasesName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetEquilibriumPhasesName(int id, int num, char *name, int l1);
/**
Returns a string containing error messages related to the last call to a PhreeqcRM method to
the character argument (@a errstr).

@param id               The instance @a id returned from @ref RM_Create.
@param errstr           The error string related to the last call to a PhreeqcRM method.
@param l                Maximum number of characters that can be written to the argument (@a errstr).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@par C Example:
@htmlonly
<CODE>
<PRE>
if (status != IRM_OK)
{
  l = RM_GetErrorStringLength(id);
  errstr = (char *) malloc((size_t) (l * sizeof(char) + 1));
  RM_GetErrorString(id, errstr, l+1);
  fprintf(stderr,"%s", errstr);
  free(errstr);
  RM_Destroy(id);
  exit(1);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetErrorString(int id, char *errstr, int l);
/**
Returns the length of the string that contains error messages related to the last call to a PhreeqcRM method.
@param id               The instance @a id returned from @ref RM_Create.
@retval int             Length of the error message string (for C, equivalent to strlen, does not include terminating \0).
@par C Example:
@htmlonly
<CODE>
<PRE>
if (status != IRM_OK)
{
  l = RM_GetErrorStringLength(id);
  errstr = (char *) malloc((size_t) (l * sizeof(char) + 1));
  RM_GetErrorString(id, errstr, l+1);
  fprintf(stderr,"%s", errstr);
  free(errstr);
  RM_Destroy(id);
  exit(1);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker..
 */
IRM_DLL_EXPORT int RM_GetErrorStringLength(int id);
/**
Retrieves an item from the exchange name list.
@ref RM_FindComponents must be called before @ref RM_GetExchangeName.
The exchange names vector is the same length as the exchange species names vector
and provides the corresponding exchange site (for example, X corresponing to NaX).
This method may be useful when generating selected output definitions related to exchangers.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the exchange name to be retrieved. (0 based index.)
@param name             The exchange name associated with exchange species @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetExchangeSpeciesCount, @ref RM_GetExchangeSpeciesName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetExchangeSpeciesCount(id); i++)
{
Utilities::strcpy_safe(line, MAX_LENGTH, "");
status = RM_GetExchangeSpeciesName(id, i, line1, 100);
status = RM_GetExchangeName(id, i, line2, 100);
sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetExchangeName(int id, int num, char *name, int l1);
/**
Returns the number of exchange species in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetExchangeSpeciesCount.
This method may be useful when generating selected output definitions related to exchangers.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of exchange species in the initial-phreeqc module.
@see
@ref RM_FindComponents,
@ref RM_GetExchangeSpeciesName, @ref RM_GetExchangeName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetExchangeSpeciesCount(id); i++)
{
Utilities::strcpy_safe(line, MAX_LENGTH), "");
status = RM_GetExchangeSpeciesName(id, i, line1, 100);
status = RM_GetExchangeName(id, i, line2, 100);
sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int RM_GetExchangeSpeciesCount(int id);
/**
Retrieves an item from the exchange species list.
The list of exchange species (such as "NaX") is derived from the list of components
(@ref RM_FindComponents) and the list of all exchange names (such as "X")
that are included in EXCHANGE definitions in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetExchangeSpeciesName.
This method may be useful when generating selected output definitions related to exchangers.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the exchange species to be retrieved. (0 based index.)
@param name             The exchange species name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetExchangeSpeciesCount, @ref RM_GetExchangeName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetExchangeSpeciesCount(id); i++)
{
Utilities::strcpy_safe(line, MAX_LENGTH, "");
status = RM_GetExchangeSpeciesName(id, i, line1, 100);
status = RM_GetExchangeName(id, i, line2, 100);
sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetExchangeSpeciesName(int id, int num, char *name, int l1);
/**
Returns the reaction-module file prefix to the character argument (@a prefix).
@param id               The instance @a id returned from @ref RM_Create.
@param prefix           Character string where the prefix is written.
@param l                Maximum number of characters that can be written to the argument (@a prefix).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetFilePrefix.
@par C Example:
@htmlonly
<CODE>
<PRE>
char str[100], str1[200];
status = RM_GetFilePrefix(id, str, 100);
Utilities::strcpy_safe(str1, MAX_LENGTH, "File prefix: ");
Utilities::strcat_safe(str1, MAX_LENGTH, str);
Utilities::strcat_safe(str1, MAX_LENGTH, "\n");
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetFilePrefix(int id, char *prefix, int l);
/**
Returns the number of gas phase components in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetGasComponentsCount.
This method may be useful when generating selected output definitions related to
gas phases.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of gas phase components in the initial-phreeqc module.
@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -gases\n");
for (i = 0; i < RM_GetGasComponentsCount(id); i++)
{
status = RM_GetGasComponentsName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int        RM_GetGasComponentsCount(int id);
/**
Retrieves an item from the gas components list.
The list includes all gas components included in any GAS_PHASE definitions in
the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetGasComponentsName.
This method may be useful when generating selected output definitions related to gas phases.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the gas component name to be retrieved. (0 based index.)
@param name             The gas component name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -gases\n");
for (i = 0; i < RM_GetGasComponentsCount(id); i++)
{
status = RM_GetGasComponentsName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetGasComponentsName(int id, int num, char *name, int l1);
/**
Transfer moles of gas components from each reaction cell
to the vector given in the argument list (@a gas_moles).

@param id               The instance @a id returned from @ref RM_Create.
@param gas_moles        Vector to receive the moles of gas components.
Dimension of the vector must be @a ngas_comps times @a nxyz,
where, @a ngas_comps is the result of @ref RM_GetGasComponentsCount,
and @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If a gas component is not defined for a cell, the number of moles is set to -1.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see
@ref RM_FindComponents, 
@ref RM_GetGasComponentsCount, 
@ref RM_GetGasCompPressures,
@ref RM_GetGasCompPhi,
@ref RM_GetGasPhaseVolume,
@ref RM_SetGasCompMoles,
@ref RM_SetGasPhaseVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
ngas_comps = RM_GetGasComponentsCount();
gas_moles = (double *) malloc((size_t) (ngas_comps * nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetGasCompMoles(id, gas_moles);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetGasCompMoles(int id, double* gas_moles);
/**
Transfer pressures of gas components from each reaction cell
to the vector given in the argument list (@a gas_pressure).

@param id               The instance @a id returned from @ref RM_Create.
@param gas_pressure        Vector to receive the pressures of gas components.
Dimension of the vector must be @a ngas_comps times @a nxyz,
where, @a ngas_comps is the result of @ref RM_GetGasComponentsCount,
and @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If a gas component is not defined for a cell, the pressure is set to -1.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsCount,
@ref RM_GetGasCompMoles,
@ref RM_GetGasCompPhi,
@ref RM_GetGasPhaseVolume,
@ref RM_SetGasCompMoles,
@ref RM_SetGasPhaseVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
ngas_comps = RM_GetGasComponentsCount();
gas_pressure = (double *) malloc((size_t) (ngas_comps * nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetGasCompPressures(id, gas_pressure);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetGasCompPressures(int id, double* gas_pressure);
/**
Transfer fugacity coefficients (phi) of gas components from each reaction cell
to the vector given in the argument list (@a gas_phi). Fugacity of a gas
component is equal to its pressure times the fugacity coefficient.

@param id               The instance @a id returned from @ref RM_Create.
@param gas_phi        Vector to receive the fugacity coefficients of gas components.
Dimension of the vector must be @a ngas_comps times @a nxyz,
where, @a ngas_comps is the result of @ref RM_GetGasComponentsCount,
and @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If a gas component is not defined for a cell, the fugacity coefficient is set to -1.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsCount,
@ref RM_GetGasCompMoles,
@ref RM_GetGasCompPressures,
@ref RM_GetGasPhaseVolume,
@ref RM_SetGasCompMoles,
@ref RM_SetGasPhaseVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
ngas_comps = RM_GetGasComponentsCount();
gas_phi = (double *) malloc((size_t) (ngas_comps * nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetGasCompPhi(id, gas_phi);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetGasCompPhi(int id, double* gas_phi);
/**
Transfer volume of gas from each reaction cell
to the vector given in the argument list (@a gas_volume).

@param id               The instance @a id returned from @ref RM_Create.
@param  gas_volume               Array to receive the gas phase volumes.
Dimension of the vector must be @a nxyz,
where,  @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If a gas phase is not defined for a cell, the volume is set to -1.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsCount,
@ref RM_GetGasCompMoles,
@ref RM_GetGasCompPressures,
@ref RM_GetGasCompPhi,
@ref RM_SetGasCompMoles,
@ref RM_SetGasPhaseVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
gas_volume = (double *) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetGasPhaseVolume(id, gas_volume);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetGasPhaseVolume(int id, double* gas_volume);
/**
Returns the gram formula weights (g/mol) for the components in the reaction-module component list.
@param id               The instance id returned from @ref RM_Create.
@param gfw              Array to receive the gram formula weights. Dimension of the array is @a ncomps,
where @a ncomps is the number of components in the component list.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents,  
@ref RM_GetComponent,
@ref RM_GetComponentCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
ncomps = RM_FindComponents(id);
components = (char **) malloc((size_t) (ncomps * sizeof(char *)));
gfw = (double *) malloc((size_t) (ncomps * sizeof(double)));
status = RM_GetGfw(id, gfw);
for (i = 0; i < ncomps; i++)
{
  components[i] = (char *) malloc((size_t) (100 * sizeof(char *)));
  status = RM_GetComponent(id, i, components[i], 100);
  sprintf(str,"%10s    %10.3f\n", components[i], gfw[i]);
  status = RM_OutputMessage(id, str);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetGfw(int id, double * gfw);
/**
Returns the number of grid cells in the user's model, which is defined in the call to @ref RM_Create.
The mapping from grid cells to reaction cells is defined by @ref RM_CreateMapping.
The number of reaction cells may be less than the number of grid cells if
there are inactive regions or symmetry in the model definition.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 Number of grid cells in the user's model, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create,  
@ref RM_CreateMapping.
@par C Example:
@htmlonly
<CODE>
<PRE>
nxyz = RM_GetGridCellCount(id);
sprintf(str1, "Number of grid cells in the user's model: %d\n", nxyz);
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int        RM_GetGridCellCount(int id);

#ifdef USE_YAML
/**
@a GetGridCellCountYAML will read a YAML file and extract the value
of GridCellCount, which can be used to construct a PhreeqcRM
instance. @ref RM_Create requires a value for the number of cells.
If a GUI or preprocessor is used to write a YAML file to initialize
PhreeqcRM, the number of cells can be written to the YAML file
and extracted with this method.
@param config_file         String containing the YAML file name.
@retval Number of grid cells specified in the YAML file; returns
zero if GridCellCount is not defined.
@see @ref RM_Create, and @ref RM_InitializeYAML.
@par C Example:
@htmlonly
<CODE>
<PRE>
int nthreads = 0;
int nxyz;
nxyz = RM_GetGridCellCountYAML("myfile.yaml");
status = RM_Create(nxyz, nthreads);
status = RM_InitializeYAML("myfile.yaml");
</PRE>
</CODE>
@endhtmlonly
@par Sequence:
Called before RM_Create.
*/
IRM_DLL_EXPORT int        RM_GetGridCellCountYAML(const char* config_file);
#endif
/**
Returns an IPhreeqc id for the @a ith IPhreeqc instance in the reaction module.


For the threaded version, there are @a nthreads + 2 IPhreeqc instances, where
@a nthreads is defined in the constructor (@ref RM_Create).
The number of threads can be determined by @ref RM_GetThreadCount.
The first @a nthreads (0 based) instances will be the workers, the
next (@a nthreads) is the InitialPhreeqc instance, and the next (@a nthreads + 1) is the Utility instance.
Getting the IPhreeqc pointer for one of these instances allows the user to use any of the IPhreeqc methods
on that instance.
For MPI, each process has exactly three IPhreeqc instances, one worker (number 0),
one InitialPhreeqc instance (number 1), and one Utility instance (number 2).

@param id               The instance @a id returned from @ref RM_Create.
@param i                The number of the IPhreeqc instance to be retrieved (0 based).
@retval                 IPhreeqc id for the @a ith IPhreeqc instance, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create, 
@ref RM_GetThreadCount.
See IPhreeqc documentation for descriptions of IPhreeqc methods.
@par C Example:
@htmlonly
<CODE>
<PRE>
// Utility pointer is worker number nthreads + 1
iphreeqc_id1 = RM_GetIPhreeqcId(id, RM_GetThreadCount(id) + 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int        RM_GetIPhreeqcId(int id, int i);
/**
Returns the number of kinetic reactions in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetKineticReactionsCount.
This method may be useful when generating selected output definitions related to
kinetic reactions.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of kinetic reactions in the initial-phreeqc module.
@see
@ref RM_FindComponents,
@ref RM_GetKineticReactionsName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -kinetics\n");
for (i = 0; i < RM_GetKineticReactionsCount(id); i++)
{
status = RM_GetKineticReactionsName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int        RM_GetKineticReactionsCount(int id);
/**
Retrieves an item from the kinetic reactions list.
The list includes all kinetic reactions included in any KINETICS definitions in
the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetKineticReactionsName.
This method may be useful when generating selected output definitions related to kinetic reactions.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the kinetic reaction name to be retrieved. (0 based index.)
@param name             The kinetic reaction name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetKineticReactionsCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -kinetics\n");
for (i = 0; i < RM_GetKineticReactionsCount(id); i++)
{
status = RM_GetKineticReactionsName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetKineticReactionsName(int id, int num, char *name, int l1);
/**
Returns the MPI task number. For the OPENMP version, the task number is always
zero and the result of @ref RM_GetMpiTasks is one. For the MPI version,
the root task number is zero, and all workers have a task number greater than zero.
The number of tasks can be obtained with @ref RM_GetMpiTasks. The number of
tasks and computer hosts are determined at run time by the mpiexec command, and the
number of reaction-module processes is defined by the communicator used in
constructing the reaction modules (@ref RM_Create).
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The MPI task number for a process, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_GetMpiTasks.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str1, "MPI task number:  %d\n", RM_GetMpiMyself(id));
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int        RM_GetMpiMyself(int id);
/**
Returns the number of MPI processes (tasks) assigned to the reaction module.
For the OPENMP version, the number of tasks is always
one (although there may be multiple threads, @ref RM_GetThreadCount),
and the task number returned by @ref RM_GetMpiMyself is zero. For the MPI version, the number of
tasks and computer hosts are determined at run time by the mpiexec command. An MPI communicator
is used in constructing reaction modules for MPI. The communicator may define a subset of the
total number of MPI processes.
The root task number is zero, and all workers have a task number greater than zero.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of MPI  processes assigned to the reaction module,
negative is failure (See @ref RM_DecodeError).
@see    
@ref RM_Create,                
@ref RM_GetMpiMyself.
@par C Example:
@htmlonly
<CODE>
<PRE>
mpi_tasks = RM_GetMpiTasks(id);
sprintf(str1, "Number of MPI processes: %d\n", mpi_tasks);
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int        RM_GetMpiTasks(int id);
/**
Returns the user number for the @a nth selected-output definition.
Definitions are sorted by user number. Phreeqc allows multiple selected-output
definitions, each of which is assigned a nonnegative integer identifier by the
user. The number of definitions can be obtained by @ref RM_GetSelectedOutputCount.
To cycle through all of the definitions, RM_GetNthSelectedOutputUserNumber
can be used to identify the user number for each selected-output definition
in sequence. @ref RM_SetCurrentSelectedOutputUserNumber is then used to select
that user number for selected-output processing.
@param id               The instance @a id returned from @ref RM_Create.
@param n                The sequence number of the selected-output definition for which the user number will be returned.
C, 0 based.
@retval                 The user number of the @a nth selected-output definition, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  fprintf(stderr, "Selected output sequence number: %d\n", isel);
  fprintf(stderr, "Selected output user number:     %d\n", n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT int        RM_GetNthSelectedOutputUserNumber(int id, int n);
/**
Transfer current porosities to the array given in the argument list (@a porosity).
Porosity is not changed by PhreeqcRM; the values are either the default values
or the values set by the last call to @ref RM_SetPorosity.
@param id                The instance @a id returned from @ref RM_Create.
@param porosity           Array to receive the porosities. Dimension of the array
must be @a nxyz, where @a nxyz is the number of user grid cells
(@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
porosity = (double*)malloc(nxyz*sizeof(double));
status = RM_GetPorosity(id, porosity);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetPorosity(int id, double* porosity);
/**
Transfer current pressures to the array given in the argument list (@a pressure).
Pressure is not usually calculated by PhreeqcRM; the values are either the default values
or the values set by the last call to @ref RM_SetPressure. Pressures can be calculated
by PhreeqcRM if a fixed-volume GAS_PHASE is used.
@param id              The instance @a id returned from @ref RM_Create.
@param pressure        Array to receive the porosities. Dimension of the array must be
@a nxyz, where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
Values for inactive cells are set to 1e30.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
pressure = (double*)malloc(nxyz*sizeof(double));
status = RM_GetPressure(id, pressure);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetPressure(int id, double* pressure);

/**
Returns a vector of saturations (@a sat) as calculated by the reaction module.
This method always returns solution_volume/(rv * porosity); the method 
@ref RM_SetSaturationUser has no effect on the values returned.
Reactions will change the volume of solution in a cell.
The transport code must decide whether to ignore or account for this change in solution 
volume due to reactions. Following reactions, the cell saturation is calculated as solution 
volume (@ref RM_GetSolutionVolume) divided by the product of representative volume 
(@ref RM_SetRepresentativeVolume) and the porosity (@ref RM_SetPorosity). The cell 
saturation returned by @a RM_GetSaturationCalculated may be less than or greater 
than the saturation set by the transport code (@ref RM_SetSaturationUser and may be 
greater than or less than 1.0, even in fully saturated simulations. Only the following 
databases distributed with PhreeqcRM have molar volume information needed to accurately 
calculate solution volume and saturation: phreeqc.dat, Amm.dat, and pitzer.dat.

@param id               The instance @a id returned from @ref RM_Create.
@param sat_calc              Vector to receive the saturations. Dimension of the array is set to @a nxyz,
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@see                    
@ref RM_GetSolutionVolume, 
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser.
@par C Example:
@htmlonly
<CODE>
<PRE>
sat_calc = (double *) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSaturationCalculated(id, sat_calc);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT               RM_GetSaturationCalculated(int id, double *sat_calc);
/**
Deprecated equivalent of RM_GetSaturationCalculated.
*/
IRM_DLL_EXPORT IRM_RESULT               RM_GetSaturation(int id, double* sat_calc);

/**
Populates an array with values from the current selected-output definition. @ref RM_SetCurrentSelectedOutputUserNumber
determines which of the selected-output definitions is used to populate the array.
@param id               The instance @a id returned from @ref RM_Create.
@param so               An array to contain the selected-output values. Size of the array is @a nxyz x @a col,
where @a nxyz is the number of grid cells in the user's model (@ref RM_GetGridCellCount), and @a col is the number of
columns in the selected-output definition (@ref RM_GetSelectedOutputColumnCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT        RM_GetSelectedOutput(int id, double *so);

/**
Returns the number of columns in the current selected-output definition. @ref RM_SetCurrentSelectedOutputUserNumber
determines which of the selected-output definitions is used.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 Number of columns in the current selected-output definition, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT int        RM_GetSelectedOutputColumnCount(int id);
/**
Returns the number of selected-output definitions. @ref RM_SetCurrentSelectedOutputUserNumber
determines which of the selected-output definitions is used.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 Number of selected-output definitions, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT int        RM_GetSelectedOutputCount(int id);
/**
Returns a selected output heading. The number of headings is determined by @ref RM_GetSelectedOutputColumnCount.
@ref RM_SetCurrentSelectedOutputUserNumber
determines which of the selected-output definitions is used.
@param id               The instance @a id returned from @ref RM_Create.
@param icol             The sequence number of the heading to be retrieved. C, 0 based.
@param heading          A string buffer to receive the heading.
@param length           The maximum number of characters that can be written to the string buffer.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
char heading[100];
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  for (j = 0; j < col; j++)
  {
	status = RM_GetSelectedOutputHeading(id, j, heading, 100);
	fprintf(stderr, "          %2d %10s\n", j, heading);
  }
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT        RM_GetSelectedOutputHeading(int id, int icol, char * heading, int length);
/**
Returns the number of rows in the current selected-output definition. However, the method
is included only for convenience; the number of rows is always equal to the number of
grid cells in the user's model, and is equal to @ref RM_GetGridCellCount.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 Number of rows in the current selected-output definition, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Print results
  for (i = 0; i < RM_GetSelectedOutputRowCount(id)/2; i++)
  {
    fprintf(stderr, "Cell number %d\n", i);
    fprintf(stderr, "     Selected output: \n");
    for (j = 0; j < col; j++)
    {
      status = RM_GetSelectedOutputHeading(id, j, heading, 100);
      fprintf(stderr, "          %2d %10s: %10.4f\n", j, heading, selected_out[j*nxyz + i]);
    }
  }
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT int        RM_GetSelectedOutputRowCount(int id);
/**
Returns the number of phases in the initial-phreeqc module for which saturation indices can be calculated.
@ref RM_FindComponents must be called before @ref RM_GetSICount.
This method may be useful when generating selected output definitions related to
saturation indices.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of phases in the initial-phreeqc module for which saturation indices
could be calculated.
@see
@ref RM_FindComponents,
@ref RM_GetSIName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -saturation_indices\n");
for (i = 0; i < RM_GetSICount(id); i++)
{
status = RM_GetSIName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int        RM_GetSICount(int id);
/**
Retrieves an item from the list of all phases for which saturation indices can be calculated.
The list includes all phases that contain only elements included in the components in
the initial-phreeqc module.
The list assumes that all components are present to be able to calculate the entire list of SIs;
it may be that one or more components are missing in any specific cell.
@ref RM_FindComponents must be called before @ref RM_GetSIName.
This method may be useful when generating selected output definitions related to saturation indices.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the saturation-index-phase name to be retrieved. (0 based index.)
@param name             The saturation-index-phase name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSICount.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -saturation_indices\n");
for (i = 0; i < RM_GetSICount(id); i++)
{
status = RM_GetSIName(id, i, line1, 100);
sprintf(line, "%4s%20s\n", "    ", line1);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSIName(int id, int num, char *name, int l1);

 /**
 Returns the number of solid solution components in the initial-phreeqc module.
 @ref RM_FindComponents must be called before @ref RM_GetSolidSolutionComponentsCount.
 This method may be useful when generating selected output definitions related to solid solutions.
 @param id               The instance @a id returned from @ref RM_Create.
 @retval                 The number of solid solution components in the initial-phreeqc module.
 @see
 @ref RM_FindComponents,
 @ref RM_GetSolidSolutionComponentsName, @ref RM_GetSolidSolutionName.
 @par C Example:
 @htmlonly
 <CODE>
 <PRE>
 Utilities::strcat_safe(input, MAX_LENGTH, "  -solid_solutions\n");
 for (i = 0; i < RM_GetSolidSolutionComponentsCount(id); i++)
 {
 status = RM_GetSolidSolutionComponentsName(id, i, line1, 100);
 status = RM_GetSolidSolutionName(id, i, line2, 100);
 sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
 Utilities::strcat_safe(input, MAX_LENGTH, line);
 }
 </PRE>
 </CODE>
 @endhtmlonly
 @par MPI:
 Called by root.
 */
IRM_DLL_EXPORT int        RM_GetSolidSolutionComponentsCount(int id);
/**
Retrieves an item from the solid solution components list.
The list includes all solid solution components included in any SOLID_SOLUTIONS definitions in
the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetSolidSolutionComponentsName.
This method may be useful when generating selected output definitions related to solid solutions.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the solid solution components name to be retrieved. (0 based index.)
@param name             The solid solution compnent name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSolidSolutionComponentsCount, @ref RM_GetSolidSolutionName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -solid_solutions\n");
for (i = 0; i < RM_GetSolidSolutionComponentsCount(id); i++)
{
status = RM_GetSolidSolutionComponentsName(id, i, line1, 100);
status = RM_GetSolidSolutionName(id, i, line2, 100);
sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSolidSolutionComponentsName(int id, int num, char *name, int l1);
/**
Retrieves an item from the solid solution names list.
The list includes solid solution names included in SOLID_SOLUTIONS definitions in
the initial-phreeqc module.
The solid solution names vector is the same length as the solid solution components vector
and provides the corresponding name of solid solution containing the component.
@ref RM_FindComponents must be called before @ref RM_GetSolidSolutionName.
This method may be useful when generating selected output definitions related to solid solutions.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the solid solution name to be retrieved. (0 based index.)
@param name             The solid solution name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSolidSolutionComponentsCount, @ref RM_GetSolidSolutionComponentsName.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcat_safe(input, MAX_LENGTH, "  -solid_solutions\n");
for (i = 0; i < RM_GetSolidSolutionComponentsCount(id); i++)
{
status = RM_GetSolidSolutionComponentsName(id, i, line1, 100);
status = RM_GetSolidSolutionName(id, i, line2, 100);
sprintf(line, "%4s%20s%3s%20s\n", "    ", line1, " # ", line2);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSolidSolutionName(int id, int num, char *name, int l1);

/**
Transfer solution volumes from the reaction cells to the array given in the argument list (@a vol).
Solution volumes are those calculated by the reaction module.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume:
phreeqc.dat, Amm.dat, and pitzer.dat.

@param id                   The instance @a id returned from @ref RM_Create.
@param vol                  Array to receive the solution volumes. Dimension of the array is (@a nxyz),
where @a nxyz is the number of user grid cells. Values for inactive cells are set to 1e30.
@retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetSaturationCalculated.

@par C Example:
@htmlonly
<CODE>
<PRE>
volume = (double *) malloc((size_t) (nxyz * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSolutionVolume(id, volume);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSolutionVolume(int id, double *vol);
/**
Transfer concentrations of aqueous species to the array argument (@a species_conc)
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous
species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.
Solution volumes used to calculate mol/L are calculated by the reaction module.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume:
phreeqc.dat, Amm.dat, and pitzer.dat.

@param id               The instance @a id returned from @ref RM_Create.
@param species_conc     Array to receive the aqueous species concentrations.
Dimension of the array is (@a nxyz, @a nspecies),
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount),
and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
Concentrations are moles per liter.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesCount, 
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ,   
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double *) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesConcentrations(id, species_c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesConcentrations(int id, double *species_conc);
/**
The number of aqueous species used in the reaction module.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous
species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.

@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      The number of aqueous species, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations,
@ref RM_GetSpeciesD25,
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName, 
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ, 
@ref RM_SpeciesConcentrations2Module, 
@ref RM_SetSpeciesSaveOn.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double *) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesConcentrations(id, species_c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int RM_GetSpeciesCount(int id);
/**
Transfers diffusion coefficients at 25C to the array argument (@a diffc).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally in the database file.
Databases distributed with the reaction module that have diffusion coefficients defined are
phreeqc.dat, Amm.dat, and pitzer.dat.

@param id               The instance @a id returned from @ref RM_Create.
@param diffc            Array to receive the diffusion coefficients at 25 C, m^2/s.
Dimension of the array is @a nspecies,
where @a nspecies is is the number of aqueous species (@ref RM_GetSpeciesCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,  
@ref RM_GetSpeciesSaveOn, 
@ref RM_GetSpeciesZ, 
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
diffc = (double *) malloc((size_t) (nspecies * sizeof(double)));
status = RM_GetSpeciesD25(id, diffc);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesD25(int id, double *diffc);
/**
Transfer aqueous-species log10 activity coefficients to the array argument (@a species_log10gammas)
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous
species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.

@param id                   The instance @a id returned from @ref RM_Create.
@param species_log10gammas  Array to receive the aqueous species log10 activity coefficients.
Dimension of the array is (@a nxyz, @a nspecies),
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount),
and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSpeciesConcentrations,
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ,
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_log10gammas = (double *) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesLog10Gammas(id, species_log10gammas);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesLog10Gammas(int id, double *species_log10gammas);
/**
Transfer aqueous-species log10 molalities to the array argument (@a species_log10molalities)
To use this method @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous
species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.

@param id                   The instance @a id returned from @ref RM_Create.
@param species_log10molalities  Array to receive the aqueous species log10 molalities.
Dimension of the array is (@a nxyz, @a nspecies),
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount),
and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSpeciesConcentrations,
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25,
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ,
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_log10molalities = (double *) malloc((size_t) (nxyz * nspecies * sizeof(double)));
status = RM_RunCells(id);
status = RM_GetSpeciesLog10Molalities(id, species_log10molalities);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesLog10Molalities(int id, double* species_log10molalities);
/**
Transfers the name of the @a ith aqueous species to the character argument (@a name).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous
species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.

@param id               The instance @a id returned from @ref RM_Create.
@param i                Sequence number of the species in the species list. C, 0 based.
@param name             Character array to receive the species name.
@param length           Maximum length of string that can be stored in the character array.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ,
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module.  

@par C Example:
@htmlonly
<CODE>
<PRE>
char name[100];
...
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
for (i = 0; i < nspecies; i++)
{
  status = RM_GetSpeciesName(id, i, name, 100);
  fprintf(stderr, "%s\n", name);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesName(int id, int i, char * name, int length);
/**
Returns the value of the species-save property.
By default, concentrations of aqueous species are not saved. Setting the species-save property to true allows
aqueous species concentrations to be retrieved
with @ref RM_GetSpeciesConcentrations, and solution compositions to be set with
@ref RM_SpeciesConcentrations2Module.

@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      0, species are not saved; 1, species are saved; negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesZ, 
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module. 

@par C Example:
@htmlonly
<CODE>
<PRE>
save_on = RM_GetSpeciesSaveOn(id);
if (save_on .ne. 0)
{
  fprintf(stderr, "Reaction module is saving species concentrations\n");
}
else
{
  fprintf(stderr, "Reaction module is not saving species concentrations\n");
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT int RM_GetSpeciesSaveOn(int id);
/**
Transfers the charge of each aqueous species to the array argument (@a  z).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.

@param id               The instance @a id returned from @ref RM_Create.
@param z                Array that receives the charge for each aqueous species.
Dimension of the array is @a nspecies,
where @a nspecies is is the number of aqueous species (@ref RM_GetSpeciesCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName, 
@ref RM_GetSpeciesSaveOn,
@ref RM_SetSpeciesSaveOn,
@ref RM_SpeciesConcentrations2Module. 

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
z = (double *) malloc((size_t) (nspecies * sizeof(double)));
status = RM_GetSpeciesZ(id, z);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetSpeciesZ(int id, double *z);

/**
Returns an array with the starting cell numbers from the range of cell numbers assigned to each worker.
@param id               The instance @a id returned from @ref RM_Create.
@param sc               Array to receive the starting cell numbers. Dimension of the array is 
                        the number of threads (OpenMP) or the number of processes (MPI).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create, 
@ref RM_GetEndCell, 
@ref RM_GetMpiTasks, 
@ref RM_GetThreadCount.
@par C Example:
@htmlonly
<CODE>
<PRE>
n = RM_GetThreadCount(id) * RM_GetMpiTasks(id);
sc = (int *) malloc((size_t) (n * sizeof(int)));
status = RM_GetStartCell(id, sc);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_GetStartCell(int id, int *sc);
/**
Returns an array of temperatures (@a temperature) from the reaction module.
Reactions do not change the temperature, so the temperatures are either the
temperatures at initialization, or the values set with the last call to
@ref RM_SetTemperature.
@param id               The instance @a id returned from @ref RM_Create.
@param temperature      Allocatable array to receive the temperatures.
Dimension of the array must be @a nxyz, where @a nxyz is the number of
user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are
set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_SetTemperature.
@par C Example:
@htmlonly
<CODE>
<PRE>
temperature = (double*)malloc(nxyz*sizeof(double));
status = RM_GetTemperature(id, temperature);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetTemperature(int id, double* temperature);
/**
Retrieves the surface name (such as "Hfo") that corresponds with
the surface species name.
The lists of surface species names and surface names are the same length.
@ref RM_FindComponents must be called before @ref RM_GetSurfaceName.
This method may be useful when generating selected output definitions related to surfaces.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the surface name to be retrieved. (0 based index.)
@param name             The surface name associated with surface species @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSurfaceSpeciesCount, @ref RM_GetSurfaceSpeciesName, @ref RM_GetSurfaceType.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetSurfaceSpeciesCount(id); i++)
{
status = RM_GetSurfaceSpeciesName(id, i, line1, 100);
status = RM_GetSurfaceType(id, i, line2, 100);
status = RM_GetSurfaceName(id, i, line3, 100);
sprintf(line, "%4s%20s%3s%20s%20s\n", "    ", line1, " # ", line2, line3);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSurfaceName(int id, int num, char *name, int l1);
/**
Returns the number of surface species (such as "Hfo_wOH") in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetSurfaceSpeciesCount.
This method may be useful when generating selected output definitions related to surfaces.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of surface species in the initial-phreeqc module.
@see
@ref RM_FindComponents,
@ref RM_GetSurfaceSpeciesName, @ref RM_GetSurfaceType, @ref RM_GetSurfaceName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetSurfaceSpeciesCount(id); i++)
{
status = RM_GetSurfaceSpeciesName(id, i, line1, 100);
status = RM_GetSurfaceType(id, i, line2, 100);
status = RM_GetSurfaceName(id, i, line3, 100);
sprintf(line, "%4s%20s%3s%20s%20s\n", "    ", line1, " # ", line2, line3);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT int        RM_GetSurfaceSpeciesCount(int id);
/**
Retrieves an item from the surface species list.
The list of surface species (for example, "Hfo_wOH") is derived from the list of components
(@ref RM_FindComponents) and the list of all surface types (such as "Hfo_w")
that are included in SURFACE definitions in the initial-phreeqc module.
@ref RM_FindComponents must be called before @ref RM_GetSurfaceSpeciesName.
This method may be useful when generating selected output definitions related to surfaces.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the surface type to be retrieved. (0 based index.)
@param name             The surface species name at number @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSurfaceSpeciesCount, @ref RM_GetSurfaceType, @ref RM_GetSurfaceName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetSurfaceSpeciesCount(id); i++)
{
status = RM_GetSurfaceSpeciesName(id, i, line1, 100);
status = RM_GetSurfaceType(id, i, line2, 100);
status = RM_GetSurfaceName(id, i, line3, 100);
sprintf(line, "%4s%20s%3s%20s%20s\n", "    ", line1, " # ", line2, line3);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSurfaceSpeciesName(int id, int num, char *name, int l1);
/**
Retrieves the surface site type (such as "Hfo_w") that corresponds with
the surface species name.
The lists of surface species names and surface species types are the same length.
@ref RM_FindComponents must be called before @ref RM_GetSurfaceType.
This method may be useful when generating selected output definitions related to surfaces.
@param id               The instance @a id returned from @ref RM_Create.
@param num              The number of the surface type to be retrieved. (0 based index.)
@param name             The surface type associated with surface species @a num.
@param l1               The length of the maximum number of characters for @a name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetSurfaceSpeciesCount, @ref RM_GetSurfaceSpeciesName, @ref RM_GetSurfaceName.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (i = 0; i < RM_GetSurfaceSpeciesCount(id); i++)
{
status = RM_GetSurfaceSpeciesName(id, i, line1, 100);
status = RM_GetSurfaceType(id, i, line2, 100);
status = RM_GetSurfaceName(id, i, line3, 100);
sprintf(line, "%4s%20s%3s%20s%20s\n", "    ", line1, " # ", line2, line3);
Utilities::strcat_safe(input, MAX_LENGTH, line);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetSurfaceType(int id, int num, char *name, int l1);
/**
Returns the number of threads, which is equal to the number of workers used to run in parallel with OPENMP.
For the OPENMP version, the number of threads is set implicitly or explicitly with @ref RM_Create. For the
MPI version, the number of threads is always one for each process.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The number of threads, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_GetMpiTasks.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str1, "Number of threads: %d\n", RM_GetThreadCount(id));
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; result is always 1.
 */
IRM_DLL_EXPORT int        RM_GetThreadCount(int id);
/**
Returns the current simulation time in seconds. The reaction module does not change the time value, so the
returned value is equal to the default (0.0) or the last time set by @ref RM_SetTime.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The current simulation time in seconds.
@see                    
@ref RM_GetTimeConversion, 
@ref RM_GetTimeStep, 
@ref RM_SetTime,
@ref RM_SetTimeConversion, 
@ref RM_SetTimeStep.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str, "%s%10.1f%s", "Beginning reaction calculation ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT double     RM_GetTime(int id);
/**
Returns a multiplier to convert time from seconds to another unit, as specified by the user.
The reaction module uses seconds as the time unit. The user can set a conversion
factor (@ref RM_SetTimeConversion) and retrieve it with RM_GetTimeConversion. The
reaction module only uses the conversion factor when printing the long version
of cell chemistry (@ref RM_SetPrintChemistryOn), which is rare. Default conversion factor is 1.0.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 Multiplier to convert seconds to another time unit.
@see                    
@ref RM_GetTime, 
@ref RM_GetTimeStep, 
@ref RM_SetTime, 
@ref RM_SetTimeConversion, 
@ref RM_SetTimeStep.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str, "%s%10.1f%s", "Beginning reaction calculation ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT double     RM_GetTimeConversion(int id);
/**
Returns the current simulation time step in seconds.
This is the time over which kinetic reactions are integrated in a call to @ref RM_RunCells.
The reaction module does not change the time step value, so the
returned value is equal to the default (0.0) or the last time step set by @ref RM_SetTimeStep.
@param id               The instance @a id returned from @ref RM_Create.
@retval                 The current simulation time step in seconds.
@see                    
@ref RM_GetTime, 
@ref RM_GetTimeConversion, 
@ref RM_SetTime,
@ref RM_SetTimeConversion, 
@ref RM_SetTimeStep.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str, "%s%10.1f%s", "          Time step                  ",
        RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT double     RM_GetTimeStep(int id);
/**
Transfer current viscosities to the array given in the argument list (@a viscosity).
@param id                   The instance @a id returned from @ref RM_Create.
@param viscosity            Allocated array to receive the viscosities. Dimension of
the array must be @a nxyz, where @a nxyz is the number of user grid cells
(@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
viscosity = (double*)malloc(nxyz*sizeof(double));
status = RM_GetViscosity(id, viscosity);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_GetViscosity(int id, double* viscosity);
#ifdef USE_YAML
/**
A YAML file can be used to initialize an instance of PhreeqcRM.
@param id               The instance @a id returned from @ref RM_Create.
@param yamlfile         String containing the YAML file name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
<p>
The file contains a YAML map of PhreeqcRM methods
and the arguments corresponding to the methods.
Note that the PhreeqcRM methods do not have the "RM_" prefix
and the id argument is not included.
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
@ref RM_InitializeYAML will read the YAML file and execute the specified methods with
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
are sufficient to identify which arguments are YAML scalars (single bool, int, double, string argument),
sequences (single vector argument), or maps (multiple arguments).
</p>
@htmlonly
<CODE>
<PRE>
CloseFiles(void);
CreateMapping(std::vector< int >& grid2chem);
DumpModule();
FindComponents();
InitialPhreeqc2Module(std::vector< int > initial_conditions1);
InitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
LoadDatabase(std::string database);
OpenFiles(void);
OutputMessage(std::string str);
RunCells(void);
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
@par C Example:
@htmlonly
<CODE>
<PRE>
		id = RM_Create(nxyz, MPI_COMM_WORLD)
		status = RM_InitializeYAML(id, "myfile.yaml")
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */

IRM_DLL_EXPORT IRM_RESULT RM_InitializeYAML(int id, const char* yamlfile);
#endif
/**
Fills an array (@a c) with concentrations from solutions in the InitialPhreeqc instance.
The method is used to obtain concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell. Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
If @a boundary_solution2 and fraction1 are NULL,
no mixing is used; concentrations are derived from @a boundary_solution1 only.

@param id                  The instance @a id returned from @ref RM_Create.
@param c                   Array of concentrations extracted from the InitialPhreeqc instance.
The dimension of @a c is e@a n_boundary * @a ncomp,
where @a ncomp is the number of components returned from @ref RM_FindComponents or @ref RM_GetComponentCount.
@param n_boundary          The number of boundary condition solutions that need to be filled.
@param boundary_solution1  Array of solution index numbers that refer to solutions in the InitialPhreeqc instance.
Size is @a n_boundary.
@param boundary_solution2  Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance
and are defined to mix with @a boundary_solution1.
Size is @a n_boundary. May be NULL in C.
@param fraction1           Fraction of @a boundary_solution1 that mixes with (1-@a fraction1) of @a boundary_solution2.
Size is (n_boundary). May be NULL in C.
@retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
@see                       
@ref RM_FindComponents, 
@ref RM_GetComponentCount.

@par C Example:
@htmlonly
<CODE>
<PRE>
nbound = 1;
bc1 = (int *) malloc((size_t) (nbound * sizeof(int)));
bc2 = (int *) malloc((size_t) (nbound * sizeof(int)));
bc_f1 = (double *) malloc((size_t) (nbound * sizeof(double)));
bc_conc = (double *) malloc((size_t) (ncomps * nbound * sizeof(double)));
for (i = 0; i < nbound; i++)
{
  bc1[i]          = 0;       // Solution 0 from InitialPhreeqc instance
  bc2[i]          = -1;      // no bc2 solution for mixing
  bc_f1[i]        = 1.0;     // mixing fraction for bc1
}
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_InitialPhreeqc2Concentrations(
                int id,
                double *c,
                int n_boundary,
                int *boundary_solution1,
                int *boundary_solution2,
                double *fraction1);
/**
Transfer SOLUTION definitions from the InitialPhreeqc instance to the reaction-
module workers.
@a solutions is used to select SOLUTION definitions for each
cell of the model. @a solutions is dimensioned @a nxyz, where @a nxyz is the
number of grid cells in the  user's model (@ref RM_GetGridCellCount).
@param id           The instance @a id returned from @ref RM_Create.
@param solutions    Array of SOLUTION index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values
are ignored, resulting in no transfer of a SOLUTION definition for that cell.
(Note that all cells must have a SOLUTION definition, which could be defined
by other calls to @a RM_InitialSolutions2Module, @ref RM_InitialPhreeqc2Module,
or @ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
solutions = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) solutions[i] = 1;
status = RM_InitialSolutions2Module(id, solutions);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialSolutions2Module(int id, int* solutions);
/**
Transfer EQUILIBRIUM_PHASES definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a equilibrium_phases is used to select EQUILIBRIUM_PHASES definitions for each
cell of the model. @a equilibrium_phases is dimensioned @a nxyz, where @a nxyz is
the number of grid cells in the  user's model (@ref RM_GetGridCellCount).
@param id                 The instance @a id returned from @ref RM_Create.
@param equilibrium_phases Array of EQUILIBRIUM_PHASES index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values are
ignored, resulting in no transfer of an EQUILIBRIUM_PHASES definition for that cell.
(Note that an EQUILIBRIUM_PHASES definition for a cell could be defined by other
calls to @a RM_InitialEquilibriumPhases2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT    0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
equilibrium_phases = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) equilibrium_phases[i] = 1;
status = RM_InitialEquilibriumPhases2Module(id, equilibrium_phases);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialEquilibriumPhases2Module(int id, int* equilibrium_phases);
/**
Transfer EXCHANGE definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a exchanges is used to select EXCHANGE definitions for each cell of the model.
@a exchanges is dimensioned @a nxyz, where @a nxyz is the number of grid cells
in the user's model (@ref RM_GetGridCellCount).
@param id           The instance @a id returned from @ref RM_Create.
@param exchanges    Vector of EXCHANGE index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values
are ignored, resulting in no transfer of an EXCHANGE definition for that cell.
(Note that an EXCHANGE definition for a cell could be defined by other
calls to @a RM_InitialExchanges2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
exchanges = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) exchanges[i] = 1;
status = RM_InitialExchanges2Module(id, exchanges);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialExchanges2Module(int id, int* exchanges);
/**
Transfer SURFACE definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a surfaces is used to select SURFACE definitions for each cell of the model.
@a surfaces is dimensioned @a nxyz, where @a nxyz is the number of grid cells
in the user's model (@ref RM_GetGridCellCount).
@param id          The instance @a id returned from @ref RM_Create.
@param surfaces    Array of SURFACE index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values
are ignored, resulting in no transfer of a SURFACE definition for that cell.
(Note that an SURFACE definition for a cell could be defined by other
calls to @a RM_InitialSurfaces2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
surfaces = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) surfaces[i] = 1;
status = RM_InitialSurfaces2Module(id, surfaces);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialSurfaces2Module(int id, int* surfaces);
/**
Transfer GAS_PHASE definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a gas_phases is used to select GAS_PHASE definitions for each cell of the model.
@a gas_phases is dimensioned @a nxyz, where @a nxyz is the number of grid cells
in the user's model (@ref RM_GetGridCellCount).
@param id           The instance @a id returned from @ref RM_Create.
@param gas_phases   Vector of GAS_PHASE index numbers that refer to
definitions in the InitialPhreeqc instance.Size is @a nxyz. Negative values are
ignored, resulting in no transfer of a GAS_PHASE definition for that cell.
(Note that an GAS_PHASE definition for a cell could be defined by other
calls to @a RM_InitialGasPhases2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
gas_phases = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) gas_phases[i] = 1;
status = RM_InitialGasPhases2Module(id, gas_phases);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialGasPhases2Module(int id, int* gas_phases);
/**
Transfer SOLID_SOLUTIONS definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a solid_solutions is used to select SOLID_SOLUTIONS definitions for each cell
of the model. @a solid_solutions is dimensioned @a nxyz, where @a nxyz is the
number of grid cells in the user's model (@ref RM_GetGridCellCount).
@param id              The instance @a id returned from @ref RM_Create.
@param solid_solutions Array of SOLID_SOLUTIONS index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values
are ignored, resulting in no transfer of a SOLID_SOLUTIONS definition for that cell.
(Note that an SOLID_SOLUTIONS definition for a cell could be defined by other
calls to @a RM_InitialSolidSolutions2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialKinetics2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
solid_solutions = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) solid_solutions[i] = 1;
status = RM_InitialSolidSolutions2Module(id, solid_solutions);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialSolidSolutions2Module(int id, int* solid_solutions);
/**
Transfer KINETICS definitions from the InitialPhreeqc instance to the
reaction-module workers.
@a kinetics is used to select KINETICS definitions for each cell of the model.
@a kinetics is dimensioned @a nxyz, where @a nxyz is the number of grid cells in the
user's model (@ref RM_GetGridCellCount).
@param id          The instance @a id returned from @ref RM_Create.
@param kinetics    Array of KINETICS index numbers that refer to
definitions in the InitialPhreeqc instance. Size is @a nxyz. Negative values are
ignored, resulting in no transfer of a KINETICS definition for that cell.
(Note that an KINETICS definition for a cell could be defined by other
calls to @a RM_InitialKinetics2Module, @ref RM_InitialPhreeqc2Module, or
@ref RM_InitialPhreeqcCell2Module.)
@retval IRM_RESULT  0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_InitialSolutions2Module,
@ref RM_InitialEquilibriumPhases2Module,
@ref RM_InitialExchanges2Module,
@ref RM_InitialGasPhases2Module,
@ref RM_InitialSolidSolutions2Module,
@ref RM_InitialSurfaces2Module,
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
kinetics = (double*)malloc(nxyz*sizeof(double));
for (i=0; i < nxyz; i++) kinetics[i] = 1;
status = RM_InitialKinetics2Module(id, kinetics);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
*/
IRM_DLL_EXPORT  IRM_RESULT RM_InitialKinetics2Module(int id, int* kinetics);
/**
Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers, possibly with mixing.
In its simplest form, @a initial_conditions1 is used to select initial conditions, including solutions and reactants,
for each cell of the model, without mixing.
@a Initial_conditions1 is dimensioned (@a nxyz, 7), where @a nxyz is the number of grid cells in the user's model
(@ref RM_GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
(1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, (4) SURFACE, (5) GAS_PHASE,
(6) SOLID_SOLUTIONS, and (7) KINETICS. In C, initial_solution1[3*nxyz + 99] = 2, indicates that
cell 99 (0 based) contains the SURFACE definition with user number 2 that has been defined in the
InitialPhreeqc instance (either by @ref RM_RunFile or @ref RM_RunString).
@n@n
It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing,
@a initials_conditions2 contains numbers for a second entity that mixes with the entity defined in @a initial_conditions1.
@a Fraction1 contains the mixing fraction for @a initial_conditions1, whereas (1 - @a fraction1) is the mixing fraction for
@a initial_conditions2.
In C, initial_solution1[3*nxyz + 99] = 2, initial_solution2[3*nxyz + 99] = 3,
fraction1[3*nxyz + 99] = 0.25 indicates that
cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface compositions have been defined in the
InitialPhreeqc instance. If the user number in @a initial_conditions2 is negative, no mixing occurs.
If @a initials_conditions2 and @a fraction1 are NULL,
no mixing is used, and initial conditions are derived solely from @a initials_conditions1.

@param id                  The instance @a id returned from @ref RM_Create.
@param initial_conditions1 Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance.
Size is (@a nxyz,7). The order of definitions is given above.
Negative values are ignored, resulting in no definition of that entity for that cell.
@param initial_conditions2  Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance.
Nonnegative values of @a initial_conditions2 result in mixing with the entities defined in @a initial_conditions1.
Negative values result in no mixing.
Size is (@a nxyz,7). The order of definitions is given above.
May be NULL in C; setting to NULL results in no mixing.
@param fraction1           Fraction of initial_conditions1 that mixes with (1-@a fraction1) of initial_conditions2.
Size is (nxyz,7). The order of definitions is given above.
May be NULL in C; setting to NULL results in no mixing.
@retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
@see                        
@ref RM_InitialPhreeqcCell2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
ic1 = (int *) malloc((size_t) (7 * nxyz * sizeof(int)));
ic2 = (int *) malloc((size_t) (7 * nxyz * sizeof(int)));
f1 = (double *) malloc((size_t) (7 * nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
  ic1[i]          = 1;       // Solution 1
  ic1[nxyz + i]   = -1;      // Equilibrium phases none
  ic1[2*nxyz + i] = 1;       // Exchange 1
  ic1[3*nxyz + i] = -1;      // Surface none
  ic1[4*nxyz + i] = -1;      // Gas phase none
  ic1[5*nxyz + i] = -1;      // Solid solutions none
  ic1[6*nxyz + i] = -1;      // Kinetics none

  ic2[i]          = -1;      // Solution none
  ic2[nxyz + i]   = -1;      // Equilibrium phases none
  ic2[2*nxyz + i] = -1;      // Exchange none
  ic2[3*nxyz + i] = -1;      // Surface none
  ic2[4*nxyz + i] = -1;      // Gas phase none
  ic2[5*nxyz + i] = -1;      // Solid solutions none
  ic2[6*nxyz + i] = -1;      // Kinetics none

  f1[i]          = 1.0;      // Mixing fraction ic1 Solution
  f1[nxyz + i]   = 1.0;      // Mixing fraction ic1 Equilibrium phases
  f1[2*nxyz + i] = 1.0;      // Mixing fraction ic1 Exchange 1
  f1[3*nxyz + i] = 1.0;      // Mixing fraction ic1 Surface
  f1[4*nxyz + i] = 1.0;      // Mixing fraction ic1 Gas phase
  f1[5*nxyz + i] = 1.0;      // Mixing fraction ic1 Solid solutions
  f1[6*nxyz + i] = 1.0;      // Mixing fraction ic1 Kinetics
}
status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1);
// No mixing is defined, so the following is equivalent
status = RM_InitialPhreeqc2Module(id, ic1, NULL, NULL);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_InitialPhreeqc2Module(int id,
                int *initial_conditions1,		// 7 x nxyz end-member 1
                int *initial_conditions2,		// 7 x nxyz end-member 2
                double *fraction1);			    // 7 x nxyz fraction of end-member 1

/**
Fills an array (@a species_c) with aqueous species concentrations from solutions in the InitialPhreeqc instance.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
If @a boundary_solution2 and @a fraction1 are NULL,
no mixing is used; concentrations are derived from @a boundary_solution1 only.

@param id                  The instance @a id returned from @ref RM_Create.
@param species_c           Array of aqueous concentrations extracted from the InitialPhreeqc instance.
The dimension of @a species_c is @a n_boundary * @a nspecies,
where @a nspecies is the number of aqueous species returned from @ref RM_GetSpeciesCount.
@param n_boundary          The number of boundary condition solutions that need to be filled.
@param boundary_solution1  Array of solution index numbers that refer to solutions in the InitialPhreeqc instance.
Size is @a n_boundary.
@param boundary_solution2  Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance
and are defined to mix with @a boundary_solution1.
Size is @a n_boundary. May be NULL in C.
@param fraction1           Fraction of @a boundary_solution1 that mixes with (1-@a fraction1) of @a boundary_solution2.
Size is @a n_boundary. May be NULL in C.
@retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
@see                  
@ref RM_FindComponents, 
@ref RM_GetSpeciesCount, 
@ref RM_SetSpeciesSaveOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
nbound = 1;
nspecies = RM_GetSpeciesCount(id);
bc1 = (int *) malloc((size_t) (nbound * sizeof(int)));
bc2 = (int *) malloc((size_t) (nbound * sizeof(int)));
bc_f1 = (double *) malloc((size_t) (nbound * sizeof(double)));
bc_conc = (double *) malloc((size_t) (nspecies * nbound * sizeof(double)));
for (i = 0; i < nbound; i++)
{
  bc1[i]          = 0;       // Solution 0 from InitialPhreeqc instance
  bc2[i]          = -1;      // no bc2 solution for mixing
  bc_f1[i]        = 1.0;     // mixing fraction for bc1
}
status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_InitialPhreeqc2SpeciesConcentrations(
                int id,
                double *species_c,
                int n_boundary,
                int *boundary_solution1,
                int *boundary_solution2,
                double *fraction1);
/**
A cell numbered @a n in the InitialPhreeqc instance is selected to populate a series of cells.
All reactants with the number @a n are transferred along with the solution.
If MIX @a n exists, it is used for the definition of the solution.
If @a n is negative, @a n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
All reactants for each cell in the list @a module_numbers are removed before the cell
definition is copied from the InitialPhreeqc instance to the workers.
@param id                 The instance @a id returned from @ref RM_Create.
@param n                  Cell number refers to a solution or MIX and associated reactants in the InitialPhreeqc instance.
A negative number indicates the largest solution or MIX number in the InitialPhreeqc instance will be used.
@param module_numbers     A list of cell numbers in the user's grid-cell numbering system that will be populated with
cell @a n from the InitialPhreeqc instance.
@param dim_module_numbers The number of cell numbers in the @a module_numbers list.
@retval IRM_RESULT        0 is success, negative is failure (See @ref RM_DecodeError).
@see                      @ref RM_InitialPhreeqc2Module.
@par C Example:
@htmlonly
<CODE>
<PRE>
module_cells = (int *) malloc((size_t) (2 * sizeof(int)));
module_cells[0] = 18;
module_cells[1] = 19;
// n will be the largest SOLUTION number in InitialPhreeqc instance
// copies solution and reactants to cells 18 and 19
status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_InitialPhreeqcCell2Module(int id,
                int n,		                            // InitialPhreeqc cell number
                int *module_numbers,		            // Module cell numbers
                int dim_module_numbers);			    // Number of module cell numbers

/**
Load a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
@param id               The instance @a id returned from @ref RM_Create.
@param db_name          String containing the database name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_Create.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_LoadDatabase(id, "phreeqc.dat");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */

IRM_DLL_EXPORT IRM_RESULT RM_LoadDatabase(int id, const char *db_name);
/**
Print a message to the log file.
@param id               The instance @a id returned from @ref RM_Create.
@param str              String to be printed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see              
@ref RM_ErrorMessage,       
@ref RM_OpenFiles, 
@ref RM_OutputMessage, 
@ref RM_ScreenMessage, 
@ref RM_WarningMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ",
        RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
status = RM_LogMessage(id, str);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */

IRM_DLL_EXPORT IRM_RESULT RM_LogMessage(int id, const char *str);
/**
MPI only. Workers (processes with @ref RM_GetMpiMyself > 0) must call RM_MpiWorker to be able to
respond to messages from the root to accept data, perform calculations, and
(or) return data. RM_MpiWorker contains a loop that reads a message from root, performs a
task, and waits for another message from root. @ref RM_SetConcentrations, @ref RM_RunCells, and @ref RM_GetConcentrations
are examples of methods that send a message from root to get the workers to perform a task. The workers will
respond to all methods that are designated "workers must be in the loop of RM_MpiWorker" in the
MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls
@ref RM_MpiWorkerBreak.
@n@n
(Advanced) The list of tasks that the workers perform can be extended by using @ref RM_SetMpiWorkerCallback.
It is then possible to use the MPI processes to perform other developer-defined tasks, such as transport calculations, without
exiting from the RM_MpiWorker loop. Alternatively, root calls @ref RM_MpiWorkerBreak to allow the workers to continue
past a call to RM_MpiWorker. The workers perform developer-defined calculations, and then RM_MpiWorker is called again to respond to
requests from root to perform reaction-module tasks.

@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError). RM_MpiWorker returns a value only when
@ref RM_MpiWorkerBreak is called by root.
@see                    
@ref RM_MpiWorkerBreak, 
@ref RM_SetMpiWorkerCallback, 
@ref RM_SetMpiWorkerCallbackCookie.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_MpiWorker(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by all workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_MpiWorker(int id);
/**
MPI only. This method is called by root to force workers (processes with @ref RM_GetMpiMyself > 0)
to return from a call to @ref RM_MpiWorker.
@ref RM_MpiWorker contains a loop that reads a message from root, performs a
task, and waits for another message from root. The workers respond to all methods that are designated
"workers must be in the loop of RM_MpiWorker" in the
MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls RM_MpiWorkerBreak.
@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_MpiWorker, 
@ref RM_SetMpiWorkerCallback, 
@ref RM_SetMpiWorkerCallbackCookie.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_MpiWorkerBreak(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_MpiWorkerBreak(int id);
/**
Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
based on the prefix defined by @ref RM_SetFilePrefix.
@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see   
@ref RM_CloseFiles,  
@ref RM_ErrorMessage, 
@ref RM_GetFilePrefix, 
@ref RM_LogMessage, 
@ref RM_OutputMessage,                
@ref RM_SetFilePrefix, 
@ref RM_WarningMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetFilePrefix(id, "Advect_c");
status = RM_OpenFiles(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_OpenFiles(int id);
/**
Print a message to the output file.
@param id               The instance @a id returned from @ref RM_Create.
@param str              String to be printed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_ErrorMessage, 
@ref RM_LogMessage, 
@ref RM_ScreenMessage, 
@ref RM_WarningMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str1, "Number of threads:                                %d\n", RM_GetThreadCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of MPI processes:                          %d\n", RM_GetMpiTasks(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "MPI task number:                                  %d\n", RM_GetMpiMyself(id));
status = RM_OutputMessage(id, str1);
status = RM_GetFilePrefix(id, str, 100);
sprintf(str1, "File prefix:                                      %s\n", str);
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of grid cells in the user's model:         %d\n", RM_GetGridCellCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of chemistry cells in the reaction module: %d\n", RM_GetChemistryCellCount(id));
status = RM_OutputMessage(id, str1);
sprintf(str1, "Number of components for transport:               %d\n", RM_GetComponentCount(id));
status = RM_OutputMessage(id, str1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_OutputMessage(int id, const char *str);
/**
Runs a reaction step for all of the cells in the reaction module.
Normally, tranport concentrations are transferred to the reaction cells (@ref RM_SetConcentrations) before
reaction calculations are run. The length of time over which kinetic reactions are integrated is set
by @ref RM_SetTimeStep. Other properties that may need to be updated as a result of the transport
calculations include porosity (@ref RM_SetPorosity), saturation (@ref RM_SetSaturationUser),
temperature (@ref RM_SetTemperature), and pressure (@ref RM_SetPressure).
@param id               The instance @a id returned from @ref RM_Create.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetConcentrations,  
@ref RM_SetPorosity,
@ref RM_SetPressure, 
@ref RM_SetSaturationUser, 
@ref RM_SetTemperature, 
@ref RM_SetTimeStep.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetPorosity(id, por);              // If porosity changes
status = RM_SetSaturationUser(id, sat);            // If saturation changes
status = RM_SetTemperature(id, temperature);   // If temperature changes
status = RM_SetPressure(id, pressure);         // If pressure changes
status = RM_SetConcentrations(id, c);          // Transported concentrations
status = RM_SetTimeStep(id, time_step);        // Time step for kinetic reactions
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);          // Concentrations after reaction
status = RM_GetDensityCalculated(id, density);           // Density after reaction
status = RM_GetSolutionVolume(id, volume);     // Solution volume after reaction
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_RunCells(int id);
/**
Run a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
files that modify the thermodynamic database should be run by all three sets of instances.
Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Files that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
@param id               The instance @a id returned from @ref RM_Create.
@param workers          1, the workers will run the file; 0, the workers will not run the file.
@param initial_phreeqc  1, the InitialPhreeqc instance will run the file; 0, the InitialPhreeqc will not run the file.
@param utility          1, the Utility instance will run the file; 0, the Utility instance will not run the file.
@param chem_name        Name of the file to run.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_RunString.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_RunFile(id, 1, 1, 1, "advect.pqi");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT        RM_RunFile(int id, int workers, int initial_phreeqc, int utility, const char *chem_name);
/**
Run a PHREEQC input string. The first three arguments determine which
IPhreeqc instances will run
the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
strings that modify the thermodynamic database should be run by all three sets of instances.
Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Strings that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
@param id               The instance @a id returned from @ref RM_Create.
@param workers          1, the workers will run the string; 0, the workers will not run the string.
@param initial_phreeqc  1, the InitialPhreeqc instance will run the string; 0, the InitialPhreeqc will not run the string.
@param utility          1, the Utility instance will run the string; 0, the Utility instance will not run the string.
@param input_string     String containing PHREEQC input.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_RunFile.
@par C Example:
@htmlonly
<CODE>
<PRE>
Utilities::strcpy_safe(str, MAX_LENGTH, "DELETE; -all");
status = RM_RunString(id, 1, 0, 1, str);	// workers, initial_phreeqc, utility
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_RunString(int id, int workers, int initial_phreeqc, int utility, const char * input_string);
/**
Print message to the screen.
@param id               The instance @a id returned from @ref RM_Create.
@param str              String to be printed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_ErrorMessage,  
@ref RM_LogMessage, 
@ref RM_OutputMessage, 
@ref RM_WarningMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ",
        time * RM_GetTimeConversion(id), " days\n");
status = RM_ScreenMessage(id, str);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_ScreenMessage(int id, const char *str);
/**
Select whether to include H2O in the component list.
The concentrations of H and O must be known
accurately (8 to 10 significant digits) for the numerical method of
PHREEQC to produce accurate pH and pe values.
Because most of the H and O are in the water species,
it may be more robust (require less accuracy in transport) to
transport the excess H and O (the H and O not in water) and water.
The default setting (@a true) is to include water, excess H, and excess O as components.
A setting of @a false will include total H and total O as components.
@a RM_SetComponentH2O must be called before @ref RM_FindComponents.

@param id               The instance id returned from @ref RM_Create.
@param tf               0, total H and O are included in the component list; 1, excess H, excess O, and water
are included in the component list.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetComponentH2O(id, 0);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetComponentH2O(int id, int tf);
/**
Use the vector of concentrations (@a c) to set the moles of components in each reaction cell.
The volume of water in a cell is the product of porosity (@ref RM_SetPorosity), 
saturation (@ref RM_SetSaturationUser), and reference volume (@ref RM_SetRepresentativeVolume).
The moles of each component are determined by the volume of water and per liter concentrations.
If concentration units (@ref RM_SetUnitsSolution) are mass fraction, the
density (as specified by @ref RM_SetDensityUser) is used to convert from mass fraction to per mass per liter.

@param id               The instance @a id returned from @ref RM_Create.
@param c                Array of component concentrations. Size of array is @a nxyz * @a ncomps,
where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount), and @a ncomps is the number of components as determined
by @ref RM_FindComponents or @ref RM_GetComponentCount.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetDensityUser, 
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume,
@ref RM_SetSaturationUser, 
@ref RM_SetUnitsSolution.

@par C Example:
@htmlonly
<CODE>
<PRE>
c = (double *) malloc((size_t) (ncomps * nxyz * sizeof(double)));
...
advect_c(c, bc_conc, ncomps, nxyz, nbound);
status = RM_SetPorsity(id, por);               // If porosity changes
status = RM_SetSaturationUser(id, sat);            // If saturation changes
status = RM_SetTemperature(id, temperature);   // If temperature changes
status = RM_SetPressure(id, pressure);         // If pressure changes
status = RM_SetConcentrations(id, c);          // Transported concentrations
status = RM_SetTimeStep(id, time_step);        // Time step for kinetic reactions
status = RM_SetTime(id, time);                 // Current time
status = RM_RunCells(id);
status = RM_GetConcentrations(id, c);          // Concentrations after reaction
status = RM_GetDensityCalculated(id, density);           // Density after reaction
status = RM_GetSolutionVolume(id, volume);     // Solution volume after reaction
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetConcentrations(int id, double *c);
/**
Select the current selected output by user number. The user may define multiple SELECTED_OUTPUT
data blocks for the workers. A user number is specified for each data block. The value of
the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
for selected-output operations.
@param id               The instance @a id returned from @ref RM_Create.
@param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetNthSelectedOutput,
@ref RM_SetSelectedOutputOn.
@par C Example:
@htmlonly
<CODE>
<PRE>
for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
{
  n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
  status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
  col = RM_GetSelectedOutputColumnCount(id);
  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
  status = RM_GetSelectedOutput(id, selected_out);
  // Process results here
  free(selected_out);
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetCurrentSelectedOutputUserNumber(int id, int n_user);
/**
Set the density for each reaction cell. These density values are used
when converting from transported mass fraction concentrations (@ref RM_SetUnitsSolution) to
produce per liter concentrations during a call to @ref RM_SetConcentrations.
They are also used when converting from module concentrations to transport concentrations
of mass fraction (@ref RM_GetConcentrations), if @ref RM_UseSolutionDensityVolume is set to @a false.

@param id               The instance @a id returned from @ref RM_Create.
@param density          Array of densities. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_GetConcentrations, 
@ref RM_SetConcentrations,
@ref RM_SetUnitsSolution, 
@ref RM_UseSolutionDensityVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
density = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
	density[i] = 1.0;
}
status = RM_SetDensityUser(id, density);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetDensityUser(int id, double *density);
/**
Deprecated equivalent of RM_SetDensityUser.
*/
IRM_DLL_EXPORT IRM_RESULT RM_SetDensity(int id, double* density);

/**
Set the name of the dump file. It is the name used by @ref RM_DumpModule.
@param id               The instance @a id returned from @ref RM_Create.
@param dump_name        Name of dump file.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_DumpModule.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetDumpFileName(id, "advection_c.dmp");
dump_on = 1;
append = 0;
status = RM_DumpModule(id, dump_on, append);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetDumpFileName(int id, const char *dump_name);
/**
Set the action to be taken when the reaction module encounters an error.
Options are 0, return to calling program with an error return code (default);
1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
2, attempt to exit gracefully.
@param id               The instance id returned from @ref RM_Create.
@param mode             Error handling mode: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@par C Example:
@htmlonly
<CODE>
<PRE>
id = RM_Create(nxyz, nthreads);
status = RM_SetErrorHandlerMode(id, 2);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetErrorHandlerMode(int id, int mode);
/**
Set the property that controls whether error messages are generated and displayed.
Messages include PHREEQC "ERROR" messages, and
any messages written with @ref RM_ErrorMessage.

@param id               The instance @a id returned from @ref RM_Create.
@param tf  @a 1, enable error messages; @a 0, disable error messages. Default is 1.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_ErrorMessage,
@ref RM_ScreenMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetErrorOn(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetErrorOn(int id, int tf);
/**
Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
These files are opened by @ref RM_OpenFiles.
@param id               The instance @a id returned from @ref RM_Create.
@param prefix           Prefix used when opening the output and log files.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see  
@ref RM_CloseFiles,                  
@ref RM_OpenFiles. 
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetFilePrefix(id, "Advect_c");
status = RM_OpenFiles(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetFilePrefix(int id, const char *prefix);

/**
Transfer moles of gas components from
the vector given in the argument list (@a gas_moles) to each reaction cell.

@param id               The instance @a id returned from @ref RM_Create.
@param gas_moles        Vector of moles of gas components.
Dimension of the vector must be @a ngas_comps times @a nxyz,
where, @a ngas_comps is the result of @ref RM_GetGasComponentsCount,
and @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If the number of moles is set to a negative number, the gas component will
not be defined for the GAS_PHASE of the reaction cell.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetGasComponentsCount, 
@ref RM_GetGasCompMoles, 
@ref RM_GetGasCompPressures,
@ref RM_GetGasPhaseVolume,
@ref RM_GetGasCompPhi,
@ref RM_SetGasPhaseVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
ngas_comps = RM_GetGasComponentsCount();
gas_moles = (double *) malloc((size_t) (ngas_comps * nxyz * sizeof(double)));
...
status = RM_SetGasCompMoles(id, gas_moles);
status = RM_RunCells(id)
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetGasCompMoles(int id, double* gas_moles);

/**
Transfer volumes of gas phases from
the array given in the argument list (@a gas_volume) to each reaction cell.
The gas-phase volume affects the pressures calculated for fixed-volume
gas phases. If a gas-phase volume is defined with this method 
for a GAS_PHASE in a cell, 
the gas phase is forced to be a fixed-volume gas phase.

@param id               The instance @a id returned from @ref RM_Create.
@param  gas_volume               Vector of volumes for each gas phase.
Dimension of the vector must be @a nxyz,
where, @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
If the volume is set to a negative number for a cell, the gas-phase volume for that cell is
not changed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_FindComponents,
@ref RM_GetGasComponentsCount,
@ref RM_GetGasCompMoles,
@ref RM_GetGasCompPressures,
@ref RM_GetGasPhaseVolume,
@ref RM_GetGasCompPhi,
@ref RM_SetGasCompMoles.

@par C Example:
@htmlonly
<CODE>
<PRE>
gas_volume = (double *) malloc((size_t) (nxyz * sizeof(double)));
...
status = RM_SetGasPhaseVolume(id, gas_moles);
status = RM_RunCells(id)
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetGasPhaseVolume(int id, double* gas_volume);


/**
MPI only. Defines a callback function that allows additional tasks to be done
by the workers. The method @ref RM_MpiWorker contains a loop,
where the workers receive a message (an integer),
run a function corresponding to that integer,
and then wait for another message.
RM_SetMpiWorkerCallback allows the developer to add another function
that responds to additional integer messages by calling developer-defined functions
corresponding to those integers.
@ref RM_MpiWorker calls the callback function when the message number
is not one of the PhreeqcRM message numbers.
Messages are unique integer numbers. PhreeqcRM uses integers in a range
beginning at 0. It is suggested that developers use message numbers starting
at 1000 or higher for their tasks.
The callback function calls a developer-defined function specified
by the message number and then returns to @ref RM_MpiWorker to wait for
another message.
@n@n
In C, an additional pointer can be supplied to find the data necessary to do the task.
A void pointer may be set with @ref RM_SetMpiWorkerCallbackCookie. This pointer
is passed to the callback function through a void pointer argument in addition
to the integer message argument. The void pointer may be a pointer to a struct that
contains pointers to additional data. @ref RM_SetMpiWorkerCallbackCookie
must be called by each worker before @ref RM_MpiWorker is called.
@n@n
The motivation for this method is to allow the workers to perform other
tasks, for instance, parallel transport calculations, within the structure
of @ref RM_MpiWorker. The callback function
can be used to allow the workers to receive data, perform transport calculations,
and (or) send results, without leaving the loop of @ref RM_MpiWorker. Alternatively,
it is possible for the workers to return from @ref RM_MpiWorker
by a call to @ref RM_MpiWorkerBreak by root. The workers could then call
subroutines to receive data, calculate transport, and send data,
and then resume processing PhreeqcRM messages from root with another
call to @ref RM_MpiWorker.
@param id               The instance @a id returned from @ref RM_Create.
@param fcn              A function that returns an integer and has an integer argument.
C has an additional void * argument.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_MpiWorker, 
@ref RM_MpiWorkerBreak,
@ref RM_SetMpiWorkerCallbackCookie.
@par C Example:
@htmlonly
<CODE>
<PRE>
Code executed by root:
// root calls a function that will involve the workers
status = do_something(&comm);

Code executed by workers:
status = RM_SetMpiWorkerCallback(id, worker_tasks_c);
status = RM_SetMpiWorkerCallbackCookie(id, &comm);
status = RM_MpiWorker(id);

Code executed by root and workers:
int do_something(void *cookie)
{
	MPI_Status status;
	MPI_Comm *comm = (MPI_Comm *) cookie;
	int i, method_number, mpi_myself, mpi_tasks, worker_number;
	method_number = 1000;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, *comm);
		fprintf(stderr, "I am root.\n");
		for (i = 1; i < mpi_tasks; i++)
		{
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, *comm, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, *comm);
	}
	return 0;
}

Code called by workers from method MpiWorker:
int worker_tasks_c(int *method_number, void * cookie)
{
	if (*method_number == 1000)
	{
		do_something(cookie);
	}
	return 0;
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by workers, before call to @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetMpiWorkerCallback(int id, int (*fcn)(int *x1, void *cookie));
/**
MPI and C only. Defines a void pointer that can be used by
C functions called from the callback function (@ref RM_SetMpiWorkerCallback)
to locate data for a task. The C callback function
that is registered with @ref RM_SetMpiWorkerCallback has
two arguments, an integer message to identify a task, and a void
pointer. RM_SetMpiWorkerCallbackCookie sets the value of the
void pointer that is passed to the callback function.
@param id               The instance id returned from @ref RM_Create.
@param cookie           Void pointer that can be used by subroutines called from the callback function
to locate data needed to perform a task.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_MpiWorker, 
@ref RM_MpiWorkerBreak,
@ref RM_SetMpiWorkerCallback.
@par C Example:
@htmlonly
<CODE>
<PRE>
Code executed by root:
// root calls a function that will involve the workers
status = do_something(&comm);

Code executed by workers:
status = RM_SetMpiWorkerCallback(id, worker_tasks_c);
status = RM_SetMpiWorkerCallbackCookie(id, &comm);
status = RM_MpiWorker(id);

Code executed by root and workers:
int do_something(void *cookie)
{
	MPI_Status status;
	MPI_Comm *comm = (MPI_Comm *) cookie;
	int i, method_number, mpi_myself, mpi_tasks, worker_number;
	method_number = 1000;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, *comm);
		fprintf(stderr, "I am root.\n");
		for (i = 1; i < mpi_tasks; i++)
		{
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, *comm, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, *comm);
	}
	return 0;
}

Code called by workers from method MpiWorker:
int worker_tasks_c(int *method_number, void * cookie)
{
	if (*method_number == 1000)
	{
		do_something(cookie);
	}
	return 0;
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by workers, before call to @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetMpiWorkerCallbackCookie(int id, void *cookie);
/**
Specify the current selected output by sequence number. The user may define multiple SELECTED_OUTPUT
data blocks for the workers. A user number is specified for each data block, and the blocks are
stored in user-number order. The value of
the argument @a n selects the sequence number of the SELECTED_OUTPUT definition that will be used
for selected-output operations.
@param id               The instance id returned from @ref RM_Create.
@param n           Sequence number of the SELECTED_OUTPUT data block that is to be used.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetSelectedOutputOn.
	@par C Example:
	@htmlonly
	<CODE>
	<PRE>
	for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
	{
	  status = RM_SetNthSelectedOutputUser(id, isel);
	  n_user = RM_GetCurrentSelectedOutputUserNumber(id);
	  col = RM_GetSelectedOutputColumnCount(id);
	  selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
	  status = RM_GetSelectedOutput(id, selected_out);
	  // Process results here
	  free(selected_out);
	}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
*/
IRM_DLL_EXPORT IRM_RESULT RM_SetNthSelectedOutput(int id, int n);
/**
Sets the property for partitioning solids between the saturated and unsaturated
parts of a partially saturated cell.

The option is intended to be used by saturated-only
flow codes that allow a variable water table.
The value has meaning only when saturations
less than 1.0 are encountered. The partially saturated cells
may have a small water-to-rock ratio that causes
reactions to proceed differently relative to fully saturated cells.
By setting  @a RM_SetPartitionUZSolids to true, the
amounts of solids and gases are partioned according to the saturation.
If a cell has a saturation of 0.5, then
the water interacts with only half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell varies,
solids and gases are transferred between the
saturated and unsaturated (unreactive) reservoirs of the cell.
Unsaturated-zone flow and transport codes will probably use the default (false),
which assumes all gases and solids are reactive regardless of saturation.

@param id       The instance @a id returned from @ref RM_Create.
@param tf       @a True, the fraction of solids and gases available for
reaction is equal to the saturation;
@a False (default), all solids and gases are reactive regardless of saturation.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetPartitionUZSolids(id, 0);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetPartitionUZSolids(int id, int tf);
/**
Set the porosity for each reaction cell.
The volume of water in a reaction cell is the product of the porosity, the saturation
(@ref RM_SetSaturationUser), and the representative volume (@ref RM_SetRepresentativeVolume).
@param id               The instance @a id returned from @ref RM_Create.
@param por              Array of porosities, unitless. Default is 0.1. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_GetSaturationCalculated, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser.
@par C Example:
@htmlonly
<CODE>
<PRE>
por = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) por[i] = 0.2;
status = RM_SetPorosity(id, por);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetPorosity(int id, double *por);
/**
Set the pressure for each reaction cell. Pressure effects are considered only in three of the
databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
@param id               The instance @a id returned from @ref RM_Create.
@param p                Array of pressures, in atm. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetTemperature.
@par C Example:
@htmlonly
<CODE>
<PRE>
pressure = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) pressure[i] = 2.0;
status = RM_SetPressure(id, pressure);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetPressure(int id, double *p);
/**
Enable or disable detailed output for each reaction cell.
Printing for a cell will occur only when the
printing is enabled with @ref RM_SetPrintChemistryOn and the @a cell_mask value is 1.

@param id               The instance @a id returned from @ref RM_Create.
@param cell_mask        Array of integers. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount). A value of 0 will
disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetPrintChemistryOn.

@par C Example:
@htmlonly
<CODE>
<PRE>
print_chemistry_mask = (int *) malloc((size_t) (nxyz * sizeof(int)));
for (i = 0; i < nxyz/2; i++)
{
  print_chemistry_mask[i] = 1;
  print_chemistry_mask[i + nxyz/2] = 0;
}
status = RM_SetPrintChemistryMask(id, print_chemistry_mask);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetPrintChemistryMask(int id, int *cell_mask);
/**
Setting to enable or disable printing detailed output from reaction calculations to the output file for a set of
cells defined by @ref RM_SetPrintChemistryMask. The detailed output prints all of the output
typical of a PHREEQC reaction calculation, which includes solution descriptions and the compositions of
all other reactants. The output can be several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, @ref RM_OpenFiles). For the worker instances, the output can be limited to a set of cells
(@ref RM_SetPrintChemistryMask) and, in general, the
amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
(applied by using @ref RM_RunFile or @ref RM_RunString).
Printing the detailed output for the workers is generally used only for debugging, and PhreeqcRM will run
significantly faster when printing detailed output for the workers is disabled.

@param id               The instance @a id returned from @ref RM_Create.
@param workers          0, disable detailed printing in the worker instances, 1, enable detailed printing
in the worker instances.
@param initial_phreeqc  0, disable detailed printing in the InitialPhreeqc instance, 1, enable detailed printing
in the InitialPhreeqc instances.
@param utility          0, disable detailed printing in the Utility instance, 1, enable detailed printing
in the Utility instance.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetPrintChemistryMask.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetPrintChemistryOn(id, 0, 1, 0); // workers, initial_phreeqc, utility
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetPrintChemistryOn(int id, int workers, int initial_phreeqc, int utility);

/**
Set the load-balancing algorithm.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a @ref RM_RunCells
calculation. Two algorithms are available; one uses individual times for each cell and
accounts for cells that were not run because
saturation was zero (default), and
the other assigns an average time to all cells.
The methods are similar, but limited testing indicates the default method performs better.

@param id               The instance @a id returned from @ref RM_Create.
@param method           0, indicates average times are used in rebalancing; 1 indicates individual
cell times are used in rebalancing (default).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetRebalanceFraction.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetRebalanceByCell(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetRebalanceByCell(int id, int method);
/**
Sets the fraction of cells that are transferred among threads or processes when rebalancing.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a @ref RM_RunCells
calculation. The rebalancing transfers cell calculations among threads or processes to
try to achieve an optimum balance. @a RM_SetRebalanceFraction
adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
determine the actual number of cell transfers. A value of zero eliminates
load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
distribution and avoid possible oscillations
when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
Default is 0.5.

@param id               The instance @a id returned from @ref RM_Create.
@param f                Fraction from 0.0 to 1.0.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetRebalanceByCell.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetRebalanceFraction(id, 0.5);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetRebalanceFraction(int id, double f);
/**
Set the representative volume of each reaction cell.
By default the representative volume of each reaction cell is 1 liter.
The volume of water in a reaction cell is determined by the procuct of the representative volume,
the porosity (@ref RM_SetPorosity), and the saturation (@ref RM_SetSaturationUser).
The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
within a couple orders of magnitude of 1.0.
Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes)
may cause non-convergence of the numerical method.
In these cases, a larger representative volume may help. Note
that increasing the representative volume also increases
the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
and others), which are defined as moles per representative volume.

@param id               The instance @a id returned from @ref RM_Create.
@param rv              Vector of representative volumes, in liters. Default is 1.0 liter.
Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetPorosity, 
@ref RM_SetSaturationUser.

@par C Example:
@htmlonly
<CODE>
<PRE>
double * rv;
rv = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) rv[i] = 1.0;
status = RM_SetRepresentativeVolume(id, rv);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetRepresentativeVolume(int id, double *rv);
/**
Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
The volume of water in a cell is the product of porosity (@ref RM_SetPorosity), 
saturation (@a RM_SetSaturationUser),
and representative volume (@ref RM_SetRepresentativeVolume). As a result of a reaction 
calculation, solution properties (density and volume) will change; the databases phreeqc.dat, 
Amm.dat, and pitzer.dat have the molar volume data to calculate these changes. The methods 
@ref RM_GetDensityCalculated, @ref RM_GetSolutionVolume, and @ref RM_GetSaturationCalculated
can be used to account for these changes in the succeeding transport calculation.
@a RM_SetRepresentativeVolume should be called before initial conditions are defined for 
the reaction cells.

@param id               The instance @a id returned from @ref RM_Create.
@param sat              Array of saturations, unitless. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_GetDensityCalculated, 
@ref RM_GetSaturationCalculated, 
@ref RM_GetSolutionVolume,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
sat = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++) sat[i] = 1.0;
status = RM_SetSaturationUser(id, sat);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetSaturationUser(int id, double *sat);
/**
Deprecated equivalent of RM_SetSaturationUser.
*/
IRM_DLL_EXPORT IRM_RESULT RM_SetSaturation(int id, double* sat);
/**
Set the property that controls whether messages are written to the screen.
Messages include information about rebalancing during @ref RM_RunCells, and
any messages written with @ref RM_ScreenMessage.

@param id               The instance @a id returned from @ref RM_Create.
@param tf  @a 1, enable screen messages; @a 0, disable screen messages. Default is 1.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_RunCells, 
@ref RM_ScreenMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetScreenOn(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetScreenOn(int id, int tf);
/**
Setting determines whether selected-output results are available to be retrieved
with @ref RM_GetSelectedOutput. @a 1 indicates that selected-output results
will be accumulated during @ref RM_RunCells and can be retrieved with @ref RM_GetSelectedOutput;
@a 0 indicates that selected-output results will not
be accumulated during @ref RM_RunCells.

@param id               The instance @a id returned from @ref RM_Create.
@param selected_output  0, disable selected output; 1, enable selected output.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see
@ref RM_GetCurrentSelectedOutputUserNumber,
@ref RM_GetNthSelectedOutputUserNumber,
@ref RM_GetSelectedOutput,
@ref RM_GetSelectedOutputColumnCount,
@ref RM_GetSelectedOutputCount,
@ref RM_GetSelectedOutputHeading,
@ref RM_GetSelectedOutputRowCount,
@ref RM_SetCurrentSelectedOutputUserNumber,
@ref RM_SetNthSelectedOutput.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSelectedOutputOn(id, 1);       // enable selected output
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetSelectedOutputOn(int id, int selected_output);
/**
Sets the value of the species-save property.
This method enables use of PhreeqcRM with multicomponent-diffusion transport calculations.
By default, concentrations of aqueous species are not saved. Setting the species-save property to 1 allows
aqueous species concentrations to be retrieved
with @ref RM_GetSpeciesConcentrations, and solution compositions to be set with
@ref RM_SpeciesConcentrations2Module.
RM_SetSpeciesSaveOn must be called before calls to @ref RM_FindComponents.

@param id               The instance @a id returned from @ref RM_Create.
@param save_on          0, indicates species concentrations are not saved; 1, indicates species concentrations are
saved.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName,
@ref RM_GetSpeciesSaveOn, 
@ref RM_GetSpeciesZ, 
@ref RM_SpeciesConcentrations2Module.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetSpeciesSaveOn(int id, int save_on);
/**
Set the temperature for each reaction cell. If @a RM_SetTemperature is not called,
worker solutions will have temperatures as defined by initial conditions
(@ref RM_InitialPhreeqc2Module and @ref RM_InitialPhreeqcCell2Module).

@param id               The instance @a id returned from @ref RM_Create.
@param t                Array of temperatures, in degrees C. Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref RM_GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module,
@ref RM_InitialPhreeqcCell2Module, 
@ref RM_SetPressure.

@par C Example:
@htmlonly
<CODE>
<PRE>
temperature = (double *) malloc((size_t) (nxyz * sizeof(double)));
for (i = 0; i < nxyz; i++)
{
  temperature[i] = 20.0;
}
status = RM_SetTemperature(id, temperature);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetTemperature(int id, double *t);
/**
Set current simulation time for the reaction module.
@param id               The instance @a id returned from @ref RM_Create.
@param time             Current simulation time, in seconds.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see     
@ref RM_SetTimeConversion,               
@ref RM_SetTimeStep. 
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetTime(id, time);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetTime(int id, double time);
/**
Set a factor to convert to user time units. Factor times seconds produces user time units.
@param id               The instance @a id returned from @ref RM_Create.
@param conv_factor      Factor to convert seconds to user time units.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetTime, 
@ref RM_SetTimeStep.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetTimeConversion(id, 1.0 / 86400.0); // days
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetTimeConversion(int id, double conv_factor);
/**
Set current time step for the reaction module. This is the length
of time over which kinetic reactions are integrated.
@param id               The instance @a id returned from @ref RM_Create.
@param time_step        Current time step, in seconds.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetTime, 
@ref RM_SetTimeConversion.
@par C Example:
@htmlonly
<CODE>
<PRE>
time_step = 86400.0;
status = RM_SetTimeStep(id, time_step);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetTimeStep(int id, double time_step);
/**
Sets input units for exchangers.
In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
@a RM_SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.

If a single EXCHANGE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of exchangers will be the same regardless of porosity. 
For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for exchangers: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsExchange(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsExchange(int id, int option);
/**
Set input units for gas phases.
In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
@a RM_SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single GAS_PHASE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a gas component will be the same regardless of porosity. 
For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for gas phases: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsGasPhase(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsGasPhase(int id, int option);
/**
Set input units for kinetic reactants.

In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
@a RM_SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single KINETICS definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of kinetic reactants will be the same regardless of porosity. 
For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.

Note that the volume of water in a cell in the reaction module is equal to the product of
porosity (@ref RM_SetPorosity), the saturation (@ref RM_SetSaturationUser), and representative volume (@ref
RM_SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the RATES
definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
water, often by calculating the rate of reaction per liter of water and multiplying by the volume
of water (Basic function SOLN_VOL). 

Rates that depend on surface area of solids, are not dependent
on the volume of water. However, it is important to get the correct surface area for the kinetic
reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant) 
can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of 
reactant (Basic function M) in RATES to obtain the surface area.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for kinetic reactants: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                     
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsKinetics(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsKinetics(int id, int option);
/**
Set input units for pure phase assemblages (equilibrium phases).
In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
@a RM_SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a mineral will be the same regardless of porosity. 
For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for equilibrium phases: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsPPassemblage(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsPPassemblage(int id, int option);
/**
Solution concentration units used by the transport model.
Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
PHREEQC defines solutions by the number of moles of each
element in the solution.

To convert from mg/L to moles
of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
multiplied by the solution volume,
which is the product of porosity (@ref RM_SetPorosity), saturation (@ref RM_SetSaturationUser),
and representative volume (@ref RM_SetRepresentativeVolume).
To convert from mol/L to moles
of element in the representative volume of a reaction cell, mol/L is
multiplied by the solution volume.
To convert from mass fraction to moles
of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
(@ref RM_SetDensityUser) and
multiplied by the solution volume.

To convert from moles
of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
solution volume resulting in mol/L, and then converted to mg/L.
To convert from moles
of element in a cell to mol/L,  the number of moles of an element is divided by the
solution volume resulting in mol/L.
To convert from moles
of element in a cell to mass fraction, the number of moles of an element is converted to kg and divided
by the total mass of the solution.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of porosity (@ref RM_SetPorosity),
saturation (@ref RM_SetSaturationUser), and representative volume (@ref RM_SetRepresentativeVolume),
and the mass of solution is volume times density as defined by @ref RM_SetDensityUser.
Which option is used is determined by @ref RM_UseSolutionDensityVolume.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_SetDensityUser, 
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser,
@ref RM_UseSolutionDensityVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsSolution(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsSolution(int id, int option);
/**
Set input units for solid-solution assemblages.
In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
@a RM_SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single SOLID_SOLUTION definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a solid-solution component will be the same regardless of porosity. 
For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for solid solutions: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsSSassemblage(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsSSassemblage(int id, int option);
/**
Set input units for surfaces.
In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
@a RM_SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single SURFACE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of surface sites will be the same regardless of porosity. 
For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.

@param id               The instance @a id returned from @ref RM_Create.
@param option           Units option for surfaces: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_InitialPhreeqc2Module, 
@ref RM_InitialPhreeqcCell2Module,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetUnitsSurface(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SetUnitsSurface(int id, int option);
/**
Set solution concentrations in the reaction cells
based on the vector of aqueous species concentrations (@a species_conc).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref RM_SetSpeciesSaveOn must be set to @a true.
The list of aqueous species is determined by @ref RM_FindComponents and includes all
aqueous species that can be made from the set of components.
The method determines the total concentration of a component
by summing the molarities of the individual species times the stoichiometric
coefficient of the element in each species.
Solution compositions in the reaction cells are updated with these component concentrations.

@param id               The instance @a id returned from @ref RM_Create.
@param species_conc     Array of aqueous species concentrations. Dimension of the array is (@a nxyz, @a nspecies),
where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount), and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
Concentrations are moles per liter.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see                    
@ref RM_FindComponents, 
@ref RM_GetSpeciesConcentrations, 
@ref RM_GetSpeciesCount,
@ref RM_GetSpeciesD25, 
@ref RM_GetSpeciesLog10Gammas,
@ref RM_GetSpeciesLog10Molalities,
@ref RM_GetSpeciesName, 
@ref RM_GetSpeciesSaveOn,
@ref RM_GetSpeciesZ, 
@ref RM_SetSpeciesSaveOn.

@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_SetSpeciesSaveOn(id, 1);
ncomps = RM_FindComponents(id);
nspecies = RM_GetSpeciesCount(id);
nxyz = RM_GetGridCellCount(id);
species_c = (double *) malloc((size_t) (nxyz * nspecies * sizeof(double)));
...
status = RM_SpeciesConcentrations2Module(id, species_c);
status = RM_RunCells(id);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_SpeciesConcentrations2Module(int id, double * species_conc);
/**
Save the state of the chemistry in all model cells, including SOLUTIONs,
EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
will be saved for each cell, if they have been defined in the worker IPhreeqc instances.
The distribution of cells among the workersand the chemistry of fully or partially
unsaturated cells are also saved.The state is saved in memory; use @ref RM_DumpModule to save the state
to file.PhreeqcRM can be reset to this state by using @ref RM_StateApply.
A state is identified by an integer, and multiple states can be saved.

@param id               The instance @a id returned from @ref RM_Create.
@param istate     Integer identifying the state that is saved.
@retval IRM_RESULT      0 is success, negative is failure(See @ref RM_DecodeError).
@see                    @ref RM_DumpModule,
@ref RM_StateApply, and
@ref RM_StateDelete.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_StateSave(id, 1);
...
status = RM_StateApply(id, 1);
status = RM_StateDelete(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI :
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_StateSave(int id, int istate);
/**
Reset the state of the module to a state previously saved with @ref RM_StateSave.
The chemistry of all model cells are reset, including SOLUTIONs,
EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
will be reset for each cell, if they were defined in the worker IPhreeqc instances
at the time the state was saved.
The distribution of cells among the workersand the chemistry of fully or partially
unsaturated cells are also reset to the saved state.
The state to be applied is identified by an integer.

@param id               The instance @a id returned from @ref RM_Create.
@param istate     Integer identifying the state that is to be applied.
@retval IRM_RESULT      0 is success, negative is failure(See @ref RM_DecodeError).
@see                    @ref RM_StateSave and
@ref RM_StateDelete.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_StateSave(id, 1);
...
status = RM_StateApply(id, 1);
status = RM_StateDelete(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI :
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_StateApply(int id, int istate);
/**
Delete a state previously saved with @ref RM_StateSave.

@param id               The instance @a id returned from @ref RM_Create.
@param istate     Integer identifying the state that is to be deleted.
@retval IRM_RESULT      0 is success, negative is failure(See @ref RM_DecodeError).
@see                    @ref RM_StateSave and
ref RM_StateApply.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_StateSave(id, 1);
...
status = RM_StateApply(id, 1);
status = RM_StateDelete(id, 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI :
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_StateDelete(int id, int istate);
/**
Determines the volume and density to use when converting from the reaction-module concentrations
to transport concentrations (@ref RM_GetConcentrations).
Two options are available to convert concentration units:
(1) the density and solution volume calculated by PHREEQC are used, or
(2) the specified density (@ref RM_SetDensityUser)
and solution volume are defined by the product of
saturation (@ref RM_SetSaturationUser), porosity (@ref RM_SetPorosity),
and representative volume (@ref RM_SetRepresentativeVolume).
Transport models that consider density-dependent flow will probably use the
PHREEQC-calculated density and solution volume (default),
whereas transport models that assume constant-density flow will probably use
specified values of density and solution volume.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
Density is only used when converting to transport units of mass fraction.

@param id               The instance @a id returned from @ref RM_Create.
@param tf               @a True indicates that the solution density and volume as
calculated by PHREEQC will be used to calculate concentrations.
@a False indicates that the solution density set by @ref RM_SetDensityUser and the volume determined by the
product of  @ref RM_SetSaturationUser, @ref RM_SetPorosity, and @ref RM_SetRepresentativeVolume,
will be used to calculate concentrations retrieved by @ref RM_GetConcentrations.
@see                    
@ref RM_GetConcentrations, 
@ref RM_SetDensityUser,
@ref RM_SetPorosity, 
@ref RM_SetRepresentativeVolume, 
@ref RM_SetSaturationUser.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_UseSolutionDensityVolume(id, 0);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref RM_MpiWorker.
 */
IRM_DLL_EXPORT IRM_RESULT RM_UseSolutionDensityVolume(int id, int tf);
/**
Print a warning message to the screen and the log file.
@param id               The instance @a id returned from @ref RM_Create.
@param warn_str         String to be printed.
@retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
@see   
@ref RM_ErrorMessage,          
@ref RM_LogMessage,        
@ref RM_OpenFiles, 
@ref RM_OutputMessage, 
@ref RM_ScreenMessage.
@par C Example:
@htmlonly
<CODE>
<PRE>
status = RM_WarningMessage(id, "Parameter is out of range, using default");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; only root writes to the log file.
 */
IRM_DLL_EXPORT IRM_RESULT RM_WarningMessage(int id, const char *warn_str);

#if defined(__cplusplus)
}
#endif

#endif // RM_INTERFACE_C_H
