%define BMIPhreeqcRM_v_DOCSTRING
"Default constructor for the BMIPhreeqcRM subclass of PhreeqcRM.

Definition of the number of cells and threads (or MPI communicator) 
is deferred to the method :meth:`initialize`.

Returns:
	 A BMIPhreeqcRM instance."
%enddef
//%feature("docstring", BMIPhreeqcRM_v_DOCSTRING) BMIPhreeqcRM::BMIPhreeqcRM();
%feature("docstring") BMIPhreeqcRM::BMIPhreeqcRM() BMIPhreeqcRM_v_DOCSTRING


// %define BMIPhreeqcRM_ii_DOCSTRING
// "Constructor for the BMIPhreeqcRM subclass of PhreeqcRM. 

// A BMIPhreeqcRM instance has the BMI methods plus all of the
// PhreeqcRM methods. The constructor requires two arguments:
// the number of cells in the user's model, and either (a) the
// number of threads for OpenMP parallelization, or (b) an MPI
// communicator.

// Args:
// 	ngrid (int): Number of cells in the user's model. 
// 	nthreads (int): Number of threads for parallelization with
// 		OpenMP or an MPI communicator if PhreeqcRM is compiled with
// 		MPI. For OpenMP, a value of zero causes the program to set
// 		nthreads to the number of logical processors of the
// 		computer.
		
// Returns:
// 	 A BMIPhreeqcRM instance."
// %enddef
// %feature("docstring", BMIPhreeqcRM_ii_DOCSTRING) BMIPhreeqcRM::BMIPhreeqcRM(int, MP_TYPE);

%define add_output_vars_DOCSTRING
"Method that defines additional output variables that can be
retrieved by the BMI method :meth:`get_value`. 

add_output_vars selects sets of variables that can be
retieved by the :meth:`get_value` method. Sets of variables
can be included or excluded with multiple calls to this
method. All calls must precede the final call to
FindComponents. FindComponents generates SELECTED_OUTPUT
333 and USER_PUNCH 333 data blocks that make the selected
variables accessible. Variables will only be accessible if
the system includes the given reactant; for example, no gas
variables will be created if there are no GAS_PHASEs in the
model.

Args:
	option (str): A string value, among those listed below,
		that includes or excludes variables from
		get_output_var_names, get_value, get_var_itemsize, and
		other BMI methods.

	def (str): A string value that can be 'false', 'true', or a
		list of items to be included as accessible variables. A
		value of 'false', excludes all variables of the given type;
		a value of 'true' includes all variables of the given type
		for the current system; a list specifies a subset of items
		of the given type.

Values for the the parameter option include:

'AddOutputVars': False excludes all variables; True causes
	the settings for each variable group to determine the
	variables that will be defined. Default True;

'SolutionProperties': False excludes all solution property
	variables; True includes variables pH, pe, alkalinity,
	ionic strength, water mass, charge balance, percent error,
	and specific conductance. Default True.

'SolutionTotalMolalities': False excludes all total element
	and element redox state variables; True includes all
	elements and element redox state variables for the system
	defined for the calculation; list restricts variables to
	the specified elements and redox states. Default True.

'ExchangeMolalities': False excludes all variables related
	to exchange; True includes all variables related to
	exchange; list includes variables for the specified
	exchange species. Default True.

'SurfaceMolalities': False excludes all variables related
	to surfaces; True includes all variables related to
	surfaces; list includes variables for the specified surface
	species. Default True.

'EquilibriumPhases': False excludes all variables related
	to equilibrium phases; True includes all variables related
	to equilibrium phases; list includes variables for the
	specified equilibiurm phases. Default True.

'Gases': False excludes all variables related to gases;
	True includes all variables related to gases; list includes
	variables for the specified gas components. Default True.

'KineticReactants': False excludes all variables related to
	kinetic reactants; True includes all variables related to
	kinetic reactants; list includes variables for the
	specified kinetic reactants. Default True.

'SolidSolutions': False excludes all variables related to
	solid solutions; True includes all variables related to
	solid solutions; list includes variables for the specified
	solid solutions components. Default True.

'CalculateValues': False excludes all calculate values;
	True includes all calculate values; list includes the
	specified calculate values. CALCLUATE_VALUES can be used to
	calculate geochemical quantities not available in the other
	sets of variables. Default True.

'SolutionActivities': False excludes all aqueous species;
	True includes all aqueous species; list includes only the
	specified aqueous species. Default False.

'SolutionMolalities': False excludes all aqueous species;
	True includes all aqueous species; list includes only the
	specified aqueous species. Default False.

'SaturationIndices': False excludes all saturation indices;
	True includes all saturation indices; list includes only
	the specified saturation indices. Default False."
%enddef
%feature("docstring") BMIPhreeqcRM::AddOutputVars add_output_vars_DOCSTRING


%define finalize_DOCSTRING
"Finalize closes any files open in the BMIPhreeqcRM instance."
%enddef
%feature("docstring") BMIPhreeqcRM::Finalize finalize_DOCSTRING


%define get_component_name_DOCSTRING
"Basic Model Interface method that returns the component
name--PhreeqcRM.

The BMI interface to PhreeqcRM provides most of the
functions necessary to couple a flow and transport model
with chemistry. The native PhreeqcRM methods
provide additional methods to complete the interface. 

Returns:
	str: The string 'PhreeqcRM'."
%enddef
%feature("docstring") BMIPhreeqcRM::GetComponentName get_component_name_DOCSTRING


%define get_current_time_DOCSTRING
"Basic Model Interface method that returns the current simulation time.

The method :meth:`GetTime` also returns the current
simulation time. The current time can be set by the methods
:meth:`set_value`('Time') or :meth:`SetTime`.

Returns:
	float: The current simulation time, in seconds."
%enddef
%feature("docstring") BMIPhreeqcRM::GetCurrentTime get_current_time_DOCSTRING


%define get_end_time_DOCSTRING
"Basic Model Interface method that returns
:meth:`get_current_time` plus :meth:`get_time_step`.

Returns:
	float: The end of the time step, in seconds."
%enddef
%feature("docstring") BMIPhreeqcRM::GetEndTime get_end_time_DOCSTRING

%define get_grid_rank_DOCSTRING
"Returns a rank of 1 for grid 0. 

BMIPhreeqcRM only has a 1D series of
cells; any grid or spatial information must
be found in the user's model.

Args:
	grid (int): Grid number, only grid 0 is considered.
Returns:
	(int) Rank of 1 is returned for grid 0; 0 for all other values of 
		grid."
%enddef
%feature("docstring") BMIPhreeqcRM::GetGridRank get_grid_rank_DOCSTRING


%define get_grid_size_DOCSTRING
"Returns the number of cells specified
at initialization of the BMIPhreeqcRM instance. 

Args:
	grid (int): Grid number, only grid 0 is considered.

Returns:
	(int) Number of cells in the user's modle (same value as
		:meth:`GetGridCellCount` is returned for grid 0; 0 for all other
		values of grid."
%enddef
%feature("docstring") BMIPhreeqcRM::GetGridSize get_grid_size_DOCSTRING


%define get_grid_type_DOCSTRING
"Grid type is points. 

No grid information is available in BMIPhreeqcRM; all grid
information must be found in the user's model. 

Args:
	grid (int): Grid number, only grid 0 is considered. 
	
Returns:
	(str) 'Points' is returned for grid 0; 'Undefined grid
		identifier' is returned for all other values of @a grid."
%enddef
%feature("docstring") BMIPhreeqcRM::GetGridType get_grid_type_DOCSTRING


%define get_input_item_count_DOCSTRING
"Basic Model Interface method that returns count of input
variables that can be set with :meth:`set_value`.

Returns:
	int: Count of input variables that can be set with
		:meth:`set_value`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetInputItemCount get_input_item_count_DOCSTRING


%define get_input_var_names_DOCSTRING
"Basic Model Interface method that returns a list of the
variable names that can be set with :meth:`set_value`.

Returns:
	tuple (str): A list of the variable names that can be set with 
		:meth:`set_value`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetInputVarNames get_input_var_names_DOCSTRING


%define get_output_item_count_DOCSTRING
"Basic Model Interface method that returns count of output
variables that can be retrieved with :meth:`get_value`.

Returns:
	int: Count of output variables that can be retrieved with 
		:meth:`get_value`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetOutputItemCount get_output_item_count_DOCSTRING


%define get_output_var_names_DOCSTRING
"Basic Model Interface method that returns a list of the
variable names that can be retrieved with
:meth:`get_value`.

Returns:
	tuple (str): A list of the variable names that can be retrieved with
		:meth:`get_value`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetOutputVarNames get_output_var_names_DOCSTRING


%define get_pointable_item_count_DOCSTRING
"Basic Model Interface method that returns count of 
variables that can be reference with :meth:`get_value_ptr`.

Returns:
	int: Count of output variables that can be referenced with 
		:meth:`get_value_ptr`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetPointableItemCount get_pointable_item_count_DOCSTRING


%define get_pointable_var_names_DOCSTRING
"Basic Model Interface method that returns a list of the
variable names that can be referenced by
:meth:`get_value_ptr`.

Returns:
	tuple (str): A list of the variable names that can be
		referenced with :meth:`get_value_ptr`."
%enddef
%feature("docstring") BMIPhreeqcRM::GetPointableVarNames get_pointable_var_names_DOCSTRING


%define get_time_step_DOCSTRING
"Basic Model Interface method that returns the current
simulation time step. (Same as :meth:`GetTimeStep`.) 

The time step can be set by :meth:`set_value`('TimeStep')
or :meth:`SetTimeStep`.

Returns:
	float: The current simulation time step, in seconds."
%enddef
%feature("docstring") BMIPhreeqcRM::GetTimeStep get_time_step_DOCSTRING


%define get_time_units_DOCSTRING
"Basic Model Interface method that returns the time units
of PhreeqcRM. 

All time units are seconds for PhreeqcRM.

Returns:
	string: Returns the string 'seconds'."
%enddef
%feature("docstring") BMIPhreeqcRM::GetTimeUnits get_time_units_DOCSTRING


%define get_value_DOCSTRING
"Basic Model Interface method that retrieves model variables. 

Only variables in the list provided by
:meth:`get_output_var_names` can be retrieved. 

Args:
	name (str): Name of the variable to retrieve.
	dest (type defined below): Variable in which to place results.

Variable names for the first argument (name) and variable type of the
second argument (dest).
'ComponentCount', dest: int;
'Components', dest : tuple (str);
'Concentrations', dest : DoubleVector;
'CurrentSelectedOutputUserNumber', dest : int;
'Density', dest : DoubleVector;
'ErrorString', dest : std::string;
'FilePrefix', dest : std::string;
'Gfw', dest : std::vector<double>;
'GridCellCount', dest : int;
'InputVarNames', dest : tuple (str);
'OutputVarNames', dest : tuple (str);
'Porosity', dest : DoubleVector;
'Pressure', dest : DoubleVector;
'Saturation', dest : DoubleVector;
'SelectedOutput', dest : DoubleVector;
'SelectedOutputColumnCount', dest : int;
'SelectedOutputCount', dest : int;
'SelectedOutputHeadings', dest : tuple (str);
'SelectedOutputOn', dest : Boolean;
'SelectedOutputRowCount', dest : int;
'SolutionVolume', dest : DoubleVector;
'Temperature', dest : DoubleVector;
'Time', dest : double;
'TimeStep', dest: double.

Additional output variables can be defined through the 
method :meth:`add_output_vars`."
%enddef
%feature("docstring") BMIPhreeqcRM::get_value get_value_DOCSTRING


%define get_value_ptr_DOCSTRING
"Returns a pointer that references current model values.

Unlike the values returned from :meth:`get_value`, the
reference always points to the current values of the
variable, even if the model's state has changed. 

Args:
	name (str): Name of the variable to for which a pointer is
		retrieved.

Returns:
	Pointer to data. 

The following list gives the variable names for which
pointers can be obtained:
'ComponentCount'
'Concentrations'
'DensityCalculated'
'Gfw'
'GridCellCount'
'Porosity'
'Pressure'
'SaturationCalculated'
'SelectedOutputOn'
'SolutionVolume'
'Temperature'
'Time'
'TimeStep'
'Viscosity'"
%enddef
%feature("docstring") BMIPhreeqcRM::get_value_ptr get_value_ptr_DOCSTRING


%define get_var_itemsize_DOCSTRING
"Basic Model Interface method that retrieves size of an
individual item that can be set or retrived. 

Sizes may be sizeof(int), sizeof(float), or a character
length for string variables. Only variables in the list
provided by :meth:`get_input_var_names` can be set. Only
variables in the list provided by
:meth:`get_output_var_names` can be retrieved.

Args:
	name (str): Name of the variable to retrieve item size.

Returns:
	int: Size of one element of variable."
%enddef
%feature("docstring") BMIPhreeqcRM::GetVarItemsize get_var_itemsize_DOCSTRING


%define get_var_nbytes_DOCSTRING
"Basic Model Interface method that retrieves the total
number of bytes that are set or retrieved.

Variables are set with :meth:`set_value` and retrieved with
:meth:`get_value`. Only variables in the list provided by
:meth:`get_input_var_names` can be set. Only variables in
the list provided by :meth:`get_output_var_names` can be
retrieved.

Args:
	name (str): Name of the variable to retrieve total bytes.

Returns:
	int: Total number of bytes set or retrieved for variable."
%enddef
%feature("docstring") BMIPhreeqcRM::GetVarNbytes get_var_nbytes_DOCSTRING


%define get_var_type_DOCSTRING
"Basic Model Interface method that retrieves the type of a
variable that can be set or retrieved.

Types are 'int', 'float', or 'str'.

Args:
	name (str): Name of the variable to retrieve type.

Returns:
	str: Character string of variable type."
%enddef
%feature("docstring") BMIPhreeqcRM::GetVarType get_var_type_DOCSTRING


%define get_var_units_DOCSTRING
"Basic Model Interface method that retrieves the units of a
variable that can be set or retrieved.

Args:
	name (str): Name of the variable to retrieve type.

Returns:
	str: Character string of units for variable."
%enddef
%feature("docstring") BMIPhreeqcRM::GetVarUnits get_var_units_DOCSTRING


%define initialize_DOCSTRING
"Basic Model Interface method that can be used to
initialize a PhreeqcRM instance. 

This method is equivalent to :meth:`InitializeYAML`. A YAML
file is used to initialize the instance. The file contains
a YAML map of PhreeqcRM methods and the arguments
corresponding to the method. 

For example,

- key: LoadDatabase
  database: phreeqc.dat
- key: RunFile
  workers: true
  initial_phreeqc: true
  utility: true
  chemistry_name: advect.pqi

:meth:`initialize` will read the YAML file and execute the
specified methods with the specified arguments. Using YAML
terminology, the argument(s) for a method may be a scalar,
a sequence, or a map, depending if the argument is a single
item, a single vector, or there are multiple arguments. In
the case of a map, the name associated with each argument
(for example 'chemistry_name' above) is arbitrary. The
names of the map keys for map arguments are not used in
parsing the YAML file; only the order of the arguments is
important.

The PhreeqcRM methods that can be specified in a YAML file
include:

CloseFiles(void);
CreateMapping(int list, numpy.ndarray, or tuple: grid2chem);
DumpModule();
FindComponents();
InitialPhreeqc2Module(int list, numpy.ndarray, or tuple: initial_conditions1);
InitialPhreeqc2Module(int list, numpy.ndarray, or tuple: initial_conditions1,
int list, numpy.ndarray, or tuple: initial_conditions2,
float list, numpy.ndarray, or tuple: fraction1);
InitialPhreeqcCell2Module(int n, int list, numpy.ndarray, or tuple: cell_numbers);
LoadDatabase(string database);
OpenFiles(void);
OutputMessage(string str);
RunCells(void);
RunFile(Boolean workers, Boolean initial_phreeqc, Boolean utility, string chemistry_name);
RunString(Boolean workers, Boolean initial_phreeqc, Boolean utility, string input_string);
ScreenMessage(string str);
SetComponentH2O(Boolean tf);
SetConcentrations(float list, numpy.ndarray, or tuple: c);
SetCurrentSelectedOutputUserNumber(int n_user);
SetDensityUser(float list, numpy.ndarray, or tuple: density);
SetDumpFileName(string dump_name);
SetErrorHandlerMode(int mode);
SetErrorOn(Boolean tf);
SetFilePrefix(string prefix);
SetGasCompMoles(float list, numpy.ndarray, or tuple: gas_moles);
SetGasPhaseVolume(float list, numpy.ndarray, or tuple: gas_volume);
SetPartitionUZSolids(Boolean tf);
SetPorosity(float list, numpy.ndarray, or tuple: por);
SetPressure(float list, numpy.ndarray, or tuple: p);
SetPrintChemistryMask(int list, numpy.ndarray, or tuple: cell_mask);
SetPrintChemistryOn(Boolean workers, Boolean initial_phreeqc, Boolean utility);
SetRebalanceByCell(Boolean tf);
SetRebalanceFraction(float f);
SetRepresentativeVolume(float list, numpy.ndarray, or tuple: rv);
SetSaturationUser(float list, numpy.ndarray, or tuple: sat);
SetScreenOn(Boolean tf);
SetSelectedOutputOn(Boolean tf);
SetSpeciesSaveOn(Boolean save_on);
SetTemperature(float list, numpy.ndarray, or tuple: t);
SetTime(float time);
SetTimeConversion(float conv_factor);
SetTimeStep(float time_step);
SetUnitsExchange(int option);
SetUnitsGasPhase(int option);
SetUnitsKinetics(int option);
SetUnitsPPassemblage(int option);
SetUnitsSolution(int option);
SetUnitsSSassemblage(int option);
SetUnitsSurface(int option);
SpeciesConcentrations2Module(float list, numpy.ndarray, or tuple: species_conc);
StateSave(int istate);
StateApply(int istate);
StateDelete(int istate);
UseSolutionDensityVolume(Boolean tf);
WarningMessage(string warnstr);

Args:
	config_file (str): File with YAML definitions for PhreeqcRM 
		initialization."
%enddef
%feature("docstring") BMIPhreeqcRM::Initialize initialize_DOCSTRING


%define set_value_DOCSTRING
"Basic Model Interface method that sets model variables. 

Only variables in the list provided by
:meth:`get_input_var_names` can be set. The BMI interface
to PhreeqcRM is only partial, and provides only the most
basic functions. The native PhreeqcRM methods provide a
complete interface, and it is expected that the native
methods will be used in preference to the BMI methods.

Variable names for the first argument of set_value and the
equivalent PhreeqcRM method are as follows:
'Concentrations'
'DensityUser'
'FilePrefix'
'NthSelectedOutput'
'Porosity'
'Pressure'
'SaturationUser'
'SelectedOutputOn'
'Temperature'
'Time'
'TimeStep'"
%enddef
%feature("docstring") BMIPhreeqcRM::set_value set_value_DOCSTRING


%define update_DOCSTRING
"Basic Model Interface method that runs PhreeqcRM for one time step. 

PhreeqcRM will equilibrate the solutions with all
equilibrium reactants (EQUILIBRIUM_PHASES, EXCHANGE,
GAS_PHASE, SOLID_SOLUTIONS, and SURFACE) and integrate
KINETICS reactions for the specified time step
(:meth:`SetTimeStep`). The time will be updated, and all
pointer variables will be updated."
%enddef
%feature("docstring") BMIPhreeqcRM::Update update_DOCSTRING

%define update_until_DOCSTRING
"Runs chemistry for a time step calculated
as end_time minus the current time.

Current time is found with :meth:`get_current_time`.

Args:
	end_time (float): Time at the end of the simulation
		period."
%enddef
%feature("docstring") BMIPhreeqcRM::UpdateUntil update_until_DOCSTRING
