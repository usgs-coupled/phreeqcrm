%define GetGridCellCountYAML_DOCSTRING
"GetGridCellCountYAML will read the YAML file and extract
the value of GridCellCount, which can be used to construct
a PhreeqcRM instance.

One constructor for a PhreeqcRM instance requires a value
for the number of cells. If a GUI or preprocessor is used
to write a YAML file to initialize PhreeqcRM, the number of
cells can be written to the YAML file and extracted with
this method.

Args:
	YAML_file (str): String containing the YAML file name.

Returns:
	int: Number of grid cells specified in the YAML file;
		returns zero if GridCellCount is not defined."
%enddef
%feature("docstring") PhreeqcRM::GetGridCellCountYAML GetGridCellCountYAML_DOCSTRING


%define PhreeqcRM_DOCSTRING
"Constructor for the PhreeqcRM reaction module. 

PhreeqcRM is a reaction module for reactive-transport
simulators. Based on IPhreeqc, it allows access to all
PHREEQC reaction capabilities. It contains methods for
initial and boundary conditions, running reactions on all
model cells, transfer of data to and from the module, and
parallelization by MPI or OpenMP.

If the code is compiled with the preprocessor directive
USE_OPENMP, the reaction module use OPENMP and multiple
threads. If the code is compiled with the preprocessor
directive USE_MPI, the reaction module will use MPI and
multiple processes. If neither preprocessor directive is
used, the reaction module will be serial (unparallelized).

Args:
	nxyz (int): The number of grid cells in the users model. 
	
	thread_count_or_communicator (int): If multithreaded, the
		number of threads to use in parallel segments of the code.
		If thread_count_or_communicator is <= 0, the number of
		threads is set equal to the number of processors in the
		computer. If multiprocessor, the MPI communicator to use
		within the reaction module.
	
	io(PHRQ_io): Optionally, a PHRQ_io input/output object can
		be provided to the constructor. By default a PHRQ_io object
		is constructed to handle reading and writing files."
%enddef
%feature("docstring") PhreeqcRM::PhreeqcRM PhreeqcRM_DOCSTRING


%define CloseFiles_DOCSTRING
"Close the output and log files.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::CloseFiles CloseFiles_DOCSTRING


%define Concentrations2Utility_DOCSTRING
"N sets of component concentrations are converted to
SOLUTIONs numbered 1-n in the Utility IPhreeqc.

The solutions can be reacted and manipulated with the
methods of IPhreeqc. If solution concentration units
(:meth:`SetUnitsSolution`) are per liter, one liter of
solution is created in the Utility instance; if solution
concentration units are mass fraction, one kilogram of
solution is created in the Utility instance. The motivation
for this method is the mixing of solutions in wells, where
it may be necessary to calculate solution properties (pH
for example) or react the mixture to form scale minerals.

Args:
	c (float list, numpy.ndarray, or tuple): Array of
		concentrations to be made SOLUTIONs in Utility IPhreeqc.
		Array contains n values for each component
		(:meth:`GetComponentCount`) in sequence.
	tc (float list, numpy.ndarray, or tuple): Array of
		temperatures to apply to the SOLUTIONs, in degrees C.
		Array of size n.
	p_atm (float list, numpy.ndarray, or tuple): Array of
		pressures to apply to the SOLUTIONs, in atm. Array of size
		n.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::Concentrations2Utility Concentrations2Utility_DOCSTRING


%define CreateMapping_DOCSTRING
"Provides a mapping from grid cells in the user model to
reaction cells for which chemistry needs to be run.

The mapping is used to eliminate inactive cells and to use
symmetry to decrease the number of cells for which
chemistry must be run. The array grid2chem of size nxyz
(the number of grid cells, :meth:`GetGridCellCount`) must
contain the set of all integers 0 <= i < count_chemistry,
where count_chemistry is a number less than or equal to
nxyz. Inactive cells are assigned a negative integer. The
mapping may be many-to-one to account for symmetry. Default
is a one-to-one mapping--all user grid cells are reaction
cells (equivalent to grid2chem values of
0,1,2,3,...,nxyz-1).

Args:
	grid2chem (int list, numpy.ndarray, or tuple): An array of
		integers: Nonnegative is a reaction-cell number (0 based),
		negative is an inactive cell. Array is of size nxyz (number
		of grid cells, :meth:`GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::CreateMapping CreateMapping_DOCSTRING


%define DecodeError_DOCSTRING
"If result is negative, this method prints an error message
corresponding to IRM_RESULT result. If result is
non-negative, no action is taken.

Args:
	result (IRM_RESULT): An IRM_RESULT value returned by one of
		the reaction- module methods.

IRM_RESULT definition:
typedef enum {
IRM_OK            =  0,  //Success
IRM_OUTOFMEMORY   = -1,  //Failure, Out of memory
IRM_BADVARTYPE    = -2,  //Failure, Invalid VAR type
IRM_INVALIDARG    = -3,  //Failure, Invalid argument
IRM_INVALIDROW    = -4,  //Failure, Invalid row
IRM_INVALIDCOL    = -5,  //Failure, Invalid column
IRM_BADINSTANCE   = -6,  //Failure, Invalid rm instance id
IRM_FAIL          = -7,  //Failure, Unspecified
} IRM_RESULT;"
%enddef
%feature("docstring") PhreeqcRM::DecodeError DecodeError_DOCSTRING


%define DumpModule_DOCSTRING
"Writes the contents of all workers to file in _RAW formats
(see appendix of PHREEQC version 3 manual), including
SOLUTIONs and all reactants.

Args:
	dump_on (Boolean): Signal for writing the dump file, true
		or false.
	append (Boolean): Signal to append to the contents of the
		dump file, true or false.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::DumpModule DumpModule_DOCSTRING


%define ErrorHandler_DOCSTRING
"Checks result for an error code. 

If result is negative, the result is decoded
(:meth:`DecodeError`), and printed as an error message
along with the e_string, and an exception is thrown. If the
result is nonnegative, no action is taken.

Args:
	result (IRM_RESULT): IRM_RESULT to be checked for an error.
	e_string (str): String to be printed if an error is found."
%enddef
%feature("docstring") PhreeqcRM::ErrorHandler ErrorHandler_DOCSTRING


%define ErrorMessage_DOCSTRING
"Send an error message to the screen, the output file, and
the log file.

Args:
	error_string (str): String to be printed.
	prepend (Boolean): True, prepends error_string with 'Error: ';
		false, error_string is used with no prepended text."
%enddef 
%feature("docstring") PhreeqcRM::ErrorMessage ErrorMessage_DOCSTRING

%define FindComponents_DOCSTRING 
"This method accumulates a list of elements that have been
defined in a solution or any other reactant
(EQUILIBRIUM_PHASE, KINETICS, and others), including charge
imbalance.
	
This method can be called multiple times and the list that
is created is cummulative.The list is the set of components
that needs to be transported.By default the list includes
water, excess H and excess O(the H and O not contained in
water); alternatively, the list may be set to contain total
H and total O(:meth:`SetComponentH2O`), which requires
transport results to be accurate to eight or nine
significant digits.If multicomponent diffusion(MCD) is to
be modeled, there is a capability to retrieve aqueous
species concentrations(:meth:`GetSpeciesConcentrations`)
and to set new solution concentrations after MCD by using
individual species concentrations
(:meth:`SpeciesConcentrations2Module`) .To use these methods the
save - species property needs to be turned on
(:meth:`SetSpeciesSaveOn`) .If the save - species property is on,
FindComponents will generate a list of aqueous species
(:meth:`GetSpeciesCount`, :meth:`GetSpeciesNames`), their
diffusion coefficients at 25 C(:meth:`GetSpeciesD25`),
and their charge(:meth:`GetSpeciesZ`) .

returns: 
	int : Number of components currently in the list, or
		IRM_RESULT error code(negative value, see
		:meth:`DecodeError`)."

%enddef 
%feature("docstring") PhreeqcRM::FindComponents FindComponents_DOCSTRING

%define GetBackwardMapping_DOCSTRING 
"Returns a dictionary of numpy arrays, where the reaction cell
n is mapped to an array of user grid cell numbers that are mapped
to n.

Each reaction cell number has an array with one or more grid
cell numbers.

Returns: 
	dict (numpy.ndarray) : Dict of array of ints. For each
	reaction cell n, n is mapped to an array that contains
	contains the grid-cell numbers that map to the reaction
	cell."
%enddef 
%feature("docstring") PhreeqcRM::GetBackwardMapping GetBackwardMapping_DOCSTRING


%define GetChemistryCellCount_DOCSTRING 
"Returns the number of reaction cells in the reaction module. 

The number of reaction cells is defined by the set of
non-negative integers in the mapping from grid cells 
(:meth:`CreateMapping`), or, by default, the number of grid
cells (:meth:`GetGridCellCount`) .The number of reaction
cells is less than or equal to the number of grid cells in
the user model.

Returns: 
	int: Number of reaction cells."
%enddef 
%feature("docstring") PhreeqcRM::GetChemistryCellCount GetChemistryCellCount_DOCSTRING

%define GetComponentCount_DOCSTRING 
"Returns the number of components in the reaction-module
component list.

returns: 
	int : The number of components in the reaction - module
		component list.The component list is generated by calls to
		:meth:`FindComponents`.The return value from the last call
		to :meth:`FindComponents` is equal to the return value from
		GetComponentCount."
%enddef 
%feature("docstring") PhreeqcRM::GetComponentCount GetComponentCount_DOCSTRING

%define GetComponents_DOCSTRING 
"Returns the reaction-module component list that was
generated by calls to :meth:`FindComponents`.

returns: 
	numpy ndarray (str)): Each string is a component name."
%enddef
%feature("docstring") PhreeqcRM::GetComponents GetComponents_DOCSTRING


%define GetConcentrations_DOCSTRING
"Return solution concentrations for all model cells. 

Units of concentration for c are defined by
:meth:`SetUnitsSolution`. For per liter concentration
units, solution volume is used to calculate the
concentrations for c. For mass-fraction concentration
units, the solution mass is used to calculate
concentrations for c. Two options are available for the
volume and mass of solution that are used in converting to
transport concentrations: (1) the volume and mass of
solution are calculated by PHREEQC, or (2) the volume of
solution is the product of saturation
(:meth:`SetSaturationUser`), porosity (:meth:`SetPorosity`),
and representative volume
(:meth:`SetRepresentativeVolume`), and the mass of solution
is volume times density as defined by :meth:`SetDensityUser`.
:meth:`UseSolutionDensityVolume` determines which option is
used. For option 1 , the databases that have partial molar
volume definitions needed to accurately calculate solution
volume are phreeqc.dat, Amm.dat, and pitzer.dat.

Returns:
	numpy.ndarray: Dimension of the array is ncomps times nxyz,
		where, ncomps is the result of :meth:` FindComponents` or
		:meth:`GetComponentCount`, and nxyz is the number of user
		grid cells (:meth:`GetGridCellCount`). Values for inactive
		cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetConcentrations GetConcentrations_DOCSTRING


%define GetCurrentSelectedOutputUserNumber_DOCSTRING
"Returns the user number of the current selected-output definition.

:meth:`SetCurrentSelectedOutputUserNumber` or
:meth:`SetNthSelectedOutput` specifies which of the
selected-output definitions is used.

Returns:
	int: User number of the the current selected-output definition,
		negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::GetCurrentSelectedOutputUserNumber GetCurrentSelectedOutputUserNumber_DOCSTRING


%define GetDatabaseFileName_DOCSTRING
"Returns the file name of the database. Should be called
after :meth:`LoadDatabase`.

Returns:
	string: The file name defined by :meth:`LoadDatabase`."
%enddef
%feature("docstring") PhreeqcRM::GetDatabaseFileName GetDatabaseFileName_DOCSTRING


%define GetDensityCalculated_DOCSTRING
"Returns calculated densities from the reaction-module cells. 

This method always returns the calculated densities;
:meth:`SetDensityUser` does not affect the result.

Returns:
	numpy.darray: Dimension of the array is set to nxyz, where
		nxyz is the number of user grid cells
		(:meth:`GetGridCellCount`). Values for inactive cells are
		set to 1e30. Densities are those calculated by the reaction
		module. Only the following databases distributed with
		PhreeqcRM have molar volume information needed to
		accurately calculate density: phreeqc.dat, Amm.dat, and
		pitzer.dat."
%enddef
%feature("docstring") PhreeqcRM::GetDensityCalculated GetDensityCalculated_DOCSTRING


%define GetEndCell_DOCSTRING
"Returns an array of integers that contains the largest
reaction-cell number assigned to each worker.

Each worker is assigned a range of reaction-cell numbers
that are run during a call to :meth:`RunCells`. The range
of reaction cells for a worker may vary as load rebalancing
occurs. At any point in the calculations, the first cell
and last cell to be run by a worker can be found in the
arrays returned by :meth:`GetStartCell` and
:meth:`GetEndCell`. Each method returns an array of
integers that has length of the number of threads
(:meth:`GetThreadCount`), if using OPENMP, or the number of
processes (:meth:`GetMpiTasks`), if using MPI.

Returns:
	numpy.ndarray (int): Array of integers, one for each
		worker, that gives the last reaction cell to be run by each
		worker."
%enddef
%feature("docstring") PhreeqcRM::GetEndCell GetEndCell_DOCSTRING


%define GetEquilibriumPhases_DOCSTRING
"Returns a the names of all equilibrium phases defined in the module. 

The list includes all phases included in any
EQUILIBRIUM_PHASES definitions in the initial-phreeqc
module. :meth:`FindComponents` must be called before
:meth:`GetEquilibriumPhases`. This method may be useful
when generating selected output definitions related to
equilibrium phases.

Returns:
	tuple (str): each string is a unique equilibrium
		phases name."

%enddef 
%feature("docstring") PhreeqcRM::GetEquilibriumPhases GetEquilibriumPhases_DOCSTRING

%define GetEquilibriumPhasesCount_DOCSTRING
"Returns the number of equilibrium phases found by
:meth:`FindComonents`.
	
:meth:`FindComponents` must be called before
:meth:`GetEquilibriumPhasesCount`.This method may be useful
when generating selected output definitions related to
equilibrium phases.

returns: 
	int : The number of equilibrium phases."
%enddef
%feature("docstring") PhreeqcRM::GetEquilibriumPhasesCount GetEquilibriumPhasesCount_DOCSTRING

%define GetErrorHandlerMode_DOCSTRING
"Get the setting for the action to be taken when the reaction module 
encounters an error.

Options are 0, return to calling program with an error
return code(default); 1, throw an exception, which can be
caught in C++ (for C and Fortran, the program will exit);
2, attempt to exit gracefully.

Returns:
	int: Current setting for the error handling mode: 0, 1, or 2."
%enddef
%feature("docstring") PhreeqcRM::GetErrorHandlerMode GetErrorHandlerMode_DOCSTRING


%define GetErrorString_DOCSTRING
"Returns a standard string containing error messages
related to the last call to a BMIPhreeqcRM method.

Returns:
	string: Error messages related to the last call to a
		PhreeqcRM method."
%enddef
%feature("docstring") PhreeqcRM::GetErrorString GetErrorString_DOCSTRING


%define GetExchangeNames_DOCSTRING
"Returns a tuple of exchange names (such as 'X') that
correspond with the exchange species names.

:meth:`FindComponents` must be called before
:meth:`GetExchangeNames`. The exchange names tuple is the
same length as the exchange species names tuple and
provides the corresponding exchange site. This method may
be useful when generating selected output definitions
related to exchangers.

Returns:
	tuple (str): Each string is an exchange name
		corresponding to the exchange species names; an exchange
		name may occur multiple times."
%enddef 
%feature("docstring") PhreeqcRM::GetExchangeNames GetExchangeNames_DOCSTRING

%define GetExchangeSpecies_DOCSTRING
"Returns a tuple of exchange species names (such as 'NaX').

The list of exchange species(such as 'NaX') is derived from
the list of components(:meth:`FindComponents`) and the
list of all exchange names(such as 'X') that are included
in EXCHANGE definitions in the initial phreeqc instance.
:meth:`FindComponents` must be called before
:meth:`GetExchangeSpecies`.This method may be useful when
generating selected output definitions related to
exchangers.

Returns: 
	tuple (str): Each string is a unique exchange species
		name."
%enddef
%feature("docstring") PhreeqcRM::GetExchangeSpecies GetExchangeSpecies_DOCSTRING


%define GetExchangeSpeciesCount_DOCSTRING
"Returns the number of exchange species found by
:meth:`FindComponents`.

:meth:`FindComponents` must be called before
:meth:`GetExchangeSpeciesCount`. This method may be useful
when generating selected output definitions related to
exchangers.

Returns:
	int: The number of exchange species."
%enddef
%feature("docstring") PhreeqcRM::GetExchangeSpeciesCount GetExchangeSpeciesCount_DOCSTRING


%define GetFilePrefix_DOCSTRING
"Returns the file prefix for the output (.chem.txt) and log
files (.log.txt).

Returns:
	string: The file prefix as set by :meth:`SetFilePrefix`, or
		'myrun', by default."
%enddef
%feature("docstring") PhreeqcRM::GetFilePrefix GetFilePrefix_DOCSTRING


%define GetForwardMapping_DOCSTRING
"Returns an array of ints that is a mapping from grid cells
to reaction cells.

The mapping is used to eliminate cells that are inactive
and cells that are unnecessary because of symmetry from the
list of cells for which reactions must be run. The mapping
may be many-to-one to account for symmetry. The mapping is
set by :meth:`CreateMapping`, or, by default, is a
one-to-one mapping--all grid cells are reaction cells
(array contains 0,1,2,3,...,nxyz-1).

Returns:
	numpy.ndarray (int): An array of integers of size nxyz
		(number of grid cells, :meth:`GetGridCellCount`).
		Nonnegative is a reaction-cell number (0 based), negative
		is an inactive cell."
%enddef
%feature("docstring") PhreeqcRM::GetForwardMapping GetForwardMapping_DOCSTRING


%define GetGasComponents_DOCSTRING
"Returns a tuple of the names of all gas components found by 
:meth:`FindComponents`. 

:meth:`FindComponents` must be called before
:meth:`GetGasComponents`. This method may be useful when
generating selected output definitions related to gas
phases.

Returns:
	tuple (str): Each string is a unique gas component
		name."
%enddef
%feature("docstring") PhreeqcRM::GetGasComponents GetGasComponents_DOCSTRING


%define GetGasComponentsCount_DOCSTRING
"Returns the number of gas phase components found by 
:meth:`FindComponents`. 

:meth:`FindComponents` must be called before
:meth:`GetGasComponentsCount`. This method may be useful
when generating selected output definitions related to gas
phases.

Returns:
	int: The number of gas phase components."
%enddef
%feature("docstring") PhreeqcRM::GetGasComponentsCount GetGasComponentsCount_DOCSTRING


%define GetGasCompMoles_DOCSTRING
"Return moles of gas components for each model cell.

Returns:
	numpy.ndarray (float): Dimension of the array is
		ngas_comps times nxyz, where, ngas_comps is the result of
		:meth:`GetGasComponentsCount`, and nxyz is the number of
		user grid cells (:meth:`GetGridCellCount`). If a gas
		component is not defined for a cell, the number of moles is
		set to -1. Values for inactive cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetGasCompMoles GetGasCompMoles_DOCSTRING


%define GetGasCompPressures_DOCSTRING
"Return pressures of gas components for each model cell.

Returns:
	numpy.ndarray (float): Dimension of the array is set
		to ngas_comps times nxyz, where ngas_comps is the result of
		:meth:`GetGasComponentsCount`, and nxyz is the number of
		user grid cells (:meth:`GetGridCellCount`). If a gas
		component is not defined for a cell, the pressure is set to
		-1. Values for inactive cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetGasCompPressures GetGasCompPressures_DOCSTRING


%define GetGasCompPhi_DOCSTRING
"Return fugacity coefficients (phi) of each gas component
for each model cell.

Fugacity is equal to the gas component pressure times the
fugacity coefficient.

Returns: 
	numpy.ndarray (float): Dimension of the array is ngas_comps
		times nxyz, where, ngas_comps is the result of
		:meth:`GetGasComponentsCount`, and nxyz is the number of
		user grid cells (:meth:`GetGridCellCount`). If a gas
		component is not defined for a cell, the fugacity
		coefficient is set to -1. Values for inactive cells are set
		to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetGasCompPhi GetGasCompPhi_DOCSTRING


%define GetGasPhaseVolume_DOCSTRING
"Returns volume of gas phase for each model cell.

Returns:
	numpy.ndarray (float): Dimension of the array is 
		nxyz, where, nxyz is the number of user grid cells
		(:meth:`GetGridCellCount`). If a gas phase is not defined
		for a cell, the volume is set to -1. Values for inactive
		cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetGasPhaseVolume GetGasPhaseVolume_DOCSTRING


%define GetGfw_DOCSTRING
"Returns an array with the gram-formula weight of each
component.

Called after :meth:`FindComponents`. Order of weights
corresponds to the list of components from
:meth:`GetComponents`.

Returns:
	numpy.ndarray (float): Each value is a component gram-
		formula weight, g/mol."
%enddef
%feature("docstring") PhreeqcRM::GetGfw GetGfw_DOCSTRING


%define GetGridCellCount_DOCSTRING
"Returns the number of grid cells in the user model.

The mapping from grid cells to 
reaction cells is defined by :meth:`CreateMapping`. The number of reaction cells 
may be less than the number of grid cells if there are inactive regions or there 
is symmetry in the model definition.

Returns:
	int: Number of grid cells in the user model."
%enddef
%feature("docstring") PhreeqcRM::GetGridCellCount GetGridCellCount_DOCSTRING


// %define GetIPhreeqcPointer_DOCSTRING
// "Returns an IPhreeqc pointer to the ith IPhreeqc instance in the reaction 
// module. For the threaded version, there are nthreads + 2 IPhreeqc instances, 
// where nthreads is defined in the constructor (:meth:`PhreeqcRM`::PhreeqcRM). The 
// number of threads can be determined by :meth:`GetThreadCount`. The first 
// nthreads (0 based) instances will be the workers, the next (nthreads) is the 
// InitialPhreeqc instance, and the next (nthreads + 1) is the Utility instance. 
// Getting the IPhreeqc pointer for one of these instances allows the user to use 
// any of the IPhreeqc methods on that instance. For MPI, each process has exactly 
// three IPhreeqc instances, one worker (number 0), one InitialPhreeqc instance 
// (number 1), and one Utility instance (number 2).

// Args:
// 	i (int): The number of the IPhreeqc instance (0 based) to be retrieved.

// Returns:
// 	IPhreeqc pointer to the ith IPhreeqc instance (0 based) in the reaction 
// 		module."
// %enddef
// %feature("docstring") PhreeqcRM::GetIPhreeqcPointer GetIPhreeqcPointer_DOCSTRING


%define GetIthConcentration_DOCSTRING
"Returns an array of concentrations for the Ith component. 

Called after :meth:`FindComponents`. Order of components
corresponds to the list of components from
:meth:`GetComponents`.

Returns:
	numpy.ndarray (float): Concentration for each model cell 
		of the Ith component."
%enddef
%feature("docstring") PhreeqcRM::GetIthConcentration GetIthConcentration_DOCSTRING


%define GetIthSpeciesConcentration_DOCSTRING
"Returns an array of concentrations for the Ith aqueous species. 

Called after :meth:`FindComponents`. Order of species
corresponds to the list of species from
:meth:`GetSpeciesNames`. :meth:`SetSpeciesSaveOn` must be set
to True.

Returns:
	numpy.ndarray (float): Concentration for each model cell 
		of the Ith species."
%enddef
%feature("docstring") PhreeqcRM::GetIthSpeciesConcentration GetIthSpeciesConcentration_DOCSTRING


%define GetKineticReactions_DOCSTRING
"Returns a tuple of all kinetic reactions found by
:meth:`FindComponents`.

:meth:`FindComponents` must be called before
:meth:`GetKineticReactions`. This method may be useful when
generating selected output definitions related to kinetic
reactions.

Returns:
	tuple (str): Each string is a unique
		kinetic reaction name."
%enddef
%feature("docstring") PhreeqcRM::GetKineticReactions GetKineticReactions_DOCSTRING


%define GetKineticReactionsCount_DOCSTRING
"Returns the number of kinetic reactions found by
:meth:`FindComponents`.

:meth:`FindComponents` must be called before
:meth:`GetKineticReactionsCount`. This method may be useful
when generating selected output definitions related to
kinetic reactions.

Returns:
	int: The number of kinetic reactions."
%enddef
%feature("docstring") PhreeqcRM::GetKineticReactionsCount GetKineticReactionsCount_DOCSTRING


%define GetMpiMyself_DOCSTRING
"Returns the MPI process (task) number. 

For the MPI version, the root task number is zero, and all
other MPI tasks have unique task numbers greater than zero.
The number of tasks can be obtained with
:meth:`GetMpiTasks`. The number of tasks and computer hosts
is determined at run time by the mpiexec command, and the
number of reaction-module processes is defined by the
communicator used in constructing the reaction modules
(:meth:`PhreeqcRM`::PhreeqcRM). For the OPENMP version, the
task number is always zero, and the result of
:meth:`GetMpiTasks` is one.

Returns:
	int: The MPI task number for a process."
%enddef
%feature("docstring") PhreeqcRM::GetMpiMyself GetMpiMyself_DOCSTRING


%define GetMpiTasks_DOCSTRING
"Returns the number of MPI processes (tasks) assigned to
the reaction module. 

For the MPI version, the number of tasks and computer hosts
is specified at run time by the mpiexec command. The number
of MPI processes used for reaction calculations is
determined by the MPI communicator used in constructing the
reaction modules. The communicator may define as a subset of
the total number of MPI processes. The root task number is
zero, and all other MPI tasks have unique task numbers
greater than zero. For the OPENMP version, the number of
tasks is one, and the task number returned by
:meth:`GetMpiMyself` is zero.

Returns:
	int: The number of MPI processes assigned to the reaction
		module."
%enddef
%feature("docstring") PhreeqcRM::GetMpiTasks GetMpiTasks_DOCSTRING


%define GetNthSelectedOutputUserNumber_DOCSTRING
"Returns the user number for the nth selected-output
definition. 

Definitions are sorted by user number. Phreeqc allows
multiple selected-output definitions, each of which is
assigned a nonnegative integer identifier by the user. The
number of definitions can be obtained by
:meth:`GetSelectedOutputCount`. To cycle through all of the
definitions, GetNthSelectedOutputUserNumber can be used to
identify the user number for each selected-output
definition in sequence.
:meth:`SetCurrentSelectedOutputUserNumber` is then used to
select that user number for selected-output processing.

Args:
	n (int): The sequence number of the selected-output
		definition for which the user number will be returned.
		Fortran, 1 based; C, 0 based.

Returns: 
	int: The user number of the nth selected - output
		definition, negative is failure(See :meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::GetNthSelectedOutputUserNumber GetNthSelectedOutputUserNumber_DOCSTRING

%define GetPartitionUZSolids_DOCSTRING
"Returns the setting for partitioning solids between the
saturated and unsaturated parts of a partially saturated
cell.

The option is intended to be used by saturated-only flow
codes that allow a variable water table.The value has
meaning only when saturations less than 1.0 are
encountered.The partially saturated cells may have a small
water-to-rock ratio that causes reactions to proceed
slowly relative to fully saturated cells. By setting
:meth:`SetPartitionUZSolids` to true, the amounts of solids
and gases are partioned according to the saturation. If a
cell has a saturation of 0.5, then the water interacts with
only half of the solids and gases; the other half is
unreactive until the water table rises. As the saturation
in a cell varies, solids and gases are transferred between
the saturated and unsaturated (unreactive) reservoirs of
the cell. Unsaturated-zone flow and transport codes will
probably use the default (false), which assumes all gases
and solids are reactive regardless of saturation.

Returns:
	Boolean: True, the fraction of solids and gases available
		for reaction is equal to the saturation; False(default),
		all solids and gases are reactive regardless of
		saturation."
%enddef 
%feature("docstring") PhreeqcRM::GetPartitionUZSolids GetPartitionUZSolids_DOCSTRING


%define GetPorosity_DOCSTRING 
"Returns the porosity for each cell. 

By default, the porosity array is initialized with 0.1,
unitless.PhreeqcRM does not change the porosity, so the
values that are retrieved are either the default, or the
values set by the last call to :meth:`SetPorosity`.

Returns: 
	numpy.ndarray (float): The porosities in each cell,
		unitless. Size of the array is nxyz, the number of grid
		cells in the user model (:meth:`GetGridCellCount`)."
%enddef 
%feature("docstring") PhreeqcRM::GetPorosity GetPorosity_DOCSTRING

%define GetPressure_DOCSTRING 
"Returns the pressure for each cell. 

By default, the pressure array is initialized with 1 atm;
if :meth:`SetPressure` has not been called, worker
solutions will have pressures as defined the SOLUTION
definitions at initialization; if :meth:`SetPressure` has
been called, worker solutions will have the pressures as
defined by :meth:`SetPressure`. Pressure effects are
considered by three PHREEQC databases: phreeqc.dat,
Amm.dat, and pitzer.dat.

Returns:
	numpy.ndarray (float): Size of the array is nxyz, the
		number of grid cells in the user model
		(:meth:`GetGridCellCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetPressure GetPressure_DOCSTRING


%define GetPrintChemistryMask_DOCSTRING
"Return an array of print flags that enable or disable
detailed output for each cell.

Printing for a cell will occur only when the printing is
enabled with :meth:`SetPrintChemistryOn`, and the value in
the array for the cell is 1 as set by
:meth:`SetPrintChemistryMask`.

Returns:
	numpy.ndarray (int): Size of array is nxyz, where nxyz is
		the number of grid cells in the user model
		(:meth:`GetGridCellCount`). A value of 0 for a cell
		indicates printing is disabled; a value of 1 for a cell
		indicates printing is enabled."
%enddef
%feature("docstring") PhreeqcRM::GetPrintChemistryMask  GetPrintChemistryMask_DOCSTRING


%define GetPrintChemistryOn_DOCSTRING
"Returns an array of the current print flags for detailed
output for the three sets of IPhreeqc instances: the
workers, the InitialPhreeqc instance, and the Utility
instance.

Dimension of the array is 3. Printing of
detailed output from reaction calculations to the output
file is enabled when the array value is true, disabled
when false. The detailed output prints all of the output
typical of a PHREEQC reaction calculation, which includes
solution descriptions and the compositions of all other
reactants. The output can be several hundred lines per
cell, which can lead to a very large output file
(prefix.chem.txt, :meth:`OpenFiles`). For the worker
instances, the output can be limited to a set of cells
(:meth:`SetPrintChemistryMask`) and, in general, the amount
of information printed can be limited by use of options in
the PRINT data block of PHREEQC (applied by using
:meth:`RunFile` or :meth:`RunString`). Printing the
detailed output for the workers is generally used only for
debugging, and PhreeqcRM will run faster when printing
detailed output for the workers is disabled
(:meth:`SetPrintChemistryOn`).

Returns:
	numpy.ndarray (Boolean): Print flag for the workers, 
		InitialPhreeqc, and Utility IPhreeqc 
		instances, respectively."
%enddef
%feature("docstring") PhreeqcRM::GetPrintChemistryOn GetPrintChemistryOn_DOCSTRING


%define GetRebalanceByCell_DOCSTRING
"Returns the load-rebalancing method used for parallel
processing. 

PhreeqcRM attempts to rebalance the load of each thread or
process such that each thread or process takes the same
amount of time to run its part of a :meth:`RunCells`
calculation. Two algorithms are available: one accounts for
cells that were not run because saturation was zero (true),
and the other uses the average time to run all of the cells
assigned to a process or thread (false), . The methods are
similar, but preliminary results indicate the default is
better in most cases.

Returns:
	Boolean: True indicates individual cell run times are used
		in rebalancing (default); False, indicates average run
		times are used in rebalancing."
%enddef 
%feature("docstring") PhreeqcRM::GetRebalanceByCell GetRebalanceByCell_DOCSTRING

%define GetRebalanceFraction_DOCSTRING 
"Get the fraction used to determine the number of cells to
transfer among threads or processes.

PhreeqcRM attempts to rebalance the load of each thread or
process such that each thread or process takes the same
amount of time to run its part of a :meth:`RunCells`
calculation.The rebalancing transfers cell calculations
among threads or processes to try to achieve an optimum
balance. :meth:`SetRebalanceFraction` adjusts the
calculated optimum number of cell transfers by a fraction
from 0 to 1.0 to determine the number of cell transfers
that actually are made.A value of zero eliminates load
rebalancing. A value less than 1.0 is suggested to avoid
possible oscillations, where too many cells are transferred
at one iteration, requiring reverse transfers at the next
iteration.Default is 0.5.

Returns: 
	float : Fraction used in rebalance,	0.0 to 1.0. "
%enddef
%feature("docstring") PhreeqcRM::GetRebalanceFraction GetRebalanceFraction_DOCSTRING

%define GetSaturationCalculated_DOCSTRING
"Returns an array of saturations as calculated by the reaction 
module.

This method always returns solution_volume / (rv * porosity); 
the method :meth:`SetSaturationUser` has no effect on the
values returned. Reactions will change the volume of
solution in a cell. The transport code must decide whether
to ignore or account for this change in solution volume due
to reactions. Following reactions, the cell saturation is
calculated as solution volume (:meth:`GetSolutionVolume`)
divided by the product of representative volume
(:meth:`SetRepresentativeVolume`) and the porosity
(:meth:`SetPorosity`). The cell saturation returned by
GetSaturationCalculated may be less than or greater than
the saturation set by the transport code
(:meth:`SetSaturationUser`), and may be greater than or less
than 1.0, even in fully saturated simulations. Only the
following databases distributed with PhreeqcRM have molar
volume information needed to accurately calculate solution
volume and saturation: phreeqc.dat, Amm.dat, and
pitzer.dat.

Returns:
	numpy.ndarray (float): Dimension of the array is set to
		nxyz, where nxyz is the number of user grid cells
		(:meth:`GetGridCellCount`). Values for inactive cells are
		set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetSaturationCalculated GetSaturationCalculated_DOCSTRING


%define GetSelectedOutput_DOCSTRING
"Returns the array of selected-output values for the current 
selected-output definition. 

:meth:`SetCurrentSelectedOutputUserNumber` or
:meth:`SetNthSelectedOutput` specifies which of the
selected-output definitions is returned.

Returns:
	numpy.ndarray (float): Size of the array is col times nxyz,
		where col is the number of columns in the selected-output
		definition (:meth:`GetSelectedOutputColumnCount`), and nxyz
		is the number of grid cells in the user model
		(:meth:`GetGridCellCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutput GetSelectedOutput_DOCSTRING


%define GetSelectedOutputColumnCount_DOCSTRING
"Returns the number of columns in the current
selected-output definition.

:meth:`SetCurrentSelectedOutputUserNumber` or
:meth:`SetNthSelectedOutput`specifies which of the
selected-output definitions is used.

Returns:
	int: Number of columns in the current selected-output definition,
		negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutputColumnCount GetSelectedOutputColumnCount_DOCSTRING


%define GetSelectedOutputCount_DOCSTRING
"Returns the number of selected-output definitions. 

:meth:`SetCurrentSelectedOutputUserNumber` or
:meth:`SetNthSelectedOutput` specifies which of the
selected-output definitions is used.

Returns:
	int: Number of selected-output definitions, negative is
		failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutputCount GetSelectedOutputCount_DOCSTRING


// %define GetSelectedOutputHeading_DOCSTRING
// "Returns a selected-output heading. The number of headings is determined 
// by :meth:`GetSelectedOutputColumnCount`. 
// :meth:`SetCurrentSelectedOutputUserNumber` specifies which of the 
// selected-output definitions is used.

// Args:
// 	icol (int): The sequence number of the heading to be retrieved, 
// 		0 based. 
// 	heading(str): A string to receive the heading.

// Returns:
// 	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
// %enddef
// %feature("docstring") PhreeqcRM::GetSelectedOutputHeading GetSelectedOutputHeading_DOCSTRING


%define GetSelectedOutputHeadings_DOCSTRING
"Returns a list of the current selected-output headings. 

The number of headings is determined by
:meth:`GetSelectedOutputColumnCount`.
:meth:`SetCurrentSelectedOutputUserNumber` or
:meth:`SetNthSelectedOutput` are used to
specify which of the selected-output definitions is used.

Returns:
	numpy.ndarray (str): Array of the headings for the 
		current selected output."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutputHeadings GetSelectedOutputHeadings_DOCSTRING


%define GetSelectedOutputOn_DOCSTRING
"Returns the current value of the selected-output-on property. 

A value of true for this property indicates that selected
output data is requested this time step. A value of
false indicates that selected output will not be retrieved
for this time step; processing the selected output is
avoided with some time savings.

Returns: 
	Boolean: True, selected output will be requested; false,
		selected output will not be retrieved."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutputOn GetSelectedOutputOn_DOCSTRING


%define GetSelectedOutputRowCount_DOCSTRING
"Returns the number of rows in the current selected-output
definition. 

The method is included only for convenience; the number of
rows is always equal to the number of grid cells in the
user model (:meth:`GetGridCellCount`).

Returns:
	int: Number of rows in the current selected-output
		definition, negative is failure."
%enddef
%feature("docstring") PhreeqcRM::GetSelectedOutputRowCount GetSelectedOutputRowCount_DOCSTRING


%define GetSICount_DOCSTRING
"Returns the number of phases found by the method
:meth:`FindComponents` for which saturation indices 
could be calculated.

:meth:`FindComponents` must be called before
:meth:`GetSICount`. This method may be useful when
generating selected output definitions related to
saturation indices.

Returns:
	int: The number of phases in for which 
		saturation indices could be calculated."
%enddef
%feature("docstring") PhreeqcRM::GetSICount GetSICount_DOCSTRING


%define GetSINames_DOCSTRING
"Returns a tuple of the names of all
phases for which saturation indices (SIs) could be
calculated.

The list includes all phases that contain only elements
found by the method :meth:`FindComponents`. The list
assumes that all components are present to be able to
calculate the entire list of SIs; it may be that one or
more components are missing in any specific cell.
:meth:`FindComponents` must be called before
:meth:`GetSINames`.This method may be useful when
generating selected output definitions related to
saturation indices.

Returns: 
	tuple (str): Each string is a unique phase name."
%enddef 
%feature("docstring") PhreeqcRM::GetSINames GetSINames_DOCSTRING

%define GetSolidSolutionComponents_DOCSTRING
"Returns a tuple of the names of solid solution components.

The list of solid solution components includes all
components in any SOLID_SOLUTION definitions found by the
method :meth:`FindComponents`. :meth:`FindComponents` must
be called before :meth:`GetSolidSolutionComponents`. This
method may be useful when generating selected output
definitions related to solid solutions.

Returns:
	tuple (str): Each string is a unique solid solution
		component."
%enddef 
%feature("docstring") PhreeqcRM::GetSolidSolutionComponents GetSolidSolutionComponents_DOCSTRING

%define GetSolidSolutionComponentsCount_DOCSTRING
"Returns the number of solid solution components found by
the method :meth:`FindComponents`.

:meth:`FindComponents` must be called before
:meth:`GetSolidSolutionComponentsCount`. This method may be
useful when generating selected output definitions related
to solid solutions.

Returns: 
	int: The number of solid solution components in the
		initial - phreeqc module."
%enddef
%feature("docstring") PhreeqcRM::GetSolidSolutionComponentsCount GetSolidSolutionComponentsCount_DOCSTRING

%define GetSolidSolutionNames_DOCSTRING
"Returns a tuple of solid solution names that correspond
with the solid solution components.

:meth:`FindComponents` must be called before
:meth:`GetSolidSolutionNames`.The solid solution names
tuple is the same length as the solid solution components
tuple and provides the corresponding name of the solid
solution containing the component. This method may be
useful when generating selected output definitions related
to solid solutions.

Returns: 
	tuple (str): Each string is a solid solution name
		corresponding to the solid solution components tuple; a
		solid solution name may occur multiple times."
%enddef
%feature("docstring") PhreeqcRM::GetSolidSolutionNames GetSolidSolutionNames_DOCSTRING


%define GetSolutionVolume_DOCSTRING
"Return an array of the current solution volumes as
calculated by the reaction module.

Dimension of the array will be nxyz, where nxyz is the
number of user grid cells. Values for inactive cells are
set to 1e30. Only the following databases distributed with
PhreeqcRM have molar volume information needed to
accurately calculate solution volume: phreeqc.dat, Amm.dat,
and pitzer.dat.

Returns:
	numpy.ndarray (float): Current solution volumes."
%enddef
%feature("docstring") PhreeqcRM::GetSolutionVolume GetSolutionVolume_DOCSTRING


%define GetSpeciesConcentrations_DOCSTRING
"Returns an array of aqueous species concentrations.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The list of
aqueous species is determined by :meth:`FindComponents` and
includes all aqueous species that can be made from the set
of components. Solution volumes used to calculate mol/L are
calculated by the reaction module. Only the following
databases distributed with PhreeqcRM have molar volume
information needed to accurately calculate solution volume:
phreeqc.dat, Amm.dat, and pitzer.dat.

Returns:
	numpy.ndarray (float): Dimension of the array is nspecies
		times nxyz, where nspecies is the number of aqueous species
		(:meth:`GetSpeciesCount`), and nxyz is the number of grid
		cells (:meth:`GetGridCellCount`). Concentrations are moles
		per liter. Values for inactive cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesConcentrations GetSpeciesConcentrations_DOCSTRING


%define GetSpeciesCount_DOCSTRING
"Returns the number of aqueous species found by the method
:meth:`FindComponents`.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The method
:meth:`FindComponents` must be called before
GetSpeciesCount. 

Returns:
	int: The number of aqueous species."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesCount GetSpeciesCount_DOCSTRING


%define GetSpeciesD25_DOCSTRING
"Returns an array of diffusion coefficients at 25C for the
set of aqueous species.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. Diffusion
coefficients are defined in SOLUTION_SPECIES data blocks,
normally in the database file. Databases distributed with
the reaction module that have diffusion coefficients
defined are phreeqc.dat, Amm.dat, and pitzer.dat.

Returns:
	numpy.ndarray (float): Diffusion coefficients at 25 C,
		m^2/s. Dimension of the array is nspecies, where nspecies
		is the number of aqueous species
		(:meth:`GetSpeciesCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesD25 GetSpeciesD25_DOCSTRING


%define GetSpeciesLog10Gammas_DOCSTRING
"Returns an array of log10 aqueous species activity
coefficients.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The list of
aqueous species is determined by :meth:`FindComponents` and
includes all aqueous species that can be made from the set
of components.

Returns:
	numpy.ndarray (float): Log10 aqueous species activity
		coefficients. Dimension of the array is nspecies times
		nxyz, where nspecies is the number of aqueous species
		(:meth:` GetSpeciesCount`), and nxyz is the number of grid
		cells (:meth:` GetGridCellCount`). Values for inactive
		cells are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesLog10Gammas GetSpeciesLog10Gammas_DOCSTRING


%define GetSpeciesLog10Molalities_DOCSTRING
"Returns an array of log10 aqueous species molalities.

To use this method :meth:`SetSpeciesSaveOn` must be set to
true. The list of aqueous species is determined by
:meth:`FindComponents` and includes all aqueous species
that can be made from the set of components.

Returns:
	numpy.ndarray (float): Log10 aqueous species molalites.
		Dimension of the array is nspecies times nxyz, where
		nspecies is the number of aqueous species
		(:meth:`GetSpeciesCount`), and nxyz is the number of grid
		cells (:meth:`GetGridCellCount`). Values for inactive cells
		are set to 1e30."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesLog10Molalities GetSpeciesLog10Molalities_DOCSTRING


%define GetSpeciesNames_DOCSTRING
"Returns a tuple of the names of the aqueous species.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The list of
aqueous species is determined by :meth:`FindComponents` and
includes all aqueous species that can be made from the set
of components.

Returns:
	tuple (str): The names of the aqueous species. Dimension of
		the tuple is nspecies, where nspecies is the number of
		aqueous species (:meth:`GetSpeciesCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesNames GetSpeciesNames_DOCSTRING


%define GetSpeciesSaveOn_DOCSTRING
"Returns the value of the species-save property. 

By default, concentrations of aqueous species are not
saved. Setting the species-save property to true allows
aqueous species concentrations to be retrieved with
:meth:`GetSpeciesConcentrations`, and solution compositions
to be set with :meth:`SpeciesConcentrations2Module`.

Returns:
	Boolean: True indicates solution species concentrations are
		saved and can be used for multicomponent-diffusion
		calculations; False indicates that solution species
		concentrations are not saved."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesSaveOn GetSpeciesSaveOn_DOCSTRING


%define GetSpeciesStoichiometry_DOCSTRING
"Returns a dictionary with the stoichiometry of each
aqueous species.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. Each species
name is mapped to a dictionary that maps elements in the
species to stoichiometry of the element in the species.

Returns:
	dict (dict):Component names and associated stoichiometric
		coefficients for each aqueous species. Dimension of the
		dict is nspecies, where nspecies is the number of aqueous
		species (:meth:`GetSpeciesCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesStoichiometry GetSpeciesStoichiometry_DOCSTRING


%define GetSpeciesZ_DOCSTRING
"Returns an array of the charge on each aqueous species.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true.

Returns:
	numpy.ndarray (float): The charge on each aqueous species.
		Dimension of the array is nspecies, where nspecies is the
		number of aqueous species (:meth:`GetSpeciesCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetSpeciesZ GetSpeciesZ_DOCSTRING


%define GetStartCell_DOCSTRING
"Returns an array of integers that contains the smallest
reaction-cell number assigned to each worker.

Each worker is assigned a range of reaction-cell numbers
that are run during a call to :meth:`RunCells`. The range
of reaction cell numbers for a worker may vary as load
rebalancing occurs. At any point in the calculations, the
first cell and last cell to be run by a worker can be found
in the arrays returned by GetStartCell and
:meth:`GetEndCell`. Each method returns an array of
integers that has size of the number of threads
(:meth:`GetThreadCount`), if using OPENMP, or the number of
processes (:meth:`GetMpiTasks`), if using MPI.

Returns:
	numpy.ndarray (int): One integer for each worker, that
		gives the first reaction cell to be run by the worker."
%enddef
%feature("docstring") PhreeqcRM::GetStartCell GetStartCell_DOCSTRING


%define GetSurfaceNames_DOCSTRING
"Returns a tuple of surface names (such as 'Hfo') that
correspond with the surface species names.

The tuples referenced by :meth:`GetSurfaceSpecies` and
:meth:`GetSurfaceNames` are the same length.
:meth:`FindComponents` must be called before
:meth:`GetSurfaceNames`. This method may be useful when
generating selected output definitions related to surfaces.

Returns:
	tuple (str)): Each string is a surface name corresponding
		to the corresponding element of thesurface species tuple; a
		surface name may occur multiple times."
%enddef
%feature("docstring") PhreeqcRM::GetSurfaceNames GetSurfaceNames_DOCSTRING


%define GetSurfaceSpecies_DOCSTRING
"Returns a tuple of surface species names (such as
'Hfo_wOH').

The list of surface species is derived from the list of
components (:meth:`FindComponents`) and the list of all
surface site types (such as 'Hfo_w') that are found by the
method :meth:`FindComponents`. :meth:`FindComponents` must
be called before :meth:`GetSurfaceSpecies`. This method may
be useful when generating selected output definitions
related to surfaces.

Returns:
	tuple (str): Each string is a unique surface species name."
%enddef
%feature("docstring") PhreeqcRM::GetSurfaceSpecies GetSurfaceSpecies_DOCSTRING


%define GetSurfaceSpeciesCount_DOCSTRING
"Returns the number of surface species (such as 'Hfo_wOH')
found by the method :meth:`FindComponents`.

:meth:`FindComponents` must be called before
:meth:`GetSurfaceSpeciesCount`. This method may be useful
when generating selected output definitions related to
surfaces.

Returns:
	int: The number of surface species in the initial-phreeqc
		module."
%enddef
%feature("docstring") PhreeqcRM::GetSurfaceSpeciesCount GetSurfaceSpeciesCount_DOCSTRING


%define GetSurfaceTypes_DOCSTRING
"Returns a tuple of surface site types (such as 'Hfo_w')
that correspond with the surface species names (such as
'Hfo_wOH').

The tuples referenced by :meth:`GetSurfaceSpecies` and
:meth:`GetSurfaceTypes` are the same length.
:meth:`FindComponents` must be called before
:meth:`GetSurfaceTypes`. This method may be useful when
generating selected output definitions related to surfaces.

Returns:
	tuple (str)): Each string is a surface site type for the
		corresponding species in the surface species tuple; a
		surface site type may occur multiple times."
%enddef 
%feature("docstring") PhreeqcRM::GetSurfaceTypes GetSurfaceTypes_DOCSTRING

%define GetTemperature_DOCSTRING
"Returns an array of the current temperatures of the cells.

By default, the temperature array is initialized to 25 C;
if :meth:`SetTemperature` has not been called, worker
solutions will have temperatures as defined by the SOLUTION
definitions used to initialize the module; if
:meth:`SetTemperature` has been called, worker solutions
will have the temperatures as defined by
:meth:`SetTemperature`.

Returns:
	numpy.ndarray (float): Temperatures, in degrees C. Size of
		the array is nxyz, where nxyz is the number of grid cells
		in the user model (:meth:` GetGridCellCount`)."
%enddef
%feature("docstring") PhreeqcRM::GetTemperature GetTemperature_DOCSTRING


%define GetThreadCount_DOCSTRING
"Returns the number of threads, which is equal to the
number of workers used to run in parallel with OPENMP. 

For the OPENMP version, the number of threads is set
implicitly or explicitly with the constructor
(:meth:`PhreeqcRM`::PhreeqcRM). For the MPI version, the
number of threads is always one for each process.

Returns:
	int: The number of threads used for OPENMP parallel
		processing."
%enddef
%feature("docstring") PhreeqcRM::GetThreadCount GetThreadCount_DOCSTRING


%define GetTime_DOCSTRING
"Returns the current simulation time in seconds. 

Returns:
	float: The current simulation time, in seconds."
%enddef
%feature("docstring") PhreeqcRM::GetTime GetTime_DOCSTRING


%define GetTimeConversion_DOCSTRING
"Returns a multiplier to convert time from seconds to
another unit, as specified by the user. 

The reaction module uses seconds as the time unit. The user
can set a conversion factor (:meth:`SetTimeConversion`) and
retrieve it with GetTimeConversion. The reaction module
only uses the conversion factor when printing the long
version of cell chemistry (:meth:`SetPrintChemistryOn`),
which is rare. Default conversion factor is 1.0.

Returns:
	float: Multiplier to convert seconds to another 
		time unit."
%enddef
%feature("docstring") PhreeqcRM::GetTimeConversion GetTimeConversion_DOCSTRING


%define GetTimeStep_DOCSTRING
"Returns the current simulation time step in seconds. 

This is the time over which kinetic reactions are
integrated in a call to :meth:`RunCells`.

Returns:
	float: The current simulation time step, in seconds."
%enddef
%feature("docstring") PhreeqcRM::GetTimeStep GetTimeStep_DOCSTRING


%define GetUnitsExchange_DOCSTRING
"Returns the input units for exchangers. 

In PHREEQC input, exchangers are defined by moles of
exchange sites (Mp). :meth:`SetUnitsExchange` specifies how
the number of moles of exchange sites in a reaction cell
(Mc) is calculated from the input value (Mp).

Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity 
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for exchangers."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsExchange GetUnitsExchange_DOCSTRING


%define GetUnitsGasPhase_DOCSTRING
"Returns the input units for gas phases. 

In PHREEQC input, gas phases are defined by moles of
component gases (Mp). :meth:`SetUnitsGasPhase` specifies
how the number of moles of component gases in a reaction
cell (Mc) is calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for gas phases (0, 1, or 2)."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsGasPhase GetUnitsGasPhase_DOCSTRING


%define GetUnitsKinetics_DOCSTRING
"Returns the input units for kinetic reactants. 

In PHREEQC input, kinetics are defined by moles of kinetic
reactants (Mp). :meth:`SetUnitsKinetics` specifies how the
number of moles of kinetic reactants in a reaction cell
(Mc) is calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for kinetic reactants (0, 1, or 2)."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsKinetics GetUnitsKinetics_DOCSTRING


%define GetUnitsPPassemblage_DOCSTRING
"Returns the input units for pure phase assemblages (equilibrium phases). 

In PHREEQC input, equilibrium phases are defined by moles
of each phase (Mp). :meth:`SetUnitsPPassemblage` specifies
how the number of moles of phases in a reaction cell (Mc)
is calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for equilibrium phases (0, 1, or 2)."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsPPassemblage GetUnitsPPassemblage_DOCSTRING


%define GetUnitsSolution_DOCSTRING
"Returns the units of concentration used by the transport
model. Options are 1, mg/L; 2 mol / L; or 3, mass fraction,
kg/kgs.

In PHREEQC, solutions are defined by the number of moles of
each element in the solution. The units of transport
concentration are used when transport concentrations are
converted to solution moles by :meth:`SetConcentrations`
and :meth:`Concentrations2Utility`. The units of solution
concentration also are used when solution moles are
converted to transport concentrations by
:meth:`GetConcentrations`.

To convert from mg/L to moles of element in the
representative volume of a reaction cell, mg/L is converted
to mol/L and multiplied by the solution volume, which is
the product of porosity (:meth:`SetPorosity`), saturation
(:meth:` SetSaturation`), and representative volume
(:meth:` SetRepresentativeVolume`). To convert from mol/L
to moles of element in a cell , mol/L is multiplied by the
solution volume. To convert from mass fraction to moles of
element in a cell, kg /kgs is converted to mol/kgs,
multiplied by density (:meth:`SetDensityUser`) and
multiplied by the solution volume.

To convert from moles of element in the representative
volume of a reaction cell to mg/L, the number of moles of
an element is divided by the solution volume resulting in
mol/L, and then converted to mg/L. To convert from moles of
element in the representative volume of a reaction cell to
mol/L, the number of moles of an element is divided by the
solution volume resulting in mol/L. To convert from moles
of element in the representative volume of a reaction cell
to mass fraction, the number of moles of an element is
converted to kg and divided by the total mass of the
solution. Two options are available for the volume and mass
of solution that are used in converting to transport
concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the
product of porosity, saturation, and representative volume,
and the mass of solution is volume times density as defined
by :meth:`SetDensityUser`. Which option is used is
determined by :meth:`UseSolutionDensityVolume`.

Returns:
	int: Units for solute concentrations in transport model."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsSolution GetUnitsSolution_DOCSTRING


%define GetUnitsSSassemblage_DOCSTRING
"Returns the input units for solid-solution assemblages. 

In PHREEQC input, solid solutions are defined by moles of
each component (Mp). :meth:`SetUnitsSSassemblage` specifies
how the number of moles in a reaction cell (Mc) is
calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for solid solutions (0, 1, or 2)."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsSSassemblage GetUnitsSSassemblage_DOCSTRING


%define GetUnitsSurface_DOCSTRING
"Returns the input units for surfaces. 

In PHREEQC input, surfaces are defined by moles of surface
sites (Mp). :meth:`SetUnitsSurface` specifies how the
number of moles of surface sites in a reaction cell (Mc) is
calculated from the input value(Mp) .

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

Returns:
	int: Input units for surfaces (0, 1, or 2)."
%enddef
%feature("docstring") PhreeqcRM::GetUnitsSurface GetUnitsSurface_DOCSTRING


// %define GetWorkers_DOCSTRING

// "Returns a reference to the array of IPhreeqcPhast
// instances. IPhreeqcPhast inherits from IPhreeqc, and the
// array can be interpreted as a array of pointers to the
// worker, InitialPhreeqc, and Utility IPhreeqc instances. For
// OPENMP, there are nthreads workers, where nthreads is
// defined in the constructor (:meth:`PhreeqcRM`::PhreeqcRM).
// For MPI, there is a single worker. For OPENMP and MPI,
// there is one InitialPhreeqc and one Utility instance.

// Returns:
// 	Array of IPhreeqcPhast instances."
// %enddef
// %feature("docstring") PhreeqcRM::GetWorkers GetWorkers_DOCSTRING


%define InitialEquilibriumPhases2Module_DOCSTRING
"Transfer equilibrium phases definitions from the
InitialPhreeqc instance to the reaction-module workers.

Initial_conditions1 is used to select equilibrium phases
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of EQUILIBRIUM_PHASES index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		equilibrium phases for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialEquilibriumPhases2Module InitialEquilibriumPhases2Module_DOCSTRING


%define InitialExchanges2Module_DOCSTRING
"Transfer exchange definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select exchange definitions
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of EXCHANGE index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		an exchange definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialExchanges2Module InitialExchanges2Module_DOCSTRING


%define InitialGasPhases2Module_DOCSTRING
"Transfer gas phase definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select gas phase definitions
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of GAS_PHASE index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		a gas phase definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialGasPhases2Module InitialGasPhases2Module_DOCSTRING


%define InitializeYAML_DOCSTRING
"A YAML file can be used to initialize an instance of PhreeqcRM. 

Args:
	YAML_file: String containing the YAML file name.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`).

The file contains a YAML sequence of PhreeqcRM methods
and the arguments corresponding to the methods. 

For example,

- key: LoadDatabase
  database: phreeqc.dat
- key: RunFile
  workers: true
  initial_phreeqc: true
  utility: true
  chemistry_name: advect.pqi

:meth:`InitializeYAML` will read the YAML file and execute
the specified methods with the specified arguments. Using
YAML terminology, the argument(s) for a method may be a
scalar, a sequence, or a map, depending if the argument is
a single item, a single array, or there are multiple
arguments. In the case of a map, the names associated with
each argument (for example 'chemistry_name' above) is
arbitrary. The names of the map keys for map arguments are
not used in parsing the YAML file; only the order of the
arguments is important.

The class YAMLPhreeqcRM can be used to write a YAML file.
The methods defined in the YAMLPhreeqcRM class include the
following list; all but SetGridCellCount and ThreadCount
correspond to PhreeqcRM methods.

CloseFiles(void);
CreateMapping(float list, numpy.ndarray, or tuple: grid2chem);
DumpModule();
FindComponents();
InitialEquilibriumPhases2Module(int list, numpy.ndarray, or tuple: ic);
InitialExchanges2Module(int list, numpy.ndarray, or tuple: ic);
InitialGasPhases2Module(int list, numpy.ndarray, or tuple: ic);
InitialKinetics2Module(int list, numpy.ndarray, or tuple: ic);
InitialPhreeqc2Module(int list, numpy.ndarray, or tuple: initial_conditions1);
InitialPhreeqc2Module_mix(int list, numpy.ndarray, or tuple: initial_conditions1,
	int list, numpy.ndarray, or tuple: initial_conditions2,
	float list, numpy.ndarray, or tuple: fraction1);
InitialPhreeqcCell2Module(int n, 
	int list, numpy.ndarray, or tuple: cell_numbers);
InitialSolidSolutions2Module(int list, numpy.ndarray, or tuple: ic);
InitialSolutions2Module(int list, numpy.ndarray, or tuple: ic);
InitialSurfaces2Module(int list, numpy.ndarray, or tuple: ic);
LoadDatabase(string database);
OpenFiles(void);
OutputMessage(string str);
RunCells(void);
RunFile(Boolean workers, 
	Boolean initial_phreeqc, 
	Boolean utility, 
	string chemistry_name);
RunString(Boolean workers, 
	Boolean initial_phreeqc, 
	Boolean utility, 
	string input_string);
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
SetGridCellCount(int nxyz);
SetPartitionUZSolids(Boolean tf);
SetPorosity(float list, numpy.ndarray, or tuple: por);
SetPressure(float list, numpy.ndarray, or tuple: p);
SetPrintChemistryMask(int list, numpy.ndarray, or tuple: cell_mask);
SetPrintChemistryOn(Boolean workers, 
	Boolean initial_phreeqc, 
	Boolean utility);
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

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitializeYAML InitializeYAML_DOCSTRING


%define InitialKinetics2Module_DOCSTRING
"Transfer kinetics definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select kinetics definitions
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of KINETICS index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		a kinetics definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialKinetics2Module InitialKinetics2Module_DOCSTRING


%define InitialPhreeqc2Concentrations_DOCSTRING
"Returns an array of component concentrations from
solutions in the InitialPhreeqc instance. 

The method is used to obtain concentrations for boundary
conditions. If a negative value is used for a cell in
boundary_solution1, then the highest numbered solution in
the InitialPhreeqc instance will be used for that cell.

Args:
	boundary_solution1(float list, numpy.ndarray, or tuple):
		Array of solution index numbers that refer to solutions in
		the InitialPhreeqc instance. Size is n_boundary.

Returns:
	numpy.ndarray (float): The dimension of the array is
		ncomps times n_boundary, where ncomps is the number
		of components returned from :meth:`FindComponents` or
		:meth:`GetComponentCount`, and n_boundary is the dimension
		of the argument boundary_solution1."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2Concentrations InitialPhreeqc2Concentrations_DOCSTRING


%define InitialPhreeqc2Concentrations_mix_DOCSTRING
"Returns component concentrations from solutions in the
InitialPhreeqc instance.

The method is used to obtain concentrations for boundary
conditions that are mixtures of solutions. If a negative
value is used for a cell in boundary_solution1, then the
highest numbered solution in the InitialPhreeqc instance
will be used for that cell. Concentrations may be a mixture
of two solutions, boundary_solution1 and
boundary_solution2, with a mixing fraction for
boundary_solution1 of fraction1 and mixing fraction for
boundary_solution2 of (1 - fraction1). A negative value for
boundary_solution2 implies no mixing, and the associated
value for fraction1 is ignored.

Args:
	boundary_solution1 (int list, numpy.ndarray, or tuple):
		Array of solution index numbers that refer to solutions in
		the InitialPhreeqc instance. Size is designated n_boundary.
	boundary_solution2 (int list, numpy.ndarray, or tuple):
		Array of solution index numbers that that refer to
		solutions in the InitialPhreeqc instance and are defined to
		mix with boundary_solution1. Size is also n_boundary.
	fraction1 (float list, numpy.ndarray, or tuple): Fraction
		of boundary_solution1 that mixes with (1 - fraction1) of
		boundary_solution2. Size is also n_boundary.

Returns:
	numpy.ndarray (float): Array of concentrations extracted
		from the InitialPhreeqc instance. The dimension of the
		array is ncomps times n_boundary, where ncomps is the
		number of components returned from :meth:`FindComponents`
		or :meth:`GetComponentCount`, and n_boundary is the
		dimension of the arrays boundary_solution1,
		boundary_solution2, and fraction1."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2Concentrations_mix InitialPhreeqc2Concentrations_mix_DOCSTRING


%define InitialPhreeqc2Module_DOCSTRING
"Transfer solutions and reactants from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select initial conditions,
including solutions and reactants, for each cell of the
model, without mixing. Initial_conditions1 is dimensioned 7
times nxyz, where nxyz is the number of grid cells in the
user model (:meth:`GetGridCellCount`). The dimension of 7
refers to solutions and reactants in the following order:
(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3)
SURFACE, (4) GAS_PHASE, (5) SOLID_SOLUTIONS, and (6)
KINETICS. The definition initial_solution1[3*nxyz + 99] =
2, indicates that cell 99 (0 based) contains the SURFACE
definition (index 3) defined by SURFACE 2 in the
InitialPhreeqc instance (created in the InitialPhreeqc
instance either by :meth:`RunFile` or :meth:`RunString`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is 7 times
		nxyz. The order of definitions is given above. Negative
		values are ignored, resulting in no definition of that
		entity for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2Module InitialPhreeqc2Module_DOCSTRING


%define InitialPhreeqc2Module_mix_DOCSTRING
"Transfer solutions and reactants from the InitialPhreeqc
instance to the reaction-module workers, possibly with
mixing. 

In its simplest form, initial_conditions1 is used
to select initial conditions, including solutions and
reactants, for each cell of the model, without mixing.
Initial_conditions1 is dimensioned 7 times nxyz, where nxyz
is the number of grid cells in the user model
(:meth:`GetGridCellCount`). The dimension of 7 refers to
solutions and reactants in the following order: (0)
SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3)
SURFACE, (4) GAS_PHASE, (5) SOLID_SOLUTIONS, and (6)
KINETICS. The definition initial_solution1[3*nxyz + 99] =
2, indicates that cell 99 (0 based) contains the SURFACE
definition (index 3) defined by SURFACE 2 in the
InitialPhreeqc instance (either by :meth:`RunFile` or
:meth:`RunString`).

It is also possible to mix solutions and reactants to
obtain the initial conditions for cells. For mixing,
initials_conditions2 contains numbers for a second entity
that mixes with the entity defined in initial_conditions1.
Fraction1 contains the mixing fraction for
initial_conditions1, whereas (1 - fraction1) is the mixing
fraction for initial_conditions2. The definitions
initial_solution1[3*nxyz + 99] = 2,
initial_solution2[3*nxyz + 99] = 3, fraction1[3*nxyz + 99]
= 0.25 indicates that cell 99 (0 based) contains a mixture
of 0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface
compositions have been defined in the InitialPhreeqc
instance. If the user number in initial_conditions2 is
negative, no mixing occurs.

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is 7 times
		nxyz, where nxyz is the number of grid cells in the user
		model (:meth:`GetGridCellCount`). The order of definitions
		is given above. Negative values are ignored, resulting in
		no definition of that entity for that cell.
	initial_conditions2 (int list, numpy.ndarray, or tuple):
		Array of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance. Nonnegative
		values of initial_conditions2 result in mixing with the
		entities defined in initial_conditions1. Negative values
		result in no mixing. Size is 7 times nxyz. The order of
		definitions is given above.
	fraction1 (float list, numpy.ndarray, or tuple): Fraction
		of initial_conditions1 that mixes with (1 - fraction1) of
		initial_conditions2. Size is 7 times nxyz. The order of
		definitions is given above.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2Module_mix InitialPhreeqc2Module_mix_DOCSTRING


%define InitialPhreeqc2SpeciesConcentrations_DOCSTRING
"Returns an array of aqueous species concentrations from
solutions in the InitialPhreeqc instance.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The method is
used to obtain aqueous species concentrations for boundary
conditions. If a negative value is used for a cell in
boundary_solution1, then the highest numbered solution in
the InitialPhreeqc instance will be used for that cell.

Args:
	boundary_solution1 (int list, numpy.ndarray, or tuple):
	Array of solution index numbers that refer to solutions in
	the InitialPhreeqc instance.

Returns:
	numpy.ndarray (float): Array of aqueous concentrations
		extracted from the InitialPhreeqc instance. The dimension
		of the array is nspecies times n_boundary, where nspecies
		is the number of aqueous species returned from
		:meth:`GetSpeciesCount`, and n_boundary is the dimension of
		boundary_solution1."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2SpeciesConcentrations InitialPhreeqc2SpeciesConcentrations_DOCSTRING


%define InitialPhreeqc2SpeciesConcentrations_mix_DOCSTRING
"Returns an array with aqueous species concentrations from
solutions in the InitialPhreeqc instance, possibly with
mixing.

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The method is
used to obtain aqueous species concentrations for boundary
conditions. If a negative value is used for a cell in
boundary_solution1, then the highest numbered solution in
the InitialPhreeqc instance will be used for that cell.
Concentrations may be a mixture of two solutions,
boundary_solution1 and boundary_solution2, with a mixing
fraction for boundary_solution1 of fraction1 and mixing
fraction for boundary_solution2 of (1 - fraction1). A
negative value for boundary_solution2 implies no mixing,
and the associated value for fraction1 is ignored.

Args:
	boundary_solution1 (int list, numpy.ndarray, or tuple):
		Array of solution index numbers that refer to solutions in
		the InitialPhreeqc instance.
	boundary_solution2 (int list, numpy.ndarray, or tuple):
		Array of solution index numbers that refer to solutions in
		the InitialPhreeqc instance and are defined to mix with
		boundary_solution1. Size is same as boundary_solution1.
	fraction1 (float list, numpy.ndarray, or tuple): Array of
		fractions of boundary_solution1 that mix with (1 -
		fraction1) of boundary_solution2. Size is same as
		boundary_solution1.

Returns:
	numpy.ndarray (float): Array of aqueous species concentrations extracted 
		from the InitialPhreeqc instance. The dimension of species_c is nspecies 
		times n_boundary, where nspecies is the number of aqueous species returned 
		from :meth:`GetSpeciesCount`, and n_boundary is the dimension of 
		boundary_solution1."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqc2SpeciesConcentrations_mix InitialPhreeqc2SpeciesConcentrations_mix_DOCSTRING


%define InitialPhreeqcCell2Module_DOCSTRING
"A cell numbered n in the InitialPhreeqc instance is
selected to populate a series of transport cells. 

All reactants with the number n are transferred along with
the solution. If MIX n exists, it is used for the
definition of the solution. If n is negative, n is
redefined to be the largest solution or MIX number in the
InitialPhreeqc instance. All reactants for each cell in the
list cell_numbers are removed before the cell definition is
copied from the InitialPhreeqc instance to the workers.

Args:
	n (int): Number that refers to a solution or MIX and
		associated reactants in the InitialPhreeqc instance.

	cell_numbers (numpy.ndarray (int)): An array of grid-cell
		numbers (user grid-cell numbering system) that will be
		populated with cell n from the InitialPhreeqc instance.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialPhreeqcCell2Module InitialPhreeqcCell2Module_DOCSTRING


%define InitialSolidSolutions2Module_DOCSTRING
"Transfer solid solution definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select solid solution
definitions for each cell of the model. Initial_conditions1
is dimensioned nxyz, where nxyz is the number of grid cells
in the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of SOLID_SOLUTIONS index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		a kinetics definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialSolidSolutions2Module InitialSolidSolutions2Module_DOCSTRING


%define InitialSolutions2Module_DOCSTRING
"Transfer solution definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select solution definitions
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of SOLUTIONS index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		a solution definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialSolutions2Module InitialSolutions2Module_DOCSTRING


%define InitialSurfaces2Module_DOCSTRING
"Transfer surface definitions from the InitialPhreeqc
instance to the reaction-module workers.

Initial_conditions1 is used to select surface definitions
for each cell of the model. Initial_conditions1 is
dimensioned nxyz, where nxyz is the number of grid cells in
the user model (:meth:`GetGridCellCount`).

Args:
	initial_conditions1 (int list, numpy.ndarray, or tuple):
		Array of SURFACE index numbers that refer to
		definitions in the InitialPhreeqc instance. Size is nxyz.
		Negative values are ignored, resulting in no transfer of
		a surface definition for that cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::InitialSolutions2Module InitialSolutions2Module_DOCSTRING


%define LoadDatabase_DOCSTRING
"Load a database for all IPhreeqc instances--workers,
InitialPhreeqc, and Utility. 

All definitions of the reaction module are cleared
(SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the
database is read.

Args:
	database (str): String containing the database name.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::LoadDatabase LoadDatabase_DOCSTRING


%define LogMessage_DOCSTRING
"Print a message to the log file.

Args:
	string (str): String to be printed."
%enddef
%feature("docstring") PhreeqcRM::LogMessage LogMessage_DOCSTRING


// %define MpiAbort_DOCSTRING
// "MPI only. Calls MPI_Abort, which aborts MPI, and makes the reaction module 
// unusable. Should be used only on encountering an unrecoverable error."
// %enddef
// %feature("docstring") PhreeqcRM::MpiAbort MpiAbort_DOCSTRING


%define MpiWorker_DOCSTRING
"MPI only. Nonroot processes (processes with
:meth:`GetMpiMyself` > 0) must call MpiWorker to be able to
respond to messages from the root.

Once MpiWorker has been called, the non-root process can
accept data, perform calculations, and (or) return data
within the reaction module. MpiWorker contains a loop that
reads a message from root, performs a task, and waits for
another message from root. :meth:`SetConcentrations`,
:meth:`RunCells`, and :meth:`GetConcentrations` are
examples of methods that send a message from root to get
the workers to perform a task. The workers will respond to
all methods that are designated 'workers must be in the
loop of MpiWorker' in the MPI section of the method
documentation. The workers will continue to respond to
messages from root until root calls :meth:`MpiWorkerBreak`.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`). MpiWorker returns a value only when
		:meth:`MpiWorkerBreak` is called by root."
%enddef
%feature("docstring") PhreeqcRM::MpiWorker MpiWorker_DOCSTRING


%define MpiWorkerBreak_DOCSTRING
"MPI only. This method is called by root to force nonroot
processes to return from a call to :meth:`MpiWorker`.

:meth:`MpiWorker` contains a loop that reads a message from
root, performs a task, and waits for another message from
root. The workers respond to all methods that are
designated 'workers must be in the loop of MpiWorker' in
the MPI section of the method documentation. The workers
will continue to respond to messages from root until root
calls MpiWorkerBreak.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::MpiWorkerBreak MpiWorkerBreak_DOCSTRING


%define OpenFiles_DOCSTRING
"Opens the output and log files. Files are named
prefix.chem.txt and prefix.log.txt based on the prefix
defined by :meth:`SetFilePrefix`.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::OpenFiles OpenFiles_DOCSTRING


%define OutputMessage_DOCSTRING
"Print a message to the output file.

Args:
	string (str) String to be printed."
%enddef
%feature("docstring") PhreeqcRM::OutputMessage OutputMessage_DOCSTRING


%define RunCells_DOCSTRING
"Runs a reaction step for all reaction cells in the
reaction module. 

Normally, tranport concentrations are transferred to the
reaction cells (:meth:`SetConcentrations`) before reaction
calculations are run. The length of time over which kinetic
reactions are integrated is set by :meth:`SetTimeStep`.
Other properties that may need to be updated as a result of
the transport calculations include porosity
(:meth:`SetPorosity`), saturation
(:meth:`SetSaturationUser`), temperature
(:meth:`SetTemperature`), and pressure
(:meth:`SetPressure`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::RunCells RunCells_DOCSTRING


%define ReturnHandler_DOCSTRING
"Process an IRM_RESULT return code. 

If the return code is nonnegative, no action is taken. If
the return code is negative, the return code is decoded and
printed as an error message along with the second argument
(std::string). On an error, the method will return the same
return code, throw an exception, or exit the program
depending on the setting for :meth:`SetErrorHandlerMode`.

Args:
	result (IRM_RESULT): Return code to be processed.
	e_string (str): Error message to be printed in case of an error.

Returns:
	IRM_RESULT : The first argument to the method is returned."
%enddef 
%feature("docstring") PhreeqcRM::ReturnHandler ReturnHandler_DOCSTRING

%define RunFile_DOCSTRING
"Run a PHREEQC input file. 

The first three arguments determine which IPhreeqc
instances will run the file-- the workers, the
InitialPhreeqc instance, and(or) the Utility instance.
Input files that modify the thermodynamic database should
be run by all three sets of instances. Files with
SELECTED_OUTPUT definitions that will be used during the
time-stepping loop need to be run by the workers. Files
that contain initial conditions or boundary conditions
should be run by the InitialPhreeqc instance.

Args: 
	workers (Boolean): True, the workers will run the file;
		False, the workers will not run the file.

	initial_phreeqc (Boolean): True, the InitialPhreeqc
		instance will run the file; False, the InitialPhreeqc will
		not run the file.

	utility(Boolean): True, the Utility instance will run the file;
		False, the Utility instance will not run the file.

	chemistry_name (str): Name of the file to run.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::RunFile RunFile_DOCSTRING

%define RunString_DOCSTRING
"Run a PHREEQC input string. 

The first three arguments determine which IPhreeqc
instances will run the string-- the workers, the
InitialPhreeqc instance, and(or) the Utility instance. Input
strings that modify the thermodynamic database should be
run by all three sets of instances. Strings with
SELECTED_OUTPUT definitions that will be used during the
time-stepping loop need to be run by the workers.Strings
that contain initial conditions or boundary conditions
should be run by the InitialPhreeqc instance.

Args: 
	workers(Boolean) : True, the workers will run the string;
		False, the workers will not run the string.
	initial_phreeqc(Boolean) : True, the InitialPhreeqc
		instance will run the string; False, the InitialPhreeqc
		will not run the string.
	utility(Boolean) : True, the Utility instance will run the
		string; False, the Utility instance will not run the
		string. input_string (str): String containing PHREEQC
		input.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::RunString RunString_DOCSTRING


%define ScreenMessage_DOCSTRING
"Print message to the screen.

Args:
	string (str): String to be printed."
%enddef
%feature("docstring") PhreeqcRM::ScreenMessage ScreenMessage_DOCSTRING


%define SetComponentH2O_DOCSTRING
"Select whether to include H2O in the component list. 

The concentrations of H and O must be known accurately (8
to 10 significant digits) for the numerical method of
PHREEQC to produce accurate pH and pe values. Because most
of the H and O are in the water species, it may be more
robust (require less accuracy in transport) to transport
the excess H and O (the H and O not in water) and water.
The default setting (true) is to include water, excess H,
and excess O as components. A setting of false will include
total H and total O as components. SetComponentH2O must be
called before :meth:`FindComponents`.

Args:
	tf (Boolean): True (default), excess H, excess O, and water
		are included in the component list; False, total H and O
		are included in the component list.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetComponentH2O SetComponentH2O_DOCSTRING


%define SetConcentrations_DOCSTRING
"Use the array of concentrations (c) to set the moles of
components in each reaction cell. 

The volume of water in a cell is the product of porosity
(:meth:`SetPorosity`), saturation
(:meth:`SetSaturationUser`), and reference volume
(:meth:`SetRepresentativeVolume`). The moles of each
component are determined by the volume of water and per
liter concentrations. If concentration units
(:meth:`SetUnitsSolution`) are mass fraction, the density
(as specified by :meth:`SetDensityUser`) is used to convert
from mass fraction to per mass per liter.

Args:
	c (double list, numpy.ndarray, or tuple): Component
		concentrations. Size of array is ncomps times nxyz, where
		ncomps is the number of components as determined by
		:meth:`FindComponents` or :meth:`GetComponentCount` and
		nxyz is the number of grid cells in the user model
		(:meth:`GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetConcentrations SetConcentrations_DOCSTRING


%define SetCurrentSelectedOutputUserNumber_DOCSTRING
"Select the current selected output by user number. 

The user may define multiple SELECTED_OUTPUT data blocks
for the workers. A user number is specified for each data
block. The value of the argument n_user selects which of
the SELECTED_OUTPUT definitions will be used for
selected-output operations.

Args:
	n_user (int): User number of the SELECTED_OUTPUT data
		block that is to be used.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetCurrentSelectedOutputUserNumber SetCurrentSelectedOutputUserNumber_DOCSTRING


%define SetDensityUser_DOCSTRING
"Set the density for each reaction cell. 

These density values are used when converting from
transported mass-fraction concentrations
(:meth:`SetUnitsSolution`) to produce per liter
concentrations during a call to :meth:`SetConcentrations`.
They are also used when converting from reaction-cell
concentrations to transport concentrations
(:meth:`GetConcentrations`), if
:meth:`UseSolutionDensityVolume` is set to false.

Args:
	density (float list, numpy.ndarray, or tuple): Densities.
		Size of array is nxyz, where nxyz is the number of grid
		cells in the user model (:meth:` GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetDensityUser SetDensityUser_DOCSTRING


%define SetDumpFileName_DOCSTRING
"Set the name of the dump file. It is the name used by
:meth:`DumpModule`.

Args:
	dump_name (str): Name of dump file.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetDumpFileName SetDumpFileName_DOCSTRING


%define SetErrorHandlerMode_DOCSTRING
"Set the action to be taken when the reaction module encounters an error. 

Options are 0, return to calling program with an error
return code (default); 1, throw an exception, in C++, the
exception can be caught, for C and Fortran, the program
will exit; or 2, attempt to exit gracefully.

Args: 
	mode (int): Error handling mode : 0, 1, or 2.

Returns:
	IRM_RESULT: 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetErrorHandlerMode SetErrorHandlerMode_DOCSTRING

%define SetErrorOn_DOCSTRING
"Set the property that controls whether error messages are
generated and displayed.

Messages include PHREEQC 'ERROR' messages, and any messages
written with :meth:`ErrorMessage`.

Args: 
	tf(Boolean) : True, enable error messages; False, disable
		error messages. Default is true.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetErrorOn SetErrorOn_DOCSTRING


%define SetFilePrefix_DOCSTRING
"Set the prefix for the output (prefix.chem.txt) and log
(prefix.log.txt) files. 

These files are opened by :meth:`OpenFiles`.

Args:
	prefix (str): Prefix used when opening the output and log
		files.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetFilePrefix SetFilePrefix_DOCSTRING


%define SetGasCompMoles_DOCSTRING
"Transfer moles of gas components from the array given in
the argument list (gas_moles) to each reaction cell.

Args:
	gas_moles (float list, numpy.ndarray, or tuple): Moles of
		gas components. Dimension of the array is set to
		ngas_comps times nxyz, where, ngas_comps is the result of
		:meth:`GetGasComponentsCount`, and nxyz is the number of
		user grid cells (:meth:`GetGridCellCount`). If the number
		of moles is set to a negative number, the gas component
		will not be defined for the GAS_PHASE of the reaction cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetGasCompMoles SetGasCompMoles_DOCSTRING


%define SetGasPhaseVolume_DOCSTRING
"Transfer volumes of gas phases from the array given in
the argument list (gas_volume) to each reaction cell.

The gas-phase volume affects the gas-component pressures
calculated for fixed-volume gas phases. If a gas-phase
volume is defined with this methood for a GAS_PHASE in a
cell, the gas phase is forced to be a fixed-volume gas
phase.

Args:
	gas_volume (float list, numpy.ndarray, or tuple): Volumes
		for each gas phase. Dimension of the array is nxyz, where
		nxyz is the number of user grid cells
		(:meth:`GetGridCellCount`). If the volume is set to a
		negative number for a cell, the gas-phase volume for that
		cell is not changed.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetGasPhaseVolume SetGasPhaseVolume_DOCSTRING


// %define SetMpiWorkerCallbackC_DOCSTRING
// "MPI and C/C++ only. Defines a callback function that allows additional tasks 
// to be done by the workers. The method :meth:`MpiWorker` contains a loop, where 
// the workers receive a message (an integer), run a function corresponding to 
// that integer, and then wait for another message. SetMpiWorkerCallbackC allows C 
// or C++ developers to add another function that responds to additional integer 
// messages by calling developer-defined functions corresponding to those 
// integers. :meth:`MpiWorker` calls the callback function when the message number 
// is not one of the PhreeqcRM message numbers. Messages are unique integer 
// numbers. PhreeqcRM uses integers in a range beginning at 0. It is suggested 
// that developers use message numbers starting at 1000 or higher for their tasks. 
// The callback function calls a developer-defined function specified by the 
// message number and then returns to :meth:`MpiWorker` to wait for another 
// message.

// In C and C++, an additional pointer can be supplied to find the data necessary 
// to do the task. A void pointer may be set with :meth:`SetMpiWorkerCallbackCookie`. 
// This pointer is passed to the callback function through a void pointer 
// argument in addition to the integer message argument. The pointer may be to a 
// struct or class instance that provides a number of additional pointers to data. 
// :meth:`SetMpiWorkerCallbackCookie` must be called by each worker before :meth:`
// MpiWorker` is called.

// The motivation for this method is to allow the workers to perform other tasks, 
// for instance, parallel transport calculations, within the structure of 
// :meth:`MpiWorker`. The callback function can be used to allow the workers to 
// receive data, perform transport calculations, and (or) send results, without 
// leaving the loop of :meth:`MpiWorker`. Alternatively, it is possible for the 
// workers to return from :meth:`MpiWorker` by a call to :meth:`MpiWorkerBreak` by 
// root. The workers could then call subroutines to receive data, calculate 
// transport, and send data, and then resume processing PhreeqcRM messages from 
// root with another call to :meth:`MpiWorker`.

// Args:
// 	fcn (function pointer): A function that returns an integer and has an 
// 		integer argument and a void * argument.

// Returns:
// 	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
// %enddef
// %feature("docstring") PhreeqcRM::SetMpiWorkerCallbackC SetMpiWorkerCallbackC_DOCSTRING


// %define SetMpiWorkerCallbackCookie_DOCSTRING
// "MPI and C/C++ only. Defines a void pointer that can be used by C and C++ 
// functions called from the callback function (:meth:`SetMpiWorkerCallbackC`) to 
// locate data for a task. The C callback function that is registered with 
// :meth:`SetMpiWorkerCallbackC` has two arguments, an integer message to identify 
// a task, and a void pointer. SetMpiWorkerCallbackCookie sets the value of the 
// void pointer that is passed to the callback function. The void pointer may be a 
// pointer to a struct of class instance that contains additonal pointers to data.

// Args:
// 	cookie           Void pointer that can be used by subroutines called from 
// 		the callback function to locate data needed to perform a task.

// Returns:
// 	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
// %enddef
// %feature("docstring") PhreeqcRM::SetMpiWorkerCallbackCookie SetMpiWorkerCallbackCookie_DOCSTRING


%define SetNthSelectedOutput_DOCSTRING
"Specify the current selected output by sequence number. The user may define 
multiple SELECTED_OUTPUT data blocks for the workers. A user number is 
specified for each data block, and the blocks are stored in user-number order. 
The value of the argument n selects the sequence number of the SELECTED_OUTPUT 
definition that will be used for selected-output operations.

Args:
	n (int): Sequence number of the SELECTED_OUTPUT data block that is to be 
		used.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetNthSelectedOutput SetNthSelectedOutput_DOCSTRING


%define SetPartitionUZSolids_DOCSTRING
"Sets the property for partitioning solids between the
saturated and unsaturated parts of a partially saturated
cell.

The option is intended to be used by saturated-only flow
codes that allow a variable water table. The value has
meaning only when saturations less than 1.0 are
encountered. The partially saturated cells may have a small
water-to-rock ratio that causes reactions to proceed
differently relative to fully saturated cells. By setting
SetPartitionUZSolids to true, the amounts of solids and
gases are partioned according to the saturation. If a cell
has a saturation of 0.5, then the water interacts with only
half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell
varies, solids and gases are transferred between the
saturated and unsaturated ( unreactive) reservoirs of the
cell. Unsaturated-zone flow and transport codes will
probably use the default (false), which assumes all gases
and solids are reactive regardless of saturation.

Args:
	tf (Boolean): True, the fraction of solids and gases
		available for reaction is equal to the saturation; False
		(default), all solids and gases are reactive regardless of
		saturation.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetPartitionUZSolids SetPartitionUZSolids_DOCSTRING


%define SetPorosity_DOCSTRING
"Set the porosity for each reaction cell. 

The volume of water in a reaction cell is the product of
porosity, saturation (:meth:`SetSaturationUser`), and
representative volume (:meth:`SetRepresentativeVolume`).

Args:
	por (float list, numpy.ndarray, or tuple): Porosities,
		unitless. Default is 0.1. Size of array is nxyz, where
		nxyz is the number of grid cells in the user model
		(:meth:`GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetPorosity SetPorosity_DOCSTRING


%define SetPressure_DOCSTRING
"Set the pressure for each reaction cell. 

Pressure effects are considered only in three of the
databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat,
and pitzer.dat. 

Args:
	p (float list, numpy.ndarray, or tuple): Pressures, in atm.
		Size of array is nxyz, where nxyz is the number of grid
		cells in the user model (:meth:` GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetPressure SetPressure_DOCSTRING


%define SetPrintChemistryMask_DOCSTRING
"Enable or disable detailed output for each reaction cell.

Printing for a reaction cell will occur only when the
printing is enabled with :meth:`SetPrintChemistryOn` and
the cell_mask value is 1.

Args:

	cell_mask (Intfloat list, numpy.ndarray, or tuple): Array
		of integers. Size of vector is nxyz, where nxyz is the
		number of grid cells in the user model
		(:meth:`GetGridCellCount`). A value of 0 will disable
		printing detailed output for the cell; a value of 1 will
		enable printing detailed output for a cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetPrintChemistryMask SetPrintChemistryMask_DOCSTRING


%define SetPrintChemistryOn_DOCSTRING
"Set property that enables or disables printing detailed
output from reaction calculations to the output file.

Detailed output will occur for a set of cells defined by
:meth:`SetPrintChemistryMask`. The detailed output prints
all of the output typical of a PHREEQC reaction
calculation, which includes solution descriptions and the
compositions of all other reactants. The output can be
several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, :meth:`OpenFiles`). For
the worker instances, the output can be limited to a set of
cells (:meth:`SetPrintChemistryMask`) and, in general, the
amount of information printed can be limited by use of
options in the PRINT data block of PHREEQC (applied by
using :meth:`RunFile` or :meth:`RunString`). Printing the
detailed output for the workers is generally used only for
debugging, and PhreeqcRM will run significantly faster when
printing detailed output for the workers is disabled.

Args:
	workers (Boolean): True, enable detailed printing in the
		worker instances; False, disable detailed printing in the
		worker instances.
	initial_phreeqc(Boolean) : True, enable detailed printing
		in the InitialPhreeqc instance; False, disable detailed
		printing in the InitialPhreeqc instance.
	utility(Boolean) : True, enable detailed printing in the
		Utility instance; False, disable detailed printing in the
		Utility instance.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetPrintChemistryOn SetPrintChemistryOn_DOCSTRING

%define SetRebalanceByCell_DOCSTRING
"Set the load-balancing algorithm. 

PhreeqcRM attempts to rebalance the load of each thread or
process such that each thread or process takes the same
amount of time to run its part of a :meth:`RunCells`
calculation. Two algorithms are available; one uses
individual times for each cell and accounts for cells that
were not run because saturation was zero (default), and the
other assigns an average time to all cells. The methods are
similar, but limited testing indicates the default method
performs better.

Args:
	tf (Boolean): True, indicates individual cell times are
		used in rebalancing (default); False, indicates average
		times are used in rebalancing.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetRebalanceByCell SetRebalanceByCell_DOCSTRING


%define SetRebalanceFraction_DOCSTRING
"Sets the fraction of cells that are transferred among
threads or processes when rebalancing. 

PhreeqcRM attempts to rebalance the load of each thread or
process such that each thread or process takes the same
amount of time to run its part of a :meth:`RunCells`
calculation. The rebalancing transfers cell calculations
among threads or processes to try to achieve an optimum
balance. SetRebalanceFraction adjusts the calculated
optimum number of cell transfers by a fraction from 0 to
1.0 to determine the actual number of cell transfers. A
value of zero eliminates load rebalancing. A value less
than 1.0 is suggested to slow the approach to the optimum
cell distribution and avoid possible oscillations when too
many cells are transferred at one iteration, requiring
reverse transfers at the next iteration. Default is 0.5.

Args:
	f (float): Fraction from 0.0 to 1.0.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetRebalanceFraction SetRebalanceFraction_DOCSTRING


%define SetRepresentativeVolume_DOCSTRING
"Set the representative volume of each reaction cell. 

By default the representative volume of each reaction cell
is 1 liter. The volume of water in a reaction cell is
determined by the product of the representative volume, the
porosity (:meth:`SetPorosity`), and the saturation
(:meth:`SetSaturationUser`). The numerical method of
PHREEQC is more robust if the water volume for a reaction
cell is within a couple orders of magnitude of 1.0. Small
water volumes caused by small porosities and (or) small
saturations (and (or) small representative volumes) may
cause non-convergence of the numerical method. In these
cases, a larger representative volume may help. Note that
increasing the representative volume also increases the
number of moles of the reactants in the reaction cell
(minerals, surfaces, exchangers, and others).

Args:
	rv (float list, numpy.ndarray, or tuple): Representative
		volumes, in liters. Default is 1.0 liter. Size of array is
		nxyz, where nxyz is the number of grid cells in the user
		model (:meth:`GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetRepresentativeVolume SetRepresentativeVolume_DOCSTRING


%define SetSaturationUser_DOCSTRING
"Set the saturation of each reaction cell. 

Saturation is a fraction ranging from 0 to 1. The volume of
water in a cell is the product of porosity
(:meth:`SetPorosity`), saturation (SetSaturationUser), and
representative volume (:meth:`SetRepresentativeVolume`). As
a result of a reaction calculation, solution properties
(density and volume) will change; the databases
phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume
data to calculate these changes. The methods
:meth:`GetDensityCalculated`, :meth:`GetSolutionVolume`,
and :meth:`GetSaturationCalculated` can be used to account
for these changes in the succeeding transport calculation.
SetRepresentativeVolume should be called before initial
conditions are defined for the reaction cells.

Args:
	sat (float list, numpy.ndarray, or tuple): Saturations,
		unitless. Default 1.0. Size of vector is nxyz, where nxyz
		is the number of grid cells in the user model
		(:meth:`GetGridCellCount`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetSaturationUser SetSaturationUser_DOCSTRING


%define SetScreenOn_DOCSTRING
"Set the property that controls whether messages are
written to the screen.

Messages include information about rebalancing during
:meth:`RunCells`, and any messages written with
:meth:`ScreenMessage`.

Args:
	tf (Boolean): True, enable screen messages; False, disable
		screen messages.Default is true.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetScreenOn SetScreenOn_DOCSTRING

%define SetSelectedOutputOn_DOCSTRING
"Set the property that controls whether selected-output
results are available to be retrieved

If the property is set to True, selected-output results are
accumulated during :meth:`RunCells` and can be retrieved
with :meth:`GetSelectedOutput`.

Args: 
	tf(Boolean) : True, enable selected output; False, disable
		selected output.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetSelectedOutputOn SetSelectedOutputOn_DOCSTRING

%define SetSpeciesSaveOn_DOCSTRING
"Sets the value of the species-save property. 

This method enables or disables saving individual aqueous
species concentrations. By default, concentrations of
aqueous species are not saved. Setting the species save
property to true allows aqueous species concentrations to
be retrieved with :meth:`GetSpeciesConcentrations`, and
solution compositions to be set with
:meth:`SpeciesConcentrations2Module`.SetSpeciesSaveOn must
be called before calls to :meth:`FindComponents`.
Individual species concentrations are required for
multicomponent diffusion calculations.

Args: 
	save_on (Boolean): True indicates species concentrations are saved;
		False indicates species concentrations are not saved.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetSpeciesSaveOn SetSpeciesSaveOn_DOCSTRING

%define SetTemperature_DOCSTRING
"Set the temperature for each reaction cell. 

If SetTemperature is not called, worker solutions will have
temperatures as defined by initial conditions(
:meth:`InitialPhreeqc2Module`,
:meth:`InitialPhreeqcCell2Module`, or
:meth:`InitialSoluti
ons2Module`).

Args: 
	t (float list, numpy.ndarray, or tuple): Temperatures, in
		degrees C. Size of vector is nxyz, where nxyz is the number
		of grid cells in the user model( :meth:`GetGridCellCount`).

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetTemperature SetTemperature_DOCSTRING

%define SetTime_DOCSTRING
"Set current simulation time for the reaction module.

Args: 
	time(float) : Current simulation time, in seconds.

Returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetTime SetTime_DOCSTRING

%define SetTimeConversion_DOCSTRING
"Set a factor to convert from seconds to user time units.
Factor times seconds produces user time units.

args: 
	conv_factor (float): Factor to convert seconds to user time units.

returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetTimeConversion SetTimeConversion_DOCSTRING

%define SetTimeStep_DOCSTRING
"Set current time step for the reaction module. This is the
length of time over which kinetic reactions are integrated.

args: 
	time_step (float): Time step, in seconds.

returns: 
	IRM_RESULT : 0 is success, negative is failure(See
		:meth:`DecodeError`) ."
%enddef
%feature("docstring") PhreeqcRM::SetTimeStep SetTimeStep_DOCSTRING

%define SetUnitsExchange_DOCSTRING
"Sets input units for exchangers. 

In PHREEQC input, exchangers are defined by moles of
exchange sites(Mp) .SetUnitsExchange specifies how the
number of moles of exchange sites in a reaction cell(Mc) is
calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single EXCHANGE definition is used for cells with
different initial porosity, the three options scale quite
differently. For option 0, the number of moles of
exchangers will be the same regardless of porosity. For
option 1, the number of moles of exchangers will be vary
directly with porosity and inversely with rock volume. For
option 2, the number of moles of exchangers will vary
directly with rock volume and inversely with porosity.

Args:
	option (int): Units option for exchangers: 0, 1, or 2.

Returns:

	IRM_RESULT: 0 is success, negative is failure (See
		:meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsExchange SetUnitsExchange_DOCSTRING


%define SetUnitsGasPhase_DOCSTRING
"Set input units for gas phases. 

In PHREEQC input, gas phases are defined by moles of
component gases (Mp). SetUnitsGasPhase specifies how the
number of moles of component gases in a reaction cell (Mc)
is calculated from the input value(Mp) .

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single GAS_PHASE definition is used for cells with
different initial porosity, the three options scale quite
differently. For option 0, the number of moles of a gas
component will be the same regardless of porosity. For
option 1, the number of moles of a gas component will be
vary directly with porosity and inversely with rock volume.
For option 2, the number of moles of a gas component will
vary directly with rock volume and inversely with porosity.

Args:
	option (int): Units option for gas phases: 0, 1, or 2.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsGasPhase SetUnitsGasPhase_DOCSTRING


%define SetUnitsKinetics_DOCSTRING
"Set input units for kinetic reactants.

In PHREEQC input, kinetics are defined by moles of kinetic
reactants (Mp). SetUnitsKinetics specifies how the number
of moles of kinetic reactants in a reaction cell (Mc) is
calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single KINETICS definition is used for cells with
different initial porosity, the three options scale quite
differently. For option 0, the number of moles of kinetic
reactants will be the same regardless of porosity. For
option 1, the number of moles of kinetic reactants will be
vary directly with porosity and inversely with rock volume.
For option 2, the number of moles of kinetic reactants will
vary directly with rock volume and inversely with porosity.

Note that the volume of water in a cell in the reaction
module is equal to the product of porosity
(:meth:`SetPorosity`), the saturation
(:meth:`SetSaturationUser`), and representative volume
(:meth: `SetRepresentativeVolume`), which is usually less
than 1 liter. It is important to write the RATES
definitions for homogeneous (aqueous) kinetic reactions to
account for the current volume of water, often by
calculating the rate of reaction per liter of water and
multiplying by the volume of water (Basic function
SOLN_VOL).

Rates that depend on surface area of solids, are not
dependent on the volume of water. However, it is important
to get the correct surface area for the kinetic reaction.
To scale the surface area with the number of moles, the
specific area (m^2 per mole of reactant) can be defined as
a parameter (KINETICS; -parm), which is multiplied by the
number of moles of reactant (Basic function M) in RATES to
obtain the surface area.

Args:
o	ption (int): Units option for kinetic reactants: 0, 1, or 2.

Returns:
I	RM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsKinetics SetUnitsKinetics_DOCSTRING


%define SetUnitsPPassemblage_DOCSTRING
"Set input units for pure phase assemblages (equilibrium phases). 

In PHREEQC input, equilibrium phases are defined by moles
of each phase (Mp). SetUnitsPPassemblage specifies how the
number of moles of phases in a reaction cell (Mc) is
calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single EQUILIBRIUM_PHASES definition is used for cells
with different initial porosity, the three options scale
quite differently. For option 0, the number of moles of a
mineral will be the same regardless of porosity. For option
1, the number of moles of a mineral will be vary directly
with porosity and inversely with rock volume. For option 2,
the number of moles of a mineral will vary directly with
rock volume and inversely with porosity.

Args:
	option (int): Units option for equilibrium phases: 0, 1, or 2.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsPPassemblage SetUnitsPPassemblage_DOCSTRING


%define SetUnitsSolution_DOCSTRING
"Solution concentration units used by the transport model. 

Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
PHREEQC defines solutions by the number of moles of each
element in the solution.

To convert from mg/L to moles of element in the
representative volume of a reaction cell, mg/L is converted
to mol/L and multiplied by the solution volume, which is
the product of porosity (:meth:`SetPorosity`), saturation
(:meth:` SetSaturationUser`), and representative volume
(:meth:`SetRepresentativeVolume`). To convert from mol/L to
moles of element in the representative volume of a reaction
cell, mol/L is multiplied by the solution volume. To
convert from mass fraction to moles of element in the
representative volume of a reaction cell, kg /kgs is
converted to mol/kgs, multiplied by density
(:meth:`SetDensityUser`) and multiplied by the solution
volume.

To convert from moles of element in the representative
volume of a reaction cell to mg/L, the number of moles of
an element is divided by the solution volume resulting in
mol/L, and then converted to mg/L. To convert from moles of
element in a cell to mol/L, the number of moles of an
element is divided by the solution volume resulting in
mol/L. To convert from moles of element in a cell to mass
fraction, the number of moles of an element is converted to
kg and divided by the total mass of the solution. Two
options are available for the volume and mass of solution
that are used in converting to transport concentrations:
(1) the volume and mass of solution are calculated by
PHREEQC, or (2) the volume of solution is the product of
porosity (:meth:`SetPorosity`), saturation
(:meth:`SetSaturationUser`), and representative volume
(:meth:`SetRepresentativeVolume`), and the mass of solution
is volume times density as defined by
:meth:`SetDensityUser`. Which option is used is determined
by :meth:`UseSolutionDensityVolume`.

Args:
	option (int): Units option for solutions: 1, 2, or 3, default is 1, mg/L.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsSolution SetUnitsSolution_DOCSTRING


%define SetUnitsSSassemblage_DOCSTRING
"Set input units for solid-solution assemblages. 

In PHREEQC, solid solutions are defined by moles of each
component (Mp). SetUnitsSSassemblage specifies how the
number of moles of solid-solution components in a reaction
cell (Mc) is calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single SOLID_SOLUTION definition is used for cells
with different initial porosity, the three options scale
quite differently. For option 0, the number of moles of a
solid-solution component will be the same regardless of
porosity. For option 1, the number of moles of a
solid-solution component will be vary directly with
porosity and inversely with rock volume. For option 2, the
number of moles of a solid-solution component will vary
directly with rock volume and inversely with porosity.

Args:
	option (int): Units option for solid solutions: 0, 1, or 2.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsSSassemblage SetUnitsSSassemblage_DOCSTRING


%define SetUnitsSurface_DOCSTRING
"Set input units for surfaces. 

In PHREEQC input, surfaces are defined by moles of surface
sites (Mp). SetUnitsSurface specifies how the number of
moles of surface sites in a reaction cell (Mc) is
calculated from the input value(Mp).

Options are 0, Mp is mol / L of RV(default), Mc = Mp * RV,
where RV is the representative volume
(:meth:`SetRepresentativeVolume`); 1, Mp is mol / L of
water in the RV, Mc = Mp * P * RV, where P is porosity
(:meth:`SetPorosity`); or 2, Mp is mol/L of rock in the RV,
Mc = Mp*(1-P)*RV.

If a single SURFACE definition is used for cells with
different initial porosity, the three options scale quite
differently. For option 0, the number of moles of surface
sites will be the same regardless of porosity. For option
1, the number of moles of surface sites will be vary
directly with porosity and inversely with rock volume. For
option 2, the number of moles of surface sites will vary
directly with rock volume and inversely with porosity.

Args:
	option (int): Units option for surfaces: 0, 1, or 2.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SetUnitsSurface SetUnitsSurface_DOCSTRING


%define SpeciesConcentrations2Module_DOCSTRING
"Set solution concentrations in the reaction module based on
the array of aqueous species concentrations
(species_conc). 

This method is intended for use with
multicomponent-diffusion transport calculations, and
:meth:`SetSpeciesSaveOn` must be set to true. The list of
aqueous species is determined by :meth:`FindComponents` and
includes all aqueous species that can be made from the set
of components. The method determines the total
concentration of a component by summing the molarities of
the individual species times the stoichiometric coefficient
of the elements in each species. Solution compositions in
the reaction cells are updated with these component
concentrations.

Args:
	species_conc (float list, numpy.ndarray, or tuple): Aqueous
		species concentrations. Dimension of the array is nspecies
		times nxyz, where nspecies is the number of aqueous species
		(:meth:`GetSpeciesCount`), and nxyz is the number of user
		grid cells (:meth:`GetGridCellCount`). Concentrations are
		moles per liter.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::SpeciesConcentrations2Module SpeciesConcentrations2Module_DOCSTRING


%define StateSave_DOCSTRING
"Save the state of the chemistry in all model cells,
including SOLUTIONs, EQUILIBRIUM_PHASES, EXCHANGEs,
GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.


Although not generally used, MIXes, REACTIONs,
REACTION_PRESSUREs, and REACTION_TEMPERATUREs will be saved
for each cell, if they have been defined in the worker
IPhreeqc instances. The distribution of cells among the
workers and the chemistry of fully or partially unsaturated
cells are also saved. The state is saved in memory; use
:meth:`DumpModule` to save the state to file. PhreeqcRM can
be reset to this state by using :meth:`StateApply`. A state
is identified by an integer, and multiple states can be
saved.

Args:
	istate (int): Integer identifying the state that is saved. 

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::StateSave StateSave_DOCSTRING


%define StateApply_DOCSTRING
"Reset the state of the module to a state previously saved
with :meth:`StateSave`.

The chemistry of all model cells are reset, including
SOLUTIONs, EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs,
KINETICS, SOLID_SOLUTIONs, and SURFACEs. MIXes, REACTIONs,
REACTION_PRESSUREs, and REACTION_TEMPERATUREs will be reset
for each cell, if they were defined in the worker IPhreeqc
instances at the time the state was saved. The distribution
of cells among the workers and the chemistry of fully or
partially unsaturated cells are also reset to the saved
state. The state to be applied is identified by an integer.

Args:
	istate (int): Integer identifying the state that 
		is to be applied.

Returns:
	IRM_RESULT: 0 is success, negative is failure 
		(See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::StateApply StateApply_DOCSTRING


%define StateDelete_DOCSTRING
"Delete a state previously saved with :meth:`StateSave`.

Args:
	istate (int): Integer identifying the state that is to be deleted.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef
%feature("docstring") PhreeqcRM::StateDelete StateDelete_DOCSTRING


%define UseSolutionDensityVolume_DOCSTRING
"Determines the volume and density to use when converting
from the reaction-cell concentrations to transport
concentrations (:meth:`GetConcentrations`). 

Two options are available to convert concentration units:
(1) the density and solution volume calculated by PHREEQC
are used, or (2) the specified density
(:meth:`SetDensityUser`) and solution volume are determined
by the product of saturation (:meth:`SetSaturationUser`),
porosity (:meth:`SetPorosity`), and representative volume
(:meth:`SetRepresentativeVolume`). Transport models that
consider density-dependent flow will probably use the
PHREEQC-calculated density and solution volume (default),
whereas transport models that assume constant-density flow
will probably use specified values of density and solution
volume. Only the following databases distributed with
PhreeqcRM have molar-volume information needed to
accurately calculate density and solution volume:
phreeqc.dat, Amm.dat, and pitzer.dat. Density is only used
when converting to or from transport units of mass
fraction.

Args:
	tf (Boolean): True indicates that the solution density and
		volume as calculated by PHREEQC will be used to calculate
		concentrations. False indicates that the solution density
		set by :meth:`SetDensityUser` and the volume determined by
		the product of :meth:`SetSaturationUser`,
		:meth:`SetPorosity`, and :meth:`SetRepresentativeVolume`,
		will be used to calculate concentrations retrieved by
		:meth:`GetConcentrations`."
%enddef
%feature("docstring") PhreeqcRM::UseSolutionDensityVolume UseSolutionDensityVolume_DOCSTRING


%define WarningMessage_DOCSTRING
"Print a warning message to the screen and the log file.

Args:
	warnstr (str): String to be printed."
%enddef
%feature("docstring") PhreeqcRM::WarningMessage WarningMessage_DOCSTRING
