
/**
%define GetGridCellCountYAML_DOCSTRING
'@a GetGridCellCountYAML will read the YAML file and extract the value
of GridCellCount, which can be used to construct a PhreeqcRM
instance. The constructor for a PhreeqcRM instance requires a 
value for the number of cells. If a GUI or preprocessor is
used to write a YAML file to initialize PhreeqcRM, the number
of cells can be written to the YAML file and extracted with
this method.
Args:
	YAML_file(string): String containing the YAML file name.
Returns:
	int: Number of grid cells specified in the YAML file; returns zero if GridCellCount is not defined.'
%enddef
%feature PhreeqcRM::GetGridCellCountYAML GetGridCellCountYAML_DOCSTRING
*/
int IRM_DLL_EXPORT GetGridCellCountYAML(const char* YAML_file);

/**
%define PhreeqcRM_DOCSTRING
'Constructor for the PhreeqcRM reaction module. If the code is compiled with
the preprocessor directive USE_OPENMP, the reaction module use OPENMP and multiple threads.
If the code is compiled with the preprocessor directive USE_MPI, the reaction
module will use MPI and multiple processes. If neither preprocessor directive is used,
the reaction module will be serial (unparallelized).
Args:
	nxyz(int): The number of grid cells in the users model. 
	thread_count_or_communicator(int): If multithreaded, the number of threads to use in parallel segments of the code.
		If @a thread_count_or_communicator is <= 0, the number of threads is set equal to the number of processors in the computer.
		If multiprocessor, the MPI communicator to use within the reaction module.
	io(PHRQ_io): Optionally, a PHRQ_io input/output object can be provided to the constructor. By default
		a PHRQ_io object is constructed to handle reading and writing files.'
%enddef
%feature PhreeqcRM::PhreeqcRM PhreeqcRM_DOCSTRING
*/
PhreeqcRM(int nxyz, MP_TYPE thread_count_or_communicator, PHRQ_io * io=NULL);

/**
%define CloseFiles_DOCSTRING
'Close the output and log files.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::CloseFiles CloseFiles_DOCSTRING
*/
IRM_RESULT                                CloseFiles(void);
/**
%define Concentrations2Utility_DOCSTRING
'@a N sets of component concentrations are converted to SOLUTIONs numbered 1-@a n in the Utility IPhreeqc.
The solutions can be reacted and manipulated with the methods of IPhreeqc. If solution concentration units
(:meth: `SetUnitsSolution`) are per liter, one liter of solution is created in the Utility instance; if solution
concentration units are mass fraction, one kilogram of solution is created in the Utility instance.
The motivation for this
method is the mixing of solutions in wells, where it may be necessary to calculate solution properties
(pH for example) or react the mixture to form scale minerals.
The code fragment below makes a mixture of
concentrations and then calculates the pH of the mixture.
Args:
	c(DoubleVector): Vector of concentrations to be made SOLUTIONs in Utility IPhreeqc.
		Vector contains @a n values for each component (:meth: `GetComponentCount`) in sequence.
	tc(DoubleVector): Vector of temperatures to apply to the SOLUTIONs, in degrees C. Vector of size @a n.
	p_atm(DoubleVector): Vector of pressures to apply to the SOLUTIONs, in atm. Vector of size @a n.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::Concentrations2Utility Concentrations2Utility_DOCSTRING
*/
IPhreeqc * Concentrations2Utility(std::vector< double > &c,
std::vector< double > tc, std::vector< double > p_atm);
/**
%define CreateMapping_DOCSTRING
'Provides a mapping from grid cells in the user's model to reaction cells for which chemistry needs to be run.
The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells
for which chemistry must be run. 
The array @a grid2chem of size @a nxyz (the number of grid cells, :meth: `GetGridCellCount`) 
must contain the set of all integers 0 <= @a i < @a count_chemistry, 
where @a count_chemistry is a number less than or equal to @a nxyz.
Inactive cells are assigned a negative integer. 
The mapping may be many-to-one to account for symmetry.
Default is a one-to-one mapping--all user grid cells are reaction cells
(equivalent to @a grid2chem values of 0,1,2,3,...,nxyz-1).
Args:
	grid2chem(DoubleVector): A vector of integers: Nonnegative is a reaction-cell number (0 based),
		negative is an inactive cell. Vector is of size @a nxyz (number of grid cells, :meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::CreateMapping CreateMapping_DOCSTRING
*/
IRM_RESULT                                CreateMapping(std::vector< int > &grid2chem);
/**
%define DecodeError_DOCSTRING
'If @a result is negative, this method prints an error message corresponding to IRM_RESULT @a result.
If @a result is non-negative, no action is taken.
Args:
	result(IRM_RESULT): An IRM_RESULT value returned by one of the reaction-module methods.
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
@endhtmlonly'
%enddef
%feature PhreeqcRM::DecodeError DecodeError_DOCSTRING
*/
void                                      DecodeError(int result);
/**
%define DumpModule_DOCSTRING
'Writes the contents of all workers to file in _RAW formats (see appendix of PHREEQC version 3 manual),
including SOLUTIONs and all reactants.
Args:
	dump_on(Boolean): Signal for writing the dump file, true or false.
	append(Boolean): Signal to append to the contents of the dump file, true or false.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::DumpModule DumpModule_DOCSTRING
*/
IRM_RESULT                                DumpModule(bool dump_on, bool append = false);
/**
%define ErrorHandler_DOCSTRING
'Checks @a result for an error code. If result is negative, the result is decoded (:meth: `DecodeError`),
and printed as an error message along with the @a e_string, and an exception is thrown. If the result
is nonnegative, no action is taken.
Args:
	result(IRM_RESULT): IRM_RESULT to be checked for an error.
	e_string(string): String to be printed if an error is found.'
%enddef
%feature PhreeqcRM::ErrorHandler ErrorHandler_DOCSTRING
*/
void                                      ErrorHandler(int result, const std::string &e_string);
/**
%define ErrorMessage_DOCSTRING
'Send an error message to the screen, the output file, and the log file.
Args:
	error_string(string): String to be printed.
	prepend(Boolean): True, prepends @a error_string with "Error: "; false, @a error_string is used with no prepended text.'
%enddef
%feature PhreeqcRM::ErrorMessage ErrorMessage_DOCSTRING
*/
void                                      ErrorMessage(const std::string &error_string, bool prepend = true);
/**
%define FindComponents_DOCSTRING
'This method accumulates a list of elements. Elements are those that have been
defined in a solution or any other reactant
(EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
This method can be called multiple times and the list that is created is cummulative.
The list is the set of components that needs to be transported. By default the list
includes water, excess H and excess O (the H and O not contained in water);
alternatively, the list may be set to contain total H and total O (:meth: `SetComponentH2O`),
which requires transport results to be accurate to eight or nine significant digits.
If multicomponent diffusion (MCD) is to be modeled,
there is a capability to retrieve aqueous species concentrations
(:meth: `GetSpeciesConcentrations`) and to set new solution concentrations after
MCD by using individual species concentrations
(:meth: `SpeciesConcentrations2Module`).
To use these methods the save-species property needs to be turned on (:meth: `SetSpeciesSaveOn`).
If the save-species property is on, FindComponents will generate
a list of aqueous species (:meth: `GetSpeciesCount`, :meth: `GetSpeciesNames`), 
their diffusion coefficients at 25 C (:meth: `GetSpeciesD25`),
and their charge (:meth: `GetSpeciesZ`).
Returns:
	int: Number of components currently in the list, or IRM_RESULT error code 
		(negative value, see :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::FindComponents FindComponents_DOCSTRING
*/
int                                       FindComponents();
/**
%define GetBackwardMapping_DOCSTRING
'Returns a vector of vectors,
where the @a nth vector is a vector of grid-cell numbers
that are mapped to reaction-cell number @a n.
Each reaction-cell number has a vector of one or more grid-cell numbers.
Returns:
	IntVector: Vector of vectors of ints. For each reaction cell @a n,
		the @a nth vector in the vector of vectors contains
		the grid-cell numbers that map to the reaction cell.'
%enddef
%feature PhreeqcRM::GetBackwardMapping GetBackwardMapping_DOCSTRING
*/
const std::vector < std::vector <int> > & GetBackwardMapping(void) {return this->backward_mapping;}
/**
%define GetChemistryCellCount_DOCSTRING
'Returns the number of reaction cells in the reaction module. The number of reaction cells is defined by
the set of non-negative integers in the mapping from grid cells (:meth: `CreateMapping`), or, by default,
the number of grid cells (:meth: `GetGridCellCount`).
The number of reaction cells is less than or equal to the number of grid cells in the user's model.
Returns:
	int: Number of reaction cells.'
%enddef
%feature PhreeqcRM::GetChemistryCellCount GetChemistryCellCount_DOCSTRING
*/
int                                       GetChemistryCellCount(void) const {return this->count_chemistry;}
/**
%define GetComponentCount_DOCSTRING
'Returns the number of components in the reaction-module component list.
Returns:
	int: The number of components in the reaction-module component list. The component list is
		generated by calls to :meth: `FindComponents`.
		The return value from the last call to :meth: `FindComponents` is equal to the return value from GetComponentCount.'
%enddef
%feature PhreeqcRM::GetComponentCount GetComponentCount_DOCSTRING
*/
int                                       GetComponentCount(void) const {return (int) this->components.size();}
/**
%define GetComponents_DOCSTRING
'Returns a reference to the reaction-module component list that was generated by calls to :meth: `FindComponents`.
Returns:
	tuple of strings: A vector of strings; each string is a component name.'
%enddef
%feature PhreeqcRM::GetComponents GetComponents_DOCSTRING
*/
const std::vector< std::string > &          GetComponents(void) const {return this->components;}

/**
%define GetConcentrations_DOCSTRING
'Transfer solution concentrations from each reaction cell
to the concentration vector given in the argument list (@a c).
Units of concentration for @a c are defined by :meth: `SetUnitsSolution`.
For per liter concentration units,
solution volume is used to calculate the concentrations for @a c.
For mass-fraction concentration units, the solution mass is used to calculate concentrations for @a c.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of saturation (:meth: `SetSaturation`),
porosity (:meth: `SetPorosity`), and representative volume (:meth: `SetRepresentativeVolume`),
and the mass of solution is volume times density as defined by :meth: `SetDensity`.
:meth: `UseSolutionDensityVolume` determines which option is used.
For option 1, the databases that have partial molar volume definitions needed
to accurately calculate solution volume are
phreeqc.dat, Amm.dat, and pitzer.dat.

Args:
	c(DoubleVector): Vector to receive the concentrations.
		Dimension of the vector is set to @a ncomps times @a nxyz,
		where,  ncomps is the result of :meth: `FindComponents` or :meth: `GetComponentCount`,
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetConcentrations GetConcentrations_DOCSTRING
*/
IRM_RESULT                                GetConcentrations(std::vector< double > &c);
/**
%define GetCurrentSelectedOutputUserNumber_DOCSTRING
'Returns the user number of the current selected-output definition.
:meth: `SetCurrentSelectedOutputUserNumber` or :meth: `SetNthSelectedOutput` specifies which of the
selected-output definitions is used.
Returns:
	int: User number of the the current selected-output definition,
		negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetCurrentSelectedOutputUserNumber GetCurrentSelectedOutputUserNumber_DOCSTRING
*/
int                                       GetCurrentSelectedOutputUserNumber(void);
/**
%define GetDatabaseFileName_DOCSTRING
'Returns the file name of the database. Should be called after :meth: `LoadDatabase`.
Returns:
	string: The file name defined in :meth: `LoadDatabase`.'
%enddef
%feature PhreeqcRM::GetDatabaseFileName GetDatabaseFileName_DOCSTRING
*/
std::string                               GetDatabaseFileName(void) {return this->database_file_name;}
/**
%define GetDensity_DOCSTRING
'Transfer solution densities from the reaction-module workers to the vector given 
in the argument list (@a density). This method always returns the calculated
densities; :meth: `SetDensity` does not affect the result.
Args:
	density(DoubleVector): Vector to receive the densities. Dimension of the array is set to @a nxyz,
		where @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		Values for inactive cells are set to 1e30.
		Densities are those calculated by the reaction module.
		Only the following databases distributed with PhreeqcRM have molar volume information needed
		to accurately calculate density: phreeqc.dat, Amm.dat, and pitzer.dat.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetDensity GetDensity_DOCSTRING
*/
IRM_RESULT                                GetDensity(std::vector< double > & density);
/**
%define GetEndCell_DOCSTRING
'Returns a vector of integers that contains the largest reaction-cell number assigned to each worker.
Each worker is assigned a range of reaction-cell numbers that are run during a call to :meth: `RunCells`.
The range of reaction cells for a worker may vary as load rebalancing occurs.
At any point in the calculations, the first cell and last cell to be run by a worker can be found
in the vectors returned by :meth: `GetStartCell` and :meth: `GetEndCell`.
Each method returns a vector of integers that has length of the number of threads (:meth: `GetThreadCount`),
if using OPENMP, or the number of processes (:meth: `GetMpiTasks`), if using MPI.

Returns:
	IRM_RESULT: Vector of integers, one for each worker, that gives the last reaction cell
		to be run by each worker.'
%enddef
%feature PhreeqcRM::GetEndCell GetEndCell_DOCSTRING
*/
const std::vector < int> &                GetEndCell(void) const {return this->end_cell;}
/**
%define GetEquilibriumPhases_DOCSTRING
'Returns a reference to the vector of all equilibrium phases.
The list includes all phases included in any EQUILIBRIUM_PHASES definitions in
the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetEquilibriumPhases`.
This method may be useful when generating selected output definitions related to equilibrium phases.

Returns:
	const tuple of strings: A vector of strings; each string is a unique
		equilibrium phases name.'
%enddef
%feature PhreeqcRM::GetEquilibriumPhases GetEquilibriumPhases_DOCSTRING
*/
const std::vector< std::string > &          GetEquilibriumPhases(void) const { return this->EquilibriumPhasesList; }
/**
%define GetEquilibriumPhasesCount_DOCSTRING
'Returns the number of equilibrium phases in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetEquilibriumPhasesCount`.
This method may be useful when generating selected output definitions related to
equilibrium phases.

Returns:
	int: The number of equilibrium phases in the initial-phreeqc module.'
%enddef
%feature PhreeqcRM::GetEquilibriumPhasesCount GetEquilibriumPhasesCount_DOCSTRING
*/
int                                       GetEquilibriumPhasesCount(void) const { return (int) this->EquilibriumPhasesList.size(); }


/**
%define GetErrorHandlerMode_DOCSTRING
'Get the setting for the action to be taken when the reaction module encounters an error.
Options are 0, return to calling program with an error return code (default);
1, throw an exception, which can be caught in C++ (for C and Fortran, the program will exit);
2, attempt to exit gracefully.
Returns:
	int: Current setting for the error handling mode: 0, 1, or 2.'
%enddef
%feature PhreeqcRM::GetErrorHandlerMode GetErrorHandlerMode_DOCSTRING
*/
int                                       GetErrorHandlerMode(void) {return this->error_handler_mode;}
/**
%define GetErrorString_DOCSTRING
'Returns a standard string containing error messages related to the last call to a PhreeqcRM method.
Returns:
	string: Error messages related to the last call to a PhreeqcRM method.'
%enddef
%feature PhreeqcRM::GetErrorString GetErrorString_DOCSTRING
*/
std::string                                  GetErrorString(void);

/**
%define GetExchangeNames_DOCSTRING
'Returns a reference to the vector of exchange names (such as "X") that correspond with
the exchange species names.
:meth: `FindComponents` must be called before :meth: `GetExchangeNames`.
The exchange names vector is the same length as the exchange species names vector
and provides the corresponding exchange site.
This method may be useful when generating selected output definitions related to exchangers.

Returns:
	const tuple of strings: A vector of strings; each string is an
		exchange name corresponding to the exchange species vector; an exchange name may occur
		multiple times.'
%enddef
%feature PhreeqcRM::GetExchangeNames GetExchangeNames_DOCSTRING
*/
const std::vector< std::string > &          GetExchangeNames(void) const { return this->ExchangeNamesList; }
/**
%define GetExchangeSpecies_DOCSTRING
'Returns a reference to the vector of exchange species names (such as "NaX").
The list of exchange species (such as "NaX") is derived from the list of components
(:meth: `FindComponents`) and the list of all exchange names (such as "X")
that are included in EXCHANGE definitions in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetExchangeSpecies`.
This method may be useful when generating selected output definitions related to exchangers.

Returns:
	const tuple of strings: Vector of strings; each string is a
		unique exchange species name.'
%enddef
%feature PhreeqcRM::GetExchangeSpecies GetExchangeSpecies_DOCSTRING
*/
const std::vector< std::string > &          GetExchangeSpecies(void) const { return this->ExchangeSpeciesNamesList; }
/**
%define GetExchangeSpeciesCount_DOCSTRING
'Returns the number of exchange species in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetExchangeSpeciesCount`.
This method may be useful when generating selected output definitions related to exchangers.

Returns:
	int: The number of exchange species in the initial-phreeqc module.'
%enddef
%feature PhreeqcRM::GetExchangeSpeciesCount GetExchangeSpeciesCount_DOCSTRING
*/
int                                       GetExchangeSpeciesCount(void) const { return (int) this->ExchangeSpeciesNamesList.size(); }


/**
%define GetFilePrefix_DOCSTRING
'Returns the file prefix for the output (.chem.txt) and log files (.log.txt).
Returns:
	string: The file prefix as set by :meth: `SetFilePrefix`, or "myrun", by default.'
%enddef
%feature PhreeqcRM::GetFilePrefix GetFilePrefix_DOCSTRING
*/
std::string                               GetFilePrefix(void) {return this->file_prefix;}
/**
%define GetForwardMapping_DOCSTRING
'Returns a reference to a vector of ints that is a mapping from grid cells to
reaction cells.
The mapping is used to eliminate cells that are inactive and cells that are unnecessary because of symmetry
from the list of cells for which reactions must be run. The mapping may be many-to-one to account for symmetry.
The mapping is set by :meth: `CreateMapping`, or, by default, is a one-to-one mapping--all grid cells are
reaction cells (vector contains 0,1,2,3,...,@a nxyz-1).
Returns:
	const IntVector: A vector of integers of size @a nxyz (number of grid cells, :meth: `GetGridCellCount`).
		Nonnegative is a reaction-cell number (0 based), negative is an inactive cell.'
%enddef
%feature PhreeqcRM::GetForwardMapping GetForwardMapping_DOCSTRING
*/
const std::vector < int > &               GetForwardMapping(void) {return this->forward_mapping_root;}

/**
%define GetGasComponents_DOCSTRING
'Returns a reference to the vector of all gas components in the initial-phreeqc module.
The list includes all gas components included in any GAS_PHASE definitions in
the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetGasComponents`.
This method may be useful when generating selected output definitions related to gas phases.

Returns:
	const tuple of strings: Vector of strings; each string is a unique
		gas component name.'
%enddef
%feature PhreeqcRM::GetGasComponents GetGasComponents_DOCSTRING
*/
const std::vector< std::string > &          GetGasComponents(void) const { return this->GasComponentsList; }
/**
%define GetGasComponentsCount_DOCSTRING
'Returns the number of gas phase components in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetGasComponentsCount`.
This method may be useful when generating selected output definitions related to
gas phases.

Returns:
	int: The number of gas phase components in the initial-phreeqc module.
'
%enddef
%feature PhreeqcRM::GetGasComponentsCount GetGasComponentsCount_DOCSTRING
*/
int                                       GetGasComponentsCount(void) const { return (int) this->GasComponentsList.size(); }

/**
%define GetGasCompMoles_DOCSTRING
'Transfer moles of gas components from each reaction cell
to the vector given in the argument list (@a gas_moles).

Args:
	gas_moles(DoubleVector): Vector to receive the moles of gas components.
		Dimension of the vector is set to @a ngas_comps times @a nxyz,
		where, @a ngas_comps is the result of :meth: `GetGasComponentsCount`,
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If a gas component is not defined for a cell, the number of moles is set to -1.
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetGasCompMoles GetGasCompMoles_DOCSTRING
*/
IRM_RESULT                                GetGasCompMoles(std::vector< double >& gas_moles);

/**
%define GetGasCompPressures_DOCSTRING
'Transfer pressures of gas components from each reaction cell
to the vector given in the argument list (@a gas_pressure).

Args:
	gas_pressure(DoubleVector): Vector to receive the pressures of gas components.
		Dimension of the vector is set to @a ngas_comps times @a nxyz,
		where, @a ngas_comps is the result of :meth: `GetGasComponentsCount`,
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If a gas component is not defined for a cell, the pressure is set to -1.
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetGasCompPressures GetGasCompPressures_DOCSTRING
*/
IRM_RESULT                                GetGasCompPressures(std::vector< double >& gas_pressure);

/**
%define GetGasCompPhi_DOCSTRING
'Transfer fugacity coefficients (phi) of gas components from each reaction cell
to the vector given in the argument list (@a gas_phi). Fugacity is
equal to the gas component pressure times the fugacity coefficient.

Args:
	gas_phi(DoubleVector): Vector to receive the fugacity coefficients of gas components.
		Dimension of the vector is set to @a ngas_comps times @a nxyz,
		where, @a ngas_comps is the result of :meth: `GetGasComponentsCount`,
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If a gas component is not defined for a cell, the fugacity coefficient is set to -1.
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetGasCompPhi GetGasCompPhi_DOCSTRING
*/
IRM_RESULT                                GetGasCompPhi(std::vector< double >& gas_phi);

/**
%define GetGasPhaseVolume_DOCSTRING
'Transfer volume of gas phase from each reaction cell
to the vector given in the argument list (@a gas_volume). 

Args:
	gas_volume(DoubleVector): Vector to receive the gas phase volumes.
		Dimension of the vector is set to @a nxyz,
		where,  @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If a gas phase is not defined for a cell, the volume is set to -1.
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetGasPhaseVolume GetGasPhaseVolume_DOCSTRING
*/
IRM_RESULT                                GetGasPhaseVolume(std::vector< double >& gas_volume);

/**
%define GetGfw_DOCSTRING
'Returns a reference to a vector of doubles that contains the gram-formula weight of
each component. Called after :meth: `FindComponents`. Order of weights corresponds to the list of components from
:meth: `GetComponents`.
Returns:
	const DoubleVector: A vector of doubles; each value is a component gram-formula weight, g/mol.'
%enddef
%feature PhreeqcRM::GetGfw GetGfw_DOCSTRING
*/
const std::vector < double > &            GetGfw(void) {return this->gfw;}
/**
%define GetGridCellCount_DOCSTRING
'Returns the number of grid cells in the user's model, which is defined in
the call to the constructor for the reaction module.
The mapping from grid cells to reaction cells is defined by :meth: `CreateMapping`.
The number of reaction cells may be less than the number of grid cells if
there are inactive regions or there is symmetry in the model definition.
Returns:
	int: Number of grid cells in the user's model.'
%enddef
%feature PhreeqcRM::GetGridCellCount GetGridCellCount_DOCSTRING
*/
int                                       GetGridCellCount(void) {return this->nxyz;}
/**
%define GetIPhreeqcPointer_DOCSTRING
'Returns an IPhreeqc pointer to the @a ith IPhreeqc instance in the reaction module.
For the threaded version, there are @a nthreads + 2 IPhreeqc instances, where
@a nthreads is defined in the constructor (:meth: `PhreeqcRM`::PhreeqcRM).
The number of threads can be determined by :meth: `GetThreadCount`.
The first @a nthreads (0 based) instances will be the workers, the
next (@a nthreads) is the InitialPhreeqc instance, and the next (@a nthreads + 1) is the Utility instance.
Getting the IPhreeqc pointer for one of these instances allows the user to use any of the IPhreeqc methods
on that instance.
For MPI, each process has exactly three IPhreeqc instances, one worker (number 0),
one InitialPhreeqc instance (number 1), and one Utility instance (number 2).
Args:
	i(int): The number of the IPhreeqc instance (0 based) to be retrieved.
Returns:
	IPhreeqc pointer to the @a ith IPhreeqc instance (0 based) in the reaction module.'
%enddef
%feature PhreeqcRM::GetIPhreeqcPointer GetIPhreeqcPointer_DOCSTRING
*/
IPhreeqc *                                GetIPhreeqcPointer(int i);

/**
%define GetKineticReactions_DOCSTRING
'Returns a reference to the vector of all kinetic reactions in the initial-phreeqc module.
The list includes all kinetic reactions included in any KINETICS definitions in
the reaction model.
:meth: `FindComponents` must be called before :meth: `GetKineticReactions`.
This method may be useful when generating selected output definitions related to kinetic reactions.
Returns:
	const tuple of strings: Vector of strings; each string is a unique
		kinetic reaction name.
'
%enddef
%feature PhreeqcRM::GetKineticReactions GetKineticReactions_DOCSTRING
*/
const std::vector< std::string > &          GetKineticReactions(void) const { return this->KineticReactionsList; }
/**
%define GetKineticReactionsCount_DOCSTRING
'Returns the number of kinetic reactions in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetKineticReactionsCount`.
This method may be useful when generating selected output definitions related to
kinetic reactions.

Returns:
	int: The number of kinetic reactions in the initial-phreeqc module.
'
%enddef
%feature PhreeqcRM::GetKineticReactionsCount GetKineticReactionsCount_DOCSTRING
*/
int                                       GetKineticReactionsCount(void) const { return (int) this->KineticReactionsList.size(); }

/**
%define GetMpiMyself_DOCSTRING
'Returns the MPI process (task) number. For the MPI version,
the root task number is zero, and all MPI tasks have unique task numbers greater than zero.
The number of tasks can be obtained with :meth: `GetMpiTasks`. The number of
tasks and computer hosts is determined at run time by the mpiexec command, and the
number of reaction-module processes is defined by the communicator used in
constructing the reaction modules (:meth: `PhreeqcRM`::PhreeqcRM).
For the OPENMP version, the task number is always
zero, and the result of :meth: `GetMpiTasks` is one.
Returns:
	int: The MPI task number for a process.'
%enddef
%feature PhreeqcRM::GetMpiMyself GetMpiMyself_DOCSTRING
*/
int                                 GetMpiMyself(void) const {return this->mpi_myself;}
/**
%define GetMpiTasks_DOCSTRING
'Returns the number of MPI processes (tasks) assigned to the reaction module.
For the MPI version, the number of
tasks and computer hosts is specified at run time by the mpiexec command. The number of MPI processes
used for reaction calculations is determined by the MPI communicator
used in constructing the reaction modules. The communicator may define a subset of the
total number of MPI processes.
The root task number is zero, and all other MPI tasks have unique task numbers greater than zero.
For the OPENMP version, the number of tasks is
one, and the task number returned by :meth: `GetMpiMyself` is zero.
Returns:
	int: The number of MPI processes assigned to the reaction module.'
%enddef
%feature PhreeqcRM::GetMpiTasks GetMpiTasks_DOCSTRING
*/
int                                 GetMpiTasks(void) const {return this->mpi_tasks;}
/**
%define GetNthSelectedOutputUserNumber_DOCSTRING
'Returns the user number for the @a nth selected-output definition.
Definitions are sorted by user number. Phreeqc allows multiple selected-output
definitions, each of which is assigned a nonnegative integer identifier by the
user. The number of definitions can be obtained by :meth: `GetSelectedOutputCount`.
To cycle through all of the definitions, GetNthSelectedOutputUserNumber
can be used to identify the user number for each selected-output definition
in sequence. :meth: `SetCurrentSelectedOutputUserNumber` is then used to select
that user number for selected-output processing.
Args:
	n(int): The sequence number of the selected-output definition for which the user number will be returned.
		Fortran, 1 based; C, 0 based.
Returns:
	int: The user number of the @a nth selected-output definition, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetNthSelectedOutputUserNumber GetNthSelectedOutputUserNumber_DOCSTRING
*/
int                                       GetNthSelectedOutputUserNumber(int n);
/**
%define GetPartitionUZSolids_DOCSTRING
'Returns the setting for partitioning solids between the saturated and unsaturated
parts of a partially saturated cell. The option is intended to be used by saturated-only
flow codes that allow a variable water table.
The value has meaning only when saturations
less than 1.0 are encountered. The partially saturated cells
may have a small water-to-rock ratio that causes
reactions to proceed slowly relative to fully saturated cells.
By setting  :meth: `SetPartitionUZSolids` to true, the
amounts of solids and gases are partioned according to the saturation.
If a cell has a saturation of 0.5, then
the water interacts with only half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell varies,
solids and gases are transferred between the
saturated and unsaturated (unreactive) reservoirs of the cell.
Unsaturated-zone flow and transport codes will probably use the default (false),
which assumes all gases and solids are reactive regardless of saturation.
Returns:
	Boolean: @a True, the fraction of solids and gases available for
		reaction is equal to the saturation;
		@a False (default), all solids and gases are reactive regardless of saturation.'
%enddef
%feature PhreeqcRM::GetPartitionUZSolids GetPartitionUZSolids_DOCSTRING
*/
bool                                GetPartitionUZSolids(void) const {return this->partition_uz_solids;}
/**
%define GetPoreVolume_DOCSTRING
'Returns the current set of pore volumes as
defined by the last use of :meth: `SetPoreVolume` or the default (0.1 L).
Pore volume is used with cell volume (:meth: `SetCellVolume`) in calculating porosity.
Pore volumes may change as a function of pressure, in which case they can be updated
with :meth: `SetPoreVolume`.
Returns:
	const DoubleVector: A vector reference to the pore volumes.
		Size of vector is @a nxyz, the number of grid cells in the user's model (:meth: `GetGridCellCount`).'
%enddef
%feature PhreeqcRM::GetPoreVolume GetPoreVolume_DOCSTRING
*/
std::vector< double > &                     GetPoreVolume(void) {return this->pore_volume;}
#endif
/**
%define GetPorosity_DOCSTRING
'Returns the porosity for each cell.
By default, the porosity vector is initialized with 0.1, unitless.
PhreeqcRM does not change the porosity, so the values that are retrieved are
either the default, or the values set by the last call to :meth: `SetPorosity`. 
Returns:
	const DoubleVector: A vector reference to the porosities in each cell, unitless.
		Size of vector is @a nxyz, the number of grid cells in the user's model (:meth: `GetGridCellCount`).'
%enddef
%feature PhreeqcRM::GetPorosity GetPorosity_DOCSTRING
*/
const std::vector< double >& GetPorosity(void);
/**
%define GetPressure_DOCSTRING
'Returns the pressure for each cell.
By default, the pressure vector is initialized with 1 atm;
if :meth: `SetPressure` has not been called, worker solutions will have pressures as defined in
input files (:meth: `RunFile`) or input strings (:meth: `RunString`); if :meth: `SetPressure` has been called,
worker solutions will have the pressures as defined by :meth: `SetPressure`.
Pressure effects are considered by three PHREEQC databases: phreeqc.dat, Amm.dat, and pitzer.dat.
Returns:
	const DoubleVector: A vector reference to the pressures in each cell, in atm.
		Size of vector is @a nxyz, the number of grid cells in the user's model (:meth: `GetGridCellCount`).'
%enddef
%feature PhreeqcRM::GetPressure GetPressure_DOCSTRING
*/
const std::vector< double > &                     GetPressure(void);
/**
%define GetPrintChemistryMask _DOCSTRING
'Return a reference to the vector of print flags that enable or disable detailed output for each cell.
Printing for a cell will occur only when the
printing is enabled with :meth: `SetPrintChemistryOn`, and the value in the vector for the cell is 1.
Returns:
	IntVector: Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`). A value of 0 for a cell indicates
		printing is disabled;
		a value of 1 for a cell indicates printing is enabled.'
%enddef
%feature PhreeqcRM::GetPrintChemistryMask  GetPrintChemistryMask _DOCSTRING
*/
const std::vector< int > &                  GetPrintChemistryMask (void) {return this->print_chem_mask_root;}
/**
%define GetPrintChemistryOn_DOCSTRING
'Returns a vector reference to the current print flags for detailed output for the three sets of IPhreeqc instances:
the workers, the InitialPhreeqc instance, and the Utility instance. Dimension of the vector is 3.
Printing of detailed output from reaction calculations to the output file
is enabled when the vector value is true, disabled when false.
The detailed output prints all of the output
typical of a PHREEQC reaction calculation, which includes solution descriptions and the compositions of
all other reactants. The output can be several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, :meth: `OpenFiles`). For the worker instances,
the output can be limited to a set of cells
(:meth: `SetPrintChemistryMask`) and, in general, the
amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
(applied by using :meth: `RunFile` or :meth: `RunString`).
Printing the detailed output for the workers is generally used only for debugging,
and PhreeqcRM will run faster when printing detailed output for the workers is disabled (:meth: `SetPrintChemistryOn`).
Returns:
	const BoolVector: Print flag for the workers, InitialPhreeqc, and Utility IPhreeqc instances, respectively.'
%enddef
%feature PhreeqcRM::GetPrintChemistryOn GetPrintChemistryOn_DOCSTRING
*/
const std::vector <bool> &                GetPrintChemistryOn(void) const {return this->print_chemistry_on;}
/**
%define GetRebalanceByCell_DOCSTRING
'Get the load-rebalancing method used for parallel processing.
PhreeqcRM attempts to rebalance the load of each thread or
process such that each
thread or process takes the same amount of time to run its part of a :meth: `RunCells`
calculation. Two algorithms are available: one accounts for cells that were not run
because saturation was zero (true), and the other uses the average time
to run all of the cells assigned to a process or thread (false), .
The methods are similar, but preliminary results indicate the default is better in most cases.
Returns:
	Boolean: @a True indicates individual
		cell run times are used in rebalancing (default); @a False, indicates average run times are used in rebalancing.'
%enddef
%feature PhreeqcRM::GetRebalanceByCell GetRebalanceByCell_DOCSTRING
*/
bool                                      GetRebalanceByCell(void) const {return this->rebalance_by_cell;}
/**
%define GetRebalanceFraction_DOCSTRING
'Get the fraction used to determine the number of cells to transfer among threads or processes.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a :meth: `RunCells`
calculation. The rebalancing transfers cell calculations among threads or processes to
try to achieve an optimum balance. :meth: `SetRebalanceFraction`
adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
determine the number of cell transfers that actually are made. A value of zero eliminates
load rebalancing. A value less than 1.0 is suggested to avoid possible oscillations,
where too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
Default is 0.5.
Returns:
	float: Fraction used in rebalance, 0.0 to 1.0.'
%enddef
%feature PhreeqcRM::GetRebalanceFraction GetRebalanceFraction_DOCSTRING
*/
double                                    GetRebalanceFraction(void) const {return this->rebalance_fraction;}
/**
%define GetSaturation_DOCSTRING
'Returns a vector of saturations (@a sat) as calculated by the reaction module. 
This method always returns solution_volume/(rv * porosity); the method 
:meth: `SetSaturation` has no effect on the values returned.
Reactions will change the volume of solution in a cell.
The transport code must decide whether to ignore or account for this change in solution volume due to reactions.
Following reactions, the cell saturation is calculated as solution volume (:meth: `GetSolutionVolume`)
divided by the product of representative volume (:meth: `SetRepresentativeVolume`) and the porosity (:meth: `SetPorosity`).
The cell saturation returned by @a GetSaturation may be less than or greater than the saturation set by the transport code
(:meth: `SetSaturation`), and may be greater than or less than 1.0, even in fully saturated simulations.
Only the following databases distributed with PhreeqcRM have molar volume information needed
to accurately calculate solution volume and saturation: phreeqc.dat, Amm.dat, and pitzer.dat.

Args:
	sat(DoubleVector): Vector to receive the saturations. Dimension of the array is set to @a nxyz,
		where @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		Values for inactive cells are set to 1e30.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSaturation GetSaturation_DOCSTRING
*/
IRM_RESULT               GetSaturation(std::vector< double > & sat);
/**
%define GetSelectedOutput_DOCSTRING
'Returns the array of selected-output values for the current selected-output definition.
:meth: `SetCurrentSelectedOutputUserNumber`
specifies which of the selected-output definitions is returned to the vector (@a so).
Args:
	so(DoubleVector): A vector to contain the selected-output values.
		Size of the vector is set to @a col times @a nxyz, where @a col is the number of
		columns in the selected-output definition (:meth: `GetSelectedOutputColumnCount`),
		and @a nxyz is the number of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSelectedOutput GetSelectedOutput_DOCSTRING
*/
IRM_RESULT                                GetSelectedOutput(std::vector< double > &so);
/**
%define GetSelectedOutputColumnCount_DOCSTRING
'Returns the number of columns in the current selected-output definition.
:meth: `SetCurrentSelectedOutputUserNumber` specifies which of the selected-output definitions is used.
Returns:
	int: Number of columns in the current selected-output definition,
		negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSelectedOutputColumnCount GetSelectedOutputColumnCount_DOCSTRING
*/
int                                       GetSelectedOutputColumnCount(void);
/**
%define GetSelectedOutputCount_DOCSTRING
'Returns the number of selected-output definitions.
:meth: `SetCurrentSelectedOutputUserNumber` specifies which of the selected-output definitions is used.
Returns:
	int: Number of selected-output definitions, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSelectedOutputCount GetSelectedOutputCount_DOCSTRING
*/
int                                       GetSelectedOutputCount(void);
/**
%define GetSelectedOutputHeading_DOCSTRING
'Returns a selected-output heading.
The number of headings is determined by :meth: `GetSelectedOutputColumnCount`.
:meth: `SetCurrentSelectedOutputUserNumber` specifies which of the selected-output definitions is used.
Args:
	icol(int): The sequence number of the heading to be retrieved, 0 based.
	heading(string): A string to receive the heading.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSelectedOutputHeading GetSelectedOutputHeading_DOCSTRING
*/
IRM_RESULT                  GetSelectedOutputHeading(int icol, std::string &heading);
/**
%define GetSelectedOutputHeadings_DOCSTRING
'Returns a list of the current selected-output headings.
The number of headings is determined by :meth: `GetSelectedOutputColumnCount`.
:meth: `SetCurrentSelectedOutputUserNumber` or :meth: `BMI`_SetValue("NthSelectedOutput",i) are
used to specify which of the selected-output definitions is used.
Args:
	headings(tuple of strings): A vector of std::strings to receive the headings.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSelectedOutputHeadings GetSelectedOutputHeadings_DOCSTRING
*/
IRM_RESULT                  GetSelectedOutputHeadings(std::vector< std::string >& headings);
/**
%define GetSelectedOutputOn_DOCSTRING
'Returns the current value of the selected-output property.
A value of true for this property indicates that selected output data will be requested this time step.
A value of false indicates that selected output will not be retrieved for this time step;
processing the selected output is avoided with some time savings.
Returns:
	Boolean: @a True, selected output will be requested; @a false, selected output will not be retrieved.'
%enddef
%feature PhreeqcRM::GetSelectedOutputOn GetSelectedOutputOn_DOCSTRING
*/
bool                                      GetSelectedOutputOn(void) {return this->selected_output_on;}
/**
%define GetSelectedOutputRowCount_DOCSTRING
'Returns the number of rows in the current selected-output definition. However, the method
is included only for convenience; the number of rows is always equal to the number of
grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	int: Number of rows in the current selected-output definition, negative is failure.'
%enddef
%feature PhreeqcRM::GetSelectedOutputRowCount GetSelectedOutputRowCount_DOCSTRING
*/
int                                       GetSelectedOutputRowCount(void);

/**
%define GetSICount_DOCSTRING
'Returns the number of phases in the initial-phreeqc module for which saturation indices could be calculated.
:meth: `FindComponents` must be called before :meth: `GetSICount`.
This method may be useful when generating selected output definitions related to
saturation indices.

Returns:
	int: The number of phases in the initial-phreeqc module for which saturation indices
		could be calculated.'
%enddef
%feature PhreeqcRM::GetSICount GetSICount_DOCSTRING
*/
int                                       GetSICount(void) const { return (int) this->SINamesList.size(); }
/**
%define GetSINames_DOCSTRING
'Returns a reference to the vector of the names of all phases for which
saturation indices (SIs) could be calculated.
The list includes all phases that contain only elements included in the components in
the initial-phreeqc module.
The list assumes that all components are present to be able to calculate the entire list of SIs;
it may be that one or more components are missing in any specific cell.
:meth: `FindComponents` must be called before :meth: `GetSINames`.
This method may be useful when generating selected output definitions related to
saturation indices.

Returns:
	const tuple of strings: Vector of strings; each string is a unique
		phase name.'
%enddef
%feature PhreeqcRM::GetSINames GetSINames_DOCSTRING
*/
const std::vector< std::string > &          GetSINames(void) const { return this->SINamesList; }

/**
%define GetSolidSolutionComponents_DOCSTRING
'Returns a reference to the vector of solid solution components.
The list of solid solution components includes all components in any SOLID_SOLUTION
definitions in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetSolidSolutionComponents`.
This method may be useful when generating selected output definitions related to
solid solutions.

Returns:
	const tuple of strings: Vector of strings; each string is a
		unique solid solution component.'
%enddef
%feature PhreeqcRM::GetSolidSolutionComponents GetSolidSolutionComponents_DOCSTRING
*/
const std::vector< std::string > &          GetSolidSolutionComponents(void) const { return this->SolidSolutionComponentsList; }
/**
%define GetSolidSolutionComponentsCount_DOCSTRING
'Returns the number of solid solution components in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetSolidSolutionComponentsCount`.
This method may be useful when generating selected output definitions related to solid solutions.

Returns:
	int: The number of solid solution components in the initial-phreeqc module.
'
%enddef
%feature PhreeqcRM::GetSolidSolutionComponentsCount GetSolidSolutionComponentsCount_DOCSTRING
*/
int                                       GetSolidSolutionComponentsCount(void) const { return (int) this->SolidSolutionComponentsList.size(); }

/**
%define GetSolidSolutionNames_DOCSTRING
'Returns a reference to the vector of solid solution names that correspond with
the solid solution components.
:meth: `FindComponents` must be called before :meth: `GetSolidSolutionNames`.
The solid solution names vector is the same length as the solid solution components vector
and provides the corresponding name of solid solution containing the component.
This method may be useful when generating selected output definitions related to solid solutions.

Returns:
	const tuple of strings: Vector of strings; each string is a
		solid solution name corresponding to the solid solution components vector; a solid solution name may occur
		multiple times.'
%enddef
%feature PhreeqcRM::GetSolidSolutionNames GetSolidSolutionNames_DOCSTRING
*/
const std::vector< std::string > &          GetSolidSolutionNames(void) const { return this->SolidSolutionNamesList; }

/**
%define GetSolutionVolume_DOCSTRING
'Return a vector reference to the current solution volumes as calculated by the reaction module.
Dimension of the vector will be @a nxyz, where @a nxyz is the number of user grid cells.
Values for inactive cells are set to 1e30.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
Returns:
	DoubleVector: Vector reference to current solution volumes.'
%enddef
%feature PhreeqcRM::GetSolutionVolume GetSolutionVolume_DOCSTRING
*/
const std::vector< double > &               GetSolutionVolume(void);
/**
%define GetSpeciesConcentrations_DOCSTRING
'Returns a vector reference to aqueous species concentrations (@a species_conc).
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.
Solution volumes used to calculate mol/L are calculated by the reaction module.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.

Args:
	species_conc(DoubleVector): Vector to receive the aqueous species concentrations.
		Dimension of the vector is set to @a nspecies times @a nxyz,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`),
		and @a nxyz is the number of grid cells (:meth: `GetGridCellCount`).
		Concentrations are moles per liter.
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSpeciesConcentrations GetSpeciesConcentrations_DOCSTRING
*/
IRM_RESULT                                GetSpeciesConcentrations(std::vector< double > & species_conc);
/**
%define GetSpeciesCount_DOCSTRING
'Returns the number of aqueous species used in the reaction module.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.
Returns:
	int: The number of aqueous species.'
%enddef
%feature PhreeqcRM::GetSpeciesCount GetSpeciesCount_DOCSTRING
*/
int                                       GetSpeciesCount(void) {return (int) this->species_names.size();}
/**
%define GetSpeciesD25_DOCSTRING
'Returns a vector reference to diffusion coefficients at 25C for the set of aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally in the database file.
Databases distributed with the reaction module that have diffusion coefficients defined are
phreeqc.dat, Amm.dat, and pitzer.dat.
Returns:
	DoubleVector: Vector reference to the diffusion coefficients at 25 C, m^2/s. Dimension of the vector is @a nspecies,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`).'
%enddef
%feature PhreeqcRM::GetSpeciesD25 GetSpeciesD25_DOCSTRING
*/
const std::vector< double > &               GetSpeciesD25(void) {return this->species_d_25;}
/**
%define GetSpeciesLog10Gammas_DOCSTRING
'Returns a vector reference to log10 aqueous species activity coefficients (@a species_log10gammas).
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.

Args:
	species_log10gammas(DoubleVector): Vector to receive the log10 aqueous species activity coefficients.
		Dimension of the vector is set to @a nspecies times @a nxyz,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`),
		and @a nxyz is the number of grid cells (:meth: `GetGridCellCount`).
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSpeciesLog10Gammas GetSpeciesLog10Gammas_DOCSTRING
*/
IRM_RESULT                                GetSpeciesLog10Gammas(std::vector< double > & species_log10gammas);	

/**
%define GetSpeciesLog10Molalities_DOCSTRING
'Returns a vector reference to log10 aqueous species molalities (@a species_log10molalities).
To use this method :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.

Args:
	species_log10molalities(DoubleVector): Vector to receive the log10 aqueous species molalites.
		Dimension of the vector is set to @a nspecies times @a nxyz,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`),
		and @a nxyz is the number of grid cells (:meth: `GetGridCellCount`).
		Values for inactive cells are set to 1e30.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::GetSpeciesLog10Molalities GetSpeciesLog10Molalities_DOCSTRING
*/
IRM_RESULT                                GetSpeciesLog10Molalities(std::vector< double >& species_log10molalities);	

/**
%define GetSpeciesNames_DOCSTRING
'Returns a vector reference to the names of the aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.
Returns:
	tuple of strings: Vector of strings containing the names of the aqueous species. Dimension of the vector is @a nspecies,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`).'
%enddef
%feature PhreeqcRM::GetSpeciesNames GetSpeciesNames_DOCSTRING
*/
const std::vector< std::string > &          GetSpeciesNames(void) {return this->species_names;}
/**
%define GetSpeciesSaveOn_DOCSTRING
'Returns the value of the species-save property.
By default, concentrations of aqueous species are not saved. Setting the species-save property to true allows
aqueous species concentrations to be retrieved
with :meth: `GetSpeciesConcentrations`, and solution compositions to be set with
:meth: `SpeciesConcentrations2Module`.

Returns:
	Boolean: True indicates solution species concentrations are saved and can be used for multicomponent-diffusion calculations;
		@a False indicates that solution species concentrations are not saved.'
%enddef
%feature PhreeqcRM::GetSpeciesSaveOn GetSpeciesSaveOn_DOCSTRING
*/
bool                                      GetSpeciesSaveOn(void) {return this->species_save_on;}

/**
%define GetSpeciesStoichiometry_DOCSTRING
'Returns a vector reference to the stoichiometry of each aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.

Returns:
	Vector of cxxNameDouble instances (maps) that contain the component names and
		associated stoichiometric coefficients for each aqueous species.  Dimension of the vector is @a nspecies,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`).'
%enddef
%feature PhreeqcRM::GetSpeciesStoichiometry GetSpeciesStoichiometry_DOCSTRING
*/

const std::vector<cxxNameDouble> &        GetSpeciesStoichiometry(void) {return this->species_stoichiometry;}
/**
%define GetSpeciesZ_DOCSTRING
'Returns a vector reference to the charge on each aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
Returns:
	DoubleVector: Vector containing the charge on each aqueous species. Dimension of the vector is @a nspecies,
		where @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`).'
%enddef
%feature PhreeqcRM::GetSpeciesZ GetSpeciesZ_DOCSTRING
*/
const std::vector< double > &               GetSpeciesZ(void) {return this->species_z;}
/**
%define GetStartCell_DOCSTRING
'Returns a vector of integers that contains the smallest reaction-cell number assigned to each worker.
Each worker is assigned a range of reaction-cell numbers that are run during a call to :meth: `RunCells`.
The range of reaction cell numbers for a worker may vary as load rebalancing occurs.
At any point in the calculations, the first cell and last cell to be run by a worker can be found
in the vectors returned by @a GetStartCell and :meth: `GetEndCell`.
Each method returns a vector of integers that has size of the number of threads (:meth: `GetThreadCount`),
if using OPENMP, or the number of processes (:meth: `GetMpiTasks`), if using MPI.
Returns:
	IntVector: Vector of integers, one for each worker, that gives the first reaction cell
		to be run by each worker.'
%enddef
%feature PhreeqcRM::GetStartCell GetStartCell_DOCSTRING
*/
const std::vector < int> &                GetStartCell(void) const {return this->start_cell;}

/**
%define GetSurfaceNames_DOCSTRING
'Returns a reference to the vector of surface names (such as "Hfo") that correspond with
the surface species names. The vectors referenced by :meth: `GetSurfaceSpecies`
and :meth: `GetSurfaceNames` are the same length.
:meth: `FindComponents` must be called before :meth: `GetSurfaceNames`.
This method may be useful when generating selected output definitions related to surfaces.

Returns:
	const tuple of strings: Vector of strings; each string is a
		surface name corresponding to the surface species vector;
		a surface name may occur multiple times.'
%enddef
%feature PhreeqcRM::GetSurfaceNames GetSurfaceNames_DOCSTRING
*/
const std::vector< std::string > &          GetSurfaceNames(void) const { return this->SurfaceNamesList; }

/**
%define GetSurfaceSpecies_DOCSTRING
'Returns a reference to the vector of surface species names (such as "Hfo_wOH").
The list of surface species is derived from the list of components
(:meth: `FindComponents`) and the list of all surface site types (such as "Hfo_w")
that are included in SURFACE definitions in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetSurfaceSpecies`.
This method may be useful when generating selected output definitions related to surfaces.

Returns:
	const std::vector< std::string >&       A vector of strings; each string is a
		unique surface species name.'
%enddef
%feature PhreeqcRM::GetSurfaceSpecies GetSurfaceSpecies_DOCSTRING
*/
const std::vector< std::string > &          GetSurfaceSpecies(void) const { return this->SurfaceSpeciesNamesList; }

/**
%define GetSurfaceSpeciesCount_DOCSTRING
'Returns the number of surface species (such as "Hfo_wOH") in the initial-phreeqc module.
:meth: `FindComponents` must be called before :meth: `GetSurfaceSpeciesCount`.
This method may be useful when generating selected output definitions related to surfaces.

Returns:
	int: The number of surface species in the initial-phreeqc module.
'
%enddef
%feature PhreeqcRM::GetSurfaceSpeciesCount GetSurfaceSpeciesCount_DOCSTRING
*/
int                                       GetSurfaceSpeciesCount(void) const { return (int) this->SurfaceSpeciesNamesList.size(); }


/**
%define GetSurfaceTypes_DOCSTRING
'Returns a reference to the vector of surface site types (such as "Hfo_w") that correspond with
the surface species names.
The vectors referenced by :meth: `GetSurfaceSpecies` and
:meth: `GetSurfaceTypes` are the same length.
:meth: `FindComponents` must be called before :meth: `GetSurfaceTypes`.
This method may be useful when generating selected output definitions related to surfaces.

Returns:
	const tuple of strings: Vector of strings; each string is a
		surface site type for the corresponding species in the surface species vector;
		a surface site type may occur multiple times.'
%enddef
%feature PhreeqcRM::GetSurfaceTypes GetSurfaceTypes_DOCSTRING
*/
const std::vector< std::string > &          GetSurfaceTypes(void) const { return this->SurfaceTypesList; }

/**
%define GetTemperature_DOCSTRING
'Vector reference to the current temperatures of the cells.
By default, the temperature vector is initialized to 25 C;
if :meth: `SetTemperature` has not been called, worker solutions will have temperatures as defined in
input files (:meth: `RunFile`) or input strings (:meth: `RunString`); if :meth: `SetTemperature` has been called,
worker solutions will have the temperatures as defined by :meth: `SetTemperature`.
Returns:
	DoubleVector: Vector of temperatures, in degrees C. Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`).'
%enddef
%feature PhreeqcRM::GetTemperature GetTemperature_DOCSTRING
*/
//const std::vector< double > &               GetTemperature(void) {return this->tempc;}
const std::vector< double > &               GetTemperature(void);
/**
%define GetThreadCount_DOCSTRING
'Returns the number of threads, which is equal to the number of workers used to run in parallel with OPENMP.
For the OPENMP version, the number of threads is set implicitly or explicitly
with the constructor (:meth: `PhreeqcRM`::PhreeqcRM).
For the MPI version, the number of threads is always one for each process.
Returns:
	int: The number of threads used for OPENMP parallel processing.'
%enddef
%feature PhreeqcRM::GetThreadCount GetThreadCount_DOCSTRING
*/
int                                       GetThreadCount() {return this->nthreads;}
/**
%define GetTime_DOCSTRING
'Returns the current simulation time in seconds.
The reaction module does not change the time value, so the
returned value is equal to the default (0.0) or the last time set by :meth: `SetTime`.
Returns:
	float: The current simulation time, in seconds.'
%enddef
%feature PhreeqcRM::GetTime GetTime_DOCSTRING
*/
double                                    GetTime(void) const {return this->time;}
/**
%define GetTimeConversion_DOCSTRING
'Returns a multiplier to convert time from seconds to another unit, as specified by the user.
The reaction module uses seconds as the time unit. The user can set a conversion
factor (:meth: `SetTimeConversion`) and retrieve it with GetTimeConversion.
The reaction module only uses the conversion factor when printing the long version
of cell chemistry (:meth: `SetPrintChemistryOn`), which is rare.
Default conversion factor is 1.0.
Returns:
	float: Multiplier to convert seconds to another time unit.'
%enddef
%feature PhreeqcRM::GetTimeConversion GetTimeConversion_DOCSTRING
*/
double                                    GetTimeConversion(void) {return this->time_conversion;}
/**
%define GetTimeStep_DOCSTRING
'Returns the current simulation time step in seconds.
This is the time over which kinetic reactions are integrated in a call to :meth: `RunCells`.
The reaction module does not change the time-step value, so the
returned value is equal to the default (0.0) or the last time step set by :meth: `SetTimeStep`.
Returns:
	float: The current simulation time step, in seconds.'
%enddef
%feature PhreeqcRM::GetTimeStep GetTimeStep_DOCSTRING
*/
double                                    GetTimeStep(void) {return this->time_step;}
/**
%define GetUnitsExchange_DOCSTRING
'Returns the input units for exchangers.
In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
:meth: `SetUnitsExchange` specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for exchangers.'
%enddef
%feature PhreeqcRM::GetUnitsExchange GetUnitsExchange_DOCSTRING
*/
int                                       GetUnitsExchange(void) {return this->units_Exchange;}
/**
%define GetUnitsGasPhase_DOCSTRING
'Returns the input units for gas phases.
In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
:meth: `SetUnitsGasPhase` specifies how the number of moles of component gases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for gas phases (0, 1, or 2).'
%enddef
%feature PhreeqcRM::GetUnitsGasPhase GetUnitsGasPhase_DOCSTRING
*/
int                                       GetUnitsGasPhase(void) {return this->units_GasPhase;}
/**
%define GetUnitsKinetics_DOCSTRING
'Returns the input units for kinetic reactants.
In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
:meth: `SetUnitsKinetics` specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for kinetic reactants (0, 1, or 2).'
%enddef
%feature PhreeqcRM::GetUnitsKinetics GetUnitsKinetics_DOCSTRING
*/
int                                       GetUnitsKinetics(void) {return this->units_Kinetics;}
/**
%define GetUnitsPPassemblage_DOCSTRING
'Returns the input units for pure phase assemblages (equilibrium phases).
In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
:meth: `SetUnitsPPassemblage` specifies how the number of moles of phases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for equilibrium phases (0, 1, or 2).'
%enddef
%feature PhreeqcRM::GetUnitsPPassemblage GetUnitsPPassemblage_DOCSTRING
*/
int                                       GetUnitsPPassemblage(void) {return this->units_PPassemblage;}
/**
%define GetUnitsSolution_DOCSTRING
'Returns the units of concentration used by the transport model.
Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
In PHREEQC, solutions are defined by the number of moles of each
element in the solution. The units of transport concentration are used when
transport concentrations are converted to
solution moles by :meth: `SetConcentrations` and :meth: `Concentrations2Utility`.
The units of solution concentration also are used when solution moles are converted to
transport concentrations by
:meth: `GetConcentrations`.
@n@n
To convert from mg/L to moles
of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
multiplied by the solution volume,
which is the product of porosity (:meth: `SetPorosity`), saturation (:meth: `SetSaturation`), and
representative volume (:meth: `SetRepresentativeVolume`).
To convert from mol/L to moles
of element in a cell, mol/L is
multiplied by the solution volume.
To convert from mass fraction to moles
of element in a cell, kg/kgs is converted to mol/kgs, multiplied by density
(:meth: `SetDensity`) and
multiplied by the solution volume.
@n@n
To convert from moles
of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
solution volume resulting in mol/L, and then converted to
mg/L.
To convert from moles
of element in the representative volume of a reaction cell to mol/L,  the number of moles of an element is divided by the
solution volume resulting in mol/L.
To convert from moles
of element in the representative volume of a reaction cell to mass fraction,
the number of moles of an element is converted to kg and divided by the total mass of the solution.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of porosity, saturation, and representative volume,
and the mass of solution is volume times density as defined by :meth: `SetDensity`.
Which option is used is determined by :meth: `UseSolutionDensityVolume`.
Returns:
	int: Units for concentrations in transport.'
%enddef
%feature PhreeqcRM::GetUnitsSolution GetUnitsSolution_DOCSTRING
*/
int                                       GetUnitsSolution(void) {return this->units_Solution;}
/**
%define GetUnitsSSassemblage_DOCSTRING
'Returns the input units for solid-solution assemblages.
In PHREEQC input, solid solutions are defined by moles of each component (@a Mp).
:meth: `SetUnitsSSassemblage` specifies how the number of moles in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for solid solutions (0, 1, or 2).'
%enddef
%feature PhreeqcRM::GetUnitsSSassemblage GetUnitsSSassemblage_DOCSTRING
*/
int                                       GetUnitsSSassemblage(void) {return this->units_SSassemblage;}
/**
%define GetUnitsSurface_DOCSTRING
'Returns the input units for surfaces.
In PHREEQC input, surfaces are defined by moles of surface sites  (@a Mp).
:meth: `SetUnitsSurface` specifies how the number of moles of surface sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

Returns:
	int: Input units for solid surfaces (0, 1, or 2).'
%enddef
%feature PhreeqcRM::GetUnitsSurface GetUnitsSurface_DOCSTRING
*/
int                                       GetUnitsSurface(void) {return this->units_Surface;}
/**
%define GetWorkers_DOCSTRING
'Returns a reference to the vector of IPhreeqcPhast instances. IPhreeqcPhast
inherits from IPhreeqc, and the vector can be interpreted as a vector of pointers to the worker,
InitialPhreeqc, and Utility IPhreeqc instances. For OPENMP, there are @a nthreads
workers, where @a nthreads is defined in the constructor (:meth: `PhreeqcRM`::PhreeqcRM).
For MPI, there is a single worker. For OPENMP and MPI, there is one InitialPhreeqc and
one Utility instance.
Returns:
	Vector of IPhreeqcPhast instances.'
%enddef
%feature PhreeqcRM::GetWorkers GetWorkers_DOCSTRING
*/
const std::vector<IPhreeqcPhast *> &      GetWorkers() {return this->workers;}
#ifdef USE_YAML
/**
%define InitializeYAML_DOCSTRING
'A YAML file can be used to initialize an instance of PhreeqcRM. 
Args:
	yamlfile(string): String containing the YAML file name.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).

The file contains a YAML map of PhreeqcRM methods
and the arguments corresponding to the methods. For example,
@htmlonly
<CODE>
<PRE>
LoadDatabase: phreeqc.dat
RunFile:
workers: true
initial_phreeqc: true
utility: true
chemistry_name: advect.pqi
</PRE>
</CODE>
@endhtmlonly

:meth: `InitializeYAML` will read the YAML file and execute the specified methods with 
the specified arguments. Using YAML
terminology, the argument(s) for a method may be a scalar, a sequence, or a map, 
depending if the argument is
a single item, a single vector, or there are multiple arguments. 
In the case of a map, the names associated
with each argument (for example "chemistry_name" above) is arbitrary. 
The names of the map keys for map
arguments are not used in parsing the YAML file; only the order of 
the arguments is important.

The class YAMLPhreeqcRM can be used to write a YAML file.
The methods defined in the YAMLPhreeqcRM class include the
following list; 
all but SetGridCellCount correspond to PhreeqcRM methods.
@htmlonly
<CODE>
<PRE>
CloseFiles(void);
CreateMapping(DoubleVector grid2chem);
DumpModule();
FindComponents();
InitialPhreeqc2Module(IntVector initial_conditions1);
InitialPhreeqc2Module(sIntVector initial_conditions1, IntVector initial_conditions2, DoubleVector fraction1);
InitialPhreeqcCell2Module(int n, IntVector cell_numbers);
LoadDatabase(string database);
OpenFiles(void);
OutputMessage(string str);
RunCells(void);
RunFile(Boolean workers, Boolean initial_phreeqc, Boolean utility, string chemistry_name);
RunString(Boolean workers, Boolean initial_phreeqc, Boolean utility, string input_string);
ScreenMessage(string str);
SetComponentH2O(Boolean tf);
SetConcentrations(DoubleVector c);
SetCurrentSelectedOutputUserNumber(int n_user);
SetDensity(DoubleVector density);
SetDumpFileName(string dump_name);
SetErrorHandlerMode(int mode);
SetErrorOn(Boolean tf);
SetFilePrefix(string prefix);
SetGasCompMoles(DoubleVector gas_moles);
SetGasPhaseVolume(DoubleVector gas_volume);
SetGridCellCount(int nxyz);
SetPartitionUZSolids(Boolean tf);
SetPorosity(DoubleVector por);
SetPressure(DoubleVector p);
SetPrintChemistryMask(IntVector cell_mask);
SetPrintChemistryOn(Boolean workers, Boolean initial_phreeqc, Boolean utility);
SetRebalanceByCell(Boolean tf);
SetRebalanceFraction(float f);
SetRepresentativeVolume(DoubleVector rv);
SetSaturation(DoubleVector sat);
SetScreenOn(Boolean tf);
SetSelectedOutputOn(Boolean tf);
SetSpeciesSaveOn(Boolean save_on);
SetTemperature(DoubleVector t);
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
SpeciesConcentrations2Module(DoubleVector species_conc);
StateSave(int istate);
StateApply(int istate);
StateDelete(int istate);
UseSolutionDensityVolume(Boolean tf);
WarningMessage(string warnstr);
</PRE>
</CODE>
@endhtmlonly
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitializeYAML InitializeYAML_DOCSTRING
*/
IRM_RESULT		InitializeYAML(std::string yamlfile);
#endif
/**
%define InitialPhreeqc2Concentrations_DOCSTRING
'Fills a vector (@a destination_c) with concentrations from solutions in the InitialPhreeqc instance.
The method is used to obtain concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
Args:
	destination_c(DoubleVector): Vector to receive the concentrations.The dimension of @a destination_c is set to @a ncomps times @a n_boundary,
		where @a ncomps is the number of components returned from :meth: `FindComponents` or :meth: `GetComponentCount`, and @a n_boundary
		is the dimension of the vector @a boundary_solution1.
		boundary_solution1  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
		Size is @a n_boundary.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2Concentrations InitialPhreeqc2Concentrations_DOCSTRING
*/
IRM_RESULT                                InitialPhreeqc2Concentrations(
								std::vector < double > & destination_c,
								const std::vector < int >    & boundary_solution1);
/**
%define InitialPhreeqc2Concentrations_DOCSTRING
'Fills a vector (@a destination_c) with concentrations from solutions in the InitialPhreeqc instance.
The method is used to obtain concentrations for boundary conditions that are mixtures of solutions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell. Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
Args:
	destination_c(DoubleVector): Vector of concentrations extracted from the InitialPhreeqc instance.
		The dimension of @a destination_c is set to @a ncomps times @a n_boundary,
		where @a ncomps is the number of components returned from :meth: `FindComponents` or :meth: `GetComponentCount`, and @a n_boundary
		is the dimension of the vectors @a boundary_solution1, @a boundary_solution2, and @a fraction1.
	boundary_solution1(IntVector): Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
		Size is @a n_boundary.
	boundary_solution2(IntVector): Vector of solution index numbers that that refer to solutions in the InitialPhreeqc instance
		and are defined to mix with @a boundary_solution1.
		Size is @a n_boundary.
	fraction1(DoubleVector): Fraction of boundary_solution1 that mixes with (1 - @a fraction1) of @a boundary_solution2.
		Size is @a n_boundary.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2Concentrations InitialPhreeqc2Concentrations_DOCSTRING
*/
IRM_RESULT								  InitialPhreeqc2Concentrations(
								std::vector < double > & destination_c,
								const std::vector < int >    & boundary_solution1,
								const std::vector < int >    & boundary_solution2,
								const std::vector < double > & fraction1);
/**
%define InitialPhreeqc2Module_DOCSTRING
'Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers.
@a Initial_conditions1 is used to select initial conditions, including solutions and reactants,
for each cell of the model, without mixing.
@a Initial_conditions1 is dimensioned 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model
(:meth: `GetGridCellCount`). The dimension of 7 refers to solutions and reactants in the following order:
(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
(5) SOLID_SOLUTIONS, and (6) KINETICS.
The definition initial_solution1[3*nxyz + 99] = 2, indicates that
cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
(created in the InitialPhreeqc instance either by :meth: `RunFile` or :meth: `RunString`).
Args:
	initial_conditions1(DoubleVector): Vector of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance.
		Size is 7 times @a nxyz. The order of definitions is given above.
		Negative values are ignored, resulting in no definition of that entity for that cell.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2Module InitialPhreeqc2Module_DOCSTRING
*/
IRM_RESULT  InitialPhreeqc2Module(const std::vector < int >    & initial_conditions1);
/**
%define InitialPhreeqc2Module_DOCSTRING
'Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers, possibly with mixing.
In its simplest form, @a  initial_conditions1 is used to select initial conditions, including solutions and reactants,
for each cell of the model, without mixing.
@a Initial_conditions1 is dimensioned 7 times @a  nxyz, where @a  nxyz is the number of grid cells in the user's model
(:meth: `GetGridCellCount`). The dimension of 7 refers to solutions and reactants in the following order:
(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
(5) SOLID_SOLUTIONS, and (6) KINETICS.
The definition initial_solution1[3*nxyz + 99] = 2, indicates that
cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
(either by :meth: `RunFile` or :meth: `RunString`).
@n@n
It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing,
@a initials_conditions2 contains numbers for a second entity that mixes with the entity defined in @a initial_conditions1.
@a Fraction1 contains the mixing fraction for @a initial_conditions1,
whereas (1 - @a fraction1) is the mixing fraction for @a initial_conditions2.
The definitions initial_solution1[3*nxyz + 99] = 2, initial_solution2[3*nxyz + 99] = 3,
fraction1[3*nxyz + 99] = 0.25 indicates that
cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3,
where the surface compositions have been defined in the InitialPhreeqc instance.
If the user number in @a initial_conditions2 is negative, no mixing occurs.
Args:
	initial_conditions1(IntVector): Vector of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance.
		Size is 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model (:meth: `GetGridCellCount`).
		The order of definitions is given above.
		Negative values are ignored, resulting in no definition of that entity for that cell.
	initial_conditions2(IntVector): Vector of solution and reactant index numbers that refer to
		definitions in the InitialPhreeqc instance.
		Nonnegative values of @a initial_conditions2 result in mixing with the entities defined in @a initial_conditions1.
		Negative values result in no mixing.
		Size is 7 times @a nxyz. The order of definitions is given above.
	fraction1(DoubleVector): Fraction of @a initial_conditions1 that mixes with (1 - @a fraction1)
		of @a initial_conditions2.
		Size is 7 times @a nxyz. The order of definitions is given above.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2Module InitialPhreeqc2Module_DOCSTRING
*/
IRM_RESULT InitialPhreeqc2Module(
const std::vector < int >    & initial_conditions1,
std::vector < int >    & initial_conditions2,
std::vector < double > & fraction1);
/**
%define InitialPhreeqc2SpeciesConcentrations_DOCSTRING
'Fills a vector @a destination_c with aqueous species concentrations from solutions in the InitialPhreeqc instance.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
Args:
	destination_c(DoubleVector): Vector of aqueous concentrations extracted from the InitialPhreeqc instance.
		The dimension of @a species_c is @a nspecies times @a n_boundary,
		where @a nspecies is the number of aqueous species returned from :meth: `GetSpeciesCount`,
		and @a n_boundary is the dimension of @a boundary_solution1.
	boundary_solution1(IntVector): Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2SpeciesConcentrations InitialPhreeqc2SpeciesConcentrations_DOCSTRING
*/
IRM_RESULT                                InitialPhreeqc2SpeciesConcentrations(
								std::vector < double > & destination_c,
								std::vector < int >    & boundary_solution1);
/**
%define InitialPhreeqc2SpeciesConcentrations_DOCSTRING
'Fills a vector @a destination_c with aqueous species concentrations from solutions in the InitialPhreeqc instance.
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
Args:
	destination_c(DoubleVector): Vector of aqueous concentrations extracted from the InitialPhreeqc instance.
		The dimension of @a species_c is @a nspecies times @a n_boundary,
		where @a nspecies is the number of aqueous species returned from :meth: `GetSpeciesCount`,
		and @a n_boundary is the dimension
		of @a boundary_solution1.
	boundary_solution1(IntVector): Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
	boundary_solution2(IntVector): Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance
		and are defined to mix with @a boundary_solution1. Size is same as @a boundary_solution1.
	fraction1(DoubleVector): Vector of fractions of @a boundary_solution1 that mix with (1 - @a fraction1) of @a boundary_solution2.
		Size is same as @a boundary_solution1.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqc2SpeciesConcentrations InitialPhreeqc2SpeciesConcentrations_DOCSTRING
*/
IRM_RESULT								  InitialPhreeqc2SpeciesConcentrations(
								std::vector < double > & destination_c,
								std::vector < int >    & boundary_solution1,
								std::vector < int >    & boundary_solution2,
								std::vector < double > & fraction1);
/**
%define InitialPhreeqcCell2Module_DOCSTRING
'A cell numbered @a n in the InitialPhreeqc instance is selected to populate a series of transport cells.
All reactants with the number @a n are transferred along with the solution.
If MIX @a n exists, it is used for the definition of the solution.
If @a n is negative, @a n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
All reactants for each cell in the list @a cell_numbers are removed before the cell
definition is copied from the InitialPhreeqc instance to the workers.
Args:
	n(int): Number that refers to a solution or MIX and associated reactants in the InitialPhreeqc instance.
	cell_numbers(int): A vector of grid-cell numbers (user's grid-cell numbering system) that
		will be populated with cell @a n from the InitialPhreeqc instance.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::InitialPhreeqcCell2Module InitialPhreeqcCell2Module_DOCSTRING
*/
IRM_RESULT                                InitialPhreeqcCell2Module(int n, const std::vector< int > &cell_numbers);
/**
%define LoadDatabase_DOCSTRING
'Load a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
Args:
	database(string): String containing the database name.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::LoadDatabase LoadDatabase_DOCSTRING
*/
IRM_RESULT                                LoadDatabase(const std::string &database);
/**
%define LogMessage_DOCSTRING
'Print a message to the log file.
Args:
	str(string): String to be printed.'
%enddef
%feature PhreeqcRM::LogMessage LogMessage_DOCSTRING
*/
void                                      LogMessage(const std::string &str);
/**
%define MpiAbort_DOCSTRING
'MPI only. Calls MPI_Abort, which aborts MPI, and makes the reaction module unusable.
Should be used only on encountering an unrecoverable error.'
%enddef
%feature PhreeqcRM::MpiAbort MpiAbort_DOCSTRING
*/
int                                       MpiAbort();
/**
%define MpiWorker_DOCSTRING
'MPI only. Nonroot processes (processes with :meth: `GetMpiMyself` > 0) must call MpiWorker to be able to
respond to messages from the root to accept data, perform calculations, and
(or) return data within the reaction module.
MpiWorker contains a loop that reads a message from root, performs a
task, and waits for another message from root.
:meth: `SetConcentrations`, :meth: `RunCells`, and :meth: `GetConcentrations`
are examples of methods that send a message from root to get the workers to perform a task.
The workers will respond to all methods that are designated "workers must be in the loop of MpiWorker"
in the MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls
:meth: `MpiWorkerBreak`.
@n@n
(Advanced) The list of tasks that the workers perform can be extended by using :meth: `SetMpiWorkerCallbackC`.
It is then possible to use the MPI processes to perform other developer-defined tasks,
such as transport calculations, without exiting from the MpiWorker loop.
Alternatively, root calls :meth: `MpiWorkerBreak` to allow the workers to continue past a call to MpiWorker.
The workers perform developer-defined calculations, and then MpiWorker is called again to respond to
requests from root to perform reaction-module tasks.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).
MpiWorker returns a value only when :meth: `MpiWorkerBreak` is called by root.'
%enddef
%feature PhreeqcRM::MpiWorker MpiWorker_DOCSTRING
*/
IRM_RESULT                                MpiWorker();
/**
%define MpiWorkerBreak_DOCSTRING
'MPI only. This method is called by root to force nonroot processes (processes with :meth: `GetMpiMyself` > 0)
to return from a call to :meth: `MpiWorker`.
:meth: `MpiWorker` contains a loop that reads a message from root, performs a
task, and waits for another message from root. The workers respond to all methods that are designated
"workers must be in the loop of MpiWorker" in the
MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls MpiWorkerBreak.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::MpiWorkerBreak MpiWorkerBreak_DOCSTRING
*/
IRM_RESULT                                MpiWorkerBreak();
/**
%define OpenFiles_DOCSTRING
'Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
based on the prefix defined by :meth: `SetFilePrefix`.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::OpenFiles OpenFiles_DOCSTRING
*/
IRM_RESULT                                OpenFiles(void);
/**
%define OutputMessage_DOCSTRING
'Print a message to the output file.
Args:
	str(string) String to be printed.'
%enddef
%feature PhreeqcRM::OutputMessage OutputMessage_DOCSTRING
*/
void                                      OutputMessage(const std::string &str);
/**
%define RunCells_DOCSTRING
'Runs a reaction step for all reaction cells in the reaction module.
Normally, tranport concentrations are transferred to the reaction cells (:meth: `SetConcentrations`) before
reaction calculations are run. The length of time over which kinetic reactions are integrated is set
by :meth: `SetTimeStep`. Other properties that may need to be updated as a result of the transport
calculations include porosity (:meth: `SetPorosity`), saturation (:meth: `SetSaturation`),
temperature (:meth: `SetTemperature`), and pressure (:meth: `SetPressure`).

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::RunCells RunCells_DOCSTRING
*/
IRM_RESULT                                RunCells(void);
/**
%define ReturnHandler_DOCSTRING
'Process an IRM_RESULT return code. If the return code is nonnegative, no action is taken. If the return code is negative,
the return code is decoded and printed as an error message along with the second argument (std::string). On an error,
the method will return the same return code, throw an exception, or exit the program depending on the setting for
:meth: `SetErrorHandlerMode`.
Args:
	result(IRM_RESULT): Return code to be processed.
	e_string(string): Error message to be printed in case of an error.
Returns:
	IRM_RESULT: The first argument to the method is returned.'
%enddef
%feature PhreeqcRM::ReturnHandler ReturnHandler_DOCSTRING
*/
IRM_RESULT                                ReturnHandler(IRM_RESULT result, const std::string &e_string);
/**
%define RunFile_DOCSTRING
'Run a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
files that modify the thermodynamic database should be run by all three sets of instances.
Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Files that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
Args:
	workers(Boolean): @a True, the workers will run the file; @a False, the workers will not run the file.
	initial_phreeqc(Boolean): @a True, the InitialPhreeqc instance will run the file; @a False, the InitialPhreeqc will not run the file.
	utility(Boolean): @a True, the Utility instance will run the file; @a False, the Utility instance will not run the file.
	chemistry_name(string): Name of the file to run.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::RunFile RunFile_DOCSTRING
*/
IRM_RESULT                                RunFile(bool workers, bool initial_phreeqc, bool utility,  const std::string & chemistry_name);
/**
%define RunString_DOCSTRING
'Run a PHREEQC input string. The first three arguments determine which
IPhreeqc instances will run
the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
strings that modify the thermodynamic database should be run by all three sets of instances.
Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Strings that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
Args:
	workers(Boolean): @a True, the workers will run the string; @a False, the workers will not run the string.
	initial_phreeqc(Boolean): @a True, the InitialPhreeqc instance will run the string; @a False, the InitialPhreeqc will not run the string.
	utility(Boolean): @a True, the Utility instance will run the string; @a False, the Utility instance will not run the string.
	input_string(string): String containing PHREEQC input.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::RunString RunString_DOCSTRING
*/
IRM_RESULT                                RunString(bool workers, bool initial_phreeqc, bool utility, const std::string & input_string);
/**
%define ScreenMessage_DOCSTRING
'Print message to the screen.
Args:
	str(string): String to be printed.'
%enddef
%feature PhreeqcRM::ScreenMessage ScreenMessage_DOCSTRING
*/
void                                      ScreenMessage(const std::string &str);
/**
%define SetComponentH2O_DOCSTRING
'Select whether to include H2O in the component list.
The concentrations of H and O must be known
accurately (8 to 10 significant digits) for the numerical method of
PHREEQC to produce accurate pH and pe values.
Because most of the H and O are in the water species,
it may be more robust (require less accuracy in transport) to
transport the excess H and O (the H and O not in water) and water.
The default setting (@a true) is to include water, excess H, and excess O as components.
A setting of @a false will include total H and total O as components.
SetComponentH2O must be called before :meth: `FindComponents`.

Args:
	tf(Boolean): @a True (default), excess H, excess O, and water are included in the component list;
		@a False, total H and O are included in the component list.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetComponentH2O SetComponentH2O_DOCSTRING
*/
IRM_RESULT                                SetComponentH2O(bool tf);
/**
%define SetConcentrations_DOCSTRING
'Use the vector of concentrations (@a c) to set the moles of components in each reaction cell.
The volume of water in a cell is the product of porosity (:meth: `SetPorosity`), saturation (:meth: `SetSaturation`),
and reference volume (:meth: `SetRepresentativeVolume`).
The moles of each component are determined by the volume of water and per liter concentrations.
If concentration units (:meth: `SetUnitsSolution`) are mass fraction, the
density (as specified by :meth: `SetDensity`) is used to convert from mass fraction to per mass per liter.
Args:
	c(DoubleVector): Vector of component concentrations. Size of vector is @a ncomps times @a nxyz,
		where @a ncomps is the number of components as determined
		by :meth: `FindComponents` or :meth: `GetComponentCount` and
		@a nxyz is the number of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetConcentrations SetConcentrations_DOCSTRING
*/
IRM_RESULT                                SetConcentrations(const std::vector< double > &c);
/**
%define SetCurrentSelectedOutputUserNumber_DOCSTRING
'Select the current selected output by user number. The user may define multiple SELECTED_OUTPUT
data blocks for the workers. A user number is specified for each data block. The value of
the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
for selected-output operations.
Args:
	n_user(int): User number of the SELECTED_OUTPUT data block that is to be used.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetCurrentSelectedOutputUserNumber SetCurrentSelectedOutputUserNumber_DOCSTRING
*/
IRM_RESULT								  SetCurrentSelectedOutputUserNumber(int n_user);
/**
%define SetDensity_DOCSTRING
'Set the density for each reaction cell. These density values are used
when converting from transported mass-fraction concentrations (:meth: `SetUnitsSolution`) to
produce per liter concentrations during a call to :meth: `SetConcentrations`.
They are also used when converting from reaction-cell concentrations to transport concentrations
(:meth: `GetConcentrations`), if :meth: `UseSolutionDensityVolume` is set to @a false.
Args:
	density(DoubleVector): Vector of densities. Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetDensity SetDensity_DOCSTRING
*/
IRM_RESULT                                SetDensity(const std::vector< double > &density);
/**
%define SetDumpFileName_DOCSTRING
'Set the name of the dump file. It is the name used by :meth: `DumpModule`.
Args:
	dump_name(string): Name of dump file.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetDumpFileName SetDumpFileName_DOCSTRING
*/
IRM_RESULT                                SetDumpFileName(const std::string & dump_name);
/**
%define SetErrorHandlerMode_DOCSTRING
'Set the action to be taken when the reaction module encounters an error.
Options are 0, return to calling program with an error return code (default);
1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
2, attempt to exit gracefully.
Args:
	mode(int): Error handling mode: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetErrorHandlerMode SetErrorHandlerMode_DOCSTRING
*/
IRM_RESULT                                SetErrorHandlerMode(int mode);
/**
%define SetErrorOn_DOCSTRING
'Set the property that controls whether error messages are generated and displayed.
Messages include PHREEQC "ERROR" messages, and
any messages written with :meth: `ErrorMessage`.

Args:
	tf(Boolean): @a True, enable error messages; @a False, disable error messages. Default is true.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetErrorOn SetErrorOn_DOCSTRING
*/
IRM_RESULT                                SetErrorOn(bool tf);
/**
%define SetFilePrefix_DOCSTRING
'Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
These files are opened by :meth: `OpenFiles`.
Args:
	prefix(string): Prefix used when opening the output and log files.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetFilePrefix SetFilePrefix_DOCSTRING
*/
IRM_RESULT                                SetFilePrefix(const std::string & prefix);

/**
%define SetGasCompMoles_DOCSTRING
'Transfer moles of gas components from
the vector given in the argument list (@a gas_moles) to each reaction cell.

Args:
	gas_moles(DoubleVector): Vector of moles of gas components.
		Dimension of the vector is set to @a ngas_comps times @a nxyz,
		where, @a ngas_comps is the result of :meth: `GetGasComponentsCount`,
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If the number of moles is set to a negative number, the gas component will
		not be defined for the GAS_PHASE of the reaction cell.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetGasCompMoles SetGasCompMoles_DOCSTRING
*/
IRM_RESULT                                SetGasCompMoles(const std::vector< double >& gas_moles);
/**
%define SetGasPhaseVolume_DOCSTRING
'Transfer volumes of gas phases from
the vector given in the argument list (@a gas_volume) to each reaction cell.
The gas-phase volume affects the gas-component pressures calculated for fixed-volume
gas phases. If a gas-phase volume is defined with this methood 
for a GAS_PHASE in a cell, 
the gas phase is forced to be a fixed-volume gas phase.

Args:
	gas_volume(DoubleVector): Vector of volumes for each gas phase.
		Dimension of the vector is @a nxyz,
		where @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		If the volume is set to a negative number for a cell, the gas-phase volume for that cell is
		not changed.

Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetGasPhaseVolume SetGasPhaseVolume_DOCSTRING
*/
IRM_RESULT                                SetGasPhaseVolume(const std::vector< double >& gas_volume);

/**
%define SetMpiWorkerCallbackC_DOCSTRING
'MPI and C/C++ only. Defines a callback function that allows additional tasks to be done
by the workers. The method :meth: `MpiWorker` contains a loop,
where the workers receive a message (an integer),
run a function corresponding to that integer,
and then wait for another message.
SetMpiWorkerCallbackC allows C or C++ developers to add another function
that responds to additional integer messages by calling developer-defined functions
corresponding to those integers.
:meth: `MpiWorker` calls the callback function when the message number
is not one of the PhreeqcRM message numbers.
Messages are unique integer numbers. PhreeqcRM uses integers in a range
beginning at 0. It is suggested that developers use message numbers starting
at 1000 or higher for their tasks.
The callback function calls a developer-defined function specified
by the message number and then returns to :meth: `MpiWorker` to wait for
another message.
@n@n
In C and C++, an additional pointer can be supplied to find the data necessary to do the task.
A void pointer may be set with :meth: `SetMpiWorkerCallbackCookie`. This pointer
is passed to the callback function through a void pointer argument in addition
to the integer message argument. The pointer may be to a struct or class instance
that provides a number of additional pointers to data. :meth: `SetMpiWorkerCallbackCookie`
must be called by each worker before :meth: `MpiWorker` is called.
@n@n
The motivation for this method is to allow the workers to perform other
tasks, for instance, parallel transport calculations, within the structure
of :meth: `MpiWorker`. The callback function
can be used to allow the workers to receive data, perform transport calculations,
and (or) send results, without leaving the loop of :meth: `MpiWorker`. Alternatively,
it is possible for the workers to return from :meth: `MpiWorker`
by a call to :meth: `MpiWorkerBreak` by root. The workers could then call
subroutines to receive data, calculate transport, and send data,
and then resume processing PhreeqcRM messages from root with another
call to :meth: `MpiWorker`.
Args:
	fcn(function pointer): A function that returns an integer and has an integer argument
		and a void * argument.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetMpiWorkerCallbackC SetMpiWorkerCallbackC_DOCSTRING
*/
IRM_RESULT								  SetMpiWorkerCallbackC(int (*fcn)(int *method, void * cookie));
/**
%define SetMpiWorkerCallbackCookie_DOCSTRING
'MPI and C/C++ only. Defines a void pointer that can be used by
C and C++ functions called from the callback function (:meth: `SetMpiWorkerCallbackC`)
to locate data for a task. The C callback function
that is registered with :meth: `SetMpiWorkerCallbackC` has
two arguments, an integer message to identify a task, and a void
pointer. SetMpiWorkerCallbackCookie sets the value of the
void pointer that is passed to the callback function.
The void pointer may be a pointer to a struct of class instance that
contains additonal pointers to data.
@param cookie           Void pointer that can be used by subroutines called from the callback function
to locate data needed to perform a task.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetMpiWorkerCallbackCookie SetMpiWorkerCallbackCookie_DOCSTRING
*/
IRM_RESULT								  SetMpiWorkerCallbackCookie(void * cookie);
/**
%define SetNthSelectedOutput_DOCSTRING
'Specify the current selected output by sequence number. The user may define multiple SELECTED_OUTPUT
data blocks for the workers. A user number is specified for each data block, and the blocks are
stored in user-number order. The value of
the argument @a n selects the sequence number of the SELECTED_OUTPUT definition that will be used
for selected-output operations.
Args:
	n(int): Sequence number of the SELECTED_OUTPUT data block that is to be used.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetNthSelectedOutput SetNthSelectedOutput_DOCSTRING
*/
IRM_RESULT                                SetNthSelectedOutput(int n);
/**
%define SetPartitionUZSolids_DOCSTRING
'Sets the property for partitioning solids between the saturated and unsaturated
parts of a partially saturated cell.

The option is intended to be used by saturated-only
flow codes that allow a variable water table.
The value has meaning only when saturations
less than 1.0 are encountered. The partially saturated cells
may have a small water-to-rock ratio that causes
reactions to proceed differently relative to fully saturated cells.
By setting  @a SetPartitionUZSolids to true, the
amounts of solids and gases are partioned according to the saturation.
If a cell has a saturation of 0.5, then
the water interacts with only half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell varies,
solids and gases are transferred between the
saturated and unsaturated (unreactive) reservoirs of the cell.
Unsaturated-zone flow and transport codes will probably use the default (false),
which assumes all gases and solids are reactive regardless of saturation.
Args:
	tf(Boolean): @a True, the fraction of solids and gases available for
		reaction is equal to the saturation;
		@a False (default), all solids and gases are reactive regardless of saturation.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetPartitionUZSolids SetPartitionUZSolids_DOCSTRING
*/
IRM_RESULT                                SetPartitionUZSolids(bool tf);
/**
%define SetPorosity_DOCSTRING
'Set the porosity for each reaction cell.
The volume of water in a reaction cell is the product of porosity, saturation
(:meth: `SetSaturation`), and representative volume (:meth: `SetRepresentativeVolume`).
Args:
	por(DoubleVector): Vector of porosities, unitless. Default is 0.1.
		Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetPorosity SetPorosity_DOCSTRING
*/
IRM_RESULT                                SetPorosity(const std::vector< double > &por);

/**
%define SetPressure_DOCSTRING
'Set the pressure for each reaction cell. Pressure effects are considered only in three of the
databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
Args:
	p(DoubleVector): Vector of pressures, in atm. Size of vector is @a nxyz,
		where @a nxyz is the number of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetPressure SetPressure_DOCSTRING
*/
IRM_RESULT                                SetPressure(const std::vector< double > &p);
/**
%define SetPrintChemistryMask_DOCSTRING
'Enable or disable detailed output for each reaction cell.
Printing for a reaction cell will occur only when the
printing is enabled with :meth: `SetPrintChemistryOn` and the @a cell_mask value is 1.
Args:
	cell_mask(IntVector): Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`). A value of 0 will
		disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetPrintChemistryMask SetPrintChemistryMask_DOCSTRING
*/
IRM_RESULT                                SetPrintChemistryMask(const std::vector<int> & cell_mask);
/**
%define SetPrintChemistryOn_DOCSTRING
'Set property that enables or disables printing detailed output from reaction calculations
to the output file for a set of cells defined by :meth: `SetPrintChemistryMask`.
The detailed output prints all of the output typical of a PHREEQC reaction calculation,
which includes solution descriptions and the compositions of all other reactants.
The output can be several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, :meth: `OpenFiles`).
For the worker instances, the output can be limited to a set of cells
(:meth: `SetPrintChemistryMask`) and, in general, the
amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
(applied by using :meth: `RunFile` or :meth: `RunString`).
Printing the detailed output for the workers is generally used only for debugging,
and PhreeqcRM will run significantly faster
when printing detailed output for the workers is disabled.

Args:
	workersDoubleVector): @a True, enable detailed printing in the worker instances;
		@a False, disable detailed printing in the worker instances.
	initial_phreeqcDoubleVector): @a True, enable detailed printing in the InitialPhreeqc instance;
		@a False, disable detailed printing in the InitialPhreeqc instance.
	utilityDoubleVector): @a True, enable detailed printing in the Utility instance;
		@a False, disable detailed printing in the Utility instance.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetPrintChemistryOn SetPrintChemistryOn_DOCSTRING
*/
IRM_RESULT                                SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
/**
%define SetRebalanceByCell_DOCSTRING
'Set the load-balancing algorithm.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a :meth: `RunCells`
calculation. Two algorithms are available; one uses individual times for each cell and
accounts for cells that were not run because
saturation was zero (default), and
the other assigns an average time to all cells.
The methods are similar, but limited testing indicates the default method performs better.
Args:
	tfDoubleVector): @a True, indicates individual cell times are used in rebalancing (default);
		@a False, indicates average times are used in rebalancing.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetRebalanceByCell SetRebalanceByCell_DOCSTRING
*/
IRM_RESULT                                SetRebalanceByCell(bool tf);
/**
%define SetRebalanceFraction_DOCSTRING
'Sets the fraction of cells that are transferred among threads or processes when rebalancing.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a :meth: `RunCells`
calculation. The rebalancing transfers cell calculations among threads or processes to
try to achieve an optimum balance. @a SetRebalanceFraction
adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
determine the actual number of cell transfers. A value of zero eliminates
load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
distribution and avoid possible oscillations
when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
Default is 0.5.

Args:
	f(float): Fraction from 0.0 to 1.0.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetRebalanceFraction SetRebalanceFraction_DOCSTRING
*/
IRM_RESULT                                SetRebalanceFraction(double f);
/**
%define SetRepresentativeVolume_DOCSTRING
'Set the representative volume of each reaction cell.
By default the representative volume of each reaction cell is 1 liter.
The volume of water in a reaction cell is determined by the product of the representative volume,
the porosity (:meth: `SetPorosity`), and the saturation (:meth: `SetSaturation`).
The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
within a couple orders of magnitude of 1.0.
Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes)
may cause non-convergence of the numerical method.
In these cases, a larger representative volume may help. Note
that increasing the representative volume also increases
the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
and others), which are defined as moles per representative volume.
Args:
	rv(DoubleVector): Vector of representative volumes, in liters. Default is 1.0 liter.
		Size of array is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetRepresentativeVolume SetRepresentativeVolume_DOCSTRING
*/
IRM_RESULT                                SetRepresentativeVolume(const std::vector< double > &rv);
/**
%define SetSaturation_DOCSTRING
'Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
The volume of water in a cell is the product of porosity (:meth: `SetPorosity`), saturation (@a SetSaturation),
and representative volume (:meth: `SetRepresentativeVolume`). As a result of a reaction calculation,
solution properties (density and volume) will change;
the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume data to calculate these changes. The methods :meth: `GetDensity`,
:meth: `GetSolutionVolume`, and :meth: `GetSaturation` can be used to account
for these changes in the succeeding transport calculation.
@a SetRepresentativeVolume should be called before initial conditions are defined for the reaction cells.

Args:
	sat(DoubleVector): Vector of saturations, unitless. Default 1.0. Size of vector is @a nxyz,
		where @a nxyz is the number of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetSaturation SetSaturation_DOCSTRING
*/
IRM_RESULT                                SetSaturation(const std::vector< double > &sat);
/**
%define SetScreenOn_DOCSTRING
'Set the property that controls whether messages are written to the screen.
Messages include information about rebalancing during :meth: `RunCells`, and
any messages written with :meth: `ScreenMessage`.

Args:
	tf(Boolean): @a True, enable screen messages; @a False, disable screen messages. Default is true.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetScreenOn SetScreenOn_DOCSTRING
*/
IRM_RESULT                                SetScreenOn(bool tf);
/**
%define SetSelectedOutputOn_DOCSTRING
'Set the property that controls whether selected-output results are available to be retrieved
with :meth: `GetSelectedOutput`. @a True indicates that selected-output results
will be accumulated during :meth: `RunCells` and can be retrieved with :meth: `GetSelectedOutput`;
@a False indicates that selected-output results will not
be accumulated during :meth: `RunCells`.

Args:
	tf(Boolean): @a True, enable selected output; @a False, disable selected output.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetSelectedOutputOn SetSelectedOutputOn_DOCSTRING
*/
IRM_RESULT                                SetSelectedOutputOn(bool tf);
/**
%define SetSpeciesSaveOn_DOCSTRING
'Sets the value of the species-save property.
This method enables or disables use of PhreeqcRM with multicomponent-diffusion transport calculations.
By default, concentrations of aqueous species are not saved.
Setting the species-save property to @a true allows
aqueous species concentrations to be retrieved
with :meth: `GetSpeciesConcentrations`, and solution compositions to be set with
:meth: `SpeciesConcentrations2Module`.
@a SetSpeciesSaveOn must be called before calls to :meth: `FindComponents`.

Args:
	save_on(Boolean): @a True indicates species concentrations are saved;
		@a False indicates species concentrations are not saved.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetSpeciesSaveOn SetSpeciesSaveOn_DOCSTRING
*/
IRM_RESULT                                SetSpeciesSaveOn(bool save_on);
/**
%define SetTemperature_DOCSTRING
'Set the temperature for each reaction cell. If @a SetTemperature is not called,
worker solutions will have temperatures as defined by initial conditions
(:meth: `InitialPhreeqc2Module` and :meth: `InitialPhreeqcCell2Module`).

Args:
	t(DoubleVector): Vector of temperatures, in degrees C.
		Size of vector is @a nxyz, where @a nxyz is the number
		of grid cells in the user's model (:meth: `GetGridCellCount`).
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetTemperature SetTemperature_DOCSTRING
*/
IRM_RESULT                                SetTemperature(const std::vector< double > &t);
/**
%define SetTime_DOCSTRING
'Set current simulation time for the reaction module.
Args:
	time(float): Current simulation time, in seconds.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetTime SetTime_DOCSTRING
*/
IRM_RESULT                                SetTime(double time);
/**
%define SetTimeConversion_DOCSTRING
'Set a factor to convert from seconds to user time units. Factor times seconds produces user time units.

Args:
	conv_factor(float): Factor to convert seconds to user time units.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetTimeConversion SetTimeConversion_DOCSTRING
*/
IRM_RESULT                                SetTimeConversion(double conv_factor);
/**
%define SetTimeStep_DOCSTRING
'Set current time step for the reaction module. This is the length
of time over which kinetic reactions are integrated.

Args:
	time_step(float): Time step, in seconds.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetTimeStep SetTimeStep_DOCSTRING
*/
IRM_RESULT                                SetTimeStep(double time_step);
/**
%define SetUnitsExchange_DOCSTRING
'Sets input units for exchangers.
In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
@a SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.

If a single EXCHANGE definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of exchangers will be the same regardless of porosity. 
For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.

Args:
	option(int): Units option for exchangers: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsExchange SetUnitsExchange_DOCSTRING
*/
IRM_RESULT                                SetUnitsExchange(int option);
/**
%define SetUnitsGasPhase_DOCSTRING
'Set input units for gas phases.
In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
@a SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single GAS_PHASE definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of a gas component will be the same regardless of porosity. 
For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.

Args:
	option(int): Units option for gas phases: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsGasPhase SetUnitsGasPhase_DOCSTRING
*/
IRM_RESULT                                SetUnitsGasPhase(int option);
/**
%define SetUnitsKinetics_DOCSTRING
'Set input units for kinetic reactants.

In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
@a SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single KINETICS definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of kinetic reactants will be the same regardless of porosity. 
For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.

Note that the volume of water in a cell in the reaction module is equal to the product of
porosity (:meth: `SetPorosity`), the saturation (:meth: `SetSaturation`), and representative volume (:meth:
`SetRepresentativeVolume`), which is usually less than 1 liter. It is important to write the RATES
definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
water, often by calculating the rate of reaction per liter of water and multiplying by the volume
of water (Basic function SOLN_VOL). 

Rates that depend on surface area of solids, are not dependent
on the volume of water. However, it is important to get the correct surface area for the kinetic
reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant) 
can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of 
reactant (Basic function M) in RATES to obtain the surface area.

Args:
	option(int): Units option for kinetic reactants: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsKinetics SetUnitsKinetics_DOCSTRING
*/
IRM_RESULT                                SetUnitsKinetics(int option);
/**
%define SetUnitsPPassemblage_DOCSTRING
'Set input units for pure phase assemblages (equilibrium phases).
In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
@a SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.

If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of a mineral will be the same regardless of porosity. 
For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.

Args:
	option(int): Units option for equilibrium phases: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsPPassemblage SetUnitsPPassemblage_DOCSTRING
*/
IRM_RESULT                                SetUnitsPPassemblage(int option);
/**
%define SetUnitsSolution_DOCSTRING
'Solution concentration units used by the transport model.
Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
PHREEQC defines solutions by the number of moles of each
element in the solution.

To convert from mg/L to moles
of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
multiplied by the solution volume,
which is the product of porosity (:meth: `SetPorosity`), saturation (:meth: `SetSaturation`),
and representative volume (:meth: `SetRepresentativeVolume`).
To convert from mol/L to moles
of element in the representative volume of a reaction cell, mol/L is
multiplied by the solution volume.
To convert from mass fraction to moles
of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
(:meth: `SetDensity`) and
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
calculated by PHREEQC, or (2) the volume of solution is the product of porosity (:meth: `SetPorosity`),
saturation (:meth: `SetSaturation`), and representative volume (:meth: `SetRepresentativeVolume`),
and the mass of solution is volume times density as defined by :meth: `SetDensity`.
Which option is used is determined by :meth: `UseSolutionDensityVolume`.

Args:
	option(int): Units option for solutions: 1, 2, or 3, default is 1, mg/L.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsSolution SetUnitsSolution_DOCSTRING
*/
IRM_RESULT                                SetUnitsSolution(int option);
/**
%define SetUnitsSSassemblage_DOCSTRING
'Set input units for solid-solution assemblages.
In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
@a SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@ P)*RV.

If a single SOLID_SOLUTION definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of a solid-solution component will be the same regardless of porosity. 
For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.

Args:
	option(int): Units option for solid solutions: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsSSassemblage SetUnitsSSassemblage_DOCSTRING
*/
IRM_RESULT                                SetUnitsSSassemblage(int option);
/**
%define SetUnitsSurface_DOCSTRING
'Set input units for surfaces.
In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
@a SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (:meth: `SetRepresentativeVolume`);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (:meth: `SetPorosity`); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single SURFACE definition is used for cells with different initial porosity, 
the three options scale quite differently. 
For option 0, the number of moles of surface sites will be the same regardless of porosity. 
For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.

Args:
	option(int): Units option for surfaces: 0, 1, or 2.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SetUnitsSurface SetUnitsSurface_DOCSTRING
*/
IRM_RESULT                                SetUnitsSurface(int option);
/**
%define SpeciesConcentrations2Module_DOCSTRING
'Set solution concentrations in the reaction cells
based on the vector of aqueous species concentrations (@a species_conc).
This method is intended for use with multicomponent-diffusion transport calculations,
and :meth: `SetSpeciesSaveOn` must be set to @a true.
The list of aqueous species is determined by :meth: `FindComponents` and includes all
aqueous species that can be made from the set of components.
The method determines the total concentration of a component
by summing the molarities of the individual species times the stoichiometric
coefficient of the element in each species.
Solution compositions in the reaction cells are updated with these component concentrations.

Args:
	species_conc(DoubleVector): Vector of aqueous species concentrations. Dimension of the array is @a nspecies times @a nxyz,
		where  @a nspecies is the number of aqueous species (:meth: `GetSpeciesCount`),
		and @a nxyz is the number of user grid cells (:meth: `GetGridCellCount`).
		Concentrations are moles per liter.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::SpeciesConcentrations2Module SpeciesConcentrations2Module_DOCSTRING
*/
IRM_RESULT								  SpeciesConcentrations2Module(std::vector< double > & species_conc);

/**
%define StateSave_DOCSTRING
'Save the state of the chemistry in all model cells, including SOLUTIONs, 
EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs. 
Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs 
will be saved for each cell, if they have been defined in the worker IPhreeqc instances. 
The distribution of cells among the workers and the chemistry of fully or partially 
unsaturated cells are also saved. The state is saved in memory; use :meth: `DumpModule` to save the state
to file. PhreeqcRM can be reset to this state by using :meth: `StateApply`.
A state is identified by an integer, and multiple states can be saved. 

Args:
	istate(int): Integer identifying the state that is saved. 
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::StateSave StateSave_DOCSTRING
*/
IRM_RESULT StateSave(int istate);
/**
%define StateApply_DOCSTRING
'Reset the state of the module to a state previously saved with :meth: `StateSave`. 
The chemistry of all model cells are reset, including SOLUTIONs,
EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
will be reset for each cell, if they were defined in the worker IPhreeqc instances
at the time the state was saved.
The distribution of cells among the workers and the chemistry of fully or partially
unsaturated cells are also reset to the saved state. 
The state to be applied is identified by an integer.

Args:
	istate(int): Integer identifying the state that is to be applied.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::StateApply StateApply_DOCSTRING
*/
IRM_RESULT StateApply(int istate);
/**
%define StateDelete_DOCSTRING
'Delete a state previously saved with :meth: `StateSave`.

Args:
	istate(int): Integer identifying the state that is to be deleted.
Returns:
	IRM_RESULT: 0 is success, negative is failure (See :meth: `DecodeError`).'
%enddef
%feature PhreeqcRM::StateDelete StateDelete_DOCSTRING
*/
IRM_RESULT StateDelete(int istate);
/**
%define UseSolutionDensityVolume_DOCSTRING
'Determines the volume and density to use when converting from the reaction-cell concentrations
to transport concentrations (:meth: `GetConcentrations`).
Two options are available to convert concentration units:
(1) the density and solution volume calculated by PHREEQC are used, or
(2) the specified density (:meth: `SetDensity`)
and solution volume are determined by the product of
saturation (:meth: `SetSaturation`), porosity (:meth: `SetPorosity`),
and representative volume (:meth: `SetRepresentativeVolume`).
Transport models that consider density-dependent flow will probably use the
PHREEQC-calculated density and solution volume (default),
whereas transport models that assume constant-density flow will probably use
specified values of density and solution volume.
Only the following databases distributed with PhreeqcRM have molar-volume information
needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
Density is only used when converting to or from transport units of mass fraction.

Args:
	tf(Boolean): @a True indicates that the solution density and volume as
		calculated by PHREEQC will be used to calculate concentrations.
		@a False indicates that the solution density set by :meth: `SetDensity` and the volume determined by the
		product of  :meth: `SetSaturation`, :meth: `SetPorosity`, and :meth: `SetRepresentativeVolume`,
		will be used to calculate concentrations retrieved by :meth: `GetConcentrations`.'
%enddef
%feature PhreeqcRM::UseSolutionDensityVolume UseSolutionDensityVolume_DOCSTRING
*/
void UseSolutionDensityVolume(bool tf);
/**
%define WarningMessage_DOCSTRING
'Print a warning message to the screen and the log file.

Args:
	warnstr(string): String to be printed.'
%enddef
%feature PhreeqcRM::WarningMessage WarningMessage_DOCSTRING
*/
void                                      WarningMessage(const std::string &warnstr);
// BMI data and methods
private:
std::map<std::string, class BMI_Var> bmi_var_map;
std::vector< std::string > bmi_input_vars;
std::vector< std::string > bmi_output_vars;
public:

//void BMI_Finalize();
/**
%define GetComponentName_DOCSTRING
'Basic Model Interface method that returns the component name--PhreeqcRM. The BMI interface to PhreeqcRM is
only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
prefix) provide a complete interface, and it is expected that the native methods will be used in preference to the BMI_
methods.

Returns:
	string: The string "PhreeqcRM".'
%enddef
%feature PhreeqcRM::GetComponentName GetComponentName_DOCSTRING
*/
std::string BMI_GetComponentName() { return "PhreeqcRM"; }


/**
%define GetCurrentTime_DOCSTRING
'Basic Model Interface method that returns the current simulation time, in seconds. (Same as :meth: `GetTime`.)
The reaction module does not change the time value, so the
returned value is equal to the default (0.0) or the last time set by
:meth: `BMI`_SetValue("Time", time) or :meth: `SetTime`.
Returns:
	float: The current simulation time, in seconds.'
%enddef
%feature PhreeqcRM::GetCurrentTime GetCurrentTime_DOCSTRING
*/
double BMI_GetCurrentTime() { return this->GetTime(); }


/**
%define GetEndTime_DOCSTRING
'Basic Model Interface method that returns :meth: `BMI`_GetCurrentTime plus :meth: `BMI`_GetTimeStep, in seconds.
Returns:
	float: The end of the time step, in seconds.'
%enddef
%feature PhreeqcRM::GetEndTime GetEndTime_DOCSTRING
*/
double BMI_GetEndTime() { return this->GetTime() + this->GetTimeStep(); }


/**
%define GetInputItemCount_DOCSTRING
'Basic Model Interface method that returns count of input variables that can be set with :meth: `BMI`_SetValue.
Returns:
	int: Count of input variables that can be set with :meth: `BMI`_SetValue.
'
%enddef
%feature PhreeqcRM::GetInputItemCount GetInputItemCount_DOCSTRING
*/
int BMI_GetInputItemCount() { return (int)this->bmi_input_vars.size(); }


/**
%define GetInputVarNames_DOCSTRING
'Basic Model Interface method that returns a list of the variable names that can be set with :meth: `BMI`_SetValue.
Returns:
	tuple of strings: A list of the variable names that can be set with :meth: `BMI`_SetValue.'
%enddef
%feature PhreeqcRM::GetInputVarNames GetInputVarNames_DOCSTRING
*/
std::vector< std::string > BMI_GetInputVarNames() { return this->bmi_input_vars; }


/**
%define GetOutputItemCount_DOCSTRING
'Basic Model Interface method that returns count of output variables that can be retrieved with :meth: `BMI`_GetValue.
Returns:
	int: Count of output variables that can be retrieved with :meth: `BMI`_GetValue.'
%enddef
%feature PhreeqcRM::GetOutputItemCount GetOutputItemCount_DOCSTRING
*/
int BMI_GetOutputItemCount() { return (int)this->bmi_output_vars.size(); }


/**
%define GetOutputVarNames_DOCSTRING
'Basic Model Interface method that returns a list of the variable names that can be retrieved with :meth: `BMI`_GetValue.
Returns:
	tuple of strings: A list of the variable names that can be retrieved with :meth: `BMI`_GetValue.'
%enddef
%feature PhreeqcRM::GetOutputVarNames GetOutputVarNames_DOCSTRING
*/
std::vector< std::string > BMI_GetOutputVarNames() { return this->bmi_output_vars; };


/**
%define GetTimeStep_DOCSTRING
'Basic Model Interface method that returns the current simulation time step, in seconds. (Same as :meth: `GetTimeStep`.)
The reaction module does not change the time-step value, so the
returned value is equal to the last time step set by
:meth: `BMI`_SetValue("TimeStep", time_step) or :meth: `SetTimeStep`.
Returns:
	float: The current simulation time step, in seconds.'
%enddef
%feature PhreeqcRM::GetTimeStep GetTimeStep_DOCSTRING
*/
double BMI_GetTimeStep() { return this->GetTimeStep(); }


/**
%define GetTimeUnits_DOCSTRING
'Basic Model Interface method that returns the time units of PhreeqcRM.
All time units are seconds for PhreeqcRM.
Returns:
	string: Returns the string "seconds".'
%enddef
%feature PhreeqcRM::GetTimeUnits GetTimeUnits_DOCSTRING
*/
std::string BMI_GetTimeUnits() { return "seconds"; };


/**
%define GetValue_DOCSTRING
'Basic Model Interface method that retrieves model variables. Only variables in the list
provided by :meth: `BMI`_GetOutputVarNames can be retrieved. The BMI interface to PhreeqcRM is
only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
prefix) provide a complete interface.
Args:
	name(string): Name of the variable to retrieve.
	dest(type defined below): Variable in which to place results.

Variable names for the first argument (@a name) and variable type of the
second argument (@a dest).
@n "ComponentCount", @a dest: int;
@n "Components", @a dest: tuple of strings;
@n "Concentrations", @a dest: DoubleVector;
@n "CurrentSelectedOutputUserNumber", @a dest: int;
@n "Density", @a dest: DoubleVector;
@n "ErrorString", @a dest: std::string;
@n "FilePrefix", @a dest: std::string;
@n "Gfw", @a dest: std::vector< double >;
@n "GridCellCount", @a dest: int;
@n "InputVarNames", @a dest: tuple of strings;
@n "OutputVarNames", @a dest: tuple of strings;
@n "Porosity", @a dest: DoubleVector;
@n "Pressure", @a dest: DoubleVector;
@n "Saturation", @a dest: DoubleVector;
@n "SelectedOutput", @a dest: DoubleVector;
@n "SelectedOutputColumnCount", @a dest: int;
@n "SelectedOutputCount", @a dest: int;
@n "SelectedOutputHeadings", @a dest: tuple of strings;
@n "SelectedOutputOn", @a dest: Boolean;
@n "SelectedOutputRowCount", @a dest: int;
@n "SolutionVolume", @a dest: DoubleVector;
@n "Temperature", @a dest: DoubleVector;
@n "Time",	@a dest: double;
@n "TimeStep",	@a dest: double.'
%enddef
%feature PhreeqcRM::GetValue GetValue_DOCSTRING
*/
void BMI_GetValue(std::string name, void* dest);
void BMI_GetValue(std::string name, bool& dest);
void BMI_GetValue(std::string name, double& dest);
void BMI_GetValue(std::string name, int& dest);
void BMI_GetValue(std::string name, std::string& dest);
void BMI_GetValue(std::string name, std::vector < double >& dest);
void BMI_GetValue(std::string name, std::vector < std::string >& dest);


/**
%define GetVarItemsize_DOCSTRING
'Basic Model Interface method that retrieves size of an 
individual item that can be set or retrived.
Sizes may be sizeof(int), sizeof(double), 
or a character length for string variables. Only variables in the list
provided by :meth: `BMI`_GetInputVarNames can be set. 
Only variables in the list
provided by :meth: `BMI`_GetOutputVarNames can be retrieved. 
Args:
	name(string): Name of the variable to retrieve size.
Returns:
	int: Size of one element of variable.'
%enddef
%feature PhreeqcRM::GetVarItemsize GetVarItemsize_DOCSTRING
*/
int BMI_GetVarItemsize(std::string name);


/**
%define GetVarNbytes_DOCSTRING
'Basic Model Interface method that retrieves the total number of bytes that are set for a variable with
:meth: `BMI`_SetValue or retrieved for a variable with :meth: `BMI`_GetValue.
Only variables in the list
provided by :meth: `BMI`_GetInputVarNames can be set. 
Only variables in the list
provided by :meth: `BMI`_GetOutputVarNames can be retrieved. 
Args:
	name(string): Name of the variable to retrieve total bytes.
Returns:
	int: Total number of bytes set or retrieved for variable.'
%enddef
%feature PhreeqcRM::GetVarNbytes GetVarNbytes_DOCSTRING
*/
int BMI_GetVarNbytes(std::string name);


/**
%define GetVarType_DOCSTRING
'Basic Model Interface method that retrieves the type of a variable that can be set with
:meth: `BMI`_SetValue or retrieved with :meth: `BMI`_GetValue. Types are "int", "double", or "string".
Only variables in the list
provided by :meth: `BMI`_GetInputVarNames can be set. 
Only variables in the list
provided by :meth: `BMI`_GetOutputVarNames can be retrieved. 

Args:
	name(string): Name of the variable to retrieve type.
Returns:
	string: Character string of variable type.'
%enddef
%feature PhreeqcRM::GetVarType GetVarType_DOCSTRING
*/
std::string BMI_GetVarType(std::string name);


/**
%define GetVarUnits_DOCSTRING
'Basic Model Interface method that retrieves the units of a variable that can be set with
:meth: `BMI`_SetValue or retrieved with :meth: `BMI`_GetValue.
Only variables in the list
provided by :meth: `BMI`_GetInputVarNames can be set. 
Only variables in the list
provided by :meth: `BMI`_GetOutputVarNames can be retrieved. 

Args:
	name(string): Name of the variable to retrieve type.
Returns:
	string: Character string of units for variable.'
%enddef
%feature PhreeqcRM::GetVarUnits GetVarUnits_DOCSTRING
*/
std::string BMI_GetVarUnits(std::string name);


/**
%define Initialize_DOCSTRING
'Basic Model Interface method that can be used to initialize a PhreeqcRM instance. This method is equivalent to
:meth: `InitializeYAML`. A YAML file can be used in initialization. The file contains a YAML map of PhreeqcRM methods
and the arguments corresponding to the method. For example,
@htmlonly
<CODE>
<PRE>
LoadDatabase: phreeqc.dat
RunFile:
workers: true
initial_phreeqc: true
utility: true
chemistry_name: advect.pqi
</PRE>
</CODE>
@endhtmlonly

BMI_Initialize will read the YAML file and execute the specified methods with the specified arguments. Using YAML
terminology, the argument(s) for a method may be a scalar, a sequence, or a map, depending if the argument is
a single item, a single vector, or there are multiple arguments. In the case of a map, the name associated
with each argument (for example "chemistry_name" above) is arbitrary. The names of the map keys for map
arguments are not used in parsing the YAML file; only the order of the arguments is important.

The PhreeqcRM methods that can be specified in a YAML file include:
@htmlonly
<CODE>
<PRE>
CloseFiles(void);
CreateMapping(IntVector& grid2chem);
DumpModule();
FindComponents();
InitialPhreeqc2Module(IntVector initial_conditions1);
InitialPhreeqc2Module(IntVector initial_conditions1, IntVector initial_conditions2, DoubleVector fraction1);
InitialPhreeqcCell2Module(int n, IntVector cell_numbers);
LoadDatabase(string database);
OpenFiles(void);
OutputMessage(string str);
RunCells(void);
RunFile(Boolean workers, Boolean initial_phreeqc, Boolean utility, string chemistry_name);
RunString(Boolean workers, Boolean initial_phreeqc, Boolean utility, string input_string);
ScreenMessage(string str);
SetComponentH2O(Boolean tf);
SetConcentrations(DoubleVector c);
SetCurrentSelectedOutputUserNumber(int n_user);
SetDensity(DoubleVector density);
SetDumpFileName(string dump_name);
SetErrorHandlerMode(int mode);
SetErrorOn(Boolean tf);
SetFilePrefix(string prefix);
SetGasCompMoles(DoubleVector gas_moles);
SetGasPhaseVolume(DoubleVector gas_volume);
SetPartitionUZSolids(Boolean tf);
SetPorosity(DoubleVector por);
SetPressure(DoubleVector p);
SetPrintChemistryMask(IntVector cell_mask);
SetPrintChemistryOn(Boolean workers, Boolean initial_phreeqc, Boolean utility);
SetRebalanceByCell(Boolean tf);
SetRebalanceFraction(float f);
SetRepresentativeVolume(DoubleVector rv);
SetSaturation(DoubleVector sat);
SetScreenOn(Boolean tf);
SetSelectedOutputOn(Boolean tf);
SetSpeciesSaveOn(Boolean save_on);
SetTemperature(DoubleVector t);
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
SpeciesConcentrations2Module(DoubleVector species_conc);
StateSave(int istate);
StateApply(int istate);
StateDelete(int istate);
UseSolutionDensityVolume(Boolean tf);
WarningMessage(string warnstr);
</PRE>
</CODE>
@endhtmlonly'
Args:
	config_file(string): File with YAML definitions for PhreeqcRM initialization.
%enddef
%feature PhreeqcRM::Initialize Initialize_DOCSTRING
*/
void BMI_Initialize(std::string config_file) { InitializeYAML(config_file); };
#endif

/**
%define SetValue_DOCSTRING
'Basic Model Interface method that sets model variables. Only variables in the list
provided by :meth: `BMI`_GetInputVarNames can be set. The BMI interface to PhreeqcRM is
only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
prefix) provide a complete interface, and it is expected that the native methods will be used in preference to the BMI_
methods.

Variable names for the first argument
of BMI_SetValue and the equivalent PhreeqcRM method are as follows:
"Concentrations", :meth: `SetConcentrations`;
"Density", :meth: `SetDensity`;
"FilePrefix", :meth: `SetFilePrefix`;
"NthSelectedOutput", :meth: `SetNthSelectedOutput`;
"Porosity", :meth: `SetPorosity`;
"Pressure", :meth: `SetPressure`;
"Saturation", :meth: `SetSaturation`;
"SelectedOutputOn", :meth: `SetSelectedOutputOn`;
"Temperature", :meth: `SetTemperature`;
"Time", :meth: `SetTime`;
"TimeStep", :meth: `SetTimeStep`.
'
%enddef
%feature PhreeqcRM::SetValue SetValue_DOCSTRING
*/
void BMI_SetValue(std::string name, void* src);
void BMI_SetValue(std::string name, bool& src);
void BMI_SetValue(std::string name, double& src);
void BMI_SetValue(std::string name, int& src);
void BMI_SetValue(std::string name, const std::string& src);
void BMI_SetValue(std::string name, std::vector < double >& src);
void BMI_SetValue(std::string name, std::vector < int >& src);
void BMI_SetValue(std::string name, std::vector < std::string >& src);

/**
%define Update_DOCSTRING
'Basic Model Interface method that runs PhreeqcRM for one time step. This method is equivalent to
:meth: `RunCells`. PhreeqcRM will equilibrate the solutions with all equilibrium reactants (EQUILIBRIUM_PHASES,
EXCHANGE, GAS_PHASE, SOLID_SOLUTIONS, and SURFACE) and
integrate KINETICS reactions for the specified time step (:meth: `SetTimeStep`).'
%enddef
%feature PhreeqcRM::Update Update_DOCSTRING
*/
void BMI_Update(void) { this->RunCells(); };
