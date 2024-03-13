#ifdef USE_YAML
/*! @file YAMLPhreeqcRM.h
*	@brief C++ header file for YAMLPhreeqcRM
* 
* YAMLPhreeqcRM can be used in preprocessors or
* Graphical User Interfaces to store initialization data for BMIPhreeqcRM or
* PhreeqcRM instances. PhreeqcRM methods and data can be stored in a YAML file.
* After an instance of BMIPhreeqcRM or PhreeqcRM has been created, the method
* Initialize or InitializeYAML can be used to run the specified methods
* with the specified data to define properties and initial conditions for the
* instance.
*/
#ifndef INC_YAMLPHREEQCRM_H
#define INC_YAMLPHREEQCRM_H
#include <map>
#include <mutex>
#include <string>
#include "yaml-cpp/yaml.h"
#include "IrmResult.h"
//#pragma once

#include "irm_dll_export.h"

/**
 * @class YAMLPhreeqcRM
 *
 * @brief YAML helper class
 */
class IRM_DLL_EXPORT YAMLPhreeqcRM
{
private:
	YAML::Node YAML_doc;
	std::map<std::string, int> method_map;
	YAML::EmitterStyle::value style;
protected:
	friend class YAMLPhreeqcRMLib;
	static std::map<size_t, YAMLPhreeqcRM*> Instances;
	static std::mutex InstancesLock;
	static size_t InstancesIndex;
	size_t Index;
public:

	/**
	Constructor
	*/
	YAMLPhreeqcRM();
	~YAMLPhreeqcRM();
	/**
	Clears all definitions from the YAML document.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::ofstream ofs = std::ofstream(YAML_filename.c_str(), std::ofstream::out);
	ofs << yrm.GetYAMLDoc();
	ofs.close();
	yrm.clear();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void Clear();	
	/**
	 *  Retrieves the id of this object.  Each instance receives an id which is incremented for each instance
	 *  starting with the value zero.
	 *  @return                 The id.
	 */
	int                      GetId(void)const;
	/**
	Returns a constant reference to the YAML document.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::ofstream ofs = std::ofstream(YAML_filename.c_str(), std::ofstream::out);
	ofs << yrm.GetYAMLDoc();
	ofs.close();
	yrm.clear();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	const YAML::Node& GetYAMLDoc() { return this->YAML_doc; };
	/**
	Writes the YAML document to file.
	@param file_name        Name of file where YAML document will be written.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.WriteYAMLDoc(YAML_filename);
	yrm.clear();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void WriteYAMLDoc(std::string file_name);
	// methods
	/**
	@a YAMLAddOutputVars inserts data into the YAML document for the PhreeqcRM method @a CloseFiles.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance. 
	<p>
	@a AddOutputVars allows selection of sets of variables that can be retieved
	by the GetValue method.Sets of variables can be included or excluded with
	multiple calls to this method. All calls must precede the final call to
	FindComponents. FindComponents generates SELECTED_OUTPUT 333 and
	USER_PUNCH 333 data blocks that make the variables accessible. Variables will
	only be accessible if the system includes the given reactant; for example, no
	gas variables will be created if there are no GAS_PHASEs in the model.
	</p>
	@param option A string value, among those listed below, that includes or
	excludes variables from GetOutputVarNames, GetValue, and other
	BMI methods.
	@param def A string value that can be "false", "true", or a list of items to be included as
	accessible variables.A value of "false", excludes all variables of the given type; a
	value of "true" includes all variables of the given type for the current system; a list
	specifies a subset of items of the given type.
	<p>
	Values for the the parameter @a option:
	</p>
	@n@n @a AddOutputVars: False excludes all variables; True causes the settings for each variable group
	to determine the variables that will be defined.Default True;
	@n@n @a SolutionProperties: False excludes all solution property variables; True includes variables pH, pe,
	alkalinity, ionic strength, water mass, charge balance, percent error, and specific conductance.
	Default True.
	@n@n @a SolutionTotalMolalities: False excludes all total element and element redox state variables;
	True includes all elements and element redox state variables for the system defined for the
	calculation; list restricts variables to the specified elements and redox states.
	Default True.
	@n@n @a ExchangeMolalities: False excludes all variables related to exchange; True includes all
	variables related to exchange; list includes variables for the specified exchange species.
	Default True.
	@n@n @a SurfaceMolalities: False excludes all variables related to surfaces; True includes all
	variables related to surfaces; list includes variables for the specified surface species.
	Default True.
	@n@n @a EquilibriumPhases: False excludes all variables related to equilibrium phases; True includes all
	variables related to equilibrium phases; list includes variables for the specified
	equilibiurm phases.Default True.
	@n@n @a Gases: False excludes all variables related to gases; True includes all
	variables related to gases; list includes variables for the specified gas components.Default True.
	@n@n @a KineticReactants: False excludes all variables related to kinetic reactants; True includes all
	variables related to kinetic reactants; list includes variables for the specified kinetic
	reactants.Default True.
	@n@n @a SolidSolutions: False excludes all variables related to solid solutions; True includes all
	variables related to solid solutions; list includes variables for the specified solid solutions
	components.Default True.
	@n@n @a CalculateValues: False excludes all calculate values; True includes all
	calculate values; list includes the specified calculate values.CALCLUATE_VALUES can be
	used to calculate geochemical quantities not available in the other sets of variables.
	Default True.
	@n@n @a SolutionActivities: False excludes all aqueous species; True includes all
	aqueous species; list includes only the specified aqueous species.Default False.
	@n@n @a SolutionMolalities: False excludes all aqueous species; True includes all
	aqueous species; list includes only the specified aqueous species.Default False.
	@n@n @a SaturationIndices: False excludes all saturation indices; True includes all
	saturation indices; list includes only the specified saturation indices. Default False.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	YAMLPhreeqcRM yrm;
	yrm.YAMLAddOutputVars("SaturationIndices", "true");
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLAddOutputVars(std::string option, std::string def);
	/**
	@a YAMLCloseFiles inserts data into the YAML document for the PhreeqcRM method @a CloseFiles.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance. 
	<p>
	@a CloseFiles closes the output and log files.
	</p>
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	YAMLPhreeqcRM yrm;
	yrm.YAMLCloseFiles();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLCloseFiles(void);
	/**
	@a YAMLCreateMapping inserts data into the YAML document for the PhreeqcRM method @a CreateMapping.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a CreateMapping
	provides a mapping from grid cells in the user's model to reaction cells for which chemistry needs to be run.
	The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells
	for which chemistry must be run.
	The array @a grid2chem of size @a nxyz (the number of grid cells)
	must contain the set of all integers 0 <= @a i < @a count_chemistry,
	where @a count_chemistry is a number less than or equal to @a nxyz.
	Inactive cells are assigned a negative integer.
	The mapping may be many-to-one to account for symmetry.
	Default is a one-to-one mapping--all user grid cells are reaction cells
	(equivalent to @a grid2chem values of 0,1,2,3,...,nxyz-1).
	</p>
	@param grid2chem        A vector of integers: Nonnegative is a reaction-cell number (0 based),
	negative is an inactive cell. Vector is of size @a nxyz (number of grid cells).
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	// For demonstation, two equivalent rows by symmetry
	std::vector<int> grid2chem;
	grid2chem.resize(nxyz, -1);
	for (int i = 0; i < nxyz/2; i++)
	{
	  grid2chem[i] = i;
	  grid2chem[i + nxyz/2] = i;
	}
	yrm.YAMLCreateMapping(grid2chem);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLCreateMapping(std::vector< int >& grid2chem);
	/**
	@a YAMLDumpModule inserts data into the YAML document for the PhreeqcRM method @a DumpModule.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a DumpModule writes the contents of all workers to file in _RAW formats (see appendix of PHREEQC version 3 manual),
	including SOLUTIONs and all reactants.
	</p>
	@param dump_on          Signal for writing the dump file, true or false.
	@param append           Signal to append to the contents of the dump file, true or false.
	@see                    @ref YAMLSetDumpFileName.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	bool dump_on = true;
	bool append = false;
	yrm.YAMLSetDumpFileName("Advect_cpp.dmp");
	yrm.YAMLDumpModule(dump_on, append);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLDumpModule(bool dump_on, bool append);
	/**
	@a YAMLThreadCount inserts data into the YAML document for the PhreeqcRM method @a ThreadCount.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a ThreadCount provides the number of threads to use in OpenMP multiprocessing when used
	to initialize a BMIPhreeqcRM instance, provided the BMIPhreeqcRM instance was created
	with the default constructor--the constructor with no arguments.
	</p>
	@param nthreads Number of threads to use in
	parallelnprocessing with OpenMP.
	*/
	void YAMLThreadCount(int nthreads);
	/**
	@a YAMLFindComponents inserts data into the YAML document for the PhreeqcRM method @a FindComponents.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a FindComponents accumulates a list of elements. Elements are those that have been
	defined in a solution or any other reactant
	(EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
	This method can be called multiple times and the list that is created is cummulative.
	The list is the set of components that needs to be transported. By default the list
	includes water, excess H and excess O (the H and O not contained in water);
	alternatively, the list may be set to contain total H and total O (@ref YAMLSetComponentH2O),
	which requires transport results to be accurate to eight or nine significant digits.
	If multicomponent diffusion (MCD) is to be modeled,
	there is a capability to retrieve aqueous species concentrations 
	and to set new solution concentrations after
	MCD by using individual species concentrations
	(@ref YAMLSpeciesConcentrations2Module).
	To use these methods, the save-species property needs to be turned on (@ref YAMLSetSpeciesSaveOn).
	If the save-species property is on, FindComponents will generate
	a list of aqueous species,
	their diffusion coefficients at 25 C,
	and their charge.
	</p>
	@see
	@ref YAMLSetComponentH2O,
	@ref YAMLSetSpeciesSaveOn,
	@ref YAMLSpeciesConcentrations2Module.
	<p>
	The @a FindComponents method also generates lists of reactants--equilibrium phases,
	exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
	The lists are cumulative, including all reactants that were
	defined in the initial phreeqc instance at any time FindComponents was called.
	In addition, a list of phases is generated for which saturation indices may be calculated from the
	cumulative list of components.
	</p>
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLFindComponents();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLFindComponents();
	/**
	@a YAMLInitialSolutions2Module inserts data into the YAML document for the PhreeqcRM method 
	@a InitialSolutions2Module.
	When the YAML document is written to file it can be processed by the method 
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialSolutions2Module transfers SOLUTION definitions from the InitialPhreeqc 
	instance to the reaction-module workers.
	@a solutions is a vector of SOLUTION index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param solutions Vector of index numbers that is dimensioned @a nxyz, 
	where @a nxyz is the number of grid cells in the user's model. 
	 */
	void YAMLInitialSolutions2Module(std::vector< int > solutions);
	/**
	@a YAMLInitialEquilibriumPhases2Module inserts data into the YAML document 
	for the PhreeqcRM method @a InitialEquilibriumPhases2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialEquilibriumPhases2Module transfers EQUILIBRIUM_PHASES definitions from the 
	InitialPhreeqc instance to the reaction-module workers.
	@a equilibrium_phases is a vector of EQUILIBRIUM_PHASES index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param equilibrium_phases Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialEquilibriumPhases2Module(std::vector< int > equilibrium_phases);
	/**
	@a YAMLInitialExchanges2Module inserts data into the YAML document for the 
	PhreeqcRM method @a InitialExchanges2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialExchanges2Module transfers EXCHANGE definitions from the
	InitialPhreeqc instance to the reaction-module workers.
	@a exchanges is a vector of EXCHANGE index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param exchanges Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialExchanges2Module(std::vector< int > exchanges);
	/**
	@a YAMLInitialSurfaces2Module inserts data into the YAML document for 
	the PhreeqcRM method @a InitialSurfaces2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialSurfaces2Module transfers SURFACE definitions from the
	InitialPhreeqc instance to the reaction-module workers.
	@a surfaces is a vector of SURFACE index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param surfaces Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialSurfaces2Module(std::vector< int > surfaces);
	/**
	@a YAMLInitialGasPhases2Module inserts data into the YAML document for 
	the PhreeqcRM method @a InitialGasPhases2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialGasPhases2Module transfers GAS_PHASE definitions from the
	InitialPhreeqc instance to the reaction-module workers.
	@a gas_phases is a vector of GAS_PHASE index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param gas_phases Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialGasPhases2Module(std::vector< int > gas_phases);
	/**
	@a YAMLInitialSolidSolutions2Module inserts data into the YAML document for 
	the PhreeqcRM method @a InitialSolidSolutions2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialSolidSolutions2Module transfers SOLID_SOLUTIONS definitions from the
	InitialPhreeqc instance to the reaction-module workers.
	@a solid_solutions is a vector of SOLID_SOLUTIONS index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param solid_solutions Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialSolidSolutions2Module(std::vector< int > solid_solutions);
	/**
	@a YAMLInitialKinetics2Module inserts data into the YAML document for the 
	PhreeqcRM method @a InitialKinetics2Module.
	When the YAML document is written to file it can be processed by the method
	InitializeYAML to initialize a PhreeqcRM instance.
	<p>
	@a InitialKinetics2Module transfers KINETICS definitions from the
	InitialPhreeqc instance to the reaction-module workers.
	@a kinetics is a vector of KINETICS index numbers that refer to
	definitions in the InitialPhreeqc instance.
	</p>
	@param kinetics Vector of index numbers that is dimensioned @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	 */
	void YAMLInitialKinetics2Module(std::vector< int > kinetics);
	/**
	@a YAMLInitialPhreeqc2Module inserts data into the YAML document 
	for the PhreeqcRM method @a InitialPhreeqc2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc instance to the reaction-module workers.
	@a Initial_conditions1 is used to select initial conditions, including solutions and reactants,
	for each cell of the model, without mixing.
	@a Initial_conditions1 is dimensioned 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model. 
	The dimension of 7 refers to solutions and reactants in the following order:
	(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
	(5) SOLID_SOLUTIONS, and (6) KINETICS.
	The definition initial_solution1[3*nxyz + 99] = 2, indicates that
	cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance.
	</p>
	@param initial_conditions1 Vector of solution and reactant index numbers that refer to
	definitions in the InitialPhreeqc instance.
	Size is 7 times @a nxyz. The order of definitions is given above.
	Negative values are ignored, resulting in no definition of that entity for that cell.
	@see                        @ref YAMLInitialPhreeqcCell2Module.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<int> ic1;
	ic1.resize(nxyz*7, -1);
	for (int i = 0; i < nxyz; i++)
	{
	  ic1[i] = 1;              // Solution 1
	  ic1[nxyz + i] = -1;      // Equilibrium phases none
	  ic1[2*nxyz + i] = 1;     // Exchange 1
	  ic1[3*nxyz + i] = -1;    // Surface none
	  ic1[4*nxyz + i] = -1;    // Gas phase none
	  ic1[5*nxyz + i] = -1;    // Solid solutions none
	  ic1[6*nxyz + i] = -1;    // Kinetics none
	}
	yrm.YAMLInitialPhreeqc2Module(ic1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1);
	/**
	@a YAMLInitialPhreeqc2Module inserts data into the YAML document 
	for the PhreeqcRM method @a InitialPhreeqc2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc instance to
	the reaction-module workers, possibly with mixing.
	In its simplest form, @a  initial_conditions1 is used to select initial conditions, including solutions and reactants,
	for each cell of the model, without mixing.
	@a Initial_conditions1 is dimensioned 7 times @a  nxyz, where @a  nxyz is the number of grid cells in the user's model. 
	The dimension of 7 refers to solutions and reactants in the following order:
	(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
	(5) SOLID_SOLUTIONS, and (6) KINETICS.
	The definition initial_solution1[3*nxyz + 99] = 2, indicates that
	cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
	(either by RunFile or RunString).
	</p>
	<p>
	It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing,
	@a initials_conditions2 contains numbers for a second entity that mixes with the entity defined in @a initial_conditions1.
	@a Fraction1 contains the mixing fraction for @a initial_conditions1,
	whereas (1 - @a fraction1) is the mixing fraction for @a initial_conditions2.
	The definitions initial_conditions1[3*nxyz + 99] = 2, initial_conditions2[3*nxyz + 99] = 3,
	fraction1[3*nxyz + 99] = 0.25 indicates that
	cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3,
	where the surface compositions have been defined in the InitialPhreeqc instance.
	If the user number in @a initial_conditions2 is negative, no mixing occurs.
	@param initial_conditions1 Vector of solution and reactant index numbers that refer to
	definitions in the InitialPhreeqc instance.
	Size is 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model.
	The order of definitions is given above.
	Negative values are ignored, resulting in no definition of that entity for that cell.
	</p>
	@param initial_conditions2  Vector of solution and reactant index numbers that refer to
	definitions in the InitialPhreeqc instance.
	Nonnegative values of @a initial_conditions2 result in mixing with the entities defined in @a initial_conditions1.
	Negative values result in no mixing.
	Size is 7 times @a nxyz. The order of definitions is given above.
	@param fraction1           Fraction of @a initial_conditions1 that mixes with (1 - @a fraction1)
	of @a initial_conditions2.
	Size is 7 times @a nxyz. The order of definitions is given above.
	@see                        @ref YAMLInitialPhreeqcCell2Module.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<int> ic1, ic2;
	ic1.resize(nxyz*7, -1);
	ic2.resize(nxyz*7, -1);
	std::vector<double> f1;
	f1.resize(nxyz*7, 1.0);
	for (int i = 0; i < nxyz; i++)
	{
	  ic1[i] = 1;              // Solution 1
	  ic1[nxyz + i] = -1;      // Equilibrium phases none
	  ic1[2*nxyz + i] = 1;     // Exchange 1
	  ic1[3*nxyz + i] = -1;    // Surface none
	  ic1[4*nxyz + i] = -1;    // Gas phase none
	  ic1[5*nxyz + i] = -1;    // Solid solutions none
	  ic1[6*nxyz + i] = -1;    // Kinetics none
	}
	yrm.YAMLInitialPhreeqc2Module(ic1, ic2, f1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
	/**
	@a YAMLInitialPhreeqcCell2Module inserts data into the YAML document for the 
	PhreeqcRM method @a InitialPhreeqcCell2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a InitialPhreeqcCell2Module uses a cell numbered @a n in the InitialPhreeqc instance to 
	populate a series of transport cells.
	All reactants with the number @a n are transferred along with the solution.
	If MIX @a n exists, it is used for the definition of the solution.
	If @a n is negative, @a n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
	All reactants for each cell in the list @a cell_numbers are removed before the cell
	definition is copied from the InitialPhreeqc instance to the workers.
	</p>
	@param n                  Number that refers to a solution or MIX and associated reactants in the InitialPhreeqc instance.
	@param cell_numbers       A vector of grid-cell numbers (user's grid-cell numbering system) that
	will be populated with cell @a n from the InitialPhreeqc instance.
	@see                      @ref YAMLInitialPhreeqc2Module.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<int> module_cells;
	module_cells.push_back(18);
	module_cells.push_back(19);
	yrm.YAMLInitialPhreeqcCell2Module(-1, module_cells);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLInitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
	/**
	@a YAMLLoadDatabase inserts data into the YAML document 
	for the PhreeqcRM method @a LoadDatabase.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a LoadDatabase loads a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. 
	All definitions
	of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
	</p>
	@param database         String containing the database name.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLLoadDatabase("phreeqc.dat");
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLLoadDatabase(std::string database);
	/**
	@a YAMLLogMessage inserts data into the YAML document for 
	the PhreeqcRM method @a LogMessage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a LogMessage prints a message to the log file.
	</p>
	@param str              String to be printed.
	@see                    @ref YAMLOutputMessage, @ref YAMLScreenMessage, @ref YAMLWarningMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLLogMessage("Finished section 1 of initialization");
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLLogMessage(std::string str);
	/**
	@a YAMLOpenFiles inserts data into the YAML document for the 
	PhreeqcRM method @a OpenFiles.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a OpenFiles opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
	based on the prefix defined by @ref YAMLSetFilePrefix.
	</p>
	@see                    @ref YAMLSetFilePrefix, @ref YAMLCloseFiles,
	@ref YAMLLogMessage, @ref YAMLOutputMessage, and @ref YAMLWarningMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetFilePrefix("Advect_cpp");
	yrm.YAMLOpenFiles();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLOpenFiles(void);
	/**
	@a YAMLOutputMessage inserts data into the YAML document for the 
	PhreeqcRM method @a OutputMessage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a OutputMessage prints a message to the output file.
	</p>
	@param str              String to be printed.
	@see                    @ref YAMLLogMessage, @ref YAMLScreenMessage, @ref YAMLWarningMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLOutputMessage("Finished section 1 of initialization");
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLOutputMessage(std::string str);
	/**
	@a YAMLRunCells inserts data into the YAML document for the 
	PhreeqcRM method @a RunCells.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a RunCells runs reactions for all cells in the reaction module.
	During initialization, RunCells can be used to equilibrate 
	each solution with all reactants in a cell while
	using a time step of zero (@ref YAMLSetTimeStep) to avoid kinetic reactions.
	Other properties that may need to be initialized before RunCells is invoked
	include porosity (@ref YAMLSetPorosity),
	saturation (@ref YAMLSetSaturationUser),
	temperature (@ref YAMLSetTemperature), and pressure (@ref YAMLSetPressure).
	</p>

	@see                    @ref YAMLSetPorosity,
	@ref YAMLSetPressure, @ref YAMLSetSaturationUser, @ref YAMLSetTemperature, @ref YAMLSetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetTimeStep(0.0);
	yrm.YAMLRunCells();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLRunCells(void);
	/**
	@a YAMLRunFile inserts data into the YAML document for the 
	PhreeqcRM method @a RunFile.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a RunFile runs a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
	the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
	files that modify the thermodynamic database should be run by all three sets of instances.
	Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
	be run by the workers. Files that contain initial conditions or boundary conditions should
	be run by the InitialPhreeqc instance.
	</p>
	@param workers          @a True, the workers will run the file; @a False, the workers will not run the file.
	@param initial_phreeqc  @a True, the InitialPhreeqc instance will run the file; @a False, the InitialPhreeqc will not run the file.
	@param utility          @a True, the Utility instance will run the file; @a False, the Utility instance will not run the file.
	@param chemistry_name   Name of the file to run.
	@see                    @ref YAMLRunString.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLRunFile(true, true, true, "advect.pqi");
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLRunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
	/**
	YAMLRunString inserts data into the YAML document for the 
	PhreeqcRM method @a RunString.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a RunString runs a PHREEQC input string. The first three arguments determine which
	IPhreeqc instances will run
	the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
	strings that modify the thermodynamic database should be run by all three sets of instances.
	Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
	be run by the workers. Strings that contain initial conditions or boundary conditions should
	be run by the InitialPhreeqc instance.
	</p>
	@param workers          @a True, the workers will run the string; @a False, the workers will not run the string.
	@param initial_phreeqc  @a True, the InitialPhreeqc instance will run the string; @a False, the InitialPhreeqc will not run the string.
	@param utility          @a True, the Utility instance will run the string; @a False, the Utility instance will not run the string.
	@param input_string     String containing PHREEQC input.
	@see                    @ref YAMLRunFile.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::string input = "DELETE; -all";
	yrm.YAMLRunString(true, false, true, input.c_str());
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLRunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
	/**
	@a YAMLScreenMessage inserts data into the YAML document for the 
	PhreeqcRM method @a ScreenMessage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a ScreenMessage prints a message to the screen.
	</p>
	@param str              String to be printed.
	@see                    @ref YAMLLogMessage, @ref YAMLOutputMessage, @ref YAMLWarningMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::ostringstream strm;
	strm << "Beginning to process YAML for initial conditions";
	yrm.YAMLScreenMessage(strm.str());
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLScreenMessage(std::string str);
	/**
	@a YAMLSetComponentH2O inserts data into the YAML document for the 
	PhreeqcRM method @a SetComponentH2O.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetComponentH2O selects whether to include H2O in the component list.
	The concentrations of H and O must be known
	accurately (8 to 10 significant digits) for the numerical method of
	PHREEQC to produce accurate pH and pe values.
	Because most of the H and O are in the water species,
	it may be more robust (require less accuracy in transport) to
	transport the excess H and O (the H and O not in water) and water.
	The default setting (@a true) is to include water, excess H, and excess O as components.
	A setting of @a false will include total H and total O as components.
	YAMLSetComponentH2O must be called before @ref YAMLFindComponents.
	</p>
	@param tf               @a True (default), excess H, excess O, and water are included in the component list;
	@a False, total H and O are included in the component list.
	@see                    @ref YAMLFindComponents.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetComponentH2O(false);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetComponentH2O(bool tf);
	/**
	@a YAMLSetConcentrations inserts data into the YAML document for the 
	PhreeqcRM method @a SetConcentrations.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	The only way to use this method is to have pre-calculated PHREEQC solution concentrations,
	which is not common. Concentrations are normally initialized
	with @ref YAMLInitialPhreeqc2Module or @ref YAMLInitialPhreeqcCell2Module.
	</p>
	@param c               Vector of component concentrations. Size of vector is @a ncomps times @a nxyz,
	where @a ncomps is the number of components as determined
	by FindComponents or GetComponentCount and
	@a nxyz is the number of grid cells in the user's model.
	@see                    @ref YAMLSetDensityUser, @ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume,
	@ref YAMLSetSaturationUser, @ref YAMLSetUnitsSolution.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetConcentrations(c);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetConcentrations(std::vector< double >& c);
	/**
	@a YAMLSetCurrentSelectedOutputUserNumber inserts data into the YAML document for the 
	PhreeqcRM method @a SetCurrentSelectedOutputUserNumber.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetCurrentSelectedOutputUserNumber selects the current selected output by user number.
	The user may define multiple SELECTED_OUTPUT
	data blocks for the workers. A user number is specified for each data block. The value of
	the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
	for selected-output operations.
	</p>
	@param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
	@see
	@ref YAMLSetNthSelectedOutput,
	@ref YAMLSetSelectedOutputOn.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetCurrentSelectedOutputUserNumber(n_user);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetCurrentSelectedOutputUserNumber(int n_user);
	/**
	@a YAMLSetDensityUser inserts data into the YAML document for the 
	PhreeqcRM method @a SetDensityUser.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetDensityUser sets the density for each reaction cell. These density values are used
	when converting from transported mass-fraction concentrations (@ref YAMLSetUnitsSolution) to
	produce per liter concentrations during a call to SetConcentrations and when converting from 
	reaction-cell concentrations to transport concentrations, if UseSolutionDensityVolume is set to 
	@a false.
	</p>
	@param density          Vector of densities. Size of vector is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model.
	@see
	@ref YAMLSetUnitsSolution, @ref YAMLUseSolutionDensityVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> initial_density;
	initial_density.resize(nxyz, 1.0);
	yrm.YAMLSetDensityUser(initial_density);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetDensityUser(std::vector< double > density);
	/**
	@a YAMLSetDumpFileName inserts data into the YAML document for the 
	PhreeqcRM method @a SetDumpFileName.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetDumpFileName	sets the name of the dump file. It is the name used by the method DumpModule.
	</p>
	@param dump_name        Name of dump file.
	@see                    @ref YAMLDumpModule.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetDumpFileName("Advect_cpp.dmp");
	bool dump_on = true;
	bool append = false;
	yrm.YAMLDumpModule(dump_on, append);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetDumpFileName(std::string dump_name);
	/**
	@a YAMLSetErrorHandlerMode inserts data into the YAML document for the 
	PhreeqcRM method @a SetErrorHandlerMode.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetErrorHandlerMode sets the action to be taken when the reaction module encounters an error.
	Options are 0, return to calling program with an error return code (default);
	1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
	2, attempt to exit gracefully.
	</p>
	@param mode             Error handling mode: 0, 1, or 2.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetErrorHandlerMode(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetErrorHandlerMode(int mode);
	/**
	@a YAMLSetErrorOn inserts data into the YAML document for the 
	PhreeqcRM method @a SetErrorOn.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetErrorOn sets the property that controls whether error messages are generated and displayed.
	Messages include PHREEQC "ERROR" messages, and
	any messages written with the method ErrorMessage.
	</p>
	@param tf  @a True, enable error messages; @a False, disable error messages. Default is true.
	@see                  @ref YAMLLogMessage, @ref YAMLOutputMessage, @ref YAMLScreenMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetErrorOn(true);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetErrorOn(bool tf);
	/**
	@a YAMLSetFilePrefix inserts data into the YAML document for the 
	PhreeqcRM method @a SetFilePrefix.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetFilePrefix sets the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
	These files are opened by the method OpenFiles.
	</p>
	@param prefix           Prefix used when opening the output and log files.
	@see                    @ref YAMLOpenFiles, @ref YAMLCloseFiles.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetFilePrefix("Advect_cpp");
	yrm.YAMLOpenFiles();
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetFilePrefix(std::string prefix);
	/**
	@a YAMLSetGasCompMoles inserts data into the YAML document for the 
	PhreeqcRM method @a SetGasCompMoles.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetGasCompMoles transfers moles of gas components from
	the vector given in the argument list (@a gas_moles) to each reaction cell.
	</p>
	@param  gas_moles               Vector of moles of gas components.
	Dimension of the vector is set to @a ngas_comps times @a nxyz,
	where, @a ngas_comps is the result of GetGasComponentsCount,
	and @a nxyz is the number of user grid cells.
	If the number of moles is set to a negative number, the gas component will
	not be defined for the GAS_PHASE of the reaction cell.

	@see
	@ref YAMLFindComponents,
	@ref YAMLSetGasPhaseVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> gas_moles;
	gas_moles.resize(nxyz*ngas, 0);
	yrm.YAMLSetGasCompMoles(gas_moles);
	</PRE>
	</CODE>
	@endhtmlonly
		*/
	void YAMLSetGasCompMoles(std::vector< double > gas_moles);
	/**
	@a YAMLSetGasPhaseVolume inserts data into the YAML document for the 
	PhreeqcRM method @a SetGasPhaseVolume.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetGasPhaseVolume transfers volumes of gas phases from
	the vector given in the argument list (@a gas_volume) to each reaction cell.
	The gas-phase volume affects the gas-component pressures calculated for fixed-volume
	gas phases. If a gas-phase volume is defined with this methood
	for a GAS_PHASE in a cell,
	the gas phase is forced to be a fixed-volume gas phase.
	</p>
	@param  gas_volume               Vector of volumes for each gas phase.
	Dimension of the vector is @a nxyz,
	where @a nxyz is the number of user grid cells.
	If the volume is set to a negative number for a cell, the gas-phase volume for that cell is
	not changed.

	@see
	@ref YAMLFindComponents,
	@ref YAMLSetGasCompMoles.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> gas_volume;
	gas_volume.resize(nxyz, 0.5);
	yrm.YAMLSetGasPhaseVolume(gas_volume);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLSetGasPhaseVolume(std::vector< double > gas_volume);
	/**
	@a YAMLSetGridCellCount Inserts data into the YAML document to define the number of cells in the 
	user's model.
	Once the YAML document is written, the number of model cells can be extracted
	with the method GetGridCellCountYAML. GetGridCellCountYAML is NOT a PhreeqcRM 
	method; it is a global method and must be used BEFORE the PhreeqcRM instance
	is created. SetGridCellCount will be ignored once the PhreeqcRM instance exists.
	@param n           Number of cells for the PhreeqcRM instance. The number of cells
	can be used in the creation of the PhreeqcRM instance. The PhreeqcRM constructor
	takes two arguments. GetGridCellCountYAML
	provides the value for the first argument. If the YAML file does not contain the
	node "SetGridCellCount:", GetGridCellCountYAML will return zero.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	int nxyz = 40;
	yrm.YAMLSetGridCellCount(nxyz);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLSetGridCellCount(int n);
	/**
	@a YAMLSetNthSelectedOutput inserts data into the YAML document for the 
	PhreeqcRM method @a SetNthSelectedOutput.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetNthSelectedOutput specifies the current selected output by sequence number. 
	The user may define multiple SELECTED_OUTPUT
	data blocks for the workers. A user number is specified for each data block, and the blocks are
	stored in user-number order. The value of
	the argument @a n selects the sequence number of the SELECTED_OUTPUT definition that will be used
	for selected-output operations.
	</p>
	@param n           Sequence number of the SELECTED_OUTPUT data block that is to be used.
	@see
	@ref YAMLSetCurrentSelectedOutputUserNumber,
	@ref YAMLSetSelectedOutputOn.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetCurrentSelectedOutput(isel);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLSetNthSelectedOutput(int n);
	/**
	@a YAMLSetPartitionUZSolids inserts data into the YAML document for the 
	PhreeqcRM method @a SetPartitionUZSolids.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetPartitionUZSolids sets the property for partitioning solids between 
	the saturated and unsaturated
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
	</p>
	@param tf       @a True, the fraction of solids and gases available for
	reaction is equal to the saturation;
	@a False (default), all solids and gases are reactive regardless of saturation.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetPartitionUZSolids(false);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetPartitionUZSolids(bool tf);
	/**
	@a YAMLSetPorosity inserts data into the YAML document for the 
	PhreeqcRM method @a SetPorosity.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetPorosity sets the porosity for each reaction cell.
	The volume of water in a reaction cell is the product of porosity, saturation
	(SetSaturationUser), and representative volume (SetRepresentativeVolume).
	@param por              Vector of porosities, unitless. Default is 0.1.
	Size of vector is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model.
	@see                    @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> por(nxyz, 0.2);
	yrm.YAMLSetPorosity(por);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetPorosity(std::vector< double > por);
	/**
	@a YAMLSetPressure inserts data into the YAML document for the 
	PhreeqcRM method @a SetPressure.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetPressure sets the pressure for each reaction cell. Pressure effects are 
	considered only in three of the
	databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
	</p>
	@param p                Vector of pressures, in atm. Size of vector is @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	@see                    @ref YAMLSetTemperature.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> pressure(nxyz, 2.0);
	yrm.YAMLSetPressure(pressure);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetPressure(std::vector< double > p);
	/**
	@a YAMLSetPrintChemistryMask inserts data into the YAML document for the 
	PhreeqcRM method @a SetPrintChemistryMask.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetPrintChemistryMask enables or disables detailed output for each reaction cell.
	Printing for a reaction cell will occur only when the
	printing is enabled with SetPrintChemistryOn and the @a cell_mask value is 1.
	</p>
	@param cell_mask        Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model. A value of 0 will
	disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
	@see                    @ref YAMLSetPrintChemistryOn.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<int> print_chemistry_mask;
	print_chemistry_mask.resize(nxyz, 0);
	for (int i = 0; i < nxyz/2; i++)
	{
	  print_chemistry_mask[i] = 1;
	}
	yrm.YAMLSetPrintChemistryMask(print_chemistry_mask);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetPrintChemistryMask(std::vector< int > cell_mask);
	/**
	@a YAMLSetPrintChemistryOn inserts data into the YAML document for the 
	PhreeqcRM method @a SetPrintChemistryOn.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetPrintChemistryOn
	sets the property that enables or disables printing detailed output from reaction calculations
	to the output file for a set of cells defined by SetPrintChemistryMask.
	The detailed output prints all of the output typical of a PHREEQC reaction calculation,
	which includes solution descriptions and the compositions of all other reactants.
	The output can be several hundred lines per cell, which can lead to a very
	large output file (prefix.chem.txt opened by the method OpenFiles).
	For the worker instances, the output can be limited to a set of cells
	(method SetPrintChemistryMask) and, in general, the
	amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
	(applied by using methods RunFile or RunString).
	Printing the detailed output for the workers is generally used only for debugging,
	and PhreeqcRM will run significantly faster
	when printing detailed output for the workers is disabled.
	</p>
	@param workers          @a True, enable detailed printing in the worker instances;
	@a False, disable detailed printing in the worker instances.
	@param initial_phreeqc  @a True, enable detailed printing in the InitialPhreeqc instance;
	@a False, disable detailed printing in the InitialPhreeqc instance.
	@param utility          @a True, enable detailed printing in the Utility instance;
	@a False, disable detailed printing in the Utility instance.
	@see                    @ref YAMLOpenFiles, @ref YAMLRunFile, @ref YAMLRunString, @ref YAMLSetPrintChemistryMask.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetPrintChemistryOn(false, true, false);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
	/**
	@p YAMLSetRebalanceByCell inserts data into the YAML document for the 
	PhreeqcRM method @a SetRebalanceByCell.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetRebalanceByCell
	sets the load-balancing algorithm.
	PhreeqcRM attempts to rebalance the load of each thread or process such that each
	thread or process takes the same amount of time to run its part of a RunCells
	calculation. Two algorithms are available; one uses individual times for each cell and
	accounts for cells that were not run because
	saturation was zero (default), and
	the other assigns an average time to all cells.
	The methods are similar, but limited testing indicates the default method performs better.
	</p>
	@param tf           @a True, indicates individual cell times are used in rebalancing (default);
	@a False, indicates average times are used in rebalancing.
	@see                    @ref YAMLSetRebalanceFraction.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetRebalanceByCell(true);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetRebalanceByCell(bool tf);
	/**
	@a YAMLSetRebalanceFraction inserts data into the YAML document for the 
	PhreeqcRM method @a SetRebalanceFraction.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetRebalanceFraction sets the fraction of cells that are transferred among threads or processes when rebalancing.
	PhreeqcRM attempts to rebalance the load of each thread or process such that each
	thread or process takes the same amount of time to run its part of a RunCells
	calculation. The rebalancing transfers cell calculations among threads or processes to
	try to achieve an optimum balance. @a SetRebalanceFraction
	adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
	determine the actual number of cell transfers. A value of zero eliminates
	load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
	distribution and avoid possible oscillations
	when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
	Default is 0.5.
	</p>
	@param f                Fraction from 0.0 to 1.0.
	@see                    @ref YAMLSetRebalanceByCell.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetRebalanceFraction(0.5);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetRebalanceFraction(double f);
	/**
	@a YAMLSetRepresentativeVolume inserts data into the YAML document for the 
	PhreeqcRM method @a SetRepresentativeVolume.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetRepresentativeVolume
	sets the representative volume of each reaction cell.
	By default the representative volume of each reaction cell is 1 liter.
	The volume of water in a reaction cell is determined by the product of the representative volume,
	the porosity (SetPorosity), and the saturation (SetSaturationUser).
	The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
	within a couple orders of magnitude of 1.0.
	Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes)
	may cause non-convergence of the numerical method.
	In these cases, a larger representative volume may help. Note
	that increasing the representative volume also increases
	the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
	and others), which are defined as moles per representative volume.
	@a SetRepresentativeVolume should be called before initial conditions 
	are defined for the reaction cells.
	</p>
	@param rv              Vector of representative volumes, in liters. Default is 1.0 liter.
	Size of array is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model.
	@see                    @ref YAMLSetPorosity, @ref YAMLSetSaturationUser.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> rv(nxyz, 2.0);
	yrm.YAMLSetRepresentativeVolume(rv);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetRepresentativeVolume(std::vector< double > rv);
	/**
	@a YAMLSetSaturationUser inserts data into the YAML document for the 
	PhreeqcRM method @a SetSaturationUser.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetSaturationUser
	sets the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
	The volume of water in a cell is the product of porosity (SetPorosity), saturation 
	(SetSaturationUser), and representative volume (SetRepresentativeVolume). As a result of 
	a reaction calculation, solution properties (density and volume) will change;
	the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume data to 
	calculate these changes. 	The methods GetDensityCalculated,
	GetSolutionVolume, and GetSaturationCalculated can be used to account
	for these changes in the succeeding transport calculation.
	</p>
	@param sat              Vector of saturations, unitless. Default 1.0. Size of vector is @a nxyz,
	where @a nxyz is the number of grid cells in the user's model.
	@see
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> sat(nxyz, 1.0);
	yrm.YAMLSetSaturationUser(sat);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetSaturationUser(std::vector< double > sat);
	/**
	@a YAMLSetScreenOn inserts data into the YAML document for the 
	PhreeqcRM method @a SetScreenOn.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetScreenOn
	sets the property that controls whether messages are written to the screen.
	Messages include information about rebalancing during RunCells, and
	any messages written with ScreenMessage.
	</p>
	@param tf  @a True, enable screen messages; @a False, disable screen messages. Default is true.
	@see                    @ref YAMLRunCells, @ref YAMLScreenMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetScreenOn(true);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetScreenOn(bool tf);
	/**
	@a YAMLSetSelectedOutputOn inserts data into the YAML document for the 
	PhreeqcRM method @a SetSelectedOutputOn.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetSelectedOutputOn
	sets the property that controls whether selected-output results are available to be retrieved
	with GetSelectedOutput. @a True indicates that selected-output results
	will be accumulated during RunCells and can be retrieved with GetSelectedOutput;
	@a False indicates that selected-output results will not
	be accumulated during RunCells.
	</p>
	@param tf  @a True, enable selected output; @a False, disable selected output.
	@see
	@ref YAMLSetCurrentSelectedOutputUserNumber,
	@ref YAMLSetNthSelectedOutput,
	@ref YAMLSetSelectedOutputOn.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetSelectedOutputOn(true);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetSelectedOutputOn(bool tf);
	/**
	@a Inserts data into the YAML document for the 
	PhreeqcRM method @a SetSpeciesSaveOn.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetSpeciesSaveOn
	sets the value of the species-save property.
	This method enables or disables use of PhreeqcRM with multicomponent-diffusion transport calculations.
	By default, concentrations of aqueous species are not saved.
	Setting the species-save property to @a true allows
	aqueous species concentrations to be retrieved
	with GetSpeciesConcentrations, and solution compositions to be set with
	SpeciesConcentrations2Module.
	SetSpeciesSaveOn must be called before calls to FindComponents.
	</p>
	@param save_on          @a True indicates species concentrations are saved;
	@a False indicates species concentrations are not saved.
	@see                    @ref YAMLFindComponents,
	@ref YAMLSpeciesConcentrations2Module.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetSpeciesSaveOn(true);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetSpeciesSaveOn(bool save_on);
	/**
	@a YAMLSetTemperature inserts data into the YAML document for the 
	PhreeqcRM method @a SetTemperature.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetTemperature
	sets the temperature for each reaction cell. If SetTemperature is not called,
	worker solutions will have temperatures as defined by initial conditions
	(InitialPhreeqc2Module and InitialPhreeqcCell2Module).
	</p>
	@param t                Vector of temperatures, in degrees C.
	Size of vector is @a nxyz, where @a nxyz is the number
	of grid cells in the user's model.
	@see                    @ref YAMLInitialPhreeqc2Module,
	@ref YAMLInitialPhreeqcCell2Module, @ref YAMLSetPressure.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector<double> temperature(nxyz, 20.0);
	yrm.YAMLSetTemperature(temperature);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetTemperature(std::vector< double > t);
	/**
	@a YAMLSetTime inserts data into the YAML document for the 
	PhreeqcRM method @a SetTime.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetTime
	sets current simulation time for the reaction module.
	</p>
	@param time             Current simulation time, in seconds.
	@see                    @ref YAMLSetTimeStep, @ref YAMLSetTimeConversion.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetTime(0.0);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetTime(double time);
	/**
	@a YAMLSetTimeConversion inserts data into the YAML document for the 
	PhreeqcRM method @a SetTimeConversion.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetTimeConversion sets a factor to convert from seconds to user time units. 
	Factor times seconds 
	produces user time units, which are used in some PhreeqcRM printing.
	</p>
	@param conv_factor      Factor to convert seconds to user time units.
	@see                    @ref YAMLSetTime, @ref YAMLSetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	double time_conversion = 1.0 / 86400;
	yrm.YAMLSetTimeConversion(time_conversion);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetTimeConversion(double conv_factor);
	/**
	@a YAMLSetTimeStep inserts data into the YAML document for the 
	PhreeqcRM method @a SetTimeStep.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetTimeStep
	sets current time step for the reaction module. This is the length
	of time over which kinetic reactions are integrated.
	</p>
	@param time_step        Time step, in seconds.
	@see                    @ref YAMLSetTime, @ref YAMLSetTimeConversion.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	time_step = 86400.;
	yrm.YAMLSetTimeStep(time_step);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetTimeStep(double time_step);
	/**
	@a YAMLSetUnitsExchange inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsExchange.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsExchange
	sets input units for exchangers.
	In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
	SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
	</p>
	<p>
	If a single EXCHANGE definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of exchangers will be the same regardless of porosity.
	For option 1, the number of moles of exchangers will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.
	</p>
	@param option           Units option for exchangers: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsExchange(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsExchange(int option);
	/**
	@a YAMLSetUnitsGasPhase inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsGasPhase.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsGasPhase
	sets input units for gas phases.
	In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
	@a SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
	</p>
	<p>
	If a single GAS_PHASE definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of a gas component will be the same regardless of porosity.
	For option 1, the number of moles of a gas component will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.
	</p>
	@param option           Units option for gas phases: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsGasPhase(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsGasPhase(int option);
	/**
	@a YAMLSetUnitsKinetics inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsKinetics.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsKinetics
	sets input units for kinetic reactants.
	</p>
	<p>
	In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
	@a SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
	</p>
	<p>
	If a single KINETICS definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of kinetic reactants will be the same regardless of porosity.
	For option 1, the number of moles of kinetic reactants will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.
	</p>
	<p>
	Note that the volume of water in a cell in the reaction module is equal to the product of
	porosity (SetPorosity), the saturation (SetSaturationUser), and representative volume 
	(SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the 
	RATES 	definitions for homogeneous (aqueous) kinetic reactions to account for the current 
	volume of water, often by calculating the rate of reaction per liter of water and multiplying 
	by the volume of water (Basic function SOLN_VOL).
	</p>
	<p>
	Rates that depend on surface area of solids, are not dependent
	on the volume of water. However, it is important to get the correct surface area for the kinetic
	reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant)
	can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of
	reactant (Basic function M) in RATES to obtain the surface area.
	</p>
	@param option           Units option for kinetic reactants: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsKinetics(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsKinetics(int option);
	/**
	@a YAMLSetUnitsPPassemblage inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsPPassemblage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsPPassemblage
	sets input units for pure phase assemblages (equilibrium phases).
	In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
	@a SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
	</p>
	<p>
	If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of a mineral will be the same regardless of porosity.
	For option 1, the number of moles of a mineral will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.
	</p>
	@param option           Units option for equilibrium phases: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsPPassemblage(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsPPassemblage(int option);
	/**
	@a YAMLSetUnitsSolution inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsSolution.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsSolution
	sets solution concentration units used by the transport model.
	Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
	PHREEQC defines solutions by the number of moles of each
	element in the solution.
	</p>
	<p>
	To convert from mg/L to moles
	of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
	multiplied by the solution volume,
	which is the product of porosity (SetPorosity), saturation (SetSaturationUser),
	and representative volume (SetRepresentativeVolume).
	To convert from mol/L to moles
	of element in the representative volume of a reaction cell, mol/L is
	multiplied by the solution volume.
	To convert from mass fraction to moles
	of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
	(SetDensityUser) and
	multiplied by the solution volume.
	</p>
	<p>
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
	calculated by PHREEQC, or (2) the volume of solution is the product of porosity (SetPorosity),
	saturation (SetSaturationUser), and representative volume (SetRepresentativeVolume),
	and the mass of solution is volume times density as defined by SetDensityUser.
	Which option is used is determined by UseSolutionDensityVolume.
	</p>
	@param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
	@see                    @ref YAMLSetDensityUser, @ref YAMLSetPorosity, 
	@ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser,
	@ref YAMLUseSolutionDensityVolume.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsSolution(2);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsSolution(int option);
	/**
	@a YAMLSetUnitsSSassemblage inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsSSassemblage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsSSassemblage
	sets input units for solid-solution assemblages.
	In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
	@a SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
	</p>
	<p>
	If a single SOLID_SOLUTION definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of a solid-solution component will be the same regardless of porosity.
	For option 1, the number of moles of a solid-solution component will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.
	</p>
	@param option           Units option for solid solutions: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsSSassemblage(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsSSassemblage(int option);
	/**
	@a YAMLSetUnitsSurface inserts data into the YAML document for the 
	PhreeqcRM method @a SetUnitsSurface.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SetUnitsSurface
	sets input units for surfaces.
	In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
	@a SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
	is calculated from the input value (@a Mp).
	</p>
	<p>
	Options are
	0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
	1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (SetPorosity); or
	2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
	</p>
	<p>
	If a single SURFACE definition is used for cells with different initial porosity,
	   the three options scale quite differently.
	For option 0, the number of moles of surface sites will be the same regardless of porosity.
	For option 1, the number of moles of surface sites will vary directly with porosity and inversely with rock volume.
	For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.
	</p>
	@param option           Units option for surfaces: 0, 1, or 2.
	@see                    @ref YAMLInitialPhreeqc2Module, @ref YAMLInitialPhreeqcCell2Module,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSetUnitsSurface(1);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSetUnitsSurface(int option);
	/**
	@a YAMLSpeciesConcentrations2Module inserts data into the YAML document for the 
	PhreeqcRM method @a SpeciesConcentrations2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a SpeciesConcentrations2Module
	sets solution concentrations in the reaction cells
	based on the vector of aqueous species concentrations (@a species_conc).
	This method is intended for use with multicomponent-diffusion transport calculations,
	and SetSpeciesSaveOn must be set to @a true.
	The list of aqueous species is determined by FindComponents and includes all
	aqueous species that can be made from the set of components.
	The method determines the total concentration of a component
	by summing the molarities of the individual species times the stoichiometric
	coefficient of the element in each species.
	Solution compositions in the reaction cells are updated with these component concentrations.
	Usually, accurate concentrations will not be known to use YAMLSpeciesConcentrations2Module during
	initialization.
	</p>
	@param species_conc     Vector of aqueous species concentrations. Dimension of the array is @a nspecies times @a nxyz,
	where  @a nspecies is the number of aqueous species,
	and @a nxyz is the number of user grid cells.
	Concentrations are moles per liter.
	@see                    @ref YAMLFindComponents,
	@ref YAMLSetSpeciesSaveOn.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLSpeciesConcentrations2Module(c);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLSpeciesConcentrations2Module(std::vector< double > species_conc);
	/**
	@a YAMLStateSave inserts data into the YAML document for the 
	PhreeqcRM method @a StateSave.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a StateSave
	saves the state of the chemistry in all model cells, including SOLUTIONs,
	EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
	Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
	will be saved for each cell, if they have been defined in the worker IPhreeqc instances.
	The distribution of cells among the workers and the chemistry of fully or partially
	unsaturated cells are also saved. The state is saved in memory; use DumpModule to save the state
	to file. PhreeqcRM can be reset to this state by using StateApply.
	A state is identified by an integer, and multiple states can be saved.
	</p>
	@param istate     Integer identifying the state that is saved.
	@see                    @ref YAMLDumpModule,
	@ref YAMLStateApply, and
	@ref YAMLStateDelete.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLStateSave(1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLStateSave(int istate);
	/**
	@a YAMLStateApply inserts data into the YAML document for the 
	PhreeqcRM method @a StateApply.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a StateApply
	resets the state of the module to a state previously saved with StateSave.
	The chemistry of all model cells are reset, including SOLUTIONs,
	EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
	MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs
	will be reset for each cell, if they were defined in the worker IPhreeqc instances
	at the time the state was saved.
	The distribution of cells among the workers and the chemistry of fully or partially
	unsaturated cells are also reset to the saved state.
	The state to be applied is identified by an integer.
	</p>
	@param istate     Integer identifying the state that is to be applied.
	@see                    @ref YAMLStateSave and
	@ref YAMLStateDelete.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLStateApply(1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLStateApply(int istate);
	/**
	@a YAMLStateDelete inserts data into the YAML document for the 
	PhreeqcRM method @a StateDelete.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a StateDelete
	deletes a state previously saved with StateSave.
	</p>
	@param istate     Integer identifying the state that is to be deleted.
	@see                    @ref YAMLStateSave and
	@ref YAMLStateApply.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLStateDelete(1);
	</PRE>
	</CODE>
	@endhtmlonly
	*/
	void YAMLStateDelete(int istate);
	/**
	@a YAMLUseSolutionDensityVolume inserts data into the YAML document for the 
	PhreeqcRM method @a UseSolutionDensityVolume.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a UseSolutionDensityVolume
	determines the volume and density to use when converting from the reaction-cell concentrations
	to transport concentrations (GetConcentrations).
	Two options are available to convert concentration units:
	(1) the density and solution volume calculated by PHREEQC are used, or
	(2) the specified density (SetDensityUser)
	and solution volume are determined by the product of
	saturation (SetSaturationUser), porosity (SetPorosity),
	and representative volume (SetRepresentativeVolume).
	Transport models that consider density-dependent flow will probably use the
	PHREEQC-calculated density and solution volume (default),
	whereas transport models that assume constant-density flow will probably use
	specified values of density and solution volume.
	Only the following databases distributed with PhreeqcRM have molar-volume information
	needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
	Density is only used when converting to or from transport units of mass fraction.
	</p>
	@param tf          @a True indicates that the solution density and volume as
	calculated by PHREEQC will be used to calculate concentrations.
	@a False indicates that the solution density set by SetDensityUser and the volume determined by the
	product of  SetSaturationUser, SetPorosity, and SetRepresentativeVolume,
	will be used to calculate concentrations retrieved by GetConcentrations.
	@see                    @ref YAMLSetDensityUser,
	@ref YAMLSetPorosity, @ref YAMLSetRepresentativeVolume, @ref YAMLSetSaturationUser.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLUseSolutionDensityVolume(false);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLUseSolutionDensityVolume(bool tf);
	/**
	@a YAMLWarningMessage inserts data into the YAML document for the 
	PhreeqcRM method @a WarningMessage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	<p>
	@a WarningMessage
	prints a warning message to the screen and the log file.
	</p>
	@param warnstr          String to be printed.
	@see                    @ref YAMLOpenFiles, @ref YAMLLogMessage,
	@ref YAMLOutputMessage, @ref YAMLScreenMessage.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	yrm.YAMLWarningMessage("Need to check these definitions.");
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLWarningMessage(std::string warnstr);
	// data

};
#include "bmi.hxx"
/**
 * @class YAMLPhreeqcRMLib
 *
 * @brief Class to implement multiple instances of YAMLPhreeqcRM
 */
class YAMLPhreeqcRMLib
{
public:
	//static void CleanupYAMLPhreeqcRMInstances(void);
	static int CreateYAMLPhreeqcRM(void);
	static IRM_RESULT DestroyYAMLPhreeqcRM(int n);
	static YAMLPhreeqcRM* GetInstance(int n);
};

#endif // INC_YAMLPHREEQCRM_H
#endif
