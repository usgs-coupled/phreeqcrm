#include <map>
#include <string>
#include "yaml-cpp/yaml.h"
#pragma once
class YAMLPhreeqcRM
{
private:
	YAML::Node YAML_doc;
public:
	YAMLPhreeqcRM();
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
	void clear();
	/**
	Returns a constant reference to the YAML document to be able to write
	the document to file.
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
	const YAML::Node& GetYAML_doc() { return this->YAML_doc; };
	// methods
/**
Inserts data into the YAML document for the PhreeqcRM method CloseFiles.
When the YAML document is written to file it can be processed by the method InitializeYAML to
initialize a PhreeqcRM instance.
@par
CloseFiles closes the output and log files.
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
	Inserts data into the YAML document for the PhreeqcRM method CreateMapping.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	CreateMapping
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
	Inserts data into the YAML document for the PhreeqcRM method DumpModule.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	DumpModule writes the contents of all workers to file in _RAW formats (see appendix of PHREEQC version 3 manual),
	including SOLUTIONs and all reactants.
	@param dump_on          Signal for writing the dump file, true or false.
	@param append           Signal to append to the contents of the dump file, true or false.
	@see                    @ref YAMLSetDumpFileName.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	bool dump_on = true;
	bool append = false;
	yrm.YAMLSetDumpFileName("advection_cpp.dmp");
	yrm.YAMLDumpModule(dump_on, append);
	</PRE>
	</CODE>
	@endhtmlonly
	 */
	void YAMLDumpModule();
	/**
	Inserts data into the YAML document for the PhreeqcRM method FindComponents.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	FindComponents accumulates a list of elements. Elements are those that have been
	defined in a solution or any other reactant
	(EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
	This method can be called multiple times and the list that is created is cummulative.
	The list is the set of components that needs to be transported. By default the list
	includes water, excess H and excess O (the H and O not contained in water);
	alternatively, the list may be set to contain total H and total O (@ref YAMLSetComponentH2O),
	which requires transport results to be accurate to eight or nine significant digits.
	If multicomponent diffusion (MCD) is to be modeled,
	there is a capability to retrieve aqueous species concentrations
	(@ref YAMLGetSpeciesConcentrations) and to set new solution concentrations after
	MCD by using individual species concentrations
	(@ref YAMLSpeciesConcentrations2Module).
	To use these methods, the save-species property needs to be turned on (@ref YAMLSetSpeciesSaveOn).
	If the save-species property is on, FindComponents will generate
	a list of aqueous species,
	their diffusion coefficients at 25 C,
	and their charge.
	@ref YAMLSetComponentH2O,
	@ref YAMLSetSpeciesSaveOn,
	@ref YAMLSpeciesConcentrations2Module.
	@par The FindComponents method also generates lists of reactants--equilibrium phases,
	exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
	The lists are cumulative, including all reactants that were
	defined in the initial phreeqc instance at any time FindComponents was called.
	In addition, a list of phases is generated for which saturation indices may be calculated from the
	cumulative list of components.
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
	Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc instance to the reaction-module workers.
	@a Initial_conditions1 is used to select initial conditions, including solutions and reactants,
	for each cell of the model, without mixing.
	@a Initial_conditions1 is dimensioned 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model
	(@ref GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
	(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
	(5) SOLID_SOLUTIONS, and (6) KINETICS.
	The definition initial_solution1[3*nxyz + 99] = 2, indicates that
	cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance.
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
	Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	InitialPhreeqc2Module transfers solutions and reactants from the InitialPhreeqc instance to
	the reaction-module workers, possibly with mixing.
	In its simplest form, @a  initial_conditions1 is used to select initial conditions, including solutions and reactants,
	for each cell of the model, without mixing.
	@a Initial_conditions1 is dimensioned 7 times @a  nxyz, where @a  nxyz is the number of grid cells in the user's model
	(@ref GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
	(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
	(5) SOLID_SOLUTIONS, and (6) KINETICS.
	The definition initial_solution1[3*nxyz + 99] = 2, indicates that
	cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
	(either by @ref RunFile or @ref RunString).
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
	@param initial_conditions1 Vector of solution and reactant index numbers that refer to
	definitions in the InitialPhreeqc instance.
	Size is 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model.
	The order of definitions is given above.
	Negative values are ignored, resulting in no definition of that entity for that cell.
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
	Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqcCell2Module.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	InitialPhreeqcCell2Module uses a cell numbered @a n in the InitialPhreeqc instance to populate a series of transport cells.
	All reactants with the number @a n are transferred along with the solution.
	If MIX @a n exists, it is used for the definition of the solution.
	If @a n is negative, @a n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
	All reactants for each cell in the list @a cell_numbers are removed before the cell
	definition is copied from the InitialPhreeqc instance to the workers.
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
	Inserts data into the YAML document for the PhreeqcRM method LoadDatabase.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	LoadDatabase loads a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
	of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
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
	Inserts data into the YAML document for the PhreeqcRM method OpenFiles.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	OpenFiles opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
	based on the prefix defined by @ref YAMLSetFilePrefix.
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
	Inserts data into the YAML document for the PhreeqcRM method OutputMessage.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	OutputMessage prints a message to the output file.
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
	@par MPI:
	Called by root.
	 */
	void YAMLOutputMessage(std::string str);
	/**
	Inserts data into the YAML document for the PhreeqcRM method RunCells.
	When the YAML document is written to file it can be processed by the method InitializeYAML to
	initialize a PhreeqcRM instance.
	@par
	RunCells runs a reactions for all cells in the reaction module.
	During initialization, RunCells can be used to equilibrate each solution with all reacts in a cell while
	using a time step of zero (@ref YAMLSetTimeStep) to avoid kinetic reactions.
	Other properties that may need to be initialized before RunCells is invoked
	include porosity (@ref YAMLSetPorosity),
	saturation (@ref YAMLSetSaturation),
	temperature (@ref YAMLSetTemperature), and pressure (@ref YAMLSetPressure).

	@see                    @ref YAMLSetPorosity,
	@ref YAMLSetPressure, @ref YAMLSetSaturation, @ref YAMLSetTemperature, @ref YAMLSetTimeStep.
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
	void YAMLRunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
	void YAMLRunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
	void YAMLScreenMessage(std::string str);
	void YAMLSetComponentH2O(bool tf);
	void YAMLSetConcentrations(std::vector< double > c);
	void YAMLSetCurrentSelectedOutputUserNumber(int n_user);
	void YAMLSetDensity(std::vector< double > density);
	void YAMLSetDumpFileName(std::string dump_name);
	void YAMLSetErrorHandlerMode(int mode);
	void YAMLSetErrorOn(bool tf);
	void YAMLSetFilePrefix(std::string prefix);
	void YAMLSetGasCompMoles(std::vector< double > gas_moles);
	void YAMLSetGasPhaseVolume(std::vector< double > gas_volume);
	void YAMLSetPartitionUZSolids(bool tf);
	void YAMLSetPorosity(std::vector< double > por);
	const YAML::Node& GetYAMLDoc() { return this->YAML_doc; }
	void YAMLSetPressure(std::vector< double > p);
	void YAMLSetPrintChemistryMask(std::vector< int > cell_mask);
	void YAMLSetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
	void YAMLSetRebalanceByCell(bool tf);
	void YAMLSetRebalanceFraction(double f);
	void YAMLSetRepresentativeVolume(std::vector< double > rv);
	void YAMLSetSaturation(std::vector< double > sat);
	void YAMLSetScreenOn(bool tf);
	void YAMLSetSelectedOutputOn(bool tf);
	void YAMLSetSpeciesSaveOn(bool save_on);
	void YAMLSetTemperature(std::vector< double > t);
	void YAMLSetTime(double time);
	void YAMLSetTimeConversion(double conv_factor);
	void YAMLSetTimeStep(double time_step);
	void YAMLSetUnitsExchange(int option);
	void YAMLSetUnitsGasPhase(int option);
	void YAMLSetUnitsKinetics(int option);
	void YAMLSetUnitsPPassemblage(int option);
	void YAMLSetUnitsSolution(int option);
	void YAMLSetUnitsSSassemblage(int option);
	void YAMLSetUnitsSurface(int option);
	void YAMLSpeciesConcentrations2Module(std::vector< double > species_conc);
	void YAMLStateSave(int istate);
	void YAMLStateApply(int istate);
	void YAMLStateDelete(int istate);
	void YAMLUseSolutionDensityVolume(bool tf);
	void YAMLWarningMessage(std::string warnstr);
	// data
private:
	std::map<std::string, int> method_map;
};
