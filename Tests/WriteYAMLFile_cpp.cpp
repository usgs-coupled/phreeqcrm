#ifdef USE_YAML
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
#include <algorithm>
#include "yaml-cpp/yaml.h"
#include "YAMLPhreeqcRM.h"

void WriteYAMLFile_cpp(void)
{
	YAMLPhreeqcRM yrm;
	int nxyz = 40;
	// Set GridCellCount
	yrm.YAMLSetGridCellCount(nxyz);

	int nthreads = 3;
	// Set ThreadCount
	yrm.YAMLThreadCount(nthreads);

	// Set some properties
	yrm.YAMLSetErrorHandlerMode(1);
	yrm.YAMLSetComponentH2O(false);
	yrm.YAMLSetRebalanceFraction(0.5);
	yrm.YAMLSetRebalanceByCell(true);
	yrm.YAMLUseSolutionDensityVolume(false);
	yrm.YAMLSetPartitionUZSolids(false);
	// Open files
	yrm.YAMLSetFilePrefix("AdvectBMI_cpp");
	yrm.YAMLOpenFiles();

	// Set concentration units
	yrm.YAMLSetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
	yrm.YAMLSetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	yrm.YAMLSetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	yrm.YAMLSetUnitsSurface(1);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	yrm.YAMLSetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	yrm.YAMLSetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	yrm.YAMLSetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

	// Set conversion from seconds to user units (days) Only affects one print statement
	double time_conversion = 1.0 / 86400;
	yrm.YAMLSetTimeConversion(time_conversion);
	// Set representative volume
	std::vector<double> rv(nxyz, 1.0);
	yrm.YAMLSetRepresentativeVolume(rv);
	// Set density
	std::vector<double> density(nxyz, 1.0);
	yrm.YAMLSetDensityUser(density);
	// Set initial porosity
	std::vector<double> por(nxyz, 0.2);
	yrm.YAMLSetPorosity(por);
	// Set initial saturation
	std::vector<double> sat(nxyz, 1.0);
	yrm.YAMLSetSaturationUser(sat);
	// Set cells to print chemistry when print chemistry is turned on
	std::vector<int> print_chemistry_mask(nxyz, 0);
	for (int i = 0; i < nxyz / 2; i++)
	{
		print_chemistry_mask[i] = 1;
	}
	yrm.YAMLSetPrintChemistryMask(print_chemistry_mask);

	// Demonstation of mapping, two equivalent rows by symmetry
	std::vector<int> grid2chem;
	grid2chem.resize(nxyz, -1);
	for (int i = 0; i < nxyz / 2; i++)
	{
		grid2chem[i] = i;
		grid2chem[i + nxyz / 2] = i;
	}
	yrm.YAMLCreateMapping(grid2chem);

	// Set printing of chemistry file
	yrm.YAMLSetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
	// Load database
	yrm.YAMLLoadDatabase("phreeqc.dat");

	// Run file to define solutions and reactants for initial conditions, selected output
	bool workers = true;             // Worker instances do the reaction calculations for transport
	bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
	bool utility = true;             // Utility instance is available for processing
	yrm.YAMLRunFile(workers, initial_phreeqc, utility, "advect.pqi");

	// Clear contents of workers and utility
	initial_phreeqc = false;
	std::string input = "DELETE; -all";
	yrm.YAMLRunString(workers, initial_phreeqc, utility, input.c_str());
	// Define additional output variables
	yrm.YAMLAddOutputVars("AddOutputVars", "true");
	// Determine number of components to transport
	yrm.YAMLFindComponents();
	// set array of initial conditions
	std::vector<int> ic1, ic2;
	ic1.resize(nxyz * 7, -1);
	ic2.resize(nxyz * 7, -1);
	std::vector<double> f1;
	f1.resize(nxyz * 7, 1.0);
	for (int i = 0; i < nxyz; i++)
	{
		ic1[i] = 1;              // Solution 1
		ic1[nxyz + i] = -1;      // Equilibrium phases none
		ic1[2 * nxyz + i] = 1;     // Exchange 1
		ic1[3 * nxyz + i] = -1;    // Surface none
		ic1[4 * nxyz + i] = -1;    // Gas phase none
		ic1[5 * nxyz + i] = -1;    // Solid solutions none
		ic1[6 * nxyz + i] = -1;    // Kinetics none
	}
	yrm.YAMLInitialPhreeqc2Module(ic1, ic2, f1);
	// No mixing is defined, so the following is equivalent
	//yrm.YAMLInitialPhreeqc2Module(ic1);

	// alternative for setting initial conditions
	// cell number in first argument (-1 indicates last solution, 40 in this case)
	// in advect.pqi and any reactants with the same number--
	// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
	// will be written to cells 18 and 19 (0 based)
	std::vector<int> module_cells;
	module_cells.push_back(18);
	module_cells.push_back(19);
	yrm.YAMLInitialPhreeqcCell2Module(-1, module_cells);
	// Initial equilibration of cells
	double time_step = 0.0;  // no kinetics
	yrm.YAMLSetTimeStep(time_step);
	double time = 0.0;
	yrm.YAMLSetTime(time);
	yrm.YAMLRunCells();
	yrm.YAMLSetTimeStep(86400);

	// Write YAML file
	std::string YAML_filename = "AdvectBMI_cpp.yaml";
	yrm.WriteYAMLDoc(YAML_filename);
	yrm.Clear();
};
#endif