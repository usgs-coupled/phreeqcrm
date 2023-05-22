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

void WriteYAMLFile_cpp_test(void)
{
	YAMLPhreeqcRM yrm;
	int nxyz = 20;
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
	yrm.YAMLSetFilePrefix("AdvectBMI_cpp_test");
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
	std::vector<int> print_chemistry_mask(nxyz, 1);
	yrm.YAMLSetPrintChemistryMask(print_chemistry_mask);

	// Set printing of chemistry file
	yrm.YAMLSetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
	// Load database
	yrm.YAMLLoadDatabase("phreeqc.dat");
	yrm.YAMLRunFile(true, true, true, "all_reactants.pqi");

	// Clear contents of workers and utility
	std::string input = "DELETE; -all";
	yrm.YAMLRunString(true, false, true, input.c_str());
	// Define additional output variables
	//yrm.YAMLAddOutputVars("AddOutputVars", "true"); // default
	// AddOutputVars 
	// SolutionProperties 
	// SolutionTotalMolalities 
	// ExchangeMolalities 
	// SurfaceMolalities 
	// EquilibriumPhases 
	// Gases 
	// KineticReactants 
	// SolidSolutions
	// CalculateValues 
	// SolutionActivities 
	// SolutionMolalities 
	// SaturationIndices 

	// Determine number of components to transport
	yrm.YAMLFindComponents();

	// set arrays of initial conditions
	std::vector<int> ic(nxyz, 1);
	yrm.YAMLInitialSolutions2Module(ic);
	yrm.YAMLInitialEquilibriumPhases2Module(ic);
	yrm.YAMLInitialExchanges2Module(ic);
	yrm.YAMLInitialGasPhases2Module(ic);
	yrm.YAMLInitialKinetics2Module(ic);
	yrm.YAMLInitialSolidSolutions2Module(ic);
	yrm.YAMLInitialSurfaces2Module(ic);

	// Initial equilibration of cells
	double time_step = 0.0;  // no kinetics
	yrm.YAMLSetTimeStep(time_step);
	double time = 0.0;
	yrm.YAMLSetTime(time);
	yrm.YAMLRunCells();
	yrm.YAMLSetTimeStep(86400);

	// Write YAML file
	std::string YAML_filename = "AdvectBMI_cpp_test.yaml";
	yrm.WriteYAMLDoc(YAML_filename);
	yrm.Clear();
};
#endif