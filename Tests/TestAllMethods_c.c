#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
#include "BMI_interface_C.h"
#include "IPhreeqc.h"
//#include "yaml-cpp/yaml.h"
//#include "YAMLPhreeqcRM.h"
#define MAX_LENGTH 256
void TestAllMethods_c()
{
	//YAMLPhreeqcRM yrm;
	IRM_RESULT status;
	int nxyz = 40;
	int nthreads = 3;
	int id;
	int* i_ptr;
	// Set GridCellCount
	//yrm.YAMLSetGridCellCount(nxyz);
	//yrm.YAMLThreadCount(3);
	char YAML_filename[MAX_LENGTH];
	strcpy_safe(YAML_filename, "TestAllMethods_cpp.yaml", MAX_LENGTH);

#ifdef USE_MPI
	// MPI
	int mpi_myself;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}

	id = RM_BMI_Create(nxyz, MPI_COMM_WORLD);
	if (mpi_myself > 0)
	{
		RM_BMI_MpiWorker();
		RM_BMI_Finalize();
		return;
	}
#else
	// OpenMP or serial
	id = RM_BMI_Create(nxyz, nthreads);
#endif
	// Use YAML file to initialize
	//RM_BMI_Initialize(id, YAML_filename);   // void function
	//RM_BMI_InitializeYAML(id, YAML_filename);
	//fprintf(stderr, "Initialize\n");
	//
	// Use all BMIPhreeqcRM methods roughly in order of use
	// 

	//-------
	status = RM_BMI_GetValue_int(id, "GridCellCount", &nxyz);
	nxyz = RM_GetGridCellCount(id);
	status = RM_BMI_GetValuePtr_int(id, "GridCellCount", &i_ptr);
	fprintf(stderr, "GetValue('GridCellCount') \n");
	//-------

#ifdef SKIP
	int n = RM_BMI_GetThreadCount();
	fprintf(stderr, "GetThreadCount " << n << "\n";
	//-------
	// Inactive cells or symmetry
	std::vector<int> grid2chem(nxyz, -1);
	for (size_t i = 0; i < nxyz / 2; i++)
	{
		grid2chem[i] = i;
	}
	status = RM_BMI_CreateMapping(grid2chem);
	fprintf(stderr, "CreateMapping \n";
	//-------
	RM_BMI_LoadDatabase("phreeqc.dat");
	fprintf(stderr, "LoadDatabase\n";
	//
	// Set properties
	// 
	status = RM_BMI_SetComponentH2O(false);
	fprintf(stderr, "SetComponentH2O \n";
	//-------
	RM_BMI_SetSpeciesSaveOn(true);
	fprintf(stderr, "SetSpeciesSaveOn \n";
	//-------
	status = RM_BMI_SetErrorOn(true);
	fprintf(stderr, "SetErrorOn \n";
	//-------
	status = RM_BMI_SetErrorHandlerMode(1);
	fprintf(stderr, "SetErrorHandlerMode \n";
	//-------
	status = RM_BMI_SetDumpFileName("TestAllMethods_cpp.dump");
	fprintf(stderr, "SetDumpFileName \n";
	//-------
	status = RM_BMI_SetFilePrefix("TestAllMethods_cpp");
	RM_BMI_SetValue("FilePrefix", "TestAllMethods_cpp");
	fprintf(stderr, "SetFilePrefix \n";
	//-------
	status = RM_BMI_OpenFiles();
	fprintf(stderr, "OpenFiles \n";
	//-------
	status = RM_BMI_SetPartitionUZSolids(false);
	fprintf(stderr, "SetPartitionUZSolids \n";
	//-------
	status = RM_BMI_SetRebalanceByCell(true);
	fprintf(stderr, "SetRebalanceByCell \n";
	//-------
	status = RM_BMI_SetRebalanceFraction(0.5);
	fprintf(stderr, "SetRebalanceFraction \n";
	//-------
	status = RM_BMI_SetScreenOn(true);
	fprintf(stderr, "SetScreenOn \n";
	//-------
	status = RM_BMI_SetSelectedOutputOn(true);
	RM_BMI_SetValue("SelectedOutputOn", true);
	fprintf(stderr, "SetSelectedOutputOn \n";
	//-------
	status = RM_BMI_SetUnitsExchange(1);
	fprintf(stderr, "SetUnitsExchange \n";
	//-------
	status = RM_BMI_SetUnitsGasPhase(1);
	fprintf(stderr, "SetUnitsGasPhase \n";
	//-------
	status = RM_BMI_SetUnitsKinetics(1);
	fprintf(stderr, "SetUnitsKinetics \n";
	//-------
	status = RM_BMI_SetUnitsPPassemblage(1);
	fprintf(stderr, "SetUnitsPPassemblage \n";
	//-------
	status = RM_BMI_SetUnitsSolution(2);
	fprintf(stderr, "SetUnitsSolution \n";
	//-------
	status = RM_BMI_SetUnitsSSassemblage(1);
	fprintf(stderr, "SetUnitsSSassemblage \n";
	//-------
	status = RM_BMI_SetUnitsSurface(1);
	fprintf(stderr, "SetUnitsSurface \n";
	//-------
	RM_BMI_UseSolutionDensityVolume(false);
	fprintf(stderr, "UseSolutionDensityVolume \n";
	//-------
	double time_conversion = 1.0 / 86400.0;
	RM_BMI_SetTimeConversion(time_conversion);
	fprintf(stderr, "SetTimeConversion \n";
	//-------
	std::vector<double> v(nxyz, 1.0);
	status = RM_BMI_SetRepresentativeVolume(v);
	fprintf(stderr, "SetRepresentativeVolume \n";

	//-------Chemistry cells may be fewer than GridCellCount
	std::vector<int> vi(nxyz, 1);
	status = RM_BMI_SetPrintChemistryMask(vi);
	fprintf(stderr, "SetPrintChemistryMask \n";
	//-------
	status = RM_BMI_SetPrintChemistryOn(false, true, false);
	fprintf(stderr, "SetPrintChemistryOn \n";
	//
	// Define reactants available for initial 
	// and boundary conditions in this file
	//
	status = RM_BMI_RunFile(true, true, true, "all_reactants.pqi");
	fprintf(stderr, "RunFile \n";
	//-------
	RM_BMI_AddOutputVars("AddOutputVars", "True");
	RM_BMI_AddOutputVars("SolutionProperties", "True");
	RM_BMI_AddOutputVars("SolutionTotalMolalities", "True");
	RM_BMI_AddOutputVars("ExchangeMolalities", "True");
	RM_BMI_AddOutputVars("SurfaceMolalities", "True");
	RM_BMI_AddOutputVars("EquilibriumPhases", "True");
	RM_BMI_AddOutputVars("Gases", "True");
	RM_BMI_AddOutputVars("KineticReactants", "True");
	RM_BMI_AddOutputVars("SolidSolutions", "True");
	RM_BMI_AddOutputVars("CalculateValues", "True");
	RM_BMI_AddOutputVars("SolutionActivities", "H+ Ca+2 Na+");
	RM_BMI_AddOutputVars("SolutionMolalities", "OH- Cl-");
	RM_BMI_AddOutputVars("SaturationIndices", "Calcite Dolomite");
	fprintf(stderr, "AddOutputVars \n";
	//-------
	int ncomps = RM_BMI_FindComponents();
	fprintf(stderr, "FindComponents \n";
	//
	// Methods up to this point are useful 
	// in a YAML initialization file
	// 
	// Lists of reactants found by FindComponents follow
	// 
	int nchem = RM_BMI_GetChemistryCellCount();
	fprintf(stderr, "GetChemistryCellCount \n";
	//-------
	ncomps = RM_BMI_GetComponentCount();
	RM_BMI_GetValue("ComponentCount", ncomps);
	i_ptr = (int*) RM_BMI_GetValuePtr("ComponentCount");
	fprintf(stderr, "GetComponentCount \n";
	//-------
	std::vector<std::string> str_vector = RM_BMI_GetComponents();
	RM_BMI_GetValue("Components", str_vector);
	fprintf(stderr, "GetComponents \n";
	// Species info
	n = RM_BMI_GetSpeciesCount();
	fprintf(stderr, "GetSpeciesCount \n";
	//-------
	str_vector = RM_BMI_GetSpeciesNames();
	fprintf(stderr, "GetSpeciesNames \n";
	//-------
	v = RM_BMI_GetSpeciesD25();
	fprintf(stderr, "GetSpeciesD25 \n";
	//-------
	v = RM_BMI_GetSpeciesZ();
	fprintf(stderr, "GetSpeciesZ \n";
	// Reactant lists
	std::vector<std::string> equiuilibrium_phases = RM_BMI_GetEquilibriumPhases();
	fprintf(stderr, "GetEquilibriumPhases \n";
	//-------
	n = RM_BMI_GetEquilibriumPhasesCount();
	fprintf(stderr, "GetEquilibriumPhasesCount \n";
	//-------
	str_vector = RM_BMI_GetExchangeNames();
	fprintf(stderr, "GetExchangeNames \n";
	//-------
	str_vector = RM_BMI_GetExchangeSpecies();
	fprintf(stderr, "GetExchangeSpecies \n";
	//-------
	n = RM_BMI_GetExchangeSpeciesCount();
	fprintf(stderr, "GetExchangeSpeciesCount \n";
	//-------
	str_vector = RM_BMI_GetGasComponents();
	fprintf(stderr, "GetGasComponents \n";
	//-------
	n = RM_BMI_GetGasComponentsCount();
	fprintf(stderr, "GetGasComponentsCount \n";
	//-------
	v = RM_BMI_GetGfw();
	RM_BMI_GetValue("Gfw", v);
	double *d_ptr = (double*)RM_BMI_GetValuePtr("Gfw");
	fprintf(stderr, "GetGfw ";
	//-------
	str_vector = RM_BMI_GetKineticReactions();
	fprintf(stderr, "GetKineticReactions \n";
	//-------
	n = RM_BMI_GetKineticReactionsCount();
	fprintf(stderr, "GetKineticReactionsCount \n";
	//-------
	n = RM_BMI_GetSICount();
	fprintf(stderr, "GetSICount \n";
	//-------
	str_vector = RM_BMI_GetSINames();
	fprintf(stderr, "GetSINames \n";
	//-------
	str_vector = RM_BMI_GetSolidSolutionComponents();
	fprintf(stderr, "GetSolidSolutionComponents \n";
	//-------
	n = RM_BMI_GetSolidSolutionComponentsCount();
	fprintf(stderr, "GetSolidSolutionComponentsCount \n";
	//-------
	str_vector = RM_BMI_GetSolidSolutionNames();
	fprintf(stderr, "GetSolidSolutionNames \n";
	//-------
	str_vector = RM_BMI_GetSurfaceNames();
	fprintf(stderr, "GetSurfaceNames \n";
	//-------
	str_vector = RM_BMI_GetSurfaceSpecies();
	fprintf(stderr, "GetSurfaceSpecies \n";
	//-------
	n = RM_BMI_GetSurfaceSpeciesCount();
	fprintf(stderr, "GetSurfaceSpeciesCount \n";
	//-------
	str_vector = RM_BMI_GetSurfaceTypes();
	fprintf(stderr, "GetSurfaceTypes \n";
	//
	// Remove any reactants in workers 
	// before populating cells with reactants
	//
	std::string input = "DELETE; -all";
	RM_BMI_RunString(true, false, false, input);
	fprintf(stderr, "RunString \n";
	//-------
	// 
	// Transfer initial conditions to model
	//
	vi.clear();
	vi.resize(nxyz, 1);
	status = RM_BMI_InitialEquilibriumPhases2Module(vi);
	fprintf(stderr, "InitialEquilibriumPhases2Module \n";
	//-------
	status = RM_BMI_InitialExchanges2Module(vi);
	fprintf(stderr, "InitialExchanges2Module \n";
	//-------
	status = RM_BMI_InitialGasPhases2Module(vi);
	fprintf(stderr, "InitialGasPhases2Module \n";
	//-------
	status = RM_BMI_InitialKinetics2Module(vi);
	fprintf(stderr, "InitialKinetics2Module \n";
	//-------
	status = RM_BMI_InitialSolutions2Module(vi);
	fprintf(stderr, "InitialSolutions2Module \n";
	//-------
	status = RM_BMI_InitialSolidSolutions2Module(vi);
	fprintf(stderr, "InitialSolidSolutions2Module \n";
	//-------
	status = RM_BMI_InitialSurfaces2Module(vi);
	fprintf(stderr, "InitialSurfaces2Module \n";
	//-------
	// Alternative A.to the previous seven methods
	std::vector<int> ic(nxyz * 7, 1);
	status = RM_BMI_InitialPhreeqc2Module(ic);
	fprintf(stderr, "InitialPhreeqc2Module \n";
	//-------
	// Alternative B.to the previous seven methods, possible mixing
	std::vector<int> v1(nxyz * 7, 1);
	std::vector<int> v2(nxyz * 7, -1);
	std::vector<double> f1(nxyz * 7, 1.0);
	status = RM_BMI_InitialPhreeqc2Module(v1, v2, f1);
	fprintf(stderr, "InitialPhreeqc2Modul mix \n";
	//-------
	// Alternative C.to the previous seven methods, initialize cells 18 and 19
	std::vector<int> cells(2);
	cells[0] = 18;
	cells[1] = 19;
	status = RM_BMI_InitialPhreeqcCell2Module(1, cells);
	fprintf(stderr, "InitialPhreeqcCell2Module \n";
	//
	// Boundary conditions
	// 
	vi.clear();
	vi.resize(1, 1);
	std::vector<double> bc;
	status = RM_BMI_InitialPhreeqc2Concentrations(bc, vi);
	//-------
	std::vector<int> vi1(1, 1);
	std::vector<int> vi2(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = RM_BMI_InitialPhreeqc2Concentrations(bc, vi1, vi2, f1);
	//-------
	vi.clear();
	vi.resize(1, 1);
	status = RM_BMI_InitialPhreeqc2SpeciesConcentrations(bc, vi);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations \n";
	//-------
	std::vector<double> bc_species;
	vi1.clear();
	vi1.resize(1, 1);
	vi2.clear();
	vi2.resize(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = RM_BMI_InitialPhreeqc2SpeciesConcentrations(bc_species, vi1, vi2, f1);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations mix \n";
	//
	// Get/Set methods for time steping
	//
	double d = RM_BMI_GetTime();
	RM_BMI_GetValue("Time", d);
	d = RM_BMI_GetCurrentTime();
    d = RM_BMI_GetStartTime();
	d_ptr = (double*)RM_BMI_GetValuePtr("Time");
	fprintf(stderr, "GetTime \n";
	//-------
	status = RM_BMI_SetTime(0.0);
	RM_BMI_SetValue("Time", 0.0);
	fprintf(stderr, "SetTime \n";
	//-------
	d = RM_BMI_GetTimeStep();
	RM_BMI_GetValue("TimeStep", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("TimeStep");
	fprintf(stderr, "GetTimeStep \n";
	//-------
	status = RM_BMI_SetTimeStep(0.0);
	RM_BMI_SetValue("TimeStep", 0.0);
	fprintf(stderr, "SetTimeStep \n";
	//-------
	std::vector<double> c;
	status = RM_BMI_GetConcentrations(c);
	RM_BMI_GetValue("Concentrations", c);
	d_ptr = (double*)RM_BMI_GetValuePtr("Concentrations");
	fprintf(stderr, "GetConcentrations \n";
	//-------
	status = RM_BMI_SetConcentrations(c);
	RM_BMI_SetValue("Concentrations", c);
	fprintf(stderr, "SetConcentrations \n";
	//-------
	status = RM_BMI_GetDensityCalculated(v);
	RM_BMI_GetValue("DensityCalculated", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("DensityCalculated");
	fprintf(stderr, "GetDensityCalculated \n";
	//-------
	status = RM_BMI_SetDensityUser(v);
	RM_BMI_SetValue("DensityUser", v);
	fprintf(stderr, "SetDensityUser \n";
	//-------
	status = RM_BMI_GetGasCompMoles(v);
	fprintf(stderr, "GetGasCompMoles \n";
	//-------
	status = RM_BMI_SetGasCompMoles(v);
	fprintf(stderr, "SetGasCompMoles \n";
	//-------
	status = RM_BMI_GetGasCompPhi(v);
	fprintf(stderr, "GetGasCompPhi \n";
	//-------
	status = RM_BMI_GetGasCompPressures(v);
	fprintf(stderr, "GetGasCompPressures \n";
	//-------
	status = RM_BMI_GetGasPhaseVolume(v);
	fprintf(stderr, "GetGasPhaseVolume \n";
	//-------
	status = RM_BMI_SetGasPhaseVolume(v);
	fprintf(stderr, "SetGasPhaseVolume \n";
	//-------
	for (size_t i = 0; i < RM_BMI_GetComponentCount(); i++)
	{
		status = RM_BMI_GetIthConcentration(i, v);
		//-------
		status = RM_BMI_SetIthConcentration(i, v);
	}
	fprintf(stderr, "GetIthConcentration \n";
	fprintf(stderr, "SetIthConcentration \n";
	//-------
	for (int i = 0; i < RM_BMI_GetSpeciesCount(); i++)
	{
		status = RM_BMI_GetIthSpeciesConcentration(i, v);
		//-------
		status = RM_BMI_SetIthSpeciesConcentration(i, v);
	}
	fprintf(stderr, "GetIthSpeciesConcentration \n";
	fprintf(stderr, "SetIthSpeciesConcentration \n";
	//-------
	v = RM_BMI_GetPorosity();
	RM_BMI_GetValue("Porosity", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("Porosity");
	fprintf(stderr, "GetPorosity \n";
	//-------
	status = RM_BMI_SetPorosity(v);
	RM_BMI_SetValue("Porosity", v);
	fprintf(stderr, "SetPorosity \n";
	//-------
	v = RM_BMI_GetPressure();
	RM_BMI_GetValue("Pressure", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("Pressure");
	fprintf(stderr, "GetPressure \n";
	//-------
	status = RM_BMI_SetPressure(v);
	RM_BMI_SetValue("Pressure", v);
	fprintf(stderr, "SetPressure \n";
	//-------
	status = RM_BMI_GetSaturationCalculated(v);
	RM_BMI_GetValue("SaturationCalculated", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("SaturationCalculated");
	fprintf(stderr, "GetSaturationCalculated \n";
	//-------
	status = RM_BMI_SetSaturationUser(v);
	RM_BMI_SetValue("SaturationUser", v);
	fprintf(stderr, "SetSaturationUser \n";
	//-------
	v = RM_BMI_GetSolutionVolume();
	RM_BMI_GetValue("SolutionVolume", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("SolutionVolume");
	fprintf(stderr, "GetSolutionVolume \n";
	//-------
	status = RM_BMI_GetSpeciesConcentrations(v);
	fprintf(stderr, "GetSpeciesConcentrations \n";
	//-------
	status = RM_BMI_SpeciesConcentrations2Module(v);
	fprintf(stderr, "SpeciesConcentrations2Module \n";
	//-------
	status = RM_BMI_GetSpeciesLog10Gammas(v);
	fprintf(stderr, "GetSpeciesLog10Gammas \n";
	//-------
	status = RM_BMI_GetSpeciesLog10Molalities(v);
	fprintf(stderr, "GetSpeciesLog10Molalities \n";
	//-------
	v = RM_BMI_GetTemperature();
	RM_BMI_GetValue("Temperature", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("Temperature");
	fprintf(stderr, "GetTemperature \n";
	//-------
	status = RM_BMI_SetTemperature(v);
	RM_BMI_SetValue("Temperature", v);
	fprintf(stderr, "SetTemperature \n";
	//-------
	v = RM_BMI_GetViscosity();
	RM_BMI_GetValue("Viscosity", v);
	d_ptr = (double*)RM_BMI_GetValuePtr("Viscosity");	
	fprintf(stderr, "GetViscosity \n";
	//
	// Take a time step
	//
	RM_BMI_Update();      // void function
	fprintf(stderr, "Update\n";
	//-------
	status = RM_BMI_RunCells();
	fprintf(stderr, "RunCells\n";
	//-------
	RM_BMI_UpdateUntil(86400.0);      // void function
	fprintf(stderr, "UpdateUntil\n";
	//
	// Selected output
	//
	status = RM_BMI_SetNthSelectedOutput(0);
	RM_BMI_SetValue("NthSelectedOutput", 0);
	fprintf(stderr, "SetNthSelectedOutput \n";
	//-------
	int n_user = RM_BMI_GetCurrentSelectedOutputUserNumber();
	RM_BMI_GetValue("CurrentSelectedOutputUserNumber", n_user);
	fprintf(stderr, "GetCurrentSelectedOutputUserNumber \n";
	//-------
	n = RM_BMI_GetNthSelectedOutputUserNumber(0);
	fprintf(stderr, "GetNthSelectedOutputUserNumber \n";
	//-------
	status = RM_BMI_GetSelectedOutput(v);
	RM_BMI_GetValue("SelectedOutput", v);
	fprintf(stderr, "GetSelectedOutput \n";
	//-------
	n = RM_BMI_GetSelectedOutputColumnCount();
	RM_BMI_GetValue("SelectedOutputColumnCount", n);
	fprintf(stderr, "GetSelectedOutputColumnCount \n";
	//-------
	n = RM_BMI_GetSelectedOutputCount();
	RM_BMI_GetValue("SelectedOutputCount", n);
	fprintf(stderr, "GetSelectedOutputCount \n";
	//-------
	status = RM_BMI_GetSelectedOutputHeadings(str_vector);
	RM_BMI_GetValue("SelectedOutputHeadings", str_vector);
	fprintf(stderr, "GetSelectedOutputHeadings \n";
	//-------
	bool b = RM_BMI_GetSelectedOutputOn();
	RM_BMI_GetValue("SelectedOutputOn", b);
	bool *b_ptr = (bool*)RM_BMI_GetValuePtr("SelectedOutputOn");
	fprintf(stderr, "GetSelectedOutputOn \n";
	//-------
	n = RM_BMI_GetSelectedOutputRowCount(); 
	RM_BMI_GetValue("SelectedOutputRowCount", n);
	fprintf(stderr, "GetSelectedOutputRowCount \n";
	//-------
	status = RM_BMI_SetCurrentSelectedOutputUserNumber(333);
	fprintf(stderr, "SetCurrentSelectedOutputUserNumber \n";
	//-------
	//
	// Getters
	// 
	std::vector< std::vector<int> > back_map = RM_BMI_GetBackwardMapping();
	fprintf(stderr, "GetBackwardMapping \n";
	//-------
	std::string db_name = RM_BMI_GetDatabaseFileName();
	fprintf(stderr, "GetDatabaseFileName \n";
	//-------
	vi = RM_BMI_GetEndCell();
	fprintf(stderr, "GetEndCell\n";
	//-------
	n = RM_BMI_GetErrorHandlerMode();
	fprintf(stderr, "GetErrorHandlerMode \n";
	//-------
	std::string str = RM_BMI_GetErrorString();
	RM_BMI_GetValue("ErrorString", str);
	fprintf(stderr, "GetErrorString \n";
	//-------
	str = RM_BMI_GetFilePrefix();
	RM_BMI_GetValue("FilePrefix", str);
	fprintf(stderr, "GetFilePrefix \n";
	//-------
	vi = RM_BMI_GetForwardMapping();
	fprintf(stderr, "GetForwardMapping \n";
	//-------
	IPhreeqc* ipq = RM_BMI_GetIPhreeqcPointer(0);
	fprintf(stderr, "GetIPhreeqcPointer \n";
	//-------
	n = RM_BMI_GetMpiMyself();
	fprintf(stderr, "GetMpiMyself \n";
	//-------
	n = RM_BMI_GetMpiTasks();
	fprintf(stderr, "GetMpiTasks \n";
	//-------
	b = RM_BMI_GetPartitionUZSolids();
	fprintf(stderr, "GetPartitionUZSolids \n";
	//-------
	vi = RM_BMI_GetPrintChemistryMask();
	fprintf(stderr, "GetPrintChemistryMask \n";
	//-------
	std::vector<bool> vb = RM_BMI_GetPrintChemistryOn();
	fprintf(stderr, "GetPrintChemistryOn \n";
	//-------
	b = RM_BMI_GetRebalanceByCell();
	fprintf(stderr, "GetRebalanceByCell \n";
	//-------
	d = RM_BMI_GetRebalanceFraction();
	fprintf(stderr, "GetRebalanceFraction \n";
	//-------
	b = RM_BMI_GetSpeciesSaveOn();
	fprintf(stderr, "GetSpeciesSaveOn \n";
	//-------
	std::vector< cxxNameDouble > s = RM_BMI_GetSpeciesStoichiometry();
	fprintf(stderr, "GetSpeciesStoichiometry \n";
	//-------
	vi = RM_BMI_GetStartCell();
	fprintf(stderr, "GetStartCell \n";
	//-------
	d = RM_BMI_GetTimeConversion();
	fprintf(stderr, "GetTimeConversion \n";
	//-------
	n = RM_BMI_GetUnitsExchange();
	fprintf(stderr, "GetUnitsExchange \n";
	//-------
	n = RM_BMI_GetUnitsGasPhase();
	fprintf(stderr, "GetUnitsGasPhase \n";
	//-------
	n = RM_BMI_GetUnitsKinetics();
	fprintf(stderr, "GetUnitsKinetics \n";
	//-------
	n = RM_BMI_GetUnitsPPassemblage();
	fprintf(stderr, "GetUnitsPPassemblage \n";
	//-------
	n = RM_BMI_GetUnitsSolution();
	fprintf(stderr, "GetUnitsSolution \n";
	//-------
	n = RM_BMI_GetUnitsSSassemblage();
	fprintf(stderr, "GetUnitsSSassemblage \n";
	//-------
	n = RM_BMI_GetUnitsSurface();
	fprintf(stderr, "GetUnitsSurface \n";
	//-------
	std::vector<IPhreeqcPhast *> w = RM_BMI_GetWorkers();
	fprintf(stderr, "GetWorkers \n";
	//
	// Utilities
	//
#ifndef USE_MPI
	BMIPhreeqcRM *bmi2 = new BMIPhreeqcRM(10, 2); // Make another instance
	fprintf(stderr, "Make a new instance with new. \n";
	//-------
	status = bmi2->CloseFiles(); 
	fprintf(stderr, "CloseFiles \n";
	//-------
	bmi2->Finalize();
	delete bmi2; // delete new instance
#endif
	//-------
	vi.clear();
	vi.resize(1, 1);
	status = RM_BMI_InitialPhreeqc2Concentrations(bc, vi);
	std::vector<double> tc(1, 30.0);
	std::vector<double> p_atm(1, 1.5);
	IPhreeqc* utility_ptr = RM_BMI_Concentrations2Utility(bc, tc, p_atm);
	fprintf(stderr, "Concentrations2Utility \n";
	//-------
	RM_BMI_DecodeError(-2);	         // void function
	fprintf(stderr, "DecodeError \n";
	//-------
	status = RM_BMI_DumpModule(true);
	fprintf(stderr, "DumpModule \n";
	//-------
	RM_BMI_ErrorHandler(0, "string"); // void function
	fprintf(stderr, "OK, just a test: ErrorHandler \n";
	//-------
	RM_BMI_ErrorMessage("my error");  // void function
	fprintf(stderr, "OK, just a test: ErrorMessage \n";
	//-------
	RM_BMI_LogMessage("Log message");  // void method
	fprintf(stderr, "LogMessage \n";
	//-------
	RM_BMI_OutputMessage("Output message");  // void method
	fprintf(stderr, "OutputMessage \n";
	//-------
	RM_BMI_ScreenMessage("Screen message\n");  // void method
	fprintf(stderr, "ScreenMessage \n";
	//-------
	status = RM_BMI_StateSave(1);
	fprintf(stderr, "StateSave \n";
	//-------
	status = RM_BMI_StateApply(1);
	fprintf(stderr, "StateApply \n";
	//-------
	status = RM_BMI_StateDelete(1);
	fprintf(stderr, "StateDelete \n";
	//-------
	RM_BMI_WarningMessage("Warning message");  // void method
	fprintf(stderr, "WarningMessage \n";
	//
	// BMI Methods
	//
	str = RM_BMI_GetComponentName();
	fprintf(stderr, "GetComponentName \n";
	//-------
	d = RM_BMI_GetCurrentTime();
	fprintf(stderr, "GetCurrentTime \n";
	//-------
	d = RM_BMI_GetEndTime();
	fprintf(stderr, "GetEndTime \n";
	//-------
    n = RM_BMI_GetGridRank(0);
	fprintf(stderr, "GetGridRank \n";
	//-------
    n = RM_BMI_GetGridSize(0);
	fprintf(stderr, "GetGridSize \n";
	//-------
    std::string gtype = RM_BMI_GetGridType(0);
	fprintf(stderr, "GetGridType \n";
	//-------
	n = RM_BMI_GetInputItemCount();
	fprintf(stderr, "GetInputItemCount \n";
	//-------
	str_vector = RM_BMI_GetInputVarNames();
	fprintf(stderr, "GetInputVarNames \n";
	//-------
	n = RM_BMI_GetOutputItemCount();
	fprintf(stderr, "GetOutputItemCount \n";
	//-------
	str_vector = RM_BMI_GetOutputVarNames();
	fprintf(stderr, "GetOutputVarNames \n";
	//-------
	d = RM_BMI_GetTimeStep();
	fprintf(stderr, "GetTimeStep \n";
	//-------
	str = RM_BMI_GetTimeUnits();
	fprintf(stderr, "GetTimeUnits \n";
	//-------
	RM_BMI_GetValue("solution_saturation_index_Calcite", v);
	fprintf(stderr, "GetValue \n";
	//-------
	n = RM_BMI_GetVarItemsize("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarItemsize \n";
	//-------
	n = RM_BMI_GetVarNbytes("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarNbytes \n";
	//-------
	str = RM_BMI_GetVarType("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarType \n";
	//-------
	str = RM_BMI_GetVarUnits("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarUnits \n";
	//-------
	//RM_BMI_Initialize(YAML_filename);
	// See above
	RM_BMI_SetValue("Time", 1.0);    // void method
	fprintf(stderr, "SetValue\n";
	//-------
	RM_BMI_Update();    // void method
	fprintf(stderr, "Update\n";
	//-------	
 	RM_BMI_UpdateUntil(864000.0);      // void function
	fprintf(stderr, "UpdateUntil\n";

	fprintf(stderr, "AddOutputVars\n";
	str_vector = RM_BMI_GetOutputVarNames();
	for (size_t i = 0; i < str_vector.size(); i++)
	{
		int	itemsize = RM_BMI_GetVarItemsize(str_vector[i]);
		int	nbytes = RM_BMI_GetVarNbytes(str_vector[i]);
		std::string	vtype = RM_BMI_GetVarType(str_vector[i]);
		if (itemsize == 0) itemsize++;
		if (nbytes == 0) nbytes++;
		int dim = nbytes / itemsize;
		if (vtype == "double")
		{
			if (dim == 1)
			{
				double dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<double> dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
		else if (vtype == "int")
		{
			if (dim == 1)
			{
				int dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<int> dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		} else if (vtype == "bool")
		{
			if (dim == 1)
			{
				bool dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
		}
		else if (vtype == "std::string")
		{
			if (dim == 1)
			{
				std::string dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<std::string> dest;
				RM_BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
	}


	//-------	
	RM_BMI_MpiWorkerBreak();
	RM_BMI_Finalize();    // void method
	fprintf(stderr, "Finalize \n";
	//Should be private: status = RM_BMI_ReturnHandler();
	//TODO status = RM_BMI_MpiAbort();
	//TODO status = RM_BMI_SetMpiWorkerCallbackC();
	//TODO status = RM_BMI_SetMpiWorkerCallbackCookie();
#endif
	fprintf(stderr, "Success.\n");
	return;
}
