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
	int id = -1;
	int* i_ptr = NULL;
	int* grid2chem = NULL;
	double* v = NULL;
	int* vi = NULL;
	int i = 0, nchem = 0;
	// Set GridCellCount
	//yrm.YAMLSetGridCellCount(nxyz);
	//yrm.YAMLThreadCount(3);
	char string[MAX_LENGTH] = "", string1[MAX_LENGTH] = "";
	char YAML_filename[MAX_LENGTH] = "";
	strcpy_safe(YAML_filename, MAX_LENGTH, "TestAllMethods_c.yaml");

#ifdef USE_MPI
	// MPI
	int mpi_myself;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}

	id = BMI_Create(nxyz, MPI_COMM_WORLD);
	if (mpi_myself > 0)
	{
		BMI_MpiWorker();
		BMI_Finalize();
		return;
	}
#else
	// OpenMP or serial
	nxyz = RM_GetGridCellCountYAML(YAML_filename);
	id = BMI_Create(nxyz, nthreads);
#endif
	// Use YAML file to initialize
	BMI_Initialize(id, YAML_filename);   // void function
	RM_InitializeYAML(id, YAML_filename);
	fprintf(stderr, "Initialize\n");

	//
	// Use all BMIPhreeqcRM methods roughly in order of use
	// 
	status = BMI_GetValueInt(id, "GridCellCount", &nxyz);
	nxyz = RM_GetGridCellCount(id);
	i_ptr = BMI_GetValuePtr(id, "GridCellCount");
	fprintf(stderr, "GetValue('GridCellCount') %d %d %d\n", nxyz, *i_ptr, status);

	//-------
	int n = RM_GetThreadCount(id);
	fprintf(stderr, "GetThreadCount %d \n", n);
	//-------
	// Inactive cells or symmetry

	grid2chem = (int*)malloc(nxyz * sizeof(int));
	for (int i = 0; i < nxyz; i++)
	{
		grid2chem[i] = i;
	}
	status = RM_CreateMapping(id, grid2chem);
	fprintf(stderr, "CreateMapping %d\n", status);
	//-------
	status = RM_LoadDatabase(id, "phreeqc.dat");
	fprintf(stderr, "LoadDatabase %d\n", status);
	//
	// Set properties
	// 
	status = RM_SetComponentH2O(id, 0);
	fprintf(stderr, "SetComponentH2O %d\n", status);
	//-------
	RM_SetSpeciesSaveOn(id, 1);
	fprintf(stderr, "SetSpeciesSaveOn %d\n", status);
	//-------
	status = RM_SetErrorOn(id, 1);
	fprintf(stderr, "SetErrorOn %d\n", status);
	//-------
	status = RM_SetErrorHandlerMode(id, 1);
	fprintf(stderr, "SetErrorHandlerMode %d\n", status);
	//-------
	status = RM_SetDumpFileName(id, "TestAllMethods_cpp.dump");
	fprintf(stderr, "SetDumpFileName %d\n", status);
	//-------
	status = RM_SetFilePrefix(id, "TestAllMethods_c");
	BMI_SetValueChar(id, "FilePrefix", "TestAllMethods_c");
	BMI_GetValueChar(id, "FilePrefix", string);
	fprintf(stderr, "SetFilePrefix %s\n", string);
	//-------
	status = RM_OpenFiles(id);
	fprintf(stderr, "OpenFiles %d\n", status);
	//-------
	status = RM_SetPartitionUZSolids(id, 0);
	fprintf(stderr, "SetPartitionUZSolids %d\n", status);
	//-------
	status = RM_SetRebalanceByCell(id, 1);
	fprintf(stderr, "SetRebalanceByCell %d\n", status);
	//-------
	status = RM_SetRebalanceFraction(id, 0.5);
	fprintf(stderr, "SetRebalanceFraction %d\n", status);
	//-------
	status = RM_SetScreenOn(id, 1);
	fprintf(stderr, "SetScreenOn %d\n", status);
	//-------
	status = RM_SetSelectedOutputOn(id, 1);
	BMI_SetValueInt(id, "SelectedOutputOn", 1);
	fprintf(stderr, "SetSelectedOutputOn %d\n", status);
	//-------
	status = RM_SetUnitsExchange(id, 1);
	fprintf(stderr, "SetUnitsExchange %d\n", status);
	//-------
	status = RM_SetUnitsGasPhase(id, 1);
	fprintf(stderr, "SetUnitsGasPhase %d\n", status);
	//-------
	status = RM_SetUnitsKinetics(id, 1);
	fprintf(stderr, "SetUnitsKinetics %d\n", status);
	//-------
	status = RM_SetUnitsPPassemblage(id, 1);
	fprintf(stderr, "SetUnitsPPassemblage %d\n", status);
	//-------
	status = RM_SetUnitsSolution(id, 2);
	fprintf(stderr, "SetUnitsSolution %d\n", status);
	//-------
	status = RM_SetUnitsSSassemblage(id, 1);
	fprintf(stderr, "SetUnitsSSassemblage %d\n", status);
	//-------
	status = RM_SetUnitsSurface(id, 1);
	fprintf(stderr, "SetUnitsSurface %d\n", status);
	//-------
	status = RM_UseSolutionDensityVolume(id, 0);
	fprintf(stderr, "UseSolutionDensityVolume %d\n", status);
	//-------
	double time_conversion = 1.0 / 86400.0;
	status = RM_SetTimeConversion(id, time_conversion);
	fprintf(stderr, "SetTimeConversion %d\n", status);
	//-------
	if (v != NULL) free(v);
	v = (double*)malloc(nxyz * sizeof(double));
	for (int i = 0; i < nxyz; i++) v[i] = 1.0;
	status = RM_SetRepresentativeVolume(id, v);
	fprintf(stderr, "SetRepresentativeVolume %d\n", status);

	//-------Chemistry cells may be fewer than GridCellCount
	if (vi != NULL) free(vi);
	vi = (int*)malloc(nxyz * sizeof(int));
	for (int i = 0; i < nxyz; i++) v[i] = 1;
	status = RM_SetPrintChemistryMask(id, vi);
	fprintf(stderr, "SetPrintChemistryMask %d\n", status);
	//-------
	status = RM_SetPrintChemistryOn(id, 0, 1, 0);
	fprintf(stderr, "SetPrintChemistryOn %d\n", status);
	//
	// Define reactants available for initial 
	// and boundary conditions in this file
	//
	status = RM_RunFile(id, 1, 1, 1, "all_reactants.pqi");
	fprintf(stderr, "RunFile %d\n", status);
	//-------
	status = BMI_AddOutputVars(id, "AddOutputVars", "True");
	status = BMI_AddOutputVars(id, "SolutionProperties", "True");
	status = BMI_AddOutputVars(id, "SolutionTotalMolalities", "True");
	status = BMI_AddOutputVars(id, "ExchangeMolalities", "True");
	status = BMI_AddOutputVars(id, "SurfaceMolalities", "True");
	status = BMI_AddOutputVars(id, "EquilibriumPhases", "True");
	status = BMI_AddOutputVars(id, "Gases", "True");
	status = BMI_AddOutputVars(id, "KineticReactants", "True");
	status = BMI_AddOutputVars(id, "SolidSolutions", "True");
	status = BMI_AddOutputVars(id, "CalculateValues", "True");
	status = BMI_AddOutputVars(id, "SolutionActivities", "H+ Ca+2 Na+");
	status = BMI_AddOutputVars(id, "SolutionMolalities", "OH- Cl-");
	status = BMI_AddOutputVars(id, "SaturationIndices", "Calcite Dolomite");
	fprintf(stderr, "AddOutputVars %d\n", status);
	//-------
	int ncomps = RM_FindComponents(id);
	fprintf(stderr, "FindComponents %d\n", status);
	//
	// Methods up to this point are useful 
	// in a YAML initialization file
	// 
	// Lists of reactants found by FindComponents follow
	// 
	nchem = RM_GetChemistryCellCount(id);
	fprintf(stderr, "GetChemistryCellCount %d\n", nchem);
	//-------
	ncomps = RM_GetComponentCount(id);
	status = BMI_GetValueInt(id, "ComponentCount", &ncomps);
	i_ptr = (int*) BMI_GetValuePtr(id, "ComponentCount");
	fprintf(stderr, "GetComponentCount %d %d\n", ncomps, *i_ptr);
	//-------
	for (i = 0; i < ncomps; i++)
	{
		status = RM_GetComponent(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	//??????  BMI_GetValue("Components", str_vector);    //?????
	fprintf(stderr, "GetComponent %d\n", status);
	// Species info
	n = RM_GetSpeciesCount(id);
	fprintf(stderr, "GetSpeciesCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetSpeciesName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetSpeciesName %d\n", status);
	//-------
	status = RM_GetSpeciesD25(id, v);
	fprintf(stderr, "GetSpeciesD25 %d\n", status);
	//-------
	status = RM_GetSpeciesZ(id, v);
	fprintf(stderr, "GetSpeciesZ %d\n", status);
	// Reactant lists
	//-------
	n = RM_GetEquilibriumPhasesCount(id);
	fprintf(stderr, "GetEquilibriumPhasesCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetEquilibriumPhasesName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetEquilibriumPhases %d\n", status);
	//-------
	n = RM_GetExchangeSpeciesCount(id);
	fprintf(stderr, "GetExchangeSpeciesCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetExchangeSpeciesName(id, i, string, MAX_LENGTH);
		status = RM_GetExchangeName(id, i, string1, MAX_LENGTH);
		fprintf(stderr, "     %20s %20s\n", string, string1);
	}
	fprintf(stderr, "GetExchangeSpeciesName \n");
	fprintf(stderr, "GetExchangeName \n");
	//-------
	n = RM_GetGasComponentsCount(id);
	fprintf(stderr, "GetGasComponentsCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetGasComponentsName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetGasComponents %d\n", status);

	fprintf(stderr, "Done.\n"); //==================================================================
	return;
#ifdef SKIP
	//-------
	v = BMI_GetGfw();
	BMI_GetValue("Gfw", v);
	double *d_ptr = (double*)BMI_GetValuePtr("Gfw");
	fprintf(stderr, "GetGfw ";
	//-------
	str_vector = BMI_GetKineticReactions();
	fprintf(stderr, "GetKineticReactions \n";
	//-------
	n = BMI_GetKineticReactionsCount();
	fprintf(stderr, "GetKineticReactionsCount \n";
	//-------
	n = BMI_GetSICount();
	fprintf(stderr, "GetSICount \n";
	//-------
	str_vector = BMI_GetSINames();
	fprintf(stderr, "GetSINames \n";
	//-------
	str_vector = BMI_GetSolidSolutionComponents();
	fprintf(stderr, "GetSolidSolutionComponents \n";
	//-------
	n = BMI_GetSolidSolutionComponentsCount();
	fprintf(stderr, "GetSolidSolutionComponentsCount \n";
	//-------
	str_vector = BMI_GetSolidSolutionNames();
	fprintf(stderr, "GetSolidSolutionNames \n";
	//-------
	str_vector = BMI_GetSurfaceNames();
	fprintf(stderr, "GetSurfaceNames \n";
	//-------
	str_vector = BMI_GetSurfaceSpecies();
	fprintf(stderr, "GetSurfaceSpecies \n";
	//-------
	n = BMI_GetSurfaceSpeciesCount();
	fprintf(stderr, "GetSurfaceSpeciesCount \n";
	//-------
	str_vector = BMI_GetSurfaceTypes();
	fprintf(stderr, "GetSurfaceTypes \n";
	//
	// Remove any reactants in workers 
	// before populating cells with reactants
	//
	std::string input = "DELETE; -all";
	BMI_RunString(true, false, false, input);
	fprintf(stderr, "RunString \n";
	//-------
	// 
	// Transfer initial conditions to model
	//
	vi.clear();
	vi.resize(nxyz, 1);
	status = BMI_InitialEquilibriumPhases2Module(vi);
	fprintf(stderr, "InitialEquilibriumPhases2Module \n";
	//-------
	status = BMI_InitialExchanges2Module(vi);
	fprintf(stderr, "InitialExchanges2Module \n";
	//-------
	status = BMI_InitialGasPhases2Module(vi);
	fprintf(stderr, "InitialGasPhases2Module \n";
	//-------
	status = BMI_InitialKinetics2Module(vi);
	fprintf(stderr, "InitialKinetics2Module \n";
	//-------
	status = BMI_InitialSolutions2Module(vi);
	fprintf(stderr, "InitialSolutions2Module \n";
	//-------
	status = BMI_InitialSolidSolutions2Module(vi);
	fprintf(stderr, "InitialSolidSolutions2Module \n";
	//-------
	status = BMI_InitialSurfaces2Module(vi);
	fprintf(stderr, "InitialSurfaces2Module \n";
	//-------
	// Alternative A.to the previous seven methods
	std::vector<int> ic(nxyz * 7, 1);
	status = BMI_InitialPhreeqc2Module(ic);
	fprintf(stderr, "InitialPhreeqc2Module \n";
	//-------
	// Alternative B.to the previous seven methods, possible mixing
	std::vector<int> v1(nxyz * 7, 1);
	std::vector<int> v2(nxyz * 7, -1);
	std::vector<double> f1(nxyz * 7, 1.0);
	status = BMI_InitialPhreeqc2Module(v1, v2, f1);
	fprintf(stderr, "InitialPhreeqc2Modul mix \n";
	//-------
	// Alternative C.to the previous seven methods, initialize cells 18 and 19
	std::vector<int> cells(2);
	cells[0] = 18;
	cells[1] = 19;
	status = BMI_InitialPhreeqcCell2Module(1, cells);
	fprintf(stderr, "InitialPhreeqcCell2Module \n";
	//
	// Boundary conditions
	// 
	vi.clear();
	vi.resize(1, 1);
	std::vector<double> bc;
	status = BMI_InitialPhreeqc2Concentrations(bc, vi);
	//-------
	std::vector<int> vi1(1, 1);
	std::vector<int> vi2(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = BMI_InitialPhreeqc2Concentrations(bc, vi1, vi2, f1);
	//-------
	vi.clear();
	vi.resize(1, 1);
	status = BMI_InitialPhreeqc2SpeciesConcentrations(bc, vi);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations \n";
	//-------
	std::vector<double> bc_species;
	vi1.clear();
	vi1.resize(1, 1);
	vi2.clear();
	vi2.resize(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = BMI_InitialPhreeqc2SpeciesConcentrations(bc_species, vi1, vi2, f1);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations mix \n";
	//
	// Get/Set methods for time steping
	//
	double d = BMI_GetTime();
	BMI_GetValue("Time", d);
	d = BMI_GetCurrentTime();
    d = BMI_GetStartTime();
	d_ptr = (double*)BMI_GetValuePtr("Time");
	fprintf(stderr, "GetTime \n";
	//-------
	status = BMI_SetTime(0.0);
	BMI_SetValue("Time", 0.0);
	fprintf(stderr, "SetTime \n";
	//-------
	d = BMI_GetTimeStep();
	BMI_GetValue("TimeStep", v);
	d_ptr = (double*)BMI_GetValuePtr("TimeStep");
	fprintf(stderr, "GetTimeStep \n";
	//-------
	status = BMI_SetTimeStep(0.0);
	BMI_SetValue("TimeStep", 0.0);
	fprintf(stderr, "SetTimeStep \n";
	//-------
	std::vector<double> c;
	status = BMI_GetConcentrations(c);
	BMI_GetValue("Concentrations", c);
	d_ptr = (double*)BMI_GetValuePtr("Concentrations");
	fprintf(stderr, "GetConcentrations \n";
	//-------
	status = BMI_SetConcentrations(c);
	BMI_SetValue("Concentrations", c);
	fprintf(stderr, "SetConcentrations \n";
	//-------
	status = BMI_GetDensityCalculated(v);
	BMI_GetValue("DensityCalculated", v);
	d_ptr = (double*)BMI_GetValuePtr("DensityCalculated");
	fprintf(stderr, "GetDensityCalculated \n";
	//-------
	status = BMI_SetDensityUser(v);
	BMI_SetValue("DensityUser", v);
	fprintf(stderr, "SetDensityUser \n";
	//-------
	status = BMI_GetGasCompMoles(v);
	fprintf(stderr, "GetGasCompMoles \n";
	//-------
	status = BMI_SetGasCompMoles(v);
	fprintf(stderr, "SetGasCompMoles \n";
	//-------
	status = BMI_GetGasCompPhi(v);
	fprintf(stderr, "GetGasCompPhi \n";
	//-------
	status = BMI_GetGasCompPressures(v);
	fprintf(stderr, "GetGasCompPressures \n";
	//-------
	status = BMI_GetGasPhaseVolume(v);
	fprintf(stderr, "GetGasPhaseVolume \n";
	//-------
	status = BMI_SetGasPhaseVolume(v);
	fprintf(stderr, "SetGasPhaseVolume \n";
	//-------
	for (size_t i = 0; i < BMI_GetComponentCount(); i++)
	{
		status = BMI_GetIthConcentration(i, v);
		//-------
		status = BMI_SetIthConcentration(i, v);
	}
	fprintf(stderr, "GetIthConcentration \n";
	fprintf(stderr, "SetIthConcentration \n";
	//-------
	for (int i = 0; i < BMI_GetSpeciesCount(); i++)
	{
		status = BMI_GetIthSpeciesConcentration(i, v);
		//-------
		status = BMI_SetIthSpeciesConcentration(i, v);
	}
	fprintf(stderr, "GetIthSpeciesConcentration \n";
	fprintf(stderr, "SetIthSpeciesConcentration \n";
	//-------
	v = BMI_GetPorosity();
	BMI_GetValue("Porosity", v);
	d_ptr = (double*)BMI_GetValuePtr("Porosity");
	fprintf(stderr, "GetPorosity \n";
	//-------
	status = BMI_SetPorosity(v);
	BMI_SetValue("Porosity", v);
	fprintf(stderr, "SetPorosity \n";
	//-------
	v = BMI_GetPressure();
	BMI_GetValue("Pressure", v);
	d_ptr = (double*)BMI_GetValuePtr("Pressure");
	fprintf(stderr, "GetPressure \n";
	//-------
	status = BMI_SetPressure(v);
	BMI_SetValue("Pressure", v);
	fprintf(stderr, "SetPressure \n";
	//-------
	status = BMI_GetSaturationCalculated(v);
	BMI_GetValue("SaturationCalculated", v);
	d_ptr = (double*)BMI_GetValuePtr("SaturationCalculated");
	fprintf(stderr, "GetSaturationCalculated \n";
	//-------
	status = BMI_SetSaturationUser(v);
	BMI_SetValue("SaturationUser", v);
	fprintf(stderr, "SetSaturationUser \n";
	//-------
	v = BMI_GetSolutionVolume();
	BMI_GetValue("SolutionVolume", v);
	d_ptr = (double*)BMI_GetValuePtr("SolutionVolume");
	fprintf(stderr, "GetSolutionVolume \n";
	//-------
	status = BMI_GetSpeciesConcentrations(v);
	fprintf(stderr, "GetSpeciesConcentrations \n";
	//-------
	status = BMI_SpeciesConcentrations2Module(v);
	fprintf(stderr, "SpeciesConcentrations2Module \n";
	//-------
	status = BMI_GetSpeciesLog10Gammas(v);
	fprintf(stderr, "GetSpeciesLog10Gammas \n";
	//-------
	status = BMI_GetSpeciesLog10Molalities(v);
	fprintf(stderr, "GetSpeciesLog10Molalities \n";
	//-------
	v = BMI_GetTemperature();
	BMI_GetValue("Temperature", v);
	d_ptr = (double*)BMI_GetValuePtr("Temperature");
	fprintf(stderr, "GetTemperature \n";
	//-------
	status = BMI_SetTemperature(v);
	BMI_SetValue("Temperature", v);
	fprintf(stderr, "SetTemperature \n";
	//-------
	v = BMI_GetViscosity();
	BMI_GetValue("Viscosity", v);
	d_ptr = (double*)BMI_GetValuePtr("Viscosity");	
	fprintf(stderr, "GetViscosity \n";
	//
	// Take a time step
	//
	BMI_Update();      // void function
	fprintf(stderr, "Update\n";
	//-------
	status = BMI_RunCells();
	fprintf(stderr, "RunCells\n";
	//-------
	BMI_UpdateUntil(86400.0);      // void function
	fprintf(stderr, "UpdateUntil\n";
	//
	// Selected output
	//
	status = BMI_SetNthSelectedOutput(0);
	BMI_SetValue("NthSelectedOutput", 0);
	fprintf(stderr, "SetNthSelectedOutput \n";
	//-------
	int n_user = BMI_GetCurrentSelectedOutputUserNumber();
	BMI_GetValue("CurrentSelectedOutputUserNumber", n_user);
	fprintf(stderr, "GetCurrentSelectedOutputUserNumber \n";
	//-------
	n = BMI_GetNthSelectedOutputUserNumber(0);
	fprintf(stderr, "GetNthSelectedOutputUserNumber \n";
	//-------
	status = BMI_GetSelectedOutput(v);
	BMI_GetValue("SelectedOutput", v);
	fprintf(stderr, "GetSelectedOutput \n";
	//-------
	n = BMI_GetSelectedOutputColumnCount();
	BMI_GetValue("SelectedOutputColumnCount", n);
	fprintf(stderr, "GetSelectedOutputColumnCount \n";
	//-------
	n = BMI_GetSelectedOutputCount();
	BMI_GetValue("SelectedOutputCount", n);
	fprintf(stderr, "GetSelectedOutputCount \n";
	//-------
	status = BMI_GetSelectedOutputHeadings(str_vector);
	BMI_GetValue("SelectedOutputHeadings", str_vector);
	fprintf(stderr, "GetSelectedOutputHeadings \n";
	//-------
	bool b = BMI_GetSelectedOutputOn();
	BMI_GetValue("SelectedOutputOn", b);
	bool *b_ptr = (bool*)BMI_GetValuePtr("SelectedOutputOn");
	fprintf(stderr, "GetSelectedOutputOn \n";
	//-------
	n = BMI_GetSelectedOutputRowCount(); 
	BMI_GetValue("SelectedOutputRowCount", n);
	fprintf(stderr, "GetSelectedOutputRowCount \n";
	//-------
	status = BMI_SetCurrentSelectedOutputUserNumber(333);
	fprintf(stderr, "SetCurrentSelectedOutputUserNumber \n";
	//-------
	//
	// Getters
	// 
	std::vector< std::vector<int> > back_map = BMI_GetBackwardMapping();
	fprintf(stderr, "GetBackwardMapping \n";
	//-------
	std::string db_name = BMI_GetDatabaseFileName();
	fprintf(stderr, "GetDatabaseFileName \n";
	//-------
	vi = BMI_GetEndCell();
	fprintf(stderr, "GetEndCell\n";
	//-------
	n = BMI_GetErrorHandlerMode();
	fprintf(stderr, "GetErrorHandlerMode \n";
	//-------
	std::string str = BMI_GetErrorString();
	BMI_GetValue("ErrorString", str);
	fprintf(stderr, "GetErrorString \n";
	//-------
	str = BMI_GetFilePrefix();
	BMI_GetValue("FilePrefix", str);
	fprintf(stderr, "GetFilePrefix \n";
	//-------
	vi = BMI_GetForwardMapping();
	fprintf(stderr, "GetForwardMapping \n";
	//-------
	IPhreeqc* ipq = BMI_GetIPhreeqcPointer(0);
	fprintf(stderr, "GetIPhreeqcPointer \n";
	//-------
	n = BMI_GetMpiMyself();
	fprintf(stderr, "GetMpiMyself \n";
	//-------
	n = BMI_GetMpiTasks();
	fprintf(stderr, "GetMpiTasks \n";
	//-------
	b = BMI_GetPartitionUZSolids();
	fprintf(stderr, "GetPartitionUZSolids \n";
	//-------
	vi = BMI_GetPrintChemistryMask();
	fprintf(stderr, "GetPrintChemistryMask \n";
	//-------
	std::vector<bool> vb = BMI_GetPrintChemistryOn();
	fprintf(stderr, "GetPrintChemistryOn \n";
	//-------
	b = BMI_GetRebalanceByCell();
	fprintf(stderr, "GetRebalanceByCell \n";
	//-------
	d = BMI_GetRebalanceFraction();
	fprintf(stderr, "GetRebalanceFraction \n";
	//-------
	b = BMI_GetSpeciesSaveOn();
	fprintf(stderr, "GetSpeciesSaveOn \n";
	//-------
	std::vector< cxxNameDouble > s = BMI_GetSpeciesStoichiometry();
	fprintf(stderr, "GetSpeciesStoichiometry \n";
	//-------
	vi = BMI_GetStartCell();
	fprintf(stderr, "GetStartCell \n";
	//-------
	d = BMI_GetTimeConversion();
	fprintf(stderr, "GetTimeConversion \n";
	//-------
	n = BMI_GetUnitsExchange();
	fprintf(stderr, "GetUnitsExchange \n";
	//-------
	n = BMI_GetUnitsGasPhase();
	fprintf(stderr, "GetUnitsGasPhase \n";
	//-------
	n = BMI_GetUnitsKinetics();
	fprintf(stderr, "GetUnitsKinetics \n";
	//-------
	n = BMI_GetUnitsPPassemblage();
	fprintf(stderr, "GetUnitsPPassemblage \n";
	//-------
	n = BMI_GetUnitsSolution();
	fprintf(stderr, "GetUnitsSolution \n";
	//-------
	n = BMI_GetUnitsSSassemblage();
	fprintf(stderr, "GetUnitsSSassemblage \n";
	//-------
	n = BMI_GetUnitsSurface();
	fprintf(stderr, "GetUnitsSurface \n";
	//-------
	std::vector<IPhreeqcPhast *> w = BMI_GetWorkers();
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
	status = BMI_InitialPhreeqc2Concentrations(bc, vi);
	std::vector<double> tc(1, 30.0);
	std::vector<double> p_atm(1, 1.5);
	IPhreeqc* utility_ptr = BMI_Concentrations2Utility(bc, tc, p_atm);
	fprintf(stderr, "Concentrations2Utility \n";
	//-------
	BMI_DecodeError(-2);	         // void function
	fprintf(stderr, "DecodeError \n";
	//-------
	status = BMI_DumpModule(true);
	fprintf(stderr, "DumpModule \n";
	//-------
	BMI_ErrorHandler(0, "string"); // void function
	fprintf(stderr, "OK, just a test: ErrorHandler \n";
	//-------
	BMI_ErrorMessage("my error");  // void function
	fprintf(stderr, "OK, just a test: ErrorMessage \n";
	//-------
	BMI_LogMessage("Log message");  // void method
	fprintf(stderr, "LogMessage \n";
	//-------
	BMI_OutputMessage("Output message");  // void method
	fprintf(stderr, "OutputMessage \n";
	//-------
	BMI_ScreenMessage("Screen message\n");  // void method
	fprintf(stderr, "ScreenMessage \n";
	//-------
	status = BMI_StateSave(1);
	fprintf(stderr, "StateSave \n";
	//-------
	status = BMI_StateApply(1);
	fprintf(stderr, "StateApply \n";
	//-------
	status = BMI_StateDelete(1);
	fprintf(stderr, "StateDelete \n";
	//-------
	BMI_WarningMessage("Warning message");  // void method
	fprintf(stderr, "WarningMessage \n";
	//
	// BMI Methods
	//
	str = BMI_GetComponentName();
	fprintf(stderr, "GetComponentName \n";
	//-------
	d = BMI_GetCurrentTime();
	fprintf(stderr, "GetCurrentTime \n";
	//-------
	d = BMI_GetEndTime();
	fprintf(stderr, "GetEndTime \n";
	//-------
    n = BMI_GetGridRank(0);
	fprintf(stderr, "GetGridRank \n";
	//-------
    n = BMI_GetGridSize(0);
	fprintf(stderr, "GetGridSize \n";
	//-------
    std::string gtype = BMI_GetGridType(0);
	fprintf(stderr, "GetGridType \n";
	//-------
	n = BMI_GetInputItemCount();
	fprintf(stderr, "GetInputItemCount \n";
	//-------
	str_vector = BMI_GetInputVarNames();
	fprintf(stderr, "GetInputVarNames \n";
	//-------
	n = BMI_GetOutputItemCount();
	fprintf(stderr, "GetOutputItemCount \n";
	//-------
	str_vector = BMI_GetOutputVarNames();
	fprintf(stderr, "GetOutputVarNames \n";
	//-------
	d = BMI_GetTimeStep();
	fprintf(stderr, "GetTimeStep \n";
	//-------
	str = BMI_GetTimeUnits();
	fprintf(stderr, "GetTimeUnits \n";
	//-------
	BMI_GetValue("solution_saturation_index_Calcite", v);
	fprintf(stderr, "GetValue \n";
	//-------
	n = BMI_GetVarItemsize("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarItemsize \n";
	//-------
	n = BMI_GetVarNbytes("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarNbytes \n";
	//-------
	str = BMI_GetVarType("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarType \n";
	//-------
	str = BMI_GetVarUnits("solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarUnits \n";
	//-------
	//BMI_Initialize(YAML_filename);
	// See above
	BMI_SetValue("Time", 1.0);    // void method
	fprintf(stderr, "SetValue\n";
	//-------
	BMI_Update();    // void method
	fprintf(stderr, "Update\n";
	//-------	
 	BMI_UpdateUntil(864000.0);      // void function
	fprintf(stderr, "UpdateUntil\n";

	fprintf(stderr, "AddOutputVars\n";
	str_vector = BMI_GetOutputVarNames();
	for (size_t i = 0; i < str_vector.size(); i++)
	{
		int	itemsize = BMI_GetVarItemsize(str_vector[i]);
		int	nbytes = BMI_GetVarNbytes(str_vector[i]);
		std::string	vtype = BMI_GetVarType(str_vector[i]);
		if (itemsize == 0) itemsize++;
		if (nbytes == 0) nbytes++;
		int dim = nbytes / itemsize;
		if (vtype == "double")
		{
			if (dim == 1)
			{
				double dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<double> dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
		else if (vtype == "int")
		{
			if (dim == 1)
			{
				int dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<int> dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		} else if (vtype == "bool")
		{
			if (dim == 1)
			{
				bool dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
		}
		else if (vtype == "std::string")
		{
			if (dim == 1)
			{
				std::string dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<std::string> dest;
				BMI_GetValue(str_vector[i], dest);
				fprintf(stderr, "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
	}


	//-------	
	BMI_MpiWorkerBreak();
	BMI_Finalize();    // void method
	fprintf(stderr, "Finalize \n";
	//Should be private: status = BMI_ReturnHandler();
	//TODO status = BMI_MpiAbort();
	//TODO status = BMI_SetMpiWorkerCallbackC();
	//TODO status = BMI_SetMpiWorkerCallbackCookie();
	fprintf(stderr, "Success.\n");
	return;
#endif
}
