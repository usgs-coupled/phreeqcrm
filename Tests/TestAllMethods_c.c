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
	int* ic1 = NULL;
	int* ic2 = NULL;
	int* bc1 = NULL;
	int* bc2 = NULL;
	int* grid2chem = NULL;
	double* d_ptr = NULL;
	double* bc = NULL;
	double* v_nxyz = NULL;
	double* f1 = NULL;
	double d = 0.0;
	double* c = NULL;
	double* c_spec = NULL;
	int* vi = NULL;
	int i = 0, nchem = 0, nspec = 0, ncomps = 0;
	char string[MAX_LENGTH] = "", string1[MAX_LENGTH] = "", string2[MAX_LENGTH] = "";
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
	if (v_nxyz != NULL) free(v_nxyz);
	v_nxyz = (double*)malloc(nxyz * sizeof(double));
	for (int i = 0; i < nxyz; i++) v_nxyz[i] = 1.0;
	status = RM_SetRepresentativeVolume(id, v_nxyz);
	fprintf(stderr, "SetRepresentativeVolume %d\n", status);

	//-------Chemistry cells may be fewer than GridCellCount
	if (vi != NULL) free(vi);
	vi = (int*)malloc(nxyz * sizeof(int));
	for (int i = 0; i < nxyz; i++) v_nxyz[i] = 1;
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
	ncomps = RM_FindComponents(id);
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
	nspec = n;
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetSpeciesName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetSpeciesName %d\n", status);
	//-------
	status = RM_GetSpeciesD25(id, v_nxyz);
	fprintf(stderr, "GetSpeciesD25 %d\n", status);
	//-------
	status = RM_GetSpeciesZ(id, v_nxyz);
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
	//-------
	status = RM_GetGfw(id, v_nxyz);
	status = BMI_GetValueDouble(id, "Gfw", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "Gfw");
	//for (i = 0; i < ncomps; i++)
	//{
	//	fprintf(stderr, "     %10.2f\n", v_nxyz[i]);
	//}
	fprintf(stderr, "GetGfw %10.2f %10.2f \n ", v_nxyz[0], d_ptr[0]);
	//-------
	n = RM_GetKineticReactionsCount(id);
	fprintf(stderr, "GetKineticReactionsCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetKineticReactionsName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetKineticReactions \n");
	//-------
	n = RM_GetSICount(id);
	fprintf(stderr, "GetSICount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetSIName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetSIName \n");
	//-------
	n = RM_GetSolidSolutionComponentsCount(id);
	fprintf(stderr, "GetSolidSolutionComponentsCount %d\n", n);
	//-------	
	for (i = 0; i < n; i++)
	{
		status = RM_GetSolidSolutionComponentsName(id, i, string, MAX_LENGTH);
		status = RM_GetSolidSolutionName(id, i, string1, MAX_LENGTH);
		fprintf(stderr, "     %20s %20s\n", string, string1);
	}
	fprintf(stderr, "RM_GetSolidSolutionComponentsName \n");
	fprintf(stderr, "RM_GetSolidSolutionName \n");
	//-------
	n = RM_GetSurfaceSpeciesCount(id);
	fprintf(stderr, "GetSurfaceSpeciesCount %d\n", n);
	for (i = 0; i < n; i++)
	{
		status = RM_GetSurfaceSpeciesName(id, i, string, MAX_LENGTH);
		status = RM_GetSurfaceType(id, i, string1, MAX_LENGTH);
		status = RM_GetSurfaceName(id, i, string2, MAX_LENGTH);
		fprintf(stderr, "     %20s %20s %20s\n", string, string1, string2);
	}
	//-------
	fprintf(stderr, "GetSurfaceSpeciesName \n");
	fprintf(stderr, "GetSurfaceType \n");
	fprintf(stderr, "GetSurfaceName \n");
	//-------
	//
	// Remove any reactants in workers 
	// before populating cells with reactants
	//
	status = RM_RunString(id, 1, 0, 1, "DELETE; -all");
	fprintf(stderr, "RunString %d\n", status);
	//-------
	// 
	// Transfer initial conditions to model
	//
	if (vi != NULL) free(vi);
	vi = (int*)malloc(nxyz * sizeof(int));
	for (i = 0; i < nxyz; i++) vi[i] = 1;
	status = RM_InitialEquilibriumPhases2Module(id, vi);
	fprintf(stderr, "InitialEquilibriumPhases2Module %d\n", status);
	//-------
	status = RM_InitialExchanges2Module(id, vi);
	fprintf(stderr, "InitialExchanges2Module %d\n", status);
	//-------
	status = RM_InitialGasPhases2Module(id, vi);
	fprintf(stderr, "InitialGasPhases2Module %d\n", status);
	//-------
	status = RM_InitialKinetics2Module(id, vi);
	fprintf(stderr, "InitialKinetics2Module %d\n", status);
	//-------
	status = RM_InitialSolutions2Module(id, vi);
	fprintf(stderr, "InitialSolutions2Module %d\n", status);
	//-------
	status = RM_InitialSolidSolutions2Module(id, vi);
	fprintf(stderr, "InitialSolidSolutions2Module %d\n", status);
	//-------
	status = RM_InitialSurfaces2Module(id, vi);
	fprintf(stderr, "InitialSurfaces2Module %d\n", status);
	//-------
	// Alternative A.to the previous seven methods
	ic1 = (int*)malloc(nxyz * 7 * sizeof(double));
	ic2 = (int*)malloc(nxyz * 7 * sizeof(double));
	f1 = (double*)malloc(nxyz * 7 * sizeof(double));
	for (i = 0; i < 7 * nxyz;i++)
	{
		ic1[i] = 1;
		ic2[i] = -1;
		f1[i] = 1.0;
	}
	status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1);
	fprintf(stderr, "InitialPhreeqc2Module %d\n", status);
	//-------
	// Alternative C.to the previous seven methods, initialize cells 18 and 19
	vi[0] = 18;
	vi[1] = 19;
	status = RM_InitialPhreeqcCell2Module(id, 1, vi, 2);
	fprintf(stderr, "InitialPhreeqcCell2Module %d\n", status);
	//
	// Boundary conditions
	// 
	bc1 = (int*)malloc(1 * sizeof(int));
	bc2 = (int*)malloc(1 * sizeof(int));
	bc1[0] = 1;
	bc2[0] = -1;
	f1[0] = 1.0;
	bc = (double*)malloc(1 * ncomps * sizeof(double));
	status = RM_InitialPhreeqc2Concentrations(id, bc, 1, bc1, bc2, f1);
	//-------
	if (bc != NULL) free(bc);
	bc = (double*)malloc(1 * nspec * sizeof(int));
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc, 1, bc1, bc2, f1);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations %d\n", status);
	//
	// Get/Set methods for time steping
	//
	d = RM_GetTime(id);
	status = BMI_GetValueDouble(id, "Time", &d);
	d = BMI_GetCurrentTime(id);
    d = BMI_GetStartTime(id);
	d_ptr = (double*)BMI_GetValuePtr(id, "Time");
	fprintf(stderr, "GetTime %f\n", *d_ptr);
	//-------
	status = RM_SetTime(id, 0.0);
	status = BMI_SetValueDouble(id, "Time", 0.0);
	fprintf(stderr, "SetTime %d\n", status);
	//-------
	d = BMI_GetTimeStep(id);
	status = BMI_GetValueDouble(id, "TimeStep", &d);
	d_ptr = (double*)BMI_GetValuePtr(id, "TimeStep");
	fprintf(stderr, "GetTimeStep %f\n", *d_ptr);
	//-------
	status = RM_SetTimeStep(id, 0.0);
	status = BMI_SetValueDouble(id, "TimeStep", 0.0);
	fprintf(stderr, "SetTimeStep %d\n", status);
	//-------
	c = (double*)malloc(ncomps * nxyz * sizeof(double));
	status = RM_GetConcentrations(id, c);
	status = BMI_GetValueDouble(id, "Concentrations", c);
	d_ptr = (double*)BMI_GetValuePtr(id, "Concentrations");
	fprintf(stderr, "GetConcentrations %f\n", d_ptr[0]);
	//-------
	status = RM_SetConcentrations(id, c);
	status = BMI_SetValueDoubleArray(id, "Concentrations", c);
	fprintf(stderr, "SetConcentrations %d\n", status);
	//-------
	status = RM_GetDensityCalculated(id, v_nxyz);
	status = BMI_GetValueDouble(id, "DensityCalculated", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "DensityCalculated");
	fprintf(stderr, "GetDensityCalculated %f\n", d_ptr[0]);
	//-------
	status = RM_SetDensityUser(id, v_nxyz);
	status = BMI_SetValueDoubleArray(id,"DensityUser", v_nxyz);
	fprintf(stderr, "SetDensityUser %d\n", status);
	//-------
	status = RM_GetGasCompMoles(id, v_nxyz);
	fprintf(stderr, "GetGasCompMoles %f\n", v_nxyz[0]);
	//-------
	status = RM_SetGasCompMoles(id, v_nxyz);
	fprintf(stderr, "SetGasCompMoles %d\n", status);
	//-------
	status = RM_GetGasCompPhi(id, v_nxyz);
	fprintf(stderr, "GetGasCompPhi %f\n", v_nxyz[0]);
	//-------
	status = RM_GetGasCompPressures(id, v_nxyz);
	fprintf(stderr, "GetGasCompPressures %f\n", v_nxyz[0]);
	//-------
	status = RM_GetGasPhaseVolume(id, v_nxyz);
	fprintf(stderr, "GetGasPhaseVolume %f\n", v_nxyz[0]);
	//-------
	status = RM_SetGasPhaseVolume(id, v_nxyz);
	fprintf(stderr, "SetGasPhaseVolume %d\n", status);
	//-------
	for (size_t i = 0; i < RM_GetComponentCount(id); i++)
	{
		status = RM_GetIthConcentration(id, i, v_nxyz);
		//-------
		status = RM_SetIthConcentration(id, i, v_nxyz);
	}
	fprintf(stderr, "GetIthConcentration \n");
	fprintf(stderr, "SetIthConcentration %d\n", status);
	//-------
	for (int i = 0; i < RM_GetSpeciesCount(id); i++)
	{
		status = RM_GetIthSpeciesConcentration(id, i, v_nxyz);
		//-------
		status = RM_SetIthSpeciesConcentration(id, i, v_nxyz);
	}
	fprintf(stderr, "GetIthSpeciesConcentration \n");
	fprintf(stderr, "SetIthSpeciesConcentration %d\n", status);
	//-------
	status = RM_GetPorosity(id, v_nxyz);
	status = BMI_GetValueDouble(id, "Porosity", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "Porosity");
	fprintf(stderr, "GetPorosity %f\n", d_ptr[0]);
	//-------
	status = RM_SetPorosity(id, v_nxyz);
	status = BMI_SetValueDoubleArray(id, "Porosity", v_nxyz);
	fprintf(stderr, "SetPorosity %d\n", status);
	//-------
	status = RM_GetPressure(id, v_nxyz);
	status = BMI_GetValueDouble(id, "Pressure", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "Pressure");
	fprintf(stderr, "GetPressure %f\n", d_ptr[0]);
	//-------
	status = RM_SetPressure(id, v_nxyz);
	status = BMI_SetValueDoubleArray(id, "Pressure", v_nxyz);
	fprintf(stderr, "SetPressure %d\n", status);
	//-------
	status = RM_GetSaturationCalculated(id, v_nxyz);
	status = BMI_GetValueDouble(id, "SaturationCalculated", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "SaturationCalculated");
	fprintf(stderr, "GetSaturationCalculated %f\n", d_ptr[0]);
	//-------
	status = RM_SetSaturationUser(id, v_nxyz);
	status = BMI_SetValueDoubleArray(id, "SaturationUser", v_nxyz);
	fprintf(stderr, "SetSaturationUser %d\n", status);
	//-------
	status = RM_GetSolutionVolume(id, v_nxyz);
	status = BMI_GetValueDouble(id, "SolutionVolume", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr(id, "SolutionVolume");
	fprintf(stderr, "GetSolutionVolume %f\n", d_ptr[0]);
	//-------
	c_spec = (double*)malloc(nspec * nxyz * sizeof(double));
	status = RM_GetSpeciesConcentrations(id, c_spec);
	fprintf(stderr, "GetSpeciesConcentrations %f\n", c_spec[0]);
	//-------
	status = RM_SpeciesConcentrations2Module(id, c_spec);
	fprintf(stderr, "SpeciesConcentrations2Module %d\n", status);
	//-------
	status = RM_GetSpeciesLog10Gammas(id, v_nxyz);
	fprintf(stderr, "GetSpeciesLog10Gammas %f\n", v_nxyz[0]);
	//-------
	status = RM_GetSpeciesLog10Molalities(id, v_nxyz);
	fprintf(stderr, "GetSpeciesLog10Molalities %f\n", v_nxyz[0]);
	//-------

	fprintf(stderr, "Done.\n"); //==========================================================================================================
	return;
#ifdef SKIP
	v_nxyz = BMI_GetTemperature();
	BMI_GetValue("Temperature", v_nxyz);
	d_ptr = (double*)BMI_GetValuePtr("Temperature");
	fprintf(stderr, "GetTemperature \n";
	//-------
	status = BMI_SetTemperature(v_nxyz);
	BMI_SetValue("Temperature", v_nxyz);
	fprintf(stderr, "SetTemperature \n";
	//-------
	v_nxyz = BMI_GetViscosity();
	BMI_GetValue("Viscosity", v_nxyz);
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
	status = BMI_GetSelectedOutput(v_nxyz);
	BMI_GetValue("SelectedOutput", v_nxyz);
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
	BMI_GetValue("solution_saturation_index_Calcite", v_nxyz);
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
