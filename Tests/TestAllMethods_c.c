#ifdef USE_YAML
#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
#include "YAML_interface_C.h"
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
	int id = -1, id1 = -1;
	int* i_ptr = NULL;
	int* ic1 = NULL;
	int* ic2 = NULL;
	int* bc1 = NULL;
	int* bc2 = NULL;
	int* grid2chem = NULL;
	double* d_ptr = NULL;
	double* bc = NULL;
	double* bc_spec = NULL;
	double* v_nxyz = NULL;
	double* f1 = NULL;
	double d = 0.0, tc = 0.0, p_atm = 0.0;
	double* c = NULL;
	double* v_spec = NULL;
	double* v_spec_nxyz = NULL;
	double* v_gas = NULL;
	int* vi = NULL;
	int i = 0, nchem = 0, nspec = 0, ncomps = 0, ngas = 0, n_user = -1;
	int yid;
	int nrow = -1;
	double* so = NULL;
	char string[MAX_LENGTH] = "", string1[MAX_LENGTH] = "", string2[MAX_LENGTH] = "";
	char YAML_filename[MAX_LENGTH] = "TestAllMethods_c.yaml";

	yid = CreateYAMLPhreeqcRM();
	// Set GridCellCount
	status = YAMLSetGridCellCount(yid, nxyz);
	status = YAMLThreadCount(yid, 3);

#ifdef USE_MPI
	// MPI
	int mpi_myself;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}
	if (mpi_myself == 0)
	{
		status = WriteYAMLDoc(yid, YAML_filename);
	}
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)  //  make sure yaml is finished being written
	{
		exit(4);
	}

	id = RM_BmiCreate(nxyz, MPI_COMM_WORLD);
	if (mpi_myself > 0)
	{
		status = RM_MpiWorker(id);
		status = RM_BmiFinalize(id);
		return;
	}
#else
	// OpenMP or serial
	status = WriteYAMLDoc(yid, YAML_filename);
	nxyz = RM_GetGridCellCountYAML(YAML_filename);
	id = RM_BmiCreate(nxyz, nthreads);
#endif

	// delete YAMLPhreeqcRM by id
	DestroyYAMLPhreeqcRM(yid);

	// Use YAML file to initialize
	//RM_BmiInitialize(id, "");            // Initializes with no file
	RM_BmiInitialize(id, YAML_filename);   // Initializes with file
	RM_InitializeYAML(id, YAML_filename);  // Does not initialize BMI
	fprintf(stderr, "Initialize\n");

	//
	// Use all BMIPhreeqcRM methods roughly in order of use
	// 
	status = RM_BmiGetValueInt(id, "GridCellCount", &nxyz);
	nxyz = RM_GetGridCellCount(id);
	i_ptr = RM_BmiGetValuePtr(id, "GridCellCount");
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
	RM_BmiSetValueChar(id, "FilePrefix", "TestAllMethods_c");
	RM_BmiGetValueChar(id, "FilePrefix", string, MAX_LENGTH);
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
	RM_BmiSetValueInt(id, "SelectedOutputOn", 1);
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
	vi = (int*)malloc(nxyz * sizeof(int));
	for (int i = 0; i < nxyz; i++) vi[i] = 1;
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
	status = RM_BmiAddOutputVars(id, "AddOutputVars", "True");
	status = RM_BmiAddOutputVars(id, "SolutionProperties", "True");
	status = RM_BmiAddOutputVars(id, "SolutionTotalMolalities", "True");
	status = RM_BmiAddOutputVars(id, "ExchangeMolalities", "True");
	status = RM_BmiAddOutputVars(id, "SurfaceMolalities", "True");
	status = RM_BmiAddOutputVars(id, "EquilibriumPhases", "True");
	status = RM_BmiAddOutputVars(id, "Gases", "True");
	status = RM_BmiAddOutputVars(id, "KineticReactants", "True");
	status = RM_BmiAddOutputVars(id, "SolidSolutions", "True");
	status = RM_BmiAddOutputVars(id, "CalculateValues", "True");
	status = RM_BmiAddOutputVars(id, "SolutionActivities", "H+ Ca+2 Na+");
	status = RM_BmiAddOutputVars(id, "SolutionMolalities", "OH- Cl-");
	status = RM_BmiAddOutputVars(id, "SaturationIndices", "Calcite Dolomite");
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
	status = RM_BmiGetValueInt(id, "ComponentCount", &ncomps);
	i_ptr = (int*)RM_BmiGetValuePtr(id, "ComponentCount");
	fprintf(stderr, "GetComponentCount %d %d\n", ncomps, *i_ptr);
	//-------
	for (i = 0; i < ncomps; i++)
	{
		status = RM_GetComponent(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	// BMI version of GetComponents
	int itemsize = RM_BmiGetVarItemsize(id, "Components");
	int nbytes = RM_BmiGetVarNbytes(id, "Components");
	char* buffer = (char*)malloc(((size_t)nbytes + 1) * sizeof(char));
	status = RM_BmiGetValueChar(id, "Components", buffer, nbytes + 1);
	char** comp_list = (char**)malloc(ncomps * sizeof(char*));
	for (i = 0; i < ncomps; i++)
	{
		comp_list[i] = (char*)malloc(((size_t)itemsize + 1) * sizeof(char));
		memcpy(comp_list[i], &buffer[i * itemsize], (size_t)itemsize);
		comp_list[i][itemsize] = '\0';
	}
	fprintf(stderr, "RM_BmiGetValue(Components) %d\n", status);
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
	v_spec = (double*)malloc((size_t)nspec * sizeof(double));
	status = RM_GetSpeciesD25(id, v_spec);
	fprintf(stderr, "GetSpeciesD25 %d\n", status);
	//-------
	status = RM_GetSpeciesZ(id, v_spec);
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
	ngas = n;
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_GetGasComponentsName(id, i, string, MAX_LENGTH);
		fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetGasComponents %d\n", status);
	//-------
	status = RM_GetGfw(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "Gfw", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Gfw");
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
	ic1 = (int*)malloc((size_t)nxyz * 7 * sizeof(double));
	ic2 = (int*)malloc((size_t)nxyz * 7 * sizeof(double));
	f1 = (double*)malloc((size_t)nxyz * 7 * sizeof(double));
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
	bc1 = (int*)malloc(2 * sizeof(int));
	bc2 = (int*)malloc(2 * sizeof(int));
	bc1[0] = 1;
	bc2[0] = -1;
	f1[0] = 1.0;
	bc = (double*)malloc(1 * (size_t)ncomps * sizeof(double));
	status = RM_InitialPhreeqc2Concentrations(id, bc, 1, bc1, bc2, f1);
	//-------
	bc_spec = (double*)malloc(1 * (size_t)nspec * sizeof(double));
	status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_spec, 1, bc1, bc2, f1);
	fprintf(stderr, "InitialPhreeqc2SpeciesConcentrations %d\n", status);
	//
	// Get/Set methods for time steping
	//
	d = RM_GetTime(id);
	status = RM_BmiGetValueDouble(id, "Time", &d);
	d = RM_BmiGetCurrentTime(id);
	d = RM_BmiGetStartTime(id);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Time");
	fprintf(stderr, "GetTime %f\n", *d_ptr);
	//-------
	status = RM_SetTime(id, 0.0);
	status = RM_BmiSetValueDouble(id, "Time", 0.0);
	fprintf(stderr, "SetTime %d\n", status);
	//-------
	d = RM_BmiGetTimeStep(id);
	status = RM_BmiGetValueDouble(id, "TimeStep", &d);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "TimeStep");
	fprintf(stderr, "GetTimeStep %f\n", *d_ptr);
	//-------
	status = RM_SetTimeStep(id, 0.0);
	status = RM_BmiSetValueDouble(id, "TimeStep", 0.0);
	fprintf(stderr, "SetTimeStep %d\n", status);
	//-------
	c = (double*)malloc((size_t)ncomps * (size_t)nxyz * sizeof(double));
	status = RM_GetConcentrations(id, c);
	status = RM_BmiGetValueDouble(id, "Concentrations", c);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Concentrations");
	fprintf(stderr, "GetConcentrations %f\n", d_ptr[0]);
	//-------
	status = RM_SetConcentrations(id, c);
	status = RM_BmiSetValueDoubleArray(id, "Concentrations", c);
	fprintf(stderr, "SetConcentrations %d\n", status);
	//-------
	status = RM_GetDensityCalculated(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "DensityCalculated", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "DensityCalculated");
	fprintf(stderr, "GetDensityCalculated %f\n", d_ptr[0]);
	//-------
	status = RM_SetDensityUser(id, v_nxyz);
	status = RM_BmiSetValueDoubleArray(id, "DensityUser", v_nxyz);
	fprintf(stderr, "SetDensityUser %d\n", status);
	//-------
	v_gas = (double*)malloc((size_t)ngas * (size_t)nxyz * sizeof(double));
	status = RM_GetGasCompMoles(id, v_gas);
	fprintf(stderr, "GetGasCompMoles %f\n", v_gas[0]);
	//-------

	status = RM_SetGasCompMoles(id, v_gas);
	fprintf(stderr, "SetGasCompMoles %d\n", status);
	//-------
	status = RM_GetGasCompPhi(id, v_gas);
	fprintf(stderr, "GetGasCompPhi %f\n", v_nxyz[0]);
	//-------
	status = RM_GetGasCompPressures(id, v_gas);
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
	status = RM_BmiGetValueDouble(id, "Porosity", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Porosity");
	fprintf(stderr, "GetPorosity %f\n", d_ptr[0]);
	//-------
	status = RM_SetPorosity(id, v_nxyz);
	status = RM_BmiSetValueDoubleArray(id, "Porosity", v_nxyz);
	fprintf(stderr, "SetPorosity %d\n", status);
	//-------
	status = RM_GetPressure(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "Pressure", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Pressure");
	fprintf(stderr, "GetPressure %f\n", d_ptr[0]);
	//-------
	status = RM_SetPressure(id, v_nxyz);
	status = RM_BmiSetValueDoubleArray(id, "Pressure", v_nxyz);
	fprintf(stderr, "SetPressure %d\n", status);
	//-------
	status = RM_GetSaturationCalculated(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "SaturationCalculated", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "SaturationCalculated");
	fprintf(stderr, "GetSaturationCalculated %f\n", d_ptr[0]);
	//-------
	status = RM_SetSaturationUser(id, v_nxyz);
	status = RM_BmiSetValueDoubleArray(id, "SaturationUser", v_nxyz);
	fprintf(stderr, "SetSaturationUser %d\n", status);
	//-------
	status = RM_GetSolutionVolume(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "SolutionVolume", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "SolutionVolume");
	fprintf(stderr, "GetSolutionVolume %f\n", d_ptr[0]);
	//-------
	v_spec_nxyz = (double*)malloc((size_t)nxyz * (size_t)nspec * sizeof(double));
	status = RM_GetSpeciesConcentrations(id, v_spec_nxyz);
	fprintf(stderr, "GetSpeciesConcentrations %f\n", v_spec_nxyz[0]);
	//-------
	status = RM_SpeciesConcentrations2Module(id, v_spec_nxyz);
	fprintf(stderr, "SpeciesConcentrations2Module %d\n", status);
	//-------
	status = RM_GetSpeciesLog10Gammas(id, v_spec_nxyz);
	fprintf(stderr, "GetSpeciesLog10Gammas %f\n", v_spec_nxyz[0]);
	//-------
	status = RM_GetSpeciesLog10Molalities(id, v_spec_nxyz);
	fprintf(stderr, "GetSpeciesLog10Molalities %f\n", v_spec_nxyz[0]);
	//-------
	status = RM_GetTemperature(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "Temperature", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Temperature");
	fprintf(stderr, "GetTemperature %f\n", d_ptr[0]);
	//-------
	status = RM_SetTemperature(id, v_nxyz);
	status = RM_BmiSetValueDoubleArray(id, "Temperature", v_nxyz);
	fprintf(stderr, "SetTemperature %d\n", status);
	//-------
	status = RM_GetViscosity(id, v_nxyz);
	status = RM_BmiGetValueDouble(id, "Viscosity", v_nxyz);
	d_ptr = (double*)RM_BmiGetValuePtr(id, "Viscosity");
	fprintf(stderr, "GetViscosity %f\n", d_ptr[0]);
	//
	// Take a time step
	//
	status = RM_BmiUpdate(id);      // void function
	fprintf(stderr, "Update %d\n", status);
	//-------
	status = RM_RunCells(id);
	fprintf(stderr, "RunCells %d\n", status);
	//-------
	status = RM_BmiUpdateUntil(id, 86400.0);      // void function
	fprintf(stderr, "UpdateUntil %d\n", status);
	//
	// Selected output
	//
	status = RM_SetNthSelectedOutput(id, 0);
	status = RM_BmiSetValueInt(id, "NthSelectedOutput", 0);
	fprintf(stderr, "SetNthSelectedOutput %d\n", status);
	//-------
	n_user = RM_GetCurrentSelectedOutputUserNumber(id);
	status = RM_BmiGetValueInt(id, "CurrentSelectedOutputUserNumber", &n_user);
	fprintf(stderr, "GetCurrentSelectedOutputUserNumber %d\n", n_user);
	//-------
	n = RM_GetNthSelectedOutputUserNumber(id, 0);
	fprintf(stderr, "GetNthSelectedOutputUserNumber %d\n", n);
	//-------
	n = RM_GetSelectedOutputColumnCount(id);
	so = (double*)malloc((size_t)n * (size_t)nxyz * sizeof(double));
	status = RM_GetSelectedOutput(id, so);
	status = RM_BmiGetValueDouble(id, "SelectedOutput", so);
	fprintf(stderr, "GetSelectedOutput %f\n", so[0]);
	//-------
	n = RM_GetSelectedOutputColumnCount(id);
	status = RM_BmiGetValueInt(id, "SelectedOutputColumnCount", &n);
	fprintf(stderr, "GetSelectedOutputColumnCount %d\n", n);
	//-------
	n = RM_GetSelectedOutputCount(id);
	status = RM_BmiGetValueInt(id, "SelectedOutputCount", &n);
	fprintf(stderr, "GetSelectedOutputCount %d\n", n);
	//-------
	n = RM_GetSelectedOutputColumnCount(id);
	for (i = 0; i < n; i++)
	{
		status = RM_GetSelectedOutputHeading(id, i, string, MAX_LENGTH);
		//fprintf(stderr, "     %s\n", string);
	}
	fprintf(stderr, "GetSelectedOutputHeading %d \n", status);
	//-------
	// BMI version of GetHeadings
	itemsize = RM_BmiGetVarItemsize(id, "SelectedOutputHeadings");
	nbytes = RM_BmiGetVarNbytes(id, "SelectedOutputHeadings");
	int ncol = 0;
	status = RM_BmiGetValueInt(id, "SelectedOutputColumnCount", &ncol);
	char* so_buffer = (char*)malloc(((size_t)nbytes + 1) * sizeof(char));
	status = RM_BmiGetValueChar(id, "SelectedOutputHeadings", so_buffer, nbytes + 1);
	char** heading_list = (char**)malloc(ncol * itemsize* sizeof(char*));
	for (i = 0; i < ncol; i++)
	{
		heading_list[i] = (char*)malloc(((size_t)itemsize + 1) * sizeof(char));
		memcpy(heading_list[i], &so_buffer[i * itemsize], (size_t)itemsize);
		heading_list[i][itemsize] = '\0';
	}
	fprintf(stderr, "RM_BmiGetValueChar(SelectedOutputHeadings %d \n", status);
	//-------
	//i = RM_GetSelectedOutputOn(id);
	status = RM_BmiGetValueInt(id, "SelectedOutputOn", &i);
	i_ptr = (int*)RM_BmiGetValuePtr(id, "SelectedOutputOn");
	fprintf(stderr, "GetSelectedOutputOn %d\n", i_ptr[0]);
	//-------
	n = RM_GetSelectedOutputRowCount(id);
	status = RM_BmiGetValueInt(id, "SelectedOutputRowCount", &n);
	fprintf(stderr, "GetSelectedOutputRowCount %d\n", n);
	//-------
	status = RM_SetCurrentSelectedOutputUserNumber(id, 333);
	fprintf(stderr, "SetCurrentSelectedOutputUserNumber %d\n", status);
	//-------
	//
	// Getters
	// 
	for (i = 0; i < nchem; i++)
	{
		n = nxyz;
		status = RM_GetBackwardMapping(id, i, vi, &n);
	}
	fprintf(stderr, "GetBackwardMapping %d\n", status);
	//-------
	// Not implemented
	//status = RM_GetDatabaseFileName(id, string);
	//fprintf(stderr, "GetDatabaseFileName %s\n", string);
	//-------
	status = RM_GetEndCell(id, vi);
	fprintf(stderr, "GetEndCell %d\n", vi[0]);
	//-------
	// Not implemented
	//status = RM_GetErrorHandlerMode(id, &n);
	//fprintf(stderr, "GetErrorHandlerMode %d\n", n);
	//-------
	status = RM_GetErrorString(id, string, MAX_LENGTH);
	status = RM_BmiGetValueChar(id, "ErrorString", string, MAX_LENGTH);
	fprintf(stderr, "GetErrorString %d\n", status);
	//-------
	status = RM_GetFilePrefix(id, string, MAX_LENGTH);
	status = RM_BmiGetValueChar(id, "FilePrefix", string, MAX_LENGTH);
	fprintf(stderr, "GetFilePrefix %s\n", string);
	//-------
	// Not implemented
	//status = RM_GetForwardMapping(id, vi);
	//fprintf(stderr, "GetForwardMapping %d\n", vi[0]);
	//-------
	n = RM_GetIPhreeqcId(id, 0);
	fprintf(stderr, "GetIPhreeqcPointer %d\n", n);
	//-------
	n = RM_GetMpiMyself(id);
	fprintf(stderr, "GetMpiMyself %d\n", n);
	//-------
	n = RM_GetMpiTasks(id);
	fprintf(stderr, "GetMpiTasks %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetPartitionUZSolids(id);
	//fprintf(stderr, "GetPartitionUZSolids %d\n", n);
	//-------
	// Not implemented
	//status = RM_GetPrintChemistryMask(id, vi);
	//fprintf(stderr, "GetPrintChemistryMask %d\n", vi[0]);
	//-------
	// Not implemented
	//n = RM_GetPrintChemistryOn(id);
	//fprintf(stderr, "GetPrintChemistryOn %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetRebalanceByCell(id);
	//fprintf(stderr, "GetRebalanceByCell %d\n", n);
	//-------
	// Not implemented
	//d = RM_GetRebalanceFraction(id);
	//fprintf(stderr, "GetRebalanceFraction %f\n", d);
	//-------
	// Not implemented
	//n = RM_GetSpeciesSaveOn();
	//fprintf(stderr, "GetSpeciesSaveOn %d\n", n);
	//-------
	// Not implemented
	//RM_GetSpeciesStoichiometry();
	//fprintf(stderr, "GetSpeciesStoichiometry \n");
	//-------
	status = RM_GetStartCell(id, vi);
	fprintf(stderr, "GetStartCell %d\n", vi[0]);
	//-------
	d = RM_GetTimeConversion(id);
	fprintf(stderr, "GetTimeConversion %f\n", d);
	//-------
	// Not implemented
	//n = RM_GetUnitsExchange(id);
	//fprintf(stderr, "GetUnitsExchange %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetUnitsGasPhase(id);
	//fprintf(stderr, "GetUnitsGasPhase %d\n", n);
	//-------
	//// Not implemented
	//n = RM_GetUnitsKinetics(id);
	//fprintf(stderr, "GetUnitsKinetics %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetUnitsPPassemblage(id);
	//fprintf(stderr, "GetUnitsPPassemblage %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetUnitsSolution(id);
	//fprintf(stderr, "GetUnitsSolution %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetUnitsSSassemblage(id);
	//fprintf(stderr, "GetUnitsSSassemblage %d\n", n);
	//-------
	// Not implemented
	//n = RM_GetUnitsSurface(id);
	//fprintf(stderr, "GetUnitsSurface %d\n", n);
	//-------
	// Not implemented
	//std::vector<IPhreeqcPhast *> w = RM_BmiGetWorkers();
	//fprintf(stderr, "GetWorkers \n";
	//
	// Utilities
	//
#ifndef USE_MPI
	// Make another instance
	id1 = RM_Create(10, 1);
	fprintf(stderr, "Make a new instance. %d\n", id1);
	//-------
	status = RM_CloseFiles(id1); 
	fprintf(stderr, "CloseFiles %d\n", status);
	//-------
	status = RM_BmiFinalize(id1);
#endif
	//-------
	vi[0] = 1;
	status = RM_InitialPhreeqc2Concentrations(id, bc, 1, bc1, bc2, f1);
	tc = 30;
	p_atm = 1.5;
	n = RM_Concentrations2Utility(id, bc, 1, &tc, &p_atm);
	fprintf(stderr, "Concentrations2Utility %d\n", n);
	//-------
	status = RM_DecodeError(id, -2);	         // void function
	fprintf(stderr, "DecodeError %d\n", status);
	//-------
	status = RM_DumpModule(id, 1, 0);
	fprintf(stderr, "DumpModule %d\n", status);
	//-------
	// Not implemented
	//status = RM_ErrorHandler(id, 0, "string"); // void function
	//fprintf(stderr, "OK, just a test: ErrorHandler \n");
	//-------
	status = RM_ErrorMessage(id, "My error message");  
	fprintf(stderr, "OK, just a test: ErrorMessage %d\n", status);
	//-------
	status = RM_LogMessage(id, "My Log message"); 
	fprintf(stderr, "LogMessage  %d\n", status);
	//-------
	status = RM_OutputMessage(id, "My Output message");  // void method
	fprintf(stderr, "OutputMessage  %d\n", status);
	//-------
	status = RM_ScreenMessage(id, "My Screen message\n");  // void method
	fprintf(stderr, "ScreenMessage  %d\n", status);
	//-------
	status = RM_WarningMessage(id, "My warning message");  // void method
	fprintf(stderr, "WarningMessage %d\n", status);
	//-------
	status = RM_StateSave(id, 1);
	fprintf(stderr, "StateSave %d\n", status);
	//-------
	status = RM_StateApply(id, 1);
	fprintf(stderr, "StateApply %d\n", status);
	//-------
	status = RM_StateDelete(id, 1);
	fprintf(stderr, "StateDelete %d\n", status);
	//
	// BMI Methods
	//
	status = RM_BmiGetComponentName(id, string, MAX_LENGTH);
	fprintf(stderr, "GetComponentName %s\n", string);
	//-------
	d = RM_BmiGetCurrentTime(id);
	fprintf(stderr, "GetCurrentTime %f\n", d);
	//-------
	d = RM_BmiGetEndTime(id);
	fprintf(stderr, "GetEndTime %f\n", d);
	//-------
    n = RM_BmiGetGridRank(id, 0);
	fprintf(stderr, "GetGridRank %d\n", n);
	//-------
    n = RM_BmiGetGridSize(id, 0);
	fprintf(stderr, "GetGridSize %d\n", n);
	//-------
    status = RM_BmiGetGridType(id, 0, string, MAX_LENGTH);
	fprintf(stderr, "GetGridType %s\n", string);
	//-------
	n = RM_BmiGetInputItemCount(id);
	fprintf(stderr, "GetInputItemCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_BmiGetInputVarName(id, i, string, MAX_LENGTH);
		//fprintf(stderr, "RM_BmiGetInputVarName %s\n", string);
	}
	fprintf(stderr, "GetInputVarNames %d\n", status);
	//-------
	n = RM_BmiGetOutputItemCount(id);
	fprintf(stderr, "GetOutputItemCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_BmiGetOutputVarName(id, i, string, MAX_LENGTH);
		//fprintf(stderr, "RM_BmiGetOutputVarName %s\n", string);
	}
	fprintf(stderr, "GetOutputVarNames %d\n", n);
	//-------
	n = RM_BmiGetPointableItemCount(id);
	fprintf(stderr, "RM_BmiGetPointableItemCount %d\n", n);
	//-------
	for (i = 0; i < n; i++)
	{
		status = RM_BmiGetPointableVarName(id, i, string, MAX_LENGTH);
		//fprintf(stderr, "RM_BmiGetPointableVarName %s\n", string);
	}
	fprintf(stderr, "GetPointableVarName %d\n", n);
	//-------
	d = RM_BmiGetTimeStep(id);
	fprintf(stderr, "GetTimeStep %f\n", d);
	//-------
	status = RM_BmiGetTimeUnits(id, string, MAX_LENGTH);
	fprintf(stderr, "GetTimeUnits %s\n", string);
	//-------
	status = RM_BmiGetValueDouble(id, "solution_saturation_index_Calcite", v_nxyz);
	fprintf(stderr, "GetValue solution_saturation_index_Calcite %f\n", v_nxyz[0]);
	//-------
	n = RM_BmiGetVarItemsize(id, "solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarItemsize %d\n", n);
	//-------
	n = RM_BmiGetVarNbytes(id, "solution_saturation_index_Calcite");
	fprintf(stderr, "GetVarNbytes %d\n", n);
	//-------
	status = RM_BmiGetVarType(id, "solution_saturation_index_Calcite", string, MAX_LENGTH);
	fprintf(stderr, "GetVarType %s\n", string);
	//-------
	status = RM_BmiGetVarUnits(id, "solution_saturation_index_Calcite", string, MAX_LENGTH);
	fprintf(stderr, "GetVarUnits %s\n", string);
	//-------
	//RM_BmiInitialize(YAML_filename);
	// See above
	status = RM_BmiSetValueDouble(id, "Time", 1.0);    // void method
	fprintf(stderr, "SetValue %d\n", status);
	//-------
	status = RM_BmiUpdate(id);    
	fprintf(stderr, "Update %d\n", status);
	//-------	
 	status = RM_BmiUpdateUntil(id, 864000.0);      
	fprintf(stderr, "UpdateUntil %d\n", status);

	status = RM_BmiGetVarType(id, "SelectedOutputOn", string, MAX_LENGTH);
	fprintf(stderr, "SelectedOutputOn type: %s\n", string);
	//-------	
#ifdef USE_MPI
	status = RM_MpiWorkerBreak(id);
	fprintf(stderr, "RM_MpiWorkerBreak %d\n", status);
#endif
	status = RM_BmiFinalize(id);    // void method
	fprintf(stderr, "Finalize %d\n", status);
	//Should be private: status = RM_BmiReturnHandler();
	//TODO status = RM_BmiMpiAbort();
	//TODO status = RM_BmiSetMpiWorkerCallbackC();
	//TODO status = RM_BmiSetMpiWorkerCallbackCookie();
	free(grid2chem);
	free(v_nxyz);
	free(vi);
	free(v_spec);
	free(ic1);
	free(ic2);
	free(f1);
	free(bc1);
	free(bc2);
	free(bc);
	free(bc_spec);
	free(c);
	free(v_gas);
	free(v_spec_nxyz);
	free(so);
	for (i = 0; i < ncomps; i++)
	{
		free(comp_list[i]);
	}
	free(comp_list);
	free(buffer);
	for (i = 0; i < ncol; i++)
	{
		free(heading_list[i]); 
	}
	free(heading_list);
	free(so_buffer);
	fprintf(stderr, "Success.\n");
	return;
}
#endif // YAML
