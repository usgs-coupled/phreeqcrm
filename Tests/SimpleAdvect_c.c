#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
#include "IPhreeqc.h"

void simpleadvection_c(double* c, double* bc_conc, int ncomps, int nxyz, int dim);

size_t strcat_safe(char* dest, size_t max, const char* src);
size_t strcpy_safe(char* dest, size_t max, const char* src);
void SimpleAdvect_c()
{
	// Based on PHREEQC Example 11

	int mpi_myself = 0;
	int i;
	int nxyz;
#ifndef USE_MPI
	int nthreads;
#else
	MPI_Comm comm;
#endif
	int id;
	int status;
	double* por;
	int* print_chemistry_mask;
	int nchem;
	char str[100];
	int ncomps;
	char** components;
	int* ic1;
	int* ic2;
	double* f1;
	int nbound;
	int* bc1;
	int* bc2;
	double* bc_f1;
	double* bc_conc;
	double* c;
	double time, time_step;
	double* temperature;
	double* pressure;
	int isteps, nsteps;
	// --------------------------------------------------------------------------
	// Create PhreeqcRM
	// --------------------------------------------------------------------------
	nxyz = 20;
#ifdef USE_MPI
	// MPI
	comm = MPI_COMM_WORLD;
	id = RM_Create(nxyz, comm);
	if (MPI_Comm_rank(comm, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}
	if (mpi_myself > 0)
	{
		status = RM_MpiWorker(id);
		status = RM_Destroy(id);
		return;
	}
#else
	// OpenMP
	nthreads = 3;
	id = RM_Create(nxyz, nthreads);
#endif
	// Set properties
	status = RM_SetComponentH2O(id, 0);
	status = RM_UseSolutionDensityVolume(id, 0);
	// Open error, log, and output files
	status = RM_SetFilePrefix(id, "SimpleAdvect_c");
	status = RM_OpenFiles(id);
	// Set concentration units
	status = RM_SetUnitsSolution(id, 2);      // 1, mg/L; 2, mol/L; 3, kg/kgs
	status = RM_SetUnitsExchange(id, 1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	// Set conversion from seconds to days
	status = RM_SetTimeConversion(id, 1.0 / 86400.0);
	// Set initial porosity
	por = (double*)malloc((size_t)(nxyz * sizeof(double)));
	for (i = 0; i < nxyz; i++) por[i] = 0.2;
	status = RM_SetPorosity(id, por);
	// Set cells to print chemistry when print chemistry is turned on
	print_chemistry_mask = (int*)malloc((size_t)(nxyz * sizeof(int)));
	for (i = 0; i < nxyz; i++)
	{
		print_chemistry_mask[i] = 1;
	}
	status = RM_SetPrintChemistryMask(id, print_chemistry_mask);
	nchem = RM_GetChemistryCellCount(id);
	// --------------------------------------------------------------------------
	// Set initial conditions
	// --------------------------------------------------------------------------
	// Set printing of chemistry file
	status = RM_SetPrintChemistryOn(id, 0, 1, 0); // workers, initial_phreeqc, utility
	// Set printing of chemistry file
	status = RM_LoadDatabase(id, "phreeqc.dat");
	// Run file to define solutions and reactants for initial conditions, selected output
	// There are three types of IPhreeqc instances in PhreeqcRM
	// Argument 1 refers to the workers for doing reaction calculations for transport
	// Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
	// Argument 3 refers to the Utility instance
	status = RM_RunFile(id, 1, 1, 1, "advect.pqi");
	// Clear contents of workers and utility
	strcpy_safe(str, 100, "DELETE; -all");
	status = RM_RunString(id, 1, 0, 1, str);	// workers, initial_phreeqc, utility 
	// Determine number of components to transport
	ncomps = RM_FindComponents(id);
	// Get component information
	components = (char**)malloc((size_t)(ncomps * sizeof(char*)));
	if (components == NULL) exit(4);
	for (i = 0; i < ncomps; i++)
	{
		components[i] = (char*)malloc((size_t)(100 * sizeof(char*)));
		if (components[i] == NULL) exit(4);
		status = RM_GetComponent(id, i, components[i], 100);
		snprintf(str, sizeof(str), "%10s\n", components[i]);
		status = RM_OutputMessage(id, str);
	}
	status = RM_OutputMessage(id, "\n");
	// Set array of initial conditions
	ic1 = (int*)malloc((size_t)(7 * nxyz * sizeof(int)));
	ic2 = (int*)malloc((size_t)(7 * nxyz * sizeof(int)));
	f1 = (double*)malloc((size_t)(7 * nxyz * sizeof(double)));
	if (ic1 == NULL || ic2 == NULL || f1 == NULL) exit(4);
	for (i = 0; i < nxyz; i++)
	{
		ic1[i] = 1;       // Solution 1
		ic1[nxyz + i] = -1;      // Equilibrium phases none
		ic1[2 * nxyz + i] = 1;       // Exchange 1
		ic1[3 * nxyz + i] = -1;      // Surface none
		ic1[4 * nxyz + i] = -1;      // Gas phase none
		ic1[5 * nxyz + i] = -1;      // Solid solutions none
		ic1[6 * nxyz + i] = -1;      // Kinetics none
		ic2[i] = -1;      // Solution none
		ic2[nxyz + i] = -1;      // Equilibrium phases none
		ic2[2 * nxyz + i] = -1;      // Exchange none
		ic2[3 * nxyz + i] = -1;      // Surface none
		ic2[4 * nxyz + i] = -1;      // Gas phase none
		ic2[5 * nxyz + i] = -1;      // Solid solutions none
		ic2[6 * nxyz + i] = -1;      // Kinetics none
		f1[i] = 1.0;      // Mixing fraction ic1 Solution
		f1[nxyz + i] = 1.0;      // Mixing fraction ic1 Equilibrium phases 
		f1[2 * nxyz + i] = 1.0;      // Mixing fraction ic1 Exchange 1
		f1[3 * nxyz + i] = 1.0;      // Mixing fraction ic1 Surface 
		f1[4 * nxyz + i] = 1.0;      // Mixing fraction ic1 Gas phase 
		f1[5 * nxyz + i] = 1.0;      // Mixing fraction ic1 Solid solutions 
		f1[6 * nxyz + i] = 1.0;      // Mixing fraction ic1 Kinetics 
	}
	status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1);
	// Initial equilibration of cells
	time = 0.0;
	time_step = 0.0;
	c = (double*)malloc((size_t)(ncomps * nxyz * sizeof(double)));
	status = RM_SetTime(id, time);
	status = RM_SetTimeStep(id, time_step);
	status = RM_RunCells(id);
	status = RM_GetConcentrations(id, c);
	// --------------------------------------------------------------------------
	// Set boundary condition
	// --------------------------------------------------------------------------
	nbound = 1;
	bc1 = (int*)malloc((size_t)(nbound * sizeof(int)));
	bc2 = (int*)malloc((size_t)(nbound * sizeof(int)));
	bc_f1 = (double*)malloc((size_t)(nbound * sizeof(double)));
	bc_conc = (double*)malloc((size_t)(ncomps * nbound * sizeof(double)));
	if (bc1 == NULL || bc2 == NULL || bc_f1 == NULL || bc_conc == NULL) exit(4);
	for (i = 0; i < nbound; i++)
	{
		bc1[i] = 0;       // Solution 0 from Initial IPhreeqc instance
		bc2[i] = -1;      // no bc2 solution for mixing
		bc_f1[i] = 1.0;     // mixing fraction for bc1
	}
	status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);
	// --------------------------------------------------------------------------
	// Transient loop
	// --------------------------------------------------------------------------
	nsteps = 10;
	pressure = (double*)malloc((size_t)(nxyz * sizeof(double)));
	temperature = (double*)malloc((size_t)(nxyz * sizeof(double)));
	if (temperature == NULL || pressure == NULL) exit(4);
	for (i = 0; i < nxyz; i++)
	{
		pressure[i] = 2.0;
		temperature[i] = 20.0;
	}
	status = RM_SetPressure(id, pressure);
	status = RM_SetTemperature(id, temperature);
	time_step = 86400;
	status = RM_SetTimeStep(id, time_step);
	for (isteps = 0; isteps < nsteps; isteps++)
	{
		// Advection calculation
		snprintf(str, sizeof(str), "%s%10.1f%s", "Beginning transport calculation      ",
			RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
		status = RM_LogMessage(id, str);
		status = RM_SetScreenOn(id, 1);
		status = RM_ScreenMessage(id, str);
		snprintf(str, sizeof(str), "%s%10.1f%s", "          Time step                  ",
			RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days\n");
		status = RM_LogMessage(id, str);
		status = RM_ScreenMessage(id, str);
		// Advect one step
		simpleadvection_c(c, bc_conc, ncomps, nxyz, nbound);
		// Transfer data to PhreeqcRM for reactions
		status = RM_SetConcentrations(id, c);          // Transported concentrations
		status = RM_SetTimeStep(id, time_step);        // Time step for kinetic reactions
		time = time + time_step;
		status = RM_SetTime(id, time);                 // Current time
		// Set print flag
		if (isteps == nsteps - 1)
		{
			status = RM_SetSelectedOutputOn(id, 1);       // enable selected output
			status = RM_SetPrintChemistryOn(id, 1, 0, 0); // print at last time step, workers, initial_phreeqc, utility
		}
		else
		{
			status = RM_SetSelectedOutputOn(id, 0);       // disable selected output
			status = RM_SetPrintChemistryOn(id, 0, 0, 0); // workers, initial_phreeqc, utility
		}
		// Run cells with transported conditions
		snprintf(str, sizeof(str), "%s%10.1f%s", "Beginning reaction calculation       ", RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
		status = RM_LogMessage(id, str);
		status = RM_ScreenMessage(id, str);
		status = RM_RunCells(id);
		// Transfer data from PhreeqcRM for transport
		status = RM_GetConcentrations(id, c);          // Concentrations after reaction 
	}
	// Finalize
	status = RM_CloseFiles(id);
	status = RM_MpiWorkerBreak(id);
	status = RM_Destroy(id);
	// free space
	free(por);
	free(print_chemistry_mask);
	for (i = 0; i < ncomps; i++)
	{
		free(components[i]);
	}
	free(components);
	free(ic1);
	free(ic2);
	free(f1);
	free(bc1);
	free(bc2);
	free(bc_f1);
	free(bc_conc);
	free(c);
	free(temperature);
	free(pressure);
}
void simpleadvection_c(double* c, double* bc_conc, int ncomps, int nxyz, int dim)
{
	int i, j;
	// Advect
	for (i = nxyz - 1; i > 0; i--)
	{
		for (j = 0; j < ncomps; j++)
		{
			c[j * nxyz + i] = c[j * nxyz + i - 1];
		}
	}
	// Cell 0 gets boundary condition
	for (j = 0; j < ncomps; j++)
	{
		c[j * nxyz] = bc_conc[j * dim];
	}
}
