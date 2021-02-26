#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
void PrintCells(char** gcomps, double* gas_moles, int nxyz,
	const char* str);

    void gas_c()
	{
		int i;
		int nxyz; 
#ifndef USE_MPI
		int nthreads;
#else
		MPI_Comm comm;
#endif
		int id;
		int status;
		int nchem;
		int ncomps;
		int ngas;
		char ** gas_comps;
		int * ic1;
		int * ic2;
		double * f1;
		double * gas_moles;
		double* por;
		double* sat;
		char * errstr = NULL;

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
		status = RM_SetFilePrefix(id, "gas_c");
		// Open error, log, and output files
		status = RM_OpenFiles(id);

		// Set concentration units
		status = RM_SetUnitsSolution(id, 2);      // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = RM_SetUnitsGasPhase(id, 0);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		// Set initial porosity
		por = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) por[i] = 0.2;
		status = RM_SetPorosity(id, por);
		// Set initial saturation
		sat = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) sat[i] = 0.5;
		status = RM_SetSaturation(id, sat);

		// Set printing of chemistry file
		status = RM_SetPrintChemistryOn(id, 0, 1, 0); // workers, initial_phreeqc, utility

		nchem = RM_GetChemistryCellCount(id);

		// Set printing of chemistry file
		status = RM_LoadDatabase(id, "phreeqc.dat"); 
		if (status < 0)
		{
			status = RM_OutputMessage(id,"Unable to open database.");
		}

		// Run file to define solutions and gases for initial conditions
		status = RM_RunFile(id, 0, 1, 0, "gas.pqi");
		if (status < 0)
		{
			status = RM_OutputMessage(id, "Unable to run input file.");
		}

		// Determine number of components and gas components
		ncomps = RM_FindComponents(id);
		ngas = RM_GetGasComponentsCount(id);

		// Get gas component names
		gas_comps = (char **) malloc((size_t) (ngas * sizeof(char *)));
		for (i = 0; i < ngas; i++)
		{
			gas_comps[i] = (char *) malloc((size_t) (100 * sizeof(char *)));
			status = RM_GetGasComponentsName(id, i, gas_comps[i], 100);
		}

		// Set array of initial conditions
		ic1 = (int *) malloc((size_t) (7 * nxyz * sizeof(int)));
		ic2 = (int*) malloc((size_t)(7 * nxyz * sizeof(int)));
		f1 = (double*) malloc((size_t)(7 * nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) 
		{
			ic1[i]          = 1;       // Solution 1
			ic1[nxyz + i]   = -1;      // Equilibrium phases none
			ic1[2*nxyz + i] = -1;       // Exchange 1
			ic1[3*nxyz + i] = -1;      // Surface none
			ic1[4*nxyz + i] = i % 3 + 1;      // Gas phase none
			ic1[5*nxyz + i] = -1;      // Solid solutions none
			ic1[6*nxyz + i] = -1;      // Kinetics none
			ic2[i]          = -1;      // Solution none
			ic2[nxyz + i]   = -1;      // Equilibrium phases none
			ic2[2*nxyz + i] = -1;      // Exchange none
			ic2[3*nxyz + i] = -1;      // Surface none
			ic2[4*nxyz + i] = -1;      // Gas phase none
			ic2[5*nxyz + i] = -1;      // Solid solutions none
			ic2[6*nxyz + i] = -1;      // Kinetics none
			f1[i]          = 1.0;      // Mixing fraction ic1 Solution
			f1[nxyz + i]   = 1.0;      // Mixing fraction ic1 Equilibrium phases 
			f1[2*nxyz + i] = 1.0;      // Mixing fraction ic1 Exchange 1
			f1[3*nxyz + i] = 1.0;      // Mixing fraction ic1 Surface 
			f1[4*nxyz + i] = 1.0;      // Mixing fraction ic1 Gas phase 
			f1[5*nxyz + i] = 1.0;      // Mixing fraction ic1 Solid solutions 
			f1[6*nxyz + i] = 1.0;      // Mixing fraction ic1 Kinetics 
		}
		status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1); 

		// Get gases
		gas_moles = (double *) malloc((size_t) (ngas * nxyz * sizeof(double)));
		status = RM_GetGasPhaseMoles(id, gas_moles);
		PrintCells(gas_comps, gas_moles, nxyz, "Initial conditions");

		// multiply by 2
		for (int i = 0; i < nxyz * RM_GetGasComponentsCount(id); i++)
		{
			gas_moles[i] *= 2.0;
		}
		status = RM_SetGasPhaseMoles(id, gas_moles);
		status = RM_GetGasPhaseMoles(id, gas_moles);
		PrintCells(gas_comps, gas_moles, nxyz, "Initial conditions times 2");

		// eliminate CH4 in cell 0, and all of Cell 1 gases
		gas_moles[0] = -1.0;
		gas_moles[1] = gas_moles[nxyz + 1] = gas_moles[2 * nxyz + 1] = -1.0;
		status = RM_SetGasPhaseMoles(id, gas_moles);
		status = RM_GetGasPhaseMoles(id, gas_moles);
		PrintCells(gas_comps, gas_moles, nxyz, "Remove some components");
		// Finalize
		status = RM_MpiWorkerBreak(id);
		status = RM_CloseFiles(id);
		status = RM_Destroy(id);

		// free space
		for (i = 0; i < ngas; i++)
		{
			free(gas_comps[i]);
		}
		free(gas_comps);
		free(ic1);
		free(ic2);
		free(f1);
		free(gas_moles);
		free(por);
		free(sat);
	}
	void PrintCells(char** gcomps, double* gas_moles, int nxyz,
		const char* str)
		/* ---------------------------------------------------------------------- */
	{
		fprintf(stderr, "\n%s\n", str);
		// print cells 1,2,3
		for (int j = 0; j < 3; j++) // cell
		{
			fprintf(stderr, "Cell: %i\n", j);
			for (size_t i = 0; i < 3; i++) // component
			{
				fprintf(stderr, "  %10s  %10.4f\n", gcomps[i],
					gas_moles[i * nxyz + j]);
			}
		}
	}