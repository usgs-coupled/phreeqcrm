#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
#include "IPhreeqc.h"
#if defined(USE_MPI)
#include <mpi.h>
#endif
int worker_tasks_c(int *task_number, void * cookie);
int do_something(void *cookie);

void register_basic_callback(void *cookie);
double my_basic_callback(double x1, double x2, const char *str, void *cookie);

struct my_data
{
#ifdef USE_MPI
	MPI_Comm rm_comm;
#endif
	int phreeqcrm_id;
	double *K_ptr;
};
void advect_c(double *c, double *bc_conc, int ncomps, int nxyz, int dim);

    void advection_c()
	{
		// Based on PHREEQC Example 11

		int mpi_myself = 0;
		int i, j;
		int nxyz; 
#ifndef USE_MPI
		int nthreads;
#else
		MPI_Comm comm;
#endif
		int id;
		int status;
		double * rv;
		double * por;
		double * sat;
		int * print_chemistry_mask;
		int * grid2chem;
		int nchem;
		char str[100];
		char str1[200];
		int ncomps;
		char ** components;
		double * gfw;
		int * ic1; 
		int * ic2;
		double * f1;
		int * module_cells;
		int nbound;
		int * bc1;
		int * bc2;
		double * bc_f1;
		double * bc_conc;
		double * c;
		double time, time_step;
		double * density;
		double * volume;
		double * sat_calc;
		double * temperature;
		double * pressure;
		int isteps, nsteps;
		double * selected_out;
		int isel, n_user, col;
		char heading[100];
		double * c_well;
		double * tc;
		double * p_atm;
		double pH;
		int vtype;
		char svalue[100];
		int iphreeqc_id, iphreeqc_id1;
		int dump_on, append;
		char * errstr = NULL;
		int l;
		double * hydraulic_K;
		struct my_data some_data;
		int n;
		int * sc, * ec;

		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		nxyz = 40;
		// Bogus conductivity field for Basic callback demonstration
		hydraulic_K = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) hydraulic_K[i] = ((double) i) * 2.0;
		some_data.K_ptr = hydraulic_K;

#ifdef USE_MPI
		// MPI
		comm = MPI_COMM_WORLD;
		some_data.rm_comm = comm;
		id = RM_Create(nxyz, comm);
		some_data.phreeqcrm_id = id;
		if (MPI_Comm_rank(comm, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			status = RM_SetMpiWorkerCallback(id, worker_tasks_c);
			status = RM_SetMpiWorkerCallbackCookie(id, &some_data);
			status = RM_MpiWorker(id);
			status = RM_Destroy(id);
			return;
		}
#else
		// OpenMP
		nthreads = 3;
		id = RM_Create(nxyz, nthreads);
		some_data.phreeqcrm_id = id;
#endif
		// Set properties
		status = RM_SetErrorHandlerMode(id, 2);
		status = RM_SetComponentH2O(id, 0);
		status = RM_SetRebalanceFraction(id, 0.5);
		status = RM_SetRebalanceByCell(id, 1);
		status = RM_UseSolutionDensityVolume(id, 0);
		status = RM_SetPartitionUZSolids(id, 0);
		status = RM_SetFilePrefix(id, "Advect_c");
		// Open error, log, and output files
		status = RM_OpenFiles(id);
#ifdef USE_MPI
		// Optional callback for MPI
		status = do_something(&some_data);  // only root is calling do_something here
#endif
		// Set concentration units
		status = RM_SetUnitsSolution(id, 2);      // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = RM_SetUnitsPPassemblage(id, 1);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = RM_SetUnitsExchange(id, 1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = RM_SetUnitsSurface(id, 1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = RM_SetUnitsGasPhase(id, 1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = RM_SetUnitsSSassemblage(id, 1);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = RM_SetUnitsKinetics(id, 1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		// Set conversion from seconds to days
		status = RM_SetTimeConversion(id, 1.0 / 86400.0); 
		// Set representative volume
		rv = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) rv[i] = 1.0;
		status = RM_SetRepresentativeVolume(id, rv);
		// Set initial porosity
		por = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) por[i] = 0.2;
		status = RM_SetPorosity(id, por);
		// Set initial saturation
		sat = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) sat[i] = 1.0;
		status = RM_SetSaturation(id, sat);
		// Set cells to print chemistry when print chemistry is turned on
		print_chemistry_mask = (int *) malloc((size_t) (nxyz * sizeof(int)));
		for (i = 0; i < nxyz/2; i++) 
		{
			print_chemistry_mask[i] = 1;
			print_chemistry_mask[i + nxyz/2] = 0;
		}
		status = RM_SetPrintChemistryMask(id, print_chemistry_mask);	
		// Partitioning of uz solids
		status = RM_SetPartitionUZSolids(id, 0);
		// Demonstation of mapping, two equivalent rows by symmetry
		grid2chem = (int *) malloc((size_t) (nxyz * sizeof(int)));
		for (i = 0; i < nxyz/2; i++) 
		{
			grid2chem[i] = i;
			grid2chem[i + nxyz/2] = i;
		}
		status = RM_CreateMapping(id, grid2chem);
		if (status < 0) status = RM_DecodeError(id, status); 
		nchem = RM_GetChemistryCellCount(id);
		
		// --------------------------------------------------------------------------
		// Set initial conditions
		// --------------------------------------------------------------------------
		
		// Set printing of chemistry file
		status = RM_SetPrintChemistryOn(id, 0, 1, 0); // workers, initial_phreeqc, utility
		// Set printing of chemistry file
		status = RM_LoadDatabase(id, "phreeqc.dat"); 

		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		register_basic_callback(&some_data);

		// Demonstration of error handling if ErrorHandlerMode is 0
		if (status != IRM_OK)
		{
			l = RM_GetErrorStringLength(id);
			errstr = (char *) malloc((size_t) (l * sizeof(char) + 1));
			RM_GetErrorString(id, errstr, l+1);
			fprintf(stderr,"Beginning of error string:\n");
			fprintf(stderr,"%s", errstr);
			fprintf(stderr,"End of error string.\n");
			free(errstr);
			errstr = NULL;
			RM_Destroy(id);
			exit(1);
		}
		// Run file to define solutions and reactants for initial conditions, selected output
		// There are three types of IPhreeqc instances in PhreeqcRM
		// Argument 1 refers to the workers for doing reaction calculations for transport
		// Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
		// Argument 3 refers to the Utility instance
		status = RM_RunFile(id, 1, 1, 1, "advect.pqi");
		// Clear contents of workers and utility
		strcpy(str, "DELETE; -all");
		status = RM_RunString(id, 1, 0, 1, str);	// workers, initial_phreeqc, utility 
		// Determine number of components to transport
		ncomps = RM_FindComponents(id);
		// Print some of the reaction module information
		sprintf(str1, "Number of threads:                                %d\n", RM_GetThreadCount(id));
		status = RM_OutputMessage(id, str1);
		sprintf(str1, "Number of MPI processes:                          %d\n", RM_GetMpiTasks(id));
		status = RM_OutputMessage(id, str1);
		sprintf(str1, "MPI task number:                                  %d\n", RM_GetMpiMyself(id));
		status = RM_OutputMessage(id, str1);
		status = RM_GetFilePrefix(id, str, 100);
		sprintf(str1, "File prefix:                                      %s\n", str);
		status = RM_OutputMessage(id, str1);
		sprintf(str1, "Number of grid cells in the user's model:         %d\n", RM_GetGridCellCount(id));
		status = RM_OutputMessage(id, str1);
		sprintf(str1, "Number of chemistry cells in the reaction module: %d\n", RM_GetChemistryCellCount(id));
		status = RM_OutputMessage(id, str1);
		sprintf(str1, "Number of components for transport:               %d\n", RM_GetComponentCount(id));
		status = RM_OutputMessage(id, str1);
		// Get component information
		components = (char **) malloc((size_t) (ncomps * sizeof(char *)));
		gfw = (double *) malloc((size_t) (ncomps * sizeof(double)));
		status = RM_GetGfw(id, gfw);
		for (i = 0; i < ncomps; i++)
		{
			components[i] = (char *) malloc((size_t) (100 * sizeof(char *)));
			status = RM_GetComponent(id, i, components[i], 100);
			sprintf(str,"%10s    %10.3f\n", components[i], gfw[i]);
			status = RM_OutputMessage(id, str);
		}
		status = RM_OutputMessage(id, "\n");
		// Set array of initial conditions
		ic1 = (int *) malloc((size_t) (7 * nxyz * sizeof(int)));
		ic2 = (int *) malloc((size_t) (7 * nxyz * sizeof(int)));
		f1 = (double *) malloc((size_t) (7 * nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) 
		{
			ic1[i]          = 1;       // Solution 1
			ic1[nxyz + i]   = -1;      // Equilibrium phases none
			ic1[2*nxyz + i] = 1;       // Exchange 1
			ic1[3*nxyz + i] = -1;      // Surface none
			ic1[4*nxyz + i] = -1;      // Gas phase none
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
		// No mixing is defined, so the following is equivalent
		// status = RM_InitialPhreeqc2Module(id, ic1, NULL, NULL);

		// alternative for setting initial conditions
		// cell number in second argument (-1 indicates last solution, 40 in this case)
		// in advect.pqi and any reactants with the same number--
		// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
		// will be written to cells 18 and 19 (0 based)
		module_cells = (int *) malloc((size_t) (2 * sizeof(int)));
		module_cells[0] = 18;
		module_cells[1] = 19;
		status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2);
		// Initial equilibration of cells
		time = 0.0;
		time_step = 0.0;
		c = (double *) malloc((size_t) (ncomps * nxyz * sizeof(double)));
		status = RM_SetTime(id, time);
		status = RM_SetTimeStep(id, time_step);
		status = RM_RunCells(id); 
		status = RM_GetConcentrations(id, c);

		// --------------------------------------------------------------------------
		// Set boundary condition
		// --------------------------------------------------------------------------

		nbound = 1;
		bc1 = (int *) malloc((size_t) (nbound * sizeof(int)));
		bc2 = (int *) malloc((size_t) (nbound * sizeof(int)));
		bc_f1 = (double *) malloc((size_t) (nbound * sizeof(double)));
		bc_conc = (double *) malloc((size_t) (ncomps * nbound * sizeof(double)));
		for (i = 0; i < nbound; i++) 
		{
			bc1[i]          = 0;       // Solution 0 from Initial IPhreeqc instance
			bc2[i]          = -1;      // no bc2 solution for mixing
			bc_f1[i]        = 1.0;     // mixing fraction for bc1
		} 
		status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1);

		// --------------------------------------------------------------------------
		// Transient loop
		// --------------------------------------------------------------------------

		nsteps = 10;
		density = (double *) malloc((size_t) (nxyz * sizeof(double)));
		volume = (double *) malloc((size_t) (nxyz * sizeof(double)));
		pressure = (double *) malloc((size_t) (nxyz * sizeof(double)));
		temperature = (double *) malloc((size_t) (nxyz * sizeof(double)));
		sat_calc = (double *) malloc((size_t) (nxyz * sizeof(double)));
		for (i = 0; i < nxyz; i++) 
		{
			density[i] = 1.0;
			pressure[i] = 2.0;
			temperature[i] = 20.0;
		}
		status = RM_SetDensity(id, density);
		status = RM_SetPressure(id, pressure);      
		status = RM_SetTemperature(id, temperature); 
		time_step = 86400;
		status = RM_SetTimeStep(id, time_step);
		for (isteps = 0; isteps < nsteps; isteps++)
		{
			// Advection calculation
			sprintf(str, "%s%10.1f%s", "Beginning transport calculation      ", 
				RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
			status = RM_LogMessage(id, str);
			status = RM_SetScreenOn(id, 1);
			status = RM_ScreenMessage(id, str);
			sprintf(str, "%s%10.1f%s", "          Time step                  ", 
				RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days\n");
			status = RM_LogMessage(id, str);
			status = RM_ScreenMessage(id, str);
			advect_c(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			status = RM_SetPorosity(id, por);              // If porosity changes 
			status = RM_SetSaturation(id, sat);            // If saturation changes
			status = RM_SetTemperature(id, temperature);   // If temperature changes
			status = RM_SetPressure(id, pressure);         // If pressure changes
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
			sprintf(str, "%s%10.1f%s", "Beginning reaction calculation       ", RM_GetTime(id) * RM_GetTimeConversion(id), " days\n");
			status = RM_LogMessage(id, str);
			status = RM_ScreenMessage(id, str);
			status = RM_RunCells(id);  
			// Transfer data from PhreeqcRM for transport
			status = RM_GetConcentrations(id, c);          // Concentrations after reaction 
			status = RM_GetDensity(id, density);           // Density after reaction
			status = RM_GetSolutionVolume(id, volume);     // Solution volume after reaction
			status = RM_GetSaturation(id, sat_calc);       // Saturation after reaction
			// Print results at last time step
			if (isteps == nsteps - 1) 
			{
				fprintf(stderr, "Current distribution of cells for workers\n");
				fprintf(stderr, "Worker      First cell        Last Cell\n");
				n = RM_GetThreadCount(id) * RM_GetMpiTasks(id);
				sc = (int *) malloc((size_t) (n * sizeof(int)));
				ec = (int *) malloc((size_t) (n * sizeof(int)));
				status = RM_GetStartCell(id, sc);
				status = RM_GetEndCell(id, ec);
				for (i = 0; i < n; i++)
				{
					fprintf(stderr, "%d%s%d%s%d\n", i,"           ",sc[i], "                 ",ec[i]);
				}
				// Loop through possible multiple selected output definitions
				for (isel = 0; isel < RM_GetSelectedOutputCount(id); isel++)
				{
					n_user = RM_GetNthSelectedOutputUserNumber(id, isel);
					status = RM_SetCurrentSelectedOutputUserNumber(id, n_user);
					fprintf(stderr, "Selected output sequence number: %d\n", isel);
					fprintf(stderr, "Selected output user number:     %d\n", n_user);
					// Get double array of selected output values
					col = RM_GetSelectedOutputColumnCount(id);
					// allocate(selected_out(nxyz,col))
					selected_out = (double *) malloc((size_t) (col * nxyz * sizeof(double)));
					status = RM_GetSelectedOutput(id, selected_out);
					// Print results
					for (i = 0; i < RM_GetSelectedOutputRowCount(id)/2; i++)
					{
						fprintf(stderr, "Cell number %d\n", i);
						fprintf(stderr, "     Density: %f\n", density[i]);
						fprintf(stderr, "     Volume:  %f\n", volume[i]);
						fprintf(stderr, "     Components: \n");
						for (j = 0; j < ncomps; j++)
						{
							fprintf(stderr, "          %2d %10s: %10.4f\n", j, components[j], c[j*nxyz + i]);
						}
						fprintf(stderr, "     Selected output: \n");
						for (j = 0; j < col; j++)
						{
							status = RM_GetSelectedOutputHeading(id, j, heading, 100);  
							fprintf(stderr, "          %2d %10s: %10.4f\n", j, heading, selected_out[j*nxyz + i]);
						}
					}
					free(selected_out);
				}
			}
		}

		// --------------------------------------------------------------------------
		// Additional features and finalize
		// --------------------------------------------------------------------------

		// Use utility instance of PhreeqcRM to calculate pH of a mixture
		c_well = (double *) malloc((size_t) ((size_t) (1 * ncomps * sizeof(double))));
		for (i = 0; i < ncomps; i++)
		{
			c_well[i] = 0.5 * c[0 + nxyz*i] + 0.5 * c[9 + nxyz*i];
		}
		tc = (double *) malloc((size_t) (1 * sizeof(double)));
		p_atm = (double *) malloc((size_t) (1 * sizeof(double)));
		tc[0] = 15.0;
		p_atm[0] = 3.0;
		iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm);
		strcpy(str, "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1");
		// Alternatively, utility pointer is worker number nthreads + 1 
		iphreeqc_id1 = RM_GetIPhreeqcId(id, RM_GetThreadCount(id) + 1);
		SetOutputFileName(iphreeqc_id, "utility_c.txt");
		SetOutputFileOn(iphreeqc_id, 1);
		status = RunString(iphreeqc_id, str);
		if (status != 0) status = RM_Abort(id, status, "IPhreeqc RunString failed");
		status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5);
		status = GetSelectedOutputValue2(iphreeqc_id, 1, 0, &vtype, &pH, svalue, 100);
		// Dump results
		dump_on = 1;
		append = 0;
		status = RM_SetDumpFileName(id, "advection_c.dmp");
		status = RM_DumpModule(id, dump_on, append);   
		// Finalize
		status = RM_CloseFiles(id);
		status = RM_MpiWorkerBreak(id);
		status = RM_Destroy(id);
		// free space
		free(rv);
		free(por);
		free(sat);
		free(print_chemistry_mask);
		free(grid2chem);
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
		free(density);
		free(temperature);
		free(pressure);
		free(hydraulic_K);
	}
	void advect_c(double *c, double *bc_conc, int ncomps, int nxyz, int dim)
	{
		int i, j;
		// Advect
		for (i = nxyz/2 - 1; i > 0; i--)
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
#ifdef USE_MPI
int worker_tasks_c(int *method_number, void * cookie)
{
	if (*method_number == 1000)
	{
		do_something(cookie);
	}
	else if (*method_number == 1001)
	{
		register_basic_callback(cookie);
	}
	return 0;
}
int do_something(void *cookie)
{
	MPI_Status status;
	struct my_data *data; 
	int method_number, mpi_tasks, mpi_myself;
	int i, worker_number;

	data = (struct my_data *) cookie;

	method_number = 1000;
	MPI_Comm_size(data->rm_comm, &mpi_tasks);
	MPI_Comm_rank(data->rm_comm, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, data->rm_comm);
		fprintf(stderr, "I am Groot.\n");
		for (i = 1; i < mpi_tasks; i++)
		{
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, data->rm_comm, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, data->rm_comm);
	}
	return 0;
}
#endif
void register_basic_callback(void *cookie)
{		
	struct my_data *data; 
	int i, j, rm_id;
#ifdef USE_MPI
	int mpi_tasks, mpi_myself;
#endif
	int	method_number = 1001;
	data = (struct my_data *) cookie;

#ifdef USE_MPI
	MPI_Comm_size(data->rm_comm, &mpi_tasks);
	MPI_Comm_rank(data->rm_comm, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, data->rm_comm);
	}
#endif

	rm_id = data->phreeqcrm_id;
	for (i = 0; i < RM_GetThreadCount(rm_id) + 2; i++)
	{
		j = RM_GetIPhreeqcId(rm_id, i);
		SetBasicCallback(j, my_basic_callback, cookie);
	}
}
double my_basic_callback(double x1, double x2, const char *str, void *cookie)
{
	struct my_data *data;
	int rm_cell_number;
	int rm_id, size=4;
	int list[4];

	data = (struct my_data *) cookie;
	rm_id = data->phreeqcrm_id;

	rm_cell_number = (int) x1;
	if (rm_cell_number >= 0 && rm_cell_number < RM_GetChemistryCellCount(rm_id))
	{
		if (RM_GetBackwardMapping(rm_id, rm_cell_number, list, &size) == 0)
		{
			if (strcmp(str, "HYDRAULIC_K") == 0)
			{
				return data->K_ptr[list[0]];
			}
		}
	}
	return -999.9;
}
