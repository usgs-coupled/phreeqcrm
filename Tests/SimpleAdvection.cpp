#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"

void SimpleAdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);

int SimpleAdvection_cpp()
{
	// Based on PHREEQC Example 11
	try
	{
		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		int nxyz = 20;

#ifdef USE_MPI
		// MPI
		PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
		some_data.PhreeqcRM_ptr = &phreeqc_rm;
		MP_TYPE comm = MPI_COMM_WORLD;
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			phreeqc_rm.SetMpiWorkerCallbackC(worker_tasks_cc);
			phreeqc_rm.SetMpiWorkerCallbackCookie(&some_data);
			phreeqc_rm.MpiWorker();
			return EXIT_SUCCESS;
		}
#else
		// OpenMP
		int nthreads = 3;
		PhreeqcRM phreeqc_rm(nxyz, nthreads);
#endif
		IRM_RESULT status;
		// Set properties
		status = phreeqc_rm.SetComponentH2O(false);
		phreeqc_rm.UseSolutionDensityVolume(false);
		// Open files
		status = phreeqc_rm.SetFilePrefix("SimpleAdvection_cpp");
		phreeqc_rm.OpenFiles();

		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		// Set conversion from seconds to user units (days)
		double time_conversion = 1.0 / 86400;
		status = phreeqc_rm.SetTimeConversion(time_conversion);

		// Set initial porosity
		std::vector<double> por;
		por.resize(nxyz, 0.2);
		status = phreeqc_rm.SetPorosity(por);

		// Set initial saturation
		std::vector<double> sat(nxyz, 1.0);
		status = phreeqc_rm.SetSaturation(sat);

		// Set cells to print chemistry when print chemistry is turned on
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask);

		int nchem = phreeqc_rm.GetChemistryCellCount();

		// --------------------------------------------------------------------------
		// Set initial conditions
		// --------------------------------------------------------------------------

		// Set printing of chemistry file
		status = phreeqc_rm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
		// Load database
		status = phreeqc_rm.LoadDatabase("phreeqc.dat");

		// Run file to define solutions and reactants for initial conditions, selected output
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm.RunFile(workers, initial_phreeqc, utility, "advect.pqi");

		// Clear contents of workers and utility
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm.RunString(workers, initial_phreeqc, utility, input.c_str());

		// Determine number of components to transport
		int ncomps = phreeqc_rm.FindComponents();

		// Get component information
		const std::vector<std::string> &components = phreeqc_rm.GetComponents();
		const std::vector < double > & gfw = phreeqc_rm.GetGfw();
		for (int i = 0; i < ncomps; i++)
		{
			std::ostringstream strm;
			strm.width(10);
			strm << components[i] << "    " << gfw[i] << "\n";
			phreeqc_rm.OutputMessage(strm.str());
		}
		phreeqc_rm.OutputMessage("\n");

		// Set array of initial conditions
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
		status = phreeqc_rm.InitialPhreeqc2Module(ic1);

		// alternative for setting initial conditions
		// cell number in first argument (-1 indicates last solution, 40 in this case)
		// in advect.pqi and any reactants with the same number--
		// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
		// will be written to cells 18 and 19 (0 based)
		//std::vector<int> module_cells;
		//module_cells.push_back(18);
		//module_cells.push_back(19);
		//status = phreeqc_rm.InitialPhreeqcCell2Module(-1, module_cells);

		// Initial equilibration of cells
		double time = 0.0;
		double time_step = 0.0;
		std::vector<double> c;
		c.resize(nxyz * components.size());
		status = phreeqc_rm.SetTime(time);
		status = phreeqc_rm.SetTimeStep(time_step);
		status = phreeqc_rm.RunCells();
		status = phreeqc_rm.GetConcentrations(c);

		// --------------------------------------------------------------------------
		// Set boundary condition
		// --------------------------------------------------------------------------

		std::vector<double> bc_conc, bc_f1;
		std::vector<int> bc1, bc2;
		int nbound = 1;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
		bc2.resize(nbound, -1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, bc1, bc2, bc_f1);

		// --------------------------------------------------------------------------
		// Transient loop
		// --------------------------------------------------------------------------

		int nsteps = 10;
		std::vector<double> initial_density, temperature, pressure;
		initial_density.resize(nxyz, 1.0);
		temperature.resize(nxyz, 20.0);
		pressure.resize(nxyz, 2.0);
		phreeqc_rm.SetDensity(initial_density);
		phreeqc_rm.SetTemperature(temperature);
		phreeqc_rm.SetPressure(pressure);
		time_step = 86400.;
		status = phreeqc_rm.SetTimeStep(time_step);
		for (int steps = 0; steps < nsteps; steps++)
		{
			// Transport calculation here
			{
				std::ostringstream strm;
				strm << "Beginning transport calculation             " <<   phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion() << " days\n";
				strm << "          Time step                         " <<   phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.SetScreenOn(true);
				phreeqc_rm.ScreenMessage(strm.str());
			}
			SimpleAdvectCpp(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on);
			status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			status = phreeqc_rm.SetPorosity(por);             // If pororosity changes due to compressibility
			status = phreeqc_rm.SetSaturation(sat);           // If saturation changes
			status = phreeqc_rm.SetTemperature(temperature);  // If temperature changes
			status = phreeqc_rm.SetPressure(pressure);        // If pressure changes
			status = phreeqc_rm.SetConcentrations(c);         // Transported concentrations
			status = phreeqc_rm.SetTimeStep(time_step);		  // Time step for kinetic reactions
			time += time_step;
			status = phreeqc_rm.SetTime(time);
			// Run cells with transported conditions
			{
				std::ostringstream strm;
				strm << "Beginning reaction calculation              " << time * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.ScreenMessage(strm.str());
			}
			status = phreeqc_rm.RunCells();
			// Transfer data from PhreeqcRM for transport
			status = phreeqc_rm.GetConcentrations(c);
			std::vector<double> density;
			status = phreeqc_rm.GetDensity(density);
		}
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "Advection_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	catch (...)
	{
		std::string e_string = "Advection_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
}
void
SimpleAdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
{
	for (int i = nxyz/2 - 1 ; i > 0; i--)
	{
		for (int j = 0; j < ncomps; j++)
		{
			c[j * nxyz + i] = c[j * nxyz + i - 1];                 // component j
		}
	}
	// Cell zero gets boundary condition
	for (int j = 0; j < ncomps; j++)
	{
		c[j * nxyz] = bc_conc[j * dim];                                // component j
	}
}

