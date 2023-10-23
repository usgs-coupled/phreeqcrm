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

void simpleadvection_cpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);

int SimpleAdvect_cpp()
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
		MP_TYPE comm = MPI_COMM_WORLD;
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			phreeqc_rm.MpiWorker();
			return EXIT_SUCCESS;
		}
#else
		// OpenMP
		int nthreads = 3;
		PhreeqcRM phreeqc_rm(nxyz, nthreads);
#endif
		// Set properties
		IRM_RESULT status = phreeqc_rm.SetComponentH2O(false);
		phreeqc_rm.UseSolutionDensityVolume(false);
		// Open files
		status = phreeqc_rm.SetFilePrefix("SimpleAdvect_cpp");
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
		status = phreeqc_rm.RunFile(true, true, true, "advect.pqi");
		// Clear contents of workers and utility
		std::string input = "DELETE; -all";
		status = phreeqc_rm.RunString(true, false, true, input.c_str());
		// Determine number of components to transport
		int ncomps = phreeqc_rm.FindComponents();
		// Get component information
		const std::vector<std::string>& components = phreeqc_rm.GetComponents();
		for (int i = 0; i < ncomps; i++)
		{
			std::ostringstream strm;
			strm.width(10);
			strm << components[i] << "\n";
			phreeqc_rm.OutputMessage(strm.str());
		}
		phreeqc_rm.OutputMessage("\n");
		// Set array of initial conditions
		std::vector<int> ic1;
		ic1.resize(nxyz * 7, -1);
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
		status = phreeqc_rm.InitialPhreeqc2Module(ic1);
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
		std::vector<double> bc_conc;
		std::vector<int> bc1;
		int nbound = 1;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
		status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, bc1);
		// --------------------------------------------------------------------------
		// Transient loop
		// --------------------------------------------------------------------------
		int nsteps = 10;
		std::vector<double> temperature, pressure;
		temperature.resize(nxyz, 20.0);
		pressure.resize(nxyz, 2.0);
		phreeqc_rm.SetTemperature(temperature);
		phreeqc_rm.SetPressure(pressure);
		time_step = 86400.;
		status = phreeqc_rm.SetTimeStep(time_step);
		for (int steps = 0; steps < nsteps; steps++)
		{
			// Transport calculation here
			{
				std::ostringstream strm;
				strm << "Beginning transport calculation             " << phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion() << " days\n";
				strm << "          Time step                         " << phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.SetScreenOn(true);
				phreeqc_rm.ScreenMessage(strm.str());
			}
			simpleadvection_cpp(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on);
			status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			status = phreeqc_rm.SetConcentrations(c);         // Transported concentrations
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
		}
		// Clean up
		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "SimpleAdvect_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
#if 0
	catch (...)
	{
		std::string e_string = "SimpleAdvect_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
#endif
	return EXIT_SUCCESS;
}
void
simpleadvection_cpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
{
	for (int i = nxyz - 1; i > 0; i--)
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

