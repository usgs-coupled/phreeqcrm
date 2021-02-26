#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"

int gas_cpp()
{
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

		IRM_RESULT status;
		// Open files
		status = phreeqc_rm.SetFilePrefix("gas_cpp");
		phreeqc_rm.OpenFiles();
		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsGasPhase(0);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		// Set initial porosity
		std::vector<double> por;
		por.resize(nxyz, 0.2);
		status = phreeqc_rm.SetPorosity(por);

		// Set initial saturation
		std::vector<double> sat;
		sat.resize(nxyz, 0.5);
		status = phreeqc_rm.SetSaturation(sat);

		// Set printing of chemistry file
		status = phreeqc_rm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility

		// Load database
		std::string database = "phreeqc.dat";
		status = phreeqc_rm.LoadDatabase(database.c_str());
		if (status < 0)
		{
			std::ostringstream oss;
			oss << "Unable to open database " << database;
			phreeqc_rm.OutputMessage(oss.str());
		}

		// Run file to define solutions and gases for initial conditions
		std::string inputfile = "gas.pqi";
		status = phreeqc_rm.RunFile(false, true, false, inputfile.c_str());
		if (status < 0)
		{
			std::ostringstream oss;
			oss << "Unable to run input file " << inputfile;
			phreeqc_rm.OutputMessage(oss.str());
		}

		// Determine number of components and gas components
		int ncomps = phreeqc_rm.FindComponents();

		// Get list of gas component names
		const std::vector<std::string>& gcomps = phreeqc_rm.GetGasComponents();

		// Set array of initial conditions
		std::vector<int> ic1;
		ic1.resize((size_t)nxyz * 7, -1);
		for (size_t i = 0; i < nxyz; i++)
		{
			ic1[i] = 1;              // Solution 1
			ic1[(size_t)(4 * (size_t)nxyz + i)] = i % 3 + 1;    // Gas phase rotating 1,2,3
		}
		status = phreeqc_rm.InitialPhreeqc2Module(ic1);

		// Get gases
		std::vector<double> gas_moles;
		gas_moles.resize((size_t)nxyz * phreeqc_rm.GetGasComponentsCount());
		status = phreeqc_rm.GetGasPhaseMoles(gas_moles);
		// multiply by 2
		for (int i = 0; i < nxyz * phreeqc_rm.GetGasComponentsCount(); i++)
		{
			gas_moles[i] *= 2.0;
		}
		// set new moles
		status = phreeqc_rm.SetGasPhaseMoles(gas_moles);
		// retrieve new moles
		status = phreeqc_rm.GetGasPhaseMoles(gas_moles);
		// print cells 1,2,3
		for (size_t j = 0; j < 3; j++) // cell
		{
			std::cerr << "Cell: " << j << std::endl;
			for (size_t i = 0; i < 3; i++) // component
			{
				std::cerr << "  " << gcomps[i] << "  " << 
					gas_moles[i * (size_t)nxyz + j] << std::endl;
			}
		}
		// Clean up

		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "gas_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	catch (...)
	{
		std::string e_string = "gas_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
}
