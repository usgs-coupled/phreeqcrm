#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
void PrintCells(const std::vector<std::string>& gcomps, const std::vector<double>& gas_moles,
	const std::vector<double>& gas_p, const std::vector<double>& gas_phi, int nxyz,
	const std::string str);
int Gas_cpp()
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
		status = phreeqc_rm.SetFilePrefix("Gas_cpp");
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
		status = phreeqc_rm.SetSaturationUser(sat);

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
		std::vector<double> gas_moles, gas_p, gas_phi;
		gas_moles.resize((size_t)nxyz * phreeqc_rm.GetGasComponentsCount());
		status = phreeqc_rm.GetGasCompMoles(gas_moles);
		status = phreeqc_rm.GetGasCompPressures(gas_p);
		status = phreeqc_rm.GetGasCompPhi(gas_phi);
		PrintCells(gcomps, gas_moles, gas_p, gas_phi, nxyz, "Initial conditions");

		// multiply by 2
		for (int i = 0; i < nxyz * phreeqc_rm.GetGasComponentsCount(); i++)
		{
			gas_moles[i] *= 2.0;
		}
		status = phreeqc_rm.SetGasCompMoles(gas_moles);
		status = phreeqc_rm.GetGasCompMoles(gas_moles);
		status = phreeqc_rm.GetGasCompPressures(gas_p);
		status = phreeqc_rm.GetGasCompPhi(gas_phi);
		PrintCells(gcomps, gas_moles, gas_p, gas_phi, nxyz, "Initial conditions times 2");

		// eliminate CH4 in cell 0
		gas_moles[0] = -1.0;
		// Gas phase is removed from cell 1
		gas_moles[1] = gas_moles[nxyz + 1] = gas_moles[2 * nxyz + 1] = -1.0;
		status = phreeqc_rm.SetGasCompMoles(gas_moles);
		phreeqc_rm.RunCells();
		status = phreeqc_rm.GetGasCompMoles(gas_moles);
		status = phreeqc_rm.GetGasCompPressures(gas_p);
		status = phreeqc_rm.GetGasCompPhi(gas_phi);
		PrintCells(gcomps, gas_moles, gas_p, gas_phi, nxyz, "Remove some components");

		// add CH4 in cell 0
		gas_moles[0] = 0.02;
		// Gas phase is added to cell 1; fixed pressure by default
		gas_moles[1] = 0.01;
		gas_moles[nxyz + 1] = 0.02;
		gas_moles[2 * nxyz + 1] = 0.03;
		status = phreeqc_rm.SetGasCompMoles(gas_moles);
		// Set volume for cell 1 and convert to fixed pressure gas phase
		std::vector<double> gas_volume;
		gas_volume.resize(nxyz, -1);
		gas_volume[1] = 12.25;
		status = phreeqc_rm.SetGasPhaseVolume(gas_volume);
		phreeqc_rm.RunCells();
		status = phreeqc_rm.GetGasCompMoles(gas_moles);
		status = phreeqc_rm.GetGasCompPressures(gas_p);
		status = phreeqc_rm.GetGasCompPhi(gas_phi);
		PrintCells(gcomps, gas_moles, gas_p, gas_phi, nxyz, "Add components back");

		// Clean up

		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "Gas_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
#if 0
	catch (...)
	{
		std::string e_string = "Gas_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
#endif
	return EXIT_SUCCESS;
}
/* ---------------------------------------------------------------------- */
void
PrintCells(const std::vector<std::string>& gcomps, const std::vector<double>& gas_moles,
	const std::vector<double>& gas_p, const std::vector<double>& gas_phi, int nxyz,
	const std::string str)
	/* ---------------------------------------------------------------------- */
{
	std::cerr << "\n" << str << std::endl;
	// print cells 0,1,2
	for (int j = 0; j < 3; j++) // cell
	{
		std::cerr << "Cell: " << j << std::endl;
		std::cerr << "               Moles         P         Phi" << std::endl;
		for (int i = 0; i < 3; i++) // component
		{
			int k = i * nxyz + j;
			fprintf(stderr, "%8s  %10.4f  %10.4f  %10.4f\n",
				gcomps[i].c_str(),
				gas_moles[k],
				gas_p[k],
				gas_phi[k]);
		}
	}
}
