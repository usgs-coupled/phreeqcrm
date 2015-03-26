#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"

#if defined(USE_MPI)
#include <mpi.h>
#endif

void SpeciesAdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);

int species_cpp()
{
	// Based on PHREEQC Example 11, transporting species rather than components

	try
	{
		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		int nxyz = 40;
#ifdef USE_MPI
		// MPI
		PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
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
		// Set properties
		status = phreeqc_rm.SetErrorHandlerMode(1);        // 1 = throw exception on error
		status = phreeqc_rm.SetSpeciesSaveOn(true);
		// Open files
		status = phreeqc_rm.SetFilePrefix("Species_cpp");
		phreeqc_rm.OpenFiles();
		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(2);      // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsPPassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		status = phreeqc_rm.SetUnitsExchange(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		status = phreeqc_rm.SetUnitsSurface(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		status = phreeqc_rm.SetUnitsGasPhase(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		status = phreeqc_rm.SetUnitsSSassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		status = phreeqc_rm.SetUnitsKinetics(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
		// Set conversion from seconds to user units (days)
		double time_conversion = 1.0 / 86400;
		status = phreeqc_rm.SetTimeConversion(time_conversion);
		// Set representative volume
		std::vector<double> rv;
		rv.resize(nxyz, 1.0);
		status = phreeqc_rm.SetRepresentativeVolume(rv);
		// Set initial porosity
		std::vector<double> por;
		por.resize(nxyz, 0.2);
		status = phreeqc_rm.SetPorosity(por);
		// Set initial saturation
		std::vector<double> sat;
		sat.resize(nxyz, 1.0);
		status = phreeqc_rm.SetSaturation(sat);
		// Set cells to print chemistry when print chemistry is turned on
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 0);
		for (int i = 0; i < nxyz/2; i++)
		{
			print_chemistry_mask[i] = 1;
		}
		status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask);
		// Demonstation of mapping, two equivalent rows by symmetry
		std::vector<int> grid2chem;
		grid2chem.resize(nxyz, -1);
		for (int i = 0; i < nxyz/2; i++)
		{
			grid2chem[i] = i;
			grid2chem[i + nxyz/2] = i;
		}
		status = phreeqc_rm.CreateMapping(grid2chem);
		if (status < 0) phreeqc_rm.DecodeError(status);
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
		// Print some of the reaction module information
		{
			char str1[100];
			sprintf(str1, "Number of threads:                                %d\n", phreeqc_rm.GetThreadCount());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "Number of MPI processes:                          %d\n", phreeqc_rm.GetMpiTasks());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "MPI task number:                                  %d\n", phreeqc_rm.GetMpiMyself());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "File prefix:                                      %s\n", phreeqc_rm.GetFilePrefix().c_str());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "Number of grid cells in the user's model:         %d\n", phreeqc_rm.GetGridCellCount());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "Number of chemistry cells in the reaction module: %d\n", phreeqc_rm.GetChemistryCellCount());
			phreeqc_rm.OutputMessage(str1);
			sprintf(str1, "Number of components for transport:               %d\n", phreeqc_rm.GetComponentCount());
			phreeqc_rm.OutputMessage(str1);
		}
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
		// Determine species information
		const std::vector<std::string> &species = phreeqc_rm.GetSpeciesNames();
		const std::vector < double > & species_z = phreeqc_rm.GetSpeciesZ();
		const std::vector < double > & species_d = phreeqc_rm.GetSpeciesD25();
		bool species_on = phreeqc_rm.GetSpeciesSaveOn();
		int nspecies = phreeqc_rm.GetSpeciesCount();
		for (int i = 0; i < nspecies; i++)
		{
			std::ostringstream strm;
			strm << species[i] << "\n";
			strm << "    Charge: " << species_z[i] << std::endl;
			strm << "    Dw:     " << species_d[i] << std::endl;
			cxxNameDouble::const_iterator it = phreeqc_rm.GetSpeciesStoichiometry()[i].begin();
			for (; it != phreeqc_rm.GetSpeciesStoichiometry()[i].end(); it++)
			{
				strm << "        " << it->first << "   " << it->second << "\n";
			}
			phreeqc_rm.OutputMessage(strm.str());
		}
		phreeqc_rm.OutputMessage("\n");
		// Set array of initial conditions
		std::vector<int> ic1, ic2;
		ic1.resize(nxyz*7, -1);
		for (int i = 0; i < nxyz; i++)
		{
			ic1[i] = 1;              // Solution 1
			ic1[2*nxyz + i] = 1;     // Exchange 1
		}
		status = phreeqc_rm.InitialPhreeqc2Module(ic1);
		// Initial equilibration of cells
		double time = 0.0;
		double time_step = 0.0;
		std::vector<double> c;
		status = phreeqc_rm.SetTime(time);
		status = phreeqc_rm.SetTimeStep(time_step);
		status = phreeqc_rm.RunCells();
		status = phreeqc_rm.GetSpeciesConcentrations(c);

		// --------------------------------------------------------------------------
		// Set boundary condition
		// --------------------------------------------------------------------------

		std::vector<double> bc_conc, bc_f1;
		std::vector<int> bc1, bc2;
		int nbound = 1;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
		bc2.resize(nbound, -1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = phreeqc_rm.InitialPhreeqc2SpeciesConcentrations(bc_conc, bc1);

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
		std::vector < double > component_c;
		component_c.resize(nxyz * ncomps);
		for (int steps = 0; steps < nsteps; steps++)
		{
			// Transport calculation here
			{
				std::ostringstream strm;
				strm << "Beginning transport calculation             " <<   phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion() << " days\n";
				strm << "          Time step                         " <<   phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.ScreenMessage(strm.str());
			}
			SpeciesAdvectCpp(c, bc_conc, nspecies, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on);
			status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			status = phreeqc_rm.SetPorosity(por);                    // If porosity changes due to compressibility
			status = phreeqc_rm.SetSaturation(sat);                  // If saturation changes
			status = phreeqc_rm.SetTemperature(temperature);         // If temperature changes
			status = phreeqc_rm.SetPressure(pressure);               // If pressure changes
			status = phreeqc_rm.SpeciesConcentrations2Module(c);     // Transported concentrations
			time = time + time_step;
			status = phreeqc_rm.SetTime(time);
			// Run cells with transported conditions
			{
				std::ostringstream strm;
				strm << "Beginning reaction calculation              " << time * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.ScreenMessage(strm.str());
			}
			std::vector < int > print_mask;
			print_mask.resize(nxyz, 1);
			phreeqc_rm.SetPrintChemistryMask(print_mask);
			status = phreeqc_rm.RunCells();
			// Transfer data from PhreeqcRM for transport
			status = phreeqc_rm.GetSpeciesConcentrations(c);
			phreeqc_rm.GetConcentrations(component_c);
			std::vector<double> density;
			status = phreeqc_rm.GetDensity(density);
			const std::vector<double> &volume = phreeqc_rm.GetSolutionVolume();
			// Print results at last time step
			if (print_chemistry_on != 0)
			{
				for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
				{
					// Loop through possible multiple selected output definitions
					int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
					status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
					std::cerr << "Selected output sequence number: " << isel << "\n";
					std::cerr << "Selected output user number:     " << n_user << "\n";
					// Get double array of selected output values
					std::vector<double> so;
					int col = phreeqc_rm.GetSelectedOutputColumnCount();
					so.resize(nxyz*col, 0);
					status = phreeqc_rm.GetSelectedOutput(so);
					// Print results
					for (int i = 0; i < phreeqc_rm.GetSelectedOutputRowCount()/2; i++)
					{
						std::cerr << "Cell number " << i << "\n";
						std::cerr << "     Density: " << density[i] << "\n";
						std::cerr << "     Volume:  " << volume[i] << "\n";
						std::cerr << "     Components: " << "\n";
						for (int j = 0; j < ncomps; j++)
						{
							std::cerr << "          " << j << " " << components[j] << ": " << component_c[j*nxyz + i] << "\n";
						}
						std::vector<std::string> headings;
						headings.resize(col);
						std::cerr << "     Selected output: " << "\n";
						for (int j = 0; j < col; j++)
						{
							status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
							std::cerr << "          " << j << " " << headings[j] << ": " << so[j*nxyz + i] << "\n";
						}
					}
				}
			}
		}

		// --------------------------------------------------------------------------
		// Additional features and finalize
		// --------------------------------------------------------------------------

 		// Use utility instance of PhreeqcRM to calculate pH of a mixture
		std::vector <double> c_well;
		c_well.resize(1*ncomps, 0.0);
		for (int i = 0; i < ncomps; i++)
		{
			c_well[i] = 0.5 * component_c[0 + nxyz*i] + 0.5 * component_c[9 + nxyz*i];
		}
		std::vector<double> tc, p_atm;
		tc.resize(1, 15.0);
		p_atm.resize(1, 3.0);
		IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c_well, tc, p_atm);
		{
			std::string input;
			input = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1";
			int iphreeqc_result;
			util_ptr->SetOutputFileName("species_utility_cpp.txt");
			util_ptr->SetOutputFileOn(true);
			iphreeqc_result = util_ptr->RunString(input.c_str());
			if (iphreeqc_result != 0)
			{
				phreeqc_rm.ErrorHandler(IRM_FAIL, "IPhreeqc RunString failed");
			}
			int vtype;
			double pH;
			char svalue[100];
			util_ptr->SetCurrentSelectedOutputUserNumber(5);
			iphreeqc_result = util_ptr->GetSelectedOutputValue2(1, 0, &vtype, &pH, svalue, 100);
		}
		// Dump results
		bool dump_on = true;
		bool append = false;
		status = phreeqc_rm.SetDumpFileName("species_cpp.dmp");
		status = phreeqc_rm.DumpModule(dump_on, append);
		// Finalize
		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "Species_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
		return IRM_FAIL;
	}
	catch (...)
	{
		std::string e_string = "Species_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
}
void
SpeciesAdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
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
		c[j * nxyz] = bc_conc[j * dim];                           // component j
	}
}
