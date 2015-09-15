#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
int worker_tasks_cc(int *task_number, void * cookie);
int do_something(void *cookie);

double my_basic_callback(double x1, double x2, const char *str, void *cookie);
void register_basic_callback(void *cookie);

#if defined(USE_MPI)
#include <mpi.h>
#endif

void AdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);
double my_basic_callback(double x1, double x2, const char *str, void *cookie);

class my_data
{
public:
	PhreeqcRM *PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::vector<double> *hydraulic_K;
};

int advection_cpp()
{
	// Based on PHREEQC Example 11
	try
	{
		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		int nxyz = 40;
		std::vector<double> hydraulic_K;
		for (int i = 0; i < nxyz; i++)
		{
			hydraulic_K.push_back(i*2.0);
		}
		my_data some_data;
		some_data.hydraulic_K = &hydraulic_K;

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
		some_data.PhreeqcRM_ptr = &phreeqc_rm;
#endif
		IRM_RESULT status;
		// Set properties
		status = phreeqc_rm.SetErrorHandlerMode(1);
		status = phreeqc_rm.SetComponentH2O(false);
		status = phreeqc_rm.SetRebalanceFraction(0.5);
		status = phreeqc_rm.SetRebalanceByCell(true);
		phreeqc_rm.UseSolutionDensityVolume(false);
		phreeqc_rm.SetPartitionUZSolids(false);
		// Open files
		status = phreeqc_rm.SetFilePrefix("Advect_cpp");
		phreeqc_rm.OpenFiles();
#ifdef USE_MPI
		// Optional callback for MPI
		int istatus = do_something(&some_data);                 // Root calls do_something, workers respond
#endif
		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSurface(1);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
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
		// test getters
		const std::vector<int> & print_chemistry_mask1 = phreeqc_rm.GetPrintChemistryMask();
		const std::vector<bool> & print_on = phreeqc_rm.GetPrintChemistryOn();
		bool rebalance = phreeqc_rm.GetRebalanceByCell();
		double f_rebalance = phreeqc_rm.GetRebalanceFraction();
		bool so_on = phreeqc_rm.GetSelectedOutputOn();
		int units_exchange = phreeqc_rm.GetUnitsExchange();
		int units_gas_phase = phreeqc_rm.GetUnitsGasPhase();
		int units_kinetics = phreeqc_rm.GetUnitsKinetics();
		int units_pp_assemblage = phreeqc_rm.GetUnitsPPassemblage();
		int units_solution = phreeqc_rm.GetUnitsSolution();
		int units_ss_exchange = phreeqc_rm.GetUnitsSSassemblage();
		int units_surface = phreeqc_rm.GetUnitsSurface();
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
		status = phreeqc_rm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
		// Load database
		status = phreeqc_rm.LoadDatabase("phreeqc.dat");

		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		register_basic_callback(&some_data);

		// Demonstration of error handling if ErrorHandlerMode is 0
		if (status != IRM_OK)
		{
			std::cerr << phreeqc_rm.GetErrorString(); // retrieve error messages if needed
			throw PhreeqcRMStop();
		}
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
		// Print some of the reaction module information
		{
			std::ostringstream oss;
			oss << "Database:                                         " << phreeqc_rm.GetDatabaseFileName().c_str() << "\n";
			oss << "Number of threads:                                " << phreeqc_rm.GetThreadCount() << "\n";
			oss << "Number of MPI processes:                          " << phreeqc_rm.GetMpiTasks() << "\n";
			oss << "MPI task number:                                  " << phreeqc_rm.GetMpiMyself() << "\n";
			oss << "File prefix:                                      " << phreeqc_rm.GetFilePrefix() << "\n";
			oss << "Number of grid cells in the user's model:         " << phreeqc_rm.GetGridCellCount() << "\n";
			oss << "Number of chemistry cells in the reaction module: " << phreeqc_rm.GetChemistryCellCount() << "\n";
			oss << "Number of components for transport:               " << phreeqc_rm.GetComponentCount() << "\n";
			oss << "Partioning of UZ solids:                          " << phreeqc_rm.GetPartitionUZSolids() << "\n";
			oss << "Error handler mode:                               " << phreeqc_rm.GetErrorHandlerMode() << "\n";
			phreeqc_rm.OutputMessage(oss.str());
		}
		const std::vector<int> &f_map = phreeqc_rm.GetForwardMapping();
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
		status = phreeqc_rm.InitialPhreeqc2Module(ic1, ic2, f1);
		// No mixing is defined, so the following is equivalent
		// status = phreeqc_rm.InitialPhreeqc2Module(ic1.data());

		// alternative for setting initial conditions
		// cell number in first argument (-1 indicates last solution, 40 in this case)
		// in advect.pqi and any reactants with the same number--
		// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
		// will be written to cells 18 and 19 (0 based)
		std::vector<int> module_cells;
		module_cells.push_back(18);
		module_cells.push_back(19);
		status = phreeqc_rm.InitialPhreeqcCell2Module(-1, module_cells);
		// Get temperatures
		const std::vector<double> &  tempc = phreeqc_rm.GetTemperature();
		// get current saturation
		std::vector<double> current_sat;
		status = phreeqc_rm.GetSaturation(current_sat);
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
			AdvectCpp(c, bc_conc, ncomps, nxyz, nbound);
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
			const std::vector<double> &volume = phreeqc_rm.GetSolutionVolume();
			// Print results at last time step
			if (print_chemistry_on != 0)
			{
				{
					std::ostringstream oss;
					oss << "Current distribution of cells for workers\n";
					oss << "Worker      First cell        Last Cell\n";
					int n;
#ifdef USE_MPI
					n = phreeqc_rm.GetMpiTasks();
#else
					n = phreeqc_rm.GetThreadCount();
#endif
					for (int i = 0; i < n; i++)
					{
						oss << i << "           " << phreeqc_rm.GetStartCell()[i] << "                 "
							<< phreeqc_rm.GetEndCell()[i] << "\n";
					}
					phreeqc_rm.OutputMessage(oss.str());
				}
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
							std::cerr << "          " << j << " " << components[j] << ": " << c[j*nxyz + i] << "\n";
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
			c_well[i] = 0.5 * c[0 + nxyz*i] + 0.5 * c[9 + nxyz*i];
		}
		std::vector<double> tc, p_atm;
		tc.resize(1, 15.0);
		p_atm.resize(1, 3.0);
		IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c_well, tc, p_atm);
		input = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1";
		int iphreeqc_result;
		util_ptr->SetOutputFileName("utility_cpp.txt");
		util_ptr->SetOutputFileOn(true);
		iphreeqc_result = util_ptr->RunString(input.c_str());
		// Alternatively, utility pointer is worker nthreads + 1
		IPhreeqc * util_ptr1 = phreeqc_rm.GetIPhreeqcPointer(phreeqc_rm.GetThreadCount() + 1);
		if (iphreeqc_result != 0)
		{
			phreeqc_rm.ErrorHandler(IRM_FAIL, "IPhreeqc RunString failed");
		}
		int vtype;
		double pH;
		char svalue[100];
		util_ptr->SetCurrentSelectedOutputUserNumber(5);
		iphreeqc_result = util_ptr->GetSelectedOutputValue2(1, 0, &vtype, &pH, svalue, 100);
		// Dump results
		bool dump_on = true;
		bool append = false;
		status = phreeqc_rm.SetDumpFileName("advection_cpp.dmp");
		status = phreeqc_rm.DumpModule(dump_on, append);    // gz disabled unless compiled with #define USE_GZ
		// Get pointer to worker
		const std::vector<IPhreeqcPhast *> w = phreeqc_rm.GetWorkers();
		w[0]->AccumulateLine("Delete; -all");
		iphreeqc_result = w[0]->RunAccumulated();
		// Clean up
		status = phreeqc_rm.CloseFiles();
		status = phreeqc_rm.MpiWorkerBreak();
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
AdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
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
int units_tester()
{
	try
	{
		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		int nxyz = 3;
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
		status = phreeqc_rm.SetErrorHandlerMode(1);
		status = phreeqc_rm.SetFilePrefix("Units_InitialPhreeqc_1");
		if (phreeqc_rm.GetMpiMyself() == 0)
		{
			phreeqc_rm.OpenFiles();
		}
		// Set concentration units
		status = phreeqc_rm.SetUnitsSolution(1);      // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm.SetUnitsPPassemblage(2);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsExchange(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSurface(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsGasPhase(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsSSassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm.SetUnitsKinetics(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		// Set representative volume
		std::vector<double> rv;
		rv.resize(nxyz, 1.0);
		status = phreeqc_rm.SetRepresentativeVolume(rv);
		// Set current porosity
		std::vector<double> por;
		por.resize(nxyz, 0.2);
		status = phreeqc_rm.SetPorosity(por);
		// Set saturation
		std::vector<double> sat;
		sat.resize(nxyz, 1.0);
		status = phreeqc_rm.SetSaturation(sat);
		// Set printing of chemistry file
		status = phreeqc_rm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility

		// --------------------------------------------------------------------------
		// Set initial conditions
		// --------------------------------------------------------------------------

		// Load database
		status = phreeqc_rm.LoadDatabase("phreeqc.dat");
		// Run file to define solutions and reactants for initial conditions, selected output
		bool workers = true;
		bool initial_phreeqc = true;
		bool utility = false;
		status = phreeqc_rm.RunFile(workers, initial_phreeqc, utility, "units.pqi");
		{
			std::string input = "DELETE; -all";
			status = phreeqc_rm.RunString(true, false, true, input.c_str());
		}
		status = phreeqc_rm.SetFilePrefix("Units_InitialPhreeqc_2");
		if (phreeqc_rm.GetMpiMyself() == 0)
		{
			phreeqc_rm.OpenFiles();
		}
		// Set reference to components
		int ncomps = phreeqc_rm.FindComponents();
		const std::vector<std::string> &components = phreeqc_rm.GetComponents();
		// Set initial conditions
		std::vector < int > cell_numbers;
		cell_numbers.push_back(0);
		status = phreeqc_rm.InitialPhreeqcCell2Module(1, cell_numbers);
		cell_numbers[0] = 1;
		status = phreeqc_rm.InitialPhreeqcCell2Module(2, cell_numbers);
		cell_numbers[0] = 2;
		status = phreeqc_rm.InitialPhreeqcCell2Module(3, cell_numbers);
		// Retrieve concentrations
		std::vector<double> c;
		status = phreeqc_rm.SetFilePrefix("Units_Worker");
		if (phreeqc_rm.GetMpiMyself() == 0)
		{
			phreeqc_rm.OpenFiles();
		}
		std::vector < int > print_mask;
		print_mask.resize(3, 1);
		phreeqc_rm.SetPrintChemistryMask(print_mask);
		phreeqc_rm.SetPrintChemistryOn(true,true,true);
		status = phreeqc_rm.RunCells();
		status = phreeqc_rm.GetConcentrations(c);
		std::vector<double> so;
		status = phreeqc_rm.GetSelectedOutput(so);
		std::vector<std::string> headings;
		{
			std::string heading;
			phreeqc_rm.GetSelectedOutputHeading(0, heading);
			std::cerr << "Cell  " << heading << std::endl;
			for (int i = 0; i < nxyz; i++)
			{
				std::cerr << i << "   " << so[i] << std::endl;
			}
		}
		// Use utility instance of PhreeqcRM
		std::vector<double> tc, p_atm;
		tc.resize(nxyz, 25.0);
		p_atm.resize(nxyz, 1.0);
		IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c, tc, p_atm);
		std::string input;
		input = "RUN_CELLS; -cells 0-2";
		// Output goes to new file
		int iphreeqc_result;
		util_ptr->SetOutputFileName("Units_utility.out");
		util_ptr->SetOutputFileOn(true);
		iphreeqc_result = util_ptr->RunString(input.c_str());
		status = phreeqc_rm.MpiWorkerBreak();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "Units tester failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	catch (...)
	{
		std::string e_string = "Units tester failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
}
#ifdef USE_MPI
int worker_tasks_cc(int *task_number, void * cookie)
{
	if (*task_number == 1000)
	{
		do_something(cookie);
	}
	else if (*task_number == 1001)
	{
		register_basic_callback(cookie);
	}
	return 0;
}
int do_something(void *cookie)
{
	int method_number = 1000;
	my_data *data = (my_data *) cookie;
	int mpi_tasks, mpi_myself, worker_number;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	std::stringstream msg;
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
		fprintf(stderr, "I am root.\n");
		for (int i = 1; i < mpi_tasks; i++)
		{
			MPI_Status status;
			MPI_Recv(&worker_number, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
	}
	return 0;
}
#endif
void register_basic_callback(void *cookie)
{		
	my_data *data; 
#ifdef USE_MPI
	int mpi_tasks, mpi_myself;
#endif
	int	method_number = 1001;
	data = (my_data *) cookie;

#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	}
#endif

	const std::vector<IPhreeqcPhast *> w = data->PhreeqcRM_ptr->GetWorkers();
	for (int i = 0; i < (int) w.size(); i++)
	{
		w[i]->SetBasicCallback(my_basic_callback, cookie);
	}
}
double my_basic_callback(double x1, double x2, const char *str, void *cookie)
{
	my_data * data_ptr = (my_data *) cookie;
	PhreeqcRM *phreeqcrm_ptr = data_ptr->PhreeqcRM_ptr;
	std::string option(str);

	int rm_cell_number = (int) x1;
	if (rm_cell_number >= 0 && rm_cell_number < phreeqcrm_ptr->GetChemistryCellCount())
	{
		const std::vector < std::vector <int> > & back = phreeqcrm_ptr->GetBackwardMapping();
		if (option == "HYDRAULIC_K")
		{
			return (*data_ptr->hydraulic_K)[back[rm_cell_number][0]];
		}
	}
	return -999.9;
}