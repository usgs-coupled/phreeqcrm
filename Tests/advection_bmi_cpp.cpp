#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
#include <algorithm>
#include "yaml-cpp/yaml.h"
#include "YAMLPhreeqcRM.h"
int worker_tasks_cc(int* task_number, void* cookie);
int do_something(void* cookie);

double bmi_basic_callback(double x1, double x2, const char* str, void* cookie);
void bmi_register_basic_callback(void* cookie);

void AdvectBMICpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);
int bmi_example_selected_output(PhreeqcRM& phreeqc_rm);
double bmi_basic_callback(double x1, double x2, const char* str, void* cookie);
void BMI_testing(PhreeqcRM &phreeqc_rm);
//void GenerateYAML(int nxyz, std::string YAML_filename);
class my_data
{
public:
	PhreeqcRM* PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::vector<double>* hydraulic_K;
};

int advection_bmi_cpp()
{

	// --------------------------------------------------------------------------
	// Write file to initialize PhreeqcRM
	// Based on PHREEQC Example 11
	// --------------------------------------------------------------------------

	try
	{
		// Write file to initialize PhreeqcRM
		// emulates a file created by a GUI
		// No PhreeqcRM methods are used to create
		// this file
		int nxyz = 40;
		//std::string YAML_filename = "advect.yaml";
		//GenerateYAML(nxyz, YAML_filename);

		// Data for call_back demostration
		std::vector<double> hydraulic_K;
		for (int i = 0; i < nxyz; i++)
		{
			hydraulic_K.push_back(i * 2.0);
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
		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		bmi_register_basic_callback(&some_data);

		// Use YAML file to initialize
		phreeqc_rm.InitializeYAML("advect.yaml");

		// Get number of components
		int ncomps;
		phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
		// Example for generating automatic selected-output definitions
		bmi_example_selected_output(phreeqc_rm);

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
		const std::vector<int>& f_map = phreeqc_rm.GetForwardMapping();

		// BMI version of getting component list
		// equivalent to following statement
		// const std::vector<std::string>& components = phreeqc_rm.GetComponents();
		std::vector<std::string> components;
		{
			size_t nbytes = (size_t)phreeqc_rm.BMI_GetVarNbytes("Components");
			std::string all_comps(nbytes, ' ');
			phreeqc_rm.BMI_GetValue("Components", (void*)all_comps.data());
			size_t string_size = (size_t)phreeqc_rm.BMI_GetVarItemsize("Components");
			for (size_t i = 0; i < ncomps; i++)
			{
				std::string bmi_comp = all_comps.substr((i * string_size), string_size);
				// strip trailing spaces
				size_t end = bmi_comp.find_last_not_of(' ');
				bmi_comp = (end == std::string::npos) ? "" : bmi_comp.substr(0, end + 1);
				// store in components
				components.push_back(bmi_comp);
			}
		}

		// BMI version of getting gram formula weights of components
		// equivalent to following statement
		//const std::vector < double >& gfw = phreeqc_rm.GetGfw();
		std::vector<double> gfw(ncomps, 0);
		phreeqc_rm.BMI_GetValue("Gfw", gfw.data());

		// write component information
		for (int i = 0; i < ncomps; i++)
		{
			std::ostringstream strm;
			strm.width(10);
			strm << components[i] << "    " << gfw[i] << "\n";
			phreeqc_rm.OutputMessage(strm.str());
		}
		phreeqc_rm.OutputMessage("\n");

		// Get temperatures (unchanged by PhreeqcRM)
		//const std::vector<double>& tempc = phreeqc_rm.GetTemperature();
		std::vector<double> tempc(nxyz, 0.0);
		phreeqc_rm.BMI_GetValue("Temperature", tempc.data());

		// get current saturation (unchanged by PhreeqcRM)
		//std::vector<double> current_sat;
		//IRM_RESULT status = phreeqc_rm.GetSaturation(current_sat);		// Demonstration of error handling if ErrorHandlerMode is 0
		//if (status != IRM_OK)
		//{
		//	std::cerr << phreeqc_rm.GetErrorString(); // retrieve error messages if needed
		//	throw PhreeqcRMStop();
		//}
		IRM_RESULT status;
		std::vector<double> current_sat(nxyz, 0.0);
		phreeqc_rm.BMI_GetValue("Saturation", current_sat.data());

		std::vector<double> c;
		c.resize(nxyz * components.size());
		//status = phreeqc_rm.GetConcentrations(c);
		phreeqc_rm.BMI_GetValue("Concentrations", c.data());

		// Set density, temperature, and pressure
		std::vector<double> density(nxyz, 1.0);
		std::vector<double> temperature(nxyz, 20.0);
		std::vector<double> pressure(nxyz, 2.0);
		phreeqc_rm.BMI_SetValue("Density", density.data());
		phreeqc_rm.BMI_SetValue("Temperature", temperature.data());
		phreeqc_rm.BMI_SetValue("Pressure", pressure.data());
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
		std::vector<double> por(nxyz, 0.2);
		std::vector<double> sat(nxyz, 1.0);

		double time_step = phreeqc_rm.BMI_GetTimeStep();
		double time = phreeqc_rm.BMI_GetCurrentTime();

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
			AdvectBMICpp(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			//status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on);
			phreeqc_rm.BMI_SetValue("SelectedOutputOn",&print_selected_output_on);
			status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			//status = phreeqc_rm.SetPorosity(por);             // If pororosity changes due to compressibility
			phreeqc_rm.BMI_SetValue("Porosity", por.data());
			//status = phreeqc_rm.SetSaturation(sat);           // If saturation changes
			phreeqc_rm.BMI_SetValue("Saturation", sat.data());
			//status = phreeqc_rm.SetTemperature(temperature);  // If temperature changes
			phreeqc_rm.BMI_SetValue("Temperature", temperature.data());
			//status = phreeqc_rm.SetPressure(pressure);        // If pressure changes
			phreeqc_rm.BMI_SetValue("Pressure", pressure.data());
			//status = phreeqc_rm.SetConcentrations(c);         // Transported concentrations
			phreeqc_rm.BMI_SetValue("Concentrations", c.data());
			//status = phreeqc_rm.SetTimeStep(time_step);		  // Time step for kinetic reactions
			phreeqc_rm.BMI_SetValue("TimeStep", &time_step);
			time += time_step;
			//status = phreeqc_rm.SetTime(time);
			phreeqc_rm.BMI_SetValue("Time", &time);
			// Run cells with transported conditions
			{
				std::ostringstream strm;
				strm << "Beginning reaction calculation              " << time * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.ScreenMessage(strm.str());
			}
			// Demonstration of state
			status = phreeqc_rm.StateSave(1);
			status = phreeqc_rm.StateApply(1);
			status = phreeqc_rm.StateDelete(1);
			//status = phreeqc_rm.RunCells();
			phreeqc_rm.BMI_Update();
			// Transfer data from PhreeqcRM for transport
			//status = phreeqc_rm.GetConcentrations(c);
			phreeqc_rm.BMI_GetValue("Concentrations", c.data());
			//status = phreeqc_rm.GetDensity(density);
			phreeqc_rm.BMI_GetValue("Density", density.data());
			const std::vector<double>& volume = phreeqc_rm.GetSolutionVolume();
			// Print results at last time step
			if (print_chemistry_on != 0)
			{
				{
					std::ostringstream oss;
					oss << "Current distribution of cells for workers\n";
					oss << "Worker      First cell        Last Cell\n";
					int n;
					n = phreeqc_rm.GetThreadCount() * phreeqc_rm.GetMpiTasks();
					for (int i = 0; i < n; i++)
					{
						oss << i << "           " << phreeqc_rm.GetStartCell()[i] << "                 "
							<< phreeqc_rm.GetEndCell()[i] << "\n";
					}
					phreeqc_rm.OutputMessage(oss.str());
				}
				int selected_output_count;
				phreeqc_rm.BMI_GetValue("SelectedOutputCount", &selected_output_count);
				for (int isel = 0; isel < selected_output_count; isel++)
				{
					std::ostringstream oss;
					// Loop through possible multiple selected output definitions
					//int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
					//status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
					phreeqc_rm.BMI_SetValue("NthSelectedOutput", &isel);
					oss << "Selected output sequence number: " << isel << "\n";
					int n_user=-1;
					phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", &n_user);
					oss << "Selected output user number:     " << n_user << "\n";
					// Get double array of selected output values
					//std::vector<double> so;
					//int col = phreeqc_rm.GetSelectedOutputColumnCount();
					//status = phreeqc_rm.GetSelectedOutput(so);
					
					int col;
					phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", &col);
					int bmi_row_count;
					phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", &bmi_row_count);

					// Get selected output
					int bmi_nbytes = phreeqc_rm.BMI_GetVarNbytes("SelectedOutput");
					int bmi_itemsize = phreeqc_rm.BMI_GetVarItemsize("SelectedOutput");
					int bmi_dim = bmi_nbytes / bmi_itemsize;
					std::vector<double> so(bmi_dim, INACTIVE_CELL_VALUE);
					phreeqc_rm.BMI_GetValue("SelectedOutput", so.data());

					// Get headings
					bmi_nbytes = phreeqc_rm.BMI_GetVarNbytes("SelectedOutputHeadings");
					int bmi_string_size = phreeqc_rm.BMI_GetVarItemsize("SelectedOutputHeadings");
					std::string all_headings(bmi_nbytes, ' ');
					phreeqc_rm.BMI_GetValue("SelectedOutputHeadings", (void*)all_headings.c_str());
					std::vector<std::string> headings;
					for (int j = 0; j < col; j++)
					{
						std::string bmi_head = all_headings.substr((size_t)(j * bmi_string_size), bmi_string_size);
						size_t end = bmi_head.find_last_not_of(' ');
						bmi_head = (end == std::string::npos) ? "" : bmi_head.substr(0, end + 1);
						headings.push_back(bmi_head);
					}

					// Print results
					for (int i = 0; i < bmi_row_count / 2; i++)
					{
						oss << "Cell number " << i << "\n";
						oss << "     Density: " << density[i] << "\n";
						oss << "     Volume:  " << volume[i] << "\n";
						oss << "     Components: " << "\n";
						for (int j = 0; j < ncomps; j++)
						{
							oss << "          " << j << " " << components[j] << ": " << c[j * nxyz + i] << "\n";
						}
						//std::vector<std::string> headings;
						//headings.resize(col);

						oss << "     Selected output: " << "\n";
						for (int j = 0; j < col; j++)
						{
							//status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
							oss << "          " << j << " " << headings[j] << ": " << so[j * nxyz + i] << "\n";
						}
					}
					phreeqc_rm.OutputMessage(oss.str());
				}
			}
		}
		BMI_testing(phreeqc_rm); // Tests BMI_GetValue
		// --------------------------------------------------------------------------
		// Additional features and finalize
		// --------------------------------------------------------------------------

		// Use utility instance of PhreeqcRM to calculate pH of a mixture
		std::vector <double> c_well;
		c_well.resize(1 * ncomps, 0.0);
		for (int i = 0; i < ncomps; i++)
		{
			c_well[i] = 0.5 * c[0 + nxyz * i] + 0.5 * c[9 + nxyz * i];
		}
		std::vector<double> tc, p_atm;
		tc.resize(1, 15.0);
		p_atm.resize(1, 3.0);
		IPhreeqc* util_ptr = phreeqc_rm.Concentrations2Utility(c_well, tc, p_atm);
		std::string input = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1";
		int iphreeqc_result;
		util_ptr->SetOutputFileName("utility_cpp.txt");
		util_ptr->SetOutputFileOn(true);
		iphreeqc_result = util_ptr->RunString(input.c_str());
		// Alternatively, utility pointer is worker nthreads + 1
		IPhreeqc* util_ptr1 = phreeqc_rm.GetIPhreeqcPointer(phreeqc_rm.GetThreadCount() + 1);
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
		status = phreeqc_rm.SetDumpFileName("advection_bmi_cpp.dmp");
		status = phreeqc_rm.DumpModule(dump_on, append);    // gz disabled unless compiled with #define USE_GZ
		// Get pointer to worker
		const std::vector<IPhreeqcPhast*> w = phreeqc_rm.GetWorkers();
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
		std::string e_string = "Advection_bmi_cpp failed with an unhandled exception.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
	return EXIT_SUCCESS;
}
void
AdvectBMICpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
{
	for (int i = nxyz / 2 - 1; i > 0; i--)
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
int bmi_units_tester()
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
		status = phreeqc_rm.SetErrorOn(true);
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
		const std::vector<std::string>& components = phreeqc_rm.GetComponents();
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
		phreeqc_rm.SetPrintChemistryOn(true, true, true);
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
		IPhreeqc* util_ptr = phreeqc_rm.Concentrations2Utility(c, tc, p_atm);
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
int worker_tasks_cc(int* task_number, void* cookie)
{
	if (*task_number == 1000)
	{
		do_something(cookie);
	}
	else if (*task_number == 1001)
	{
		bmi_register_basic_callback(cookie);
	}
	return 0;
}
int do_something(void* cookie)
{
	int method_number = 1000;
	my_data* data = (my_data*)cookie;
	int mpi_tasks, mpi_myself, worker_number;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	std::stringstream msg;
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
		fprintf(stderr, "I am root.\n");
		for (int i = 1; i < mpi_tasks; i++)
		{
			MPI_Status status;
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	return 0;
}
#endif
void bmi_register_basic_callback(void* cookie)
{
	my_data* data;
#ifdef USE_MPI
	int mpi_tasks, mpi_myself;
#endif
	int	method_number = 1001;
	data = (my_data*)cookie;

#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
#endif

	const std::vector<IPhreeqcPhast*> w = data->PhreeqcRM_ptr->GetWorkers();
	for (int i = 0; i < (int)w.size(); i++)
	{
		w[i]->SetBasicCallback(bmi_basic_callback, cookie);
	}
}
double bmi_basic_callback(double x1, double x2, const char* str, void* cookie)
{
	my_data* data_ptr = (my_data*)cookie;
	PhreeqcRM* phreeqcrm_ptr = data_ptr->PhreeqcRM_ptr;
	std::string option(str);

	int rm_cell_number = (int)x1;
	if (rm_cell_number >= 0 && rm_cell_number < phreeqcrm_ptr->GetChemistryCellCount())
	{
		const std::vector < std::vector <int> >& back = phreeqcrm_ptr->GetBackwardMapping();
		if (option == "HYDRAULIC_K")
		{
			return (*data_ptr->hydraulic_K)[back[rm_cell_number][0]];
		}
	}
	return -999.9;
}
int bmi_example_selected_output(PhreeqcRM& phreeqc_rm)
{

	std::ostringstream oss;
	oss << "SELECTED_OUTPUT 2" << "\n";
	// totals
	oss << "  -totals " << "\n";
	const std::vector<std::string>& components = phreeqc_rm.GetComponents();
	for (size_t i = 0; i < phreeqc_rm.GetComponents().size(); i++)
	{
		if (components[i] != "H" &&
			components[i] != "O" &&
			components[i] != "charge" &&
			components[i] != "Charge" &&
			components[i] != "H2O")
		{
			oss << "    " << components[i] << "\n";
		}
	}

	oss << "  -molalities " << "\n";
	{
		// molalities of aqueous species 
		const std::vector<std::string>& aq_species = phreeqc_rm.GetSpeciesNames();
		for (size_t i = 0; i < phreeqc_rm.GetSpeciesNames().size(); i++)
		{
			oss << "    " << aq_species[i] << "\n";
		}
	}
	{
		// molalities of exchange species
		const std::vector<std::string>& ex_species = phreeqc_rm.GetExchangeSpecies();
		const std::vector<std::string>& ex_names = phreeqc_rm.GetExchangeNames();
		for (size_t i = 0; i < phreeqc_rm.GetExchangeSpeciesCount(); i++)
		{

			oss << "    ";
			oss.width(15);
			oss << std::left << ex_species[i];
			oss << " # " << ex_names[i] << "\n";
		}
	}
	{
		// molalities of surface species
		const std::vector<std::string>& surf_species = phreeqc_rm.GetSurfaceSpecies();
		const std::vector<std::string>& surf_types = phreeqc_rm.GetSurfaceTypes();
		const std::vector<std::string>& surf_names = phreeqc_rm.GetSurfaceNames();
		for (size_t i = 0; i < phreeqc_rm.GetSurfaceSpeciesCount(); i++)
		{
			oss << "    ";
			oss.width(15);
			oss << std::left << surf_species[i];
			oss << " # ";
			oss.width(15);
			oss << surf_types[i] << "   " << surf_names[i] << "\n";
		}
	}
	oss << "  -equilibrium_phases " << "\n";
	{
		// equilibrium phases 
		const std::vector<std::string>& eq_phases = phreeqc_rm.GetEquilibriumPhases();
		for (size_t i = 0; i < phreeqc_rm.GetEquilibriumPhasesCount(); i++)
		{
			oss << "    " << eq_phases[i] << "\n";
		}
	}
	oss << "  -gases " << "\n";
	{
		// gas components
		const std::vector<std::string>& gas_phases = phreeqc_rm.GetGasComponents();
		for (size_t i = 0; i < phreeqc_rm.GetGasComponentsCount(); i++)
		{
			oss << "    " << gas_phases[i] << "\n";
		}
	}
	oss << "  -kinetics " << "\n";
	{
		// kinetic reactions 
		const std::vector<std::string>& kin_reactions = phreeqc_rm.GetKineticReactions();
		for (size_t i = 0; i < phreeqc_rm.GetKineticReactionsCount(); i++)
		{
			oss << "    " << kin_reactions[i] << "\n";
		}
	}
	oss << "  -solid_solutions " << "\n";
	{
		// solid solutions 
		const std::vector<std::string>& ss_comps = phreeqc_rm.GetSolidSolutionComponents();
		const std::vector<std::string>& ss_names = phreeqc_rm.GetSolidSolutionNames();
		for (size_t i = 0; i < phreeqc_rm.GetSolidSolutionComponentsCount(); i++)
		{

			oss << "    ";
			oss.width(15);
			oss << std::left << ss_comps[i];
			oss << " # " << ss_names[i] << "\n";
		}
	}
	oss << "  -saturation_indices " << "\n";
	{
		// molalities of aqueous species 
		const std::vector<std::string>& si = phreeqc_rm.GetSINames();
		for (size_t i = 0; i < phreeqc_rm.GetSICount(); i++)
		{
			oss << "    " << si[i] << "\n";
		}
	}

	// Uncommenting the following line would define SELECTED_OUTPUT 2 with all species, reactants, and SIs
	// int status = phreeqc_rm.RunString(true, false, false, oss.str().c_str());
	std::cerr << oss.str();

	return(0);
}
void BMI_testing(PhreeqcRM& phreeqc_rm)
{
	std::ostringstream oss;
	//std::string BMI_GetComponentName() { return "PhreeqcRM_BMI"; }
	{
		oss << phreeqc_rm.BMI_GetComponentName() << "\n";
	}
	//double BMI_GetCurrentTime() { return this->GetTime(); }
	{
		oss << "BMI_GetCurrentTime: " << phreeqc_rm.BMI_GetCurrentTime() << "\n";
	}
	//double BMI_GetEndTime() { return this->GetTime() + this->GetTimeStep(); }
	{
		oss << "BMI_GetEndTime:     " << phreeqc_rm.BMI_GetEndTime() << "\n";
	}
	//std::string BMI_GetTimeUnits() { return "seconds"; };
	{
		oss << "BMI_GetTimeUnits:   " << phreeqc_rm.BMI_GetTimeUnits() << "\n";
	}
	//double BMI_GetTimeStep() { return this->GetTimeStep(); }
	{
		oss << "BMI_GetTimeStep:    " << phreeqc_rm.BMI_GetTimeStep() << "\n";
	}
	// std::vector<std::string> BMI_GetInputVarNames() { return this->bmi_input_vars; }
	// int BMI_GetInputItemCount() { return (int)this->bmi_input_vars.size(); }
	{
		std::vector<std::string> InputVarNames = phreeqc_rm.BMI_GetInputVarNames();
		int count = phreeqc_rm.BMI_GetInputItemCount();
		assert(InputVarNames.size() == (size_t)count);
		oss << "BMI_SetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << InputVarNames[i] << "\n";
			oss << "     Type:  " << phreeqc_rm.BMI_GetVarType(InputVarNames[i]) << "\n";
			oss << "     Units: " << phreeqc_rm.BMI_GetVarUnits(InputVarNames[i]) << "\n";
		}
	}
	// std::vector<std::string> BMI_GetOutputVarNames() { return this->bmi_output_vars; };
	// int BMI_GetOutputItemCount() { return (int)this->bmi_output_vars.size(); }
	{
		std::vector<std::string> OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
		int count = phreeqc_rm.BMI_GetOutputItemCount();
		assert(OutputVarNames.size() == (size_t)count);
		oss << "BMI_GetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << OutputVarNames[i] << "\n";
			oss << "     Type:  " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
			oss << "     Units: " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
		}
	}
	phreeqc_rm.OutputMessage(oss.str());
	// GetValue("Components")
	{
		int ncomps = -1;
		phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
		assert(ncomps == (size_t)phreeqc_rm.GetComponentCount());
		size_t nbytes = (size_t)phreeqc_rm.BMI_GetVarNbytes("Components");
		std::string all_comps(nbytes, ' ');
		phreeqc_rm.BMI_GetValue("Components", (void*) all_comps.data());
		size_t string_size = (size_t)phreeqc_rm.BMI_GetVarItemsize("Components");
		std::vector<std::string> bmi_comps;
		for (size_t i = 0; i < ncomps; i++)
		{
			std::string bmi_comp = all_comps.substr((i * string_size), string_size);
			size_t end = bmi_comp.find_last_not_of(' ');
			bmi_comp = (end == std::string::npos) ? "" : bmi_comp.substr(0, end + 1);
			bmi_comps.push_back(bmi_comp);
		}
		assert(bmi_comps == phreeqc_rm.GetComponents());
	}
	// GetValue("Concentrations")
	{
		int ncomps = -1;
		phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
		assert(ncomps == phreeqc_rm.GetComponentCount());
		int conc_nbytes = phreeqc_rm.BMI_GetVarNbytes("Concentrations");
		int conc_itemsize = phreeqc_rm.BMI_GetVarItemsize("Concentrations");
		int conc_dim = conc_nbytes / conc_itemsize;
		std::vector<double> rm_conc;
		phreeqc_rm.GetConcentrations(rm_conc);
		assert(conc_nbytes == (int)(rm_conc.size() * sizeof(double)));
		assert(conc_dim == (int)rm_conc.size());
		std::vector<double> bmi_conc(conc_dim, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Concentrations", bmi_conc.data());
		assert(rm_conc == bmi_conc);
	}
	// GetValue("Density")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", &ngrid);
		std::vector<double> bmi_density(ngrid, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Density", bmi_density.data());
		std::vector<double> rm_density;
		phreeqc_rm.GetDensity(rm_density);
		assert(bmi_density == rm_density);
	}
	// GetValue("ErrorString")
	{
		int nbytes = phreeqc_rm.BMI_GetVarNbytes("ErrorString");
		std::string bmi_err(nbytes, ' ');
		std::string rm_err = phreeqc_rm.GetErrorString();
		assert(bmi_err == rm_err);
	}
	// GetValue("Gfw")
	{
		int ncomps = -1;
		phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
		std::vector<double> bmi_gfw(ncomps, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Gfw", bmi_gfw.data());
		const std::vector<double> rm_gfw = phreeqc_rm.GetGfw();
		assert(bmi_gfw == rm_gfw);
	}
	// GetValue("GridCellCount")
	{
		int bmi_ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", &bmi_ngrid);
		int rm_ngrid = phreeqc_rm.GetGridCellCount();
		assert(bmi_ngrid == rm_ngrid);
	}
	// SetValue("NthSelectedOutput") 
	{
		// tested in "SelectedOutput"
		// tested in "SelectedOutputHeadings"
	}
	// GetValue("Pressure")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", &ngrid);
		std::vector<double> bmi_pressure(ngrid, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Pressure", bmi_pressure.data());
		const std::vector<double> rm_pressure = phreeqc_rm.GetPressure();
		assert(bmi_pressure == rm_pressure);
	}
	// GetValue("Saturation")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", &ngrid);
		std::vector<double> bmi_sat(ngrid, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Saturation", bmi_sat.data());
		std::vector<double> rm_sat;
		phreeqc_rm.GetSaturation(rm_sat);
		assert(bmi_sat == rm_sat);
	}
	// GetValue("SelectedOutput")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", &bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			phreeqc_rm.BMI_SetValue("NthSelectedOutput", &i);
			int bmi_nuser;
			phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", &bmi_nuser);
			int rm_nuser = phreeqc_rm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_col_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", &bmi_col_count);
			int rm_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			int bmi_row_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", &bmi_row_count);
			int rm_row_count = phreeqc_rm.GetSelectedOutputRowCount();
			assert(bmi_row_count == rm_row_count);

			int bmi_nbytes = phreeqc_rm.BMI_GetVarNbytes("SelectedOutput");
			int bmi_itemsize = phreeqc_rm.BMI_GetVarItemsize("SelectedOutput");
			int bmi_dim = bmi_nbytes / bmi_itemsize;
			std::vector<double> bmi_so(bmi_dim, INACTIVE_CELL_VALUE);
			phreeqc_rm.BMI_GetValue("SelectedOutput", bmi_so.data());
			std::vector<double> rm_so;
			phreeqc_rm.GetSelectedOutput(rm_so);
			assert(bmi_dim == (int)rm_so.size());
			assert(bmi_nbytes == (int)(rm_so.size() * sizeof(double)));
			assert(bmi_so == rm_so);
		}
	}
	// GetValue("SelectedOutputColumnCount")
	{
		int bmi_so_col_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", &bmi_so_col_count);
		int rm_so_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
		assert(bmi_so_col_count == rm_so_col_count);
	}
	// GetValue("SelectedOutputCount")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", &bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
	}
	// GetValue("SelectedOutputHeadings")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", &bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			phreeqc_rm.BMI_SetValue("NthSelectedOutput", &i);
			int bmi_nuser;
			phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", &bmi_nuser);
			int rm_nuser = phreeqc_rm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_nbytes = phreeqc_rm.BMI_GetVarNbytes("SelectedOutputHeadings");
			int bmi_string_size = phreeqc_rm.BMI_GetVarItemsize("SelectedOutputHeadings");
			std::string all_headings(bmi_nbytes, ' ');
			phreeqc_rm.BMI_GetValue("SelectedOutputHeadings", (void *)all_headings.c_str());

			int bmi_col_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", &bmi_col_count);
			int rm_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			for (int j = 0; j < bmi_col_count; j++)
			{
				std::string bmi_head = all_headings.substr((j * bmi_string_size), bmi_string_size);
				size_t end = bmi_head.find_last_not_of(' ');
				bmi_head = (end == std::string::npos) ? "" : bmi_head.substr(0, end + 1);
				std::string rm_head;
				phreeqc_rm.GetSelectedOutputHeading(j, rm_head);
				assert(bmi_head == rm_head);
			}
		}
	}
	// GetValue("SelectedOutputRowCount")
	{
		int bmi_so_row_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", &bmi_so_row_count);
		int rm_so_row_count = phreeqc_rm.GetSelectedOutputRowCount();
		assert(bmi_so_row_count == rm_so_row_count);
	}
	// GetValue("Temperature")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", &ngrid);
		std::vector<double> bmi_temp(ngrid, INACTIVE_CELL_VALUE);
		phreeqc_rm.BMI_GetValue("Temperature", bmi_temp.data());
		const std::vector<double> rm_temp = phreeqc_rm.GetTemperature();
		assert(bmi_temp == rm_temp);
	}
}
#ifdef SKIP
void GenerateYAML(int nxyz, std::string YAML_filename)
{
	YAMLPhreeqcRM yrm;

	// Set some properties
	yrm.YAMLSetErrorHandlerMode(1);
	yrm.YAMLSetComponentH2O(false);
	yrm.YAMLSetRebalanceFraction(0.5);
	yrm.YAMLSetRebalanceByCell(true);
	yrm.YAMLUseSolutionDensityVolume(false);
	yrm.YAMLSetPartitionUZSolids(false);
	// Open files
	yrm.YAMLSetFilePrefix("Advect_bmi_cpp");
	yrm.YAMLOpenFiles();
#ifdef SKIP
#ifdef USE_MPI
	// Optional callback for MPI
	int istatus = do_something(&some_data);                 // Root calls do_something, workers respond
#endif
#endif
	// Set concentration units
	yrm.YAMLSetUnitsSolution(2);           // 1, mg/L); 2, mol/L); 3, kg/kgs
	yrm.YAMLSetUnitsPPassemblage(1);       // 0, mol/L cell); 1, mol/L water); 2 mol/L rock
	yrm.YAMLSetUnitsExchange(1);           // 0, mol/L cell); 1, mol/L water); 2 mol/L rock
	yrm.YAMLSetUnitsSurface(1);            // 0, mol/L cell); 1, mol/L water); 2 mol/L rock
	yrm.YAMLSetUnitsGasPhase(1);           // 0, mol/L cell); 1, mol/L water); 2 mol/L rock
	yrm.YAMLSetUnitsSSassemblage(1);       // 0, mol/L cell); 1, mol/L water); 2 mol/L rock
	yrm.YAMLSetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

	// Set conversion from seconds to user units (days) Only affects one print statement
	double time_conversion = 1.0 / 86400;
	yrm.YAMLSetTimeConversion(time_conversion);
	// Set representative volume
	std::vector<double> rv(nxyz, 1.0);
	yrm.YAMLSetRepresentativeVolume(rv);
	// Set initial porosity
	std::vector<double> por(nxyz, 0.2);
	yrm.YAMLSetPorosity(por);
	// Set initial saturation
	std::vector<double> sat(nxyz, 1.0);
	yrm.YAMLSetSaturation(sat);
	// Set cells to print chemistry when print chemistry is turned on
	std::vector<int> print_chemistry_mask(nxyz, 0);
	for (int i = 0; i < nxyz / 2; i++)
	{
		print_chemistry_mask[i] = 1;
	}
	yrm.YAMLSetPrintChemistryMask(print_chemistry_mask);

	// Demonstation of mapping, two equivalent rows by symmetry
	std::vector<int> grid2chem;
	grid2chem.resize(nxyz, -1);
	for (int i = 0; i < nxyz / 2; i++)
	{
		grid2chem[i] = i;
		grid2chem[i + nxyz / 2] = i;
	}
	yrm.YAMLCreateMapping(grid2chem);

	// Set printing of chemistry file
	yrm.YAMLSetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
	// Load database
	yrm.YAMLLoadDatabase("phreeqc.dat");

	// Run file to define solutions and reactants for initial conditions, selected output
	bool workers = true;             // Worker instances do the reaction calculations for transport
	bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
	bool utility = true;             // Utility instance is available for processing
	yrm.YAMLRunFile(workers, initial_phreeqc, utility, "advect.pqi");

	// Clear contents of workers and utility
	initial_phreeqc = false;
	std::string input = "DELETE; -all";
	yrm.YAMLRunString(workers, initial_phreeqc, utility, input.c_str());
	// Determine number of components to transport
	yrm.YAMLFindComponents();
	// set array of initial conditions
	std::vector<int> ic1, ic2;
	ic1.resize(nxyz * 7, -1);
	ic2.resize(nxyz * 7, -1);
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
	yrm.YAMLInitialPhreeqc2Module(ic1, ic2, f1);
	// No mixing is defined, so the following is equivalent
	//yrm.YAMLInitialPhreeqc2Module(ic1);

	// alternative for setting initial conditions
	// cell number in first argument (-1 indicates last solution, 40 in this case)
	// in advect.pqi and any reactants with the same number--
	// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
	// will be written to cells 18 and 19 (0 based)
	std::vector<int> module_cells;
	module_cells.push_back(18);
	module_cells.push_back(19);
	yrm.YAMLInitialPhreeqcCell2Module(-1, module_cells);

	// Set density, temperature, and pressure
	std::vector<double> initial_density(nxyz, 1.0);
	std::vector<double> temperature(nxyz, 20.0);
	std::vector<double> pressure(nxyz, 2.0);
	yrm.YAMLSetDensity(initial_density);
	yrm.YAMLSetTemperature(temperature);
	yrm.YAMLSetPressure(pressure);

	// Initial equilibration of cells
	double time_step = 0.0;  // no kinetics
	yrm.YAMLSetTimeStep(time_step);
	double time = 0.0;
	yrm.YAMLSetTime(time);
	yrm.YAMLRunCells();
	yrm.YAMLSetTimeStep(86400);

	// Write YAML file
	{
		std::ostringstream oss;
		oss << yrm.GetYAMLDoc();
		std::ofstream ofs = std::ofstream(YAML_filename.c_str(), std::ofstream::out);
		ofs << oss.str();
		ofs.close();

		//phreeqc_rm.BMI_Initialize(yaml.str()); // TODO TODO TODO
		//for (YAML::const_iterator it = ref.begin();it != ref.end();++it)
		//{
		//	std::cerr << "  " << it->first << std::endl;
		//	std::cerr << "  " << it->second << std::endl;
		//	std::cerr << std::endl;
		//}
		yrm.clear();
	}
};
#endif