#ifdef USE_YAML
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

void advectionbmi_cpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);
int bmi_example_selected_output(PhreeqcRM& phreeqc_rm);
double bmi_basic_callback(double x1, double x2, const char* str, void* cookie);
void BMI_testing(PhreeqcRM& phreeqc_rm);
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

int AdvectBMI_cpp()
{

	// --------------------------------------------------------------------------
	// Write file to initialize PhreeqcRM
	// Based on PHREEQC Example 11
	// --------------------------------------------------------------------------
	try
	{
		std::string yaml_file("AdvectBMI_cpp.yaml");
		// GetGridCellCountYAML must be called BEFORE
		// the PhreeqcRM instance is created. The
		// return value can be used to create the
		// PhreeqcRM instance.
		// 
		// If the YAML file does not contain 
		// a node "SetGridCellCount:" (usually written
		// using the YAMLPhreeqcRM class and the method
		// YAMLSetGridCellCount), the return
		// value is zero.
		int nxyz = GetGridCellCountYAML(yaml_file.c_str());

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
		phreeqc_rm.BMI_SetValue("FilePrefix", "AdvectBMI_cpp");
		phreeqc_rm.OpenFiles();
		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		bmi_register_basic_callback(&some_data);
		// Use YAML file to initialize
		phreeqc_rm.BMI_Initialize(yaml_file);
		// Get number of components
		int ncomps;
		phreeqc_rm.BMI_GetValue("ComponentCount", ncomps);
		// Example for generating automatic selected-output definitions
		bmi_example_selected_output(phreeqc_rm);
		// Print some of the reaction module information
		{
			std::ostringstream oss;
			oss << "Database:                                         " << phreeqc_rm.GetDatabaseFileName().c_str() << "\n";
			oss << "Number of threads:                                " << phreeqc_rm.GetThreadCount() << "\n";
			oss << "Number of MPI processes:                          " << phreeqc_rm.GetMpiTasks() << "\n";
			oss << "MPI task number:                                  " << phreeqc_rm.GetMpiMyself() << "\n";
			std::string prefix;
			phreeqc_rm.BMI_GetValue("FilePrefix", prefix);
			oss << "File prefix:                                      " << prefix << "\n";
			int ngrid;
			phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
			oss << "Number of grid cells in the user's model:         " << ngrid << "\n";
			oss << "Number of chemistry cells in the reaction module: " << phreeqc_rm.GetChemistryCellCount() << "\n";
			int ncomps;
			phreeqc_rm.BMI_GetValue("ComponentCount", ncomps);
			oss << "Number of components for transport:               " << ncomps << "\n";
			oss << "Partioning of UZ solids:                          " << phreeqc_rm.GetPartitionUZSolids() << "\n";
			oss << "Error handler mode:                               " << phreeqc_rm.GetErrorHandlerMode() << "\n";
			phreeqc_rm.OutputMessage(oss.str());
		}
		const std::vector<int>& f_map = phreeqc_rm.GetForwardMapping();
		// Get components
		std::vector<std::string> components;
		phreeqc_rm.BMI_GetValue("Components", components);
		// Get gfw
		std::vector<double> gfw(ncomps, 0);
		phreeqc_rm.BMI_GetValue("Gfw", gfw);
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
		std::vector<double> tempc;
		phreeqc_rm.BMI_GetValue("Temperature", tempc);
		// Get initial saturation
		IRM_RESULT status;
		std::vector<double> sat;
		phreeqc_rm.BMI_GetValue("Saturation", sat);
		//Get initial porosity
		std::vector<double> por;
		phreeqc_rm.BMI_GetValue("Porosity", por);
		//Get initial solution volume
		std::vector<double> volume;
		phreeqc_rm.BMI_GetValue("SolutionVolume", volume);
		//Get initial concentrations
		std::vector<double> c;
		phreeqc_rm.BMI_GetValue("Concentrations", c);

		// Set density, temperature, and pressure
		std::vector<double> density(nxyz, 1.0);
		std::vector<double> temperature(nxyz, 20.0);
		std::vector<double> pressure(nxyz, 2.0);
		phreeqc_rm.BMI_SetValue("Density", density);
		phreeqc_rm.BMI_SetValue("Temperature", temperature);
		phreeqc_rm.BMI_SetValue("Pressure", pressure);
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
		double time = 0.0;
		phreeqc_rm.BMI_SetValue("Time", time);
		double time_step = 86400;
		phreeqc_rm.BMI_SetValue("TimeStep", time_step);

		for (int steps = 0; steps < nsteps; steps++)
		{
			{
				std::ostringstream strm;
				strm << "Beginning transport calculation             " << phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion() << " days\n";
				strm << "          Time step                         " << phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion() << " days\n";
				phreeqc_rm.LogMessage(strm.str());
				phreeqc_rm.SetScreenOn(true);
				phreeqc_rm.ScreenMessage(strm.str());
			}
			// Transport calculation here
			advectionbmi_cpp(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			phreeqc_rm.BMI_SetValue("SelectedOutputOn", print_selected_output_on);
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			phreeqc_rm.BMI_SetValue("Concentrations", c);        // Transported concentrations
			phreeqc_rm.BMI_SetValue("Porosity", por);            // If pororosity changes due to compressibility
			phreeqc_rm.BMI_SetValue("Saturation", sat);          // If saturation changes
			phreeqc_rm.BMI_SetValue("Temperature", temperature); // If temperature changes
			phreeqc_rm.BMI_SetValue("Pressure", pressure);       // If pressure changes 
			phreeqc_rm.BMI_SetValue("TimeStep", time_step);      // Time step for kinetic reactions
			time += time_step;
			phreeqc_rm.BMI_SetValue("Time", time);
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
			//Run chemistry
			phreeqc_rm.BMI_Update();
			// Transfer data from PhreeqcRM for transport
			phreeqc_rm.BMI_GetValue("Concentrations", c);
			phreeqc_rm.BMI_GetValue("Density", density);
			phreeqc_rm.BMI_GetValue("SolutionVolume", volume);
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
				phreeqc_rm.BMI_GetValue("SelectedOutputCount", selected_output_count);
				for (int isel = 0; isel < selected_output_count; isel++)
				{
					std::ostringstream oss;
					// Loop through possible multiple selected output definitions
					phreeqc_rm.BMI_SetValue("NthSelectedOutput", isel);
					oss << "Selected output sequence number: " << isel << "\n";
					int n_user;
					phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", n_user);
					oss << "Selected output user number:     " << n_user << "\n";
					// Get array of selected output values
					int col;
					phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", col);
					int bmi_row_count;
					phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", bmi_row_count);

					// Get selected output
					std::vector<double> so;
					phreeqc_rm.BMI_GetValue("SelectedOutput", so);
					// Get headings
					std::vector<std::string> headings;
					phreeqc_rm.BMI_GetValue("SelectedOutputHeadings", headings);
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
						oss << "     Selected output: " << "\n";
						for (int j = 0; j < col; j++)
						{
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
		util_ptr->SetOutputFileName("AdvectBMI_cpp_utility.txt");
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
		status = phreeqc_rm.SetDumpFileName("AdvectBMI_cpp.dmp");
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
advectionbmi_cpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
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
		status = phreeqc_rm.SetFilePrefix("AdvectBMI_cpp");
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
		//status = phreeqc_rm.SetFilePrefix("Units_InitialPhreeqc_2");
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
		status = phreeqc_rm.SetFilePrefix("AdvectBMI_cpp_units_worker");
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
		util_ptr->SetOutputFileName("AdvectBMI_cpp_units_utility.txt");
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
	// ComponentName
	oss << phreeqc_rm.BMI_GetComponentName() << "\n";
	// Time
	{
		double rm_time = 3600.0;
		phreeqc_rm.BMI_SetValue("Time", rm_time);
		double bmi_time = 0.0;
		phreeqc_rm.BMI_GetValue("Time", bmi_time);
		assert(rm_time == bmi_time);
		rm_time = phreeqc_rm.GetTime();
		assert(rm_time == bmi_time);
		bmi_time = phreeqc_rm.BMI_GetCurrentTime();
		assert(rm_time == bmi_time);
	}
	// TimeStep
	{
		double rm_timestep = 60.;
		phreeqc_rm.BMI_SetValue("TimeStep", rm_timestep);
		double bmi_timestep = 0.0;
		phreeqc_rm.BMI_GetValue("TimeStep", bmi_timestep);
		assert(rm_timestep == bmi_timestep);
		rm_timestep = phreeqc_rm.GetTimeStep();
		assert(rm_timestep == bmi_timestep);
		bmi_timestep = phreeqc_rm.BMI_GetTimeStep();
		assert(rm_timestep == bmi_timestep);
	}
	// EndTime
	oss << "BMI_GetEndTime:     " << phreeqc_rm.BMI_GetEndTime() << "\n";
	// TimeUnits
	oss << "BMI_GetTimeUnits:   " << phreeqc_rm.BMI_GetTimeUnits() << "\n";
	// TimeStep
	oss << "BMI_GetTimeStep:    " << phreeqc_rm.BMI_GetTimeStep() << "\n";
	// InputVarNames
	{
		std::vector<std::string> InputVarNames = phreeqc_rm.BMI_GetInputVarNames();
		int count = phreeqc_rm.BMI_GetInputItemCount();
		assert(InputVarNames.size() == (size_t)count);
		std::vector<std::string> InputVarNames2;
		phreeqc_rm.BMI_GetValue("InputVarNames", InputVarNames2);
		assert(InputVarNames == InputVarNames2);
		oss << "BMI_SetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << InputVarNames[i] << "\n";
			oss << "     Type:  " << phreeqc_rm.BMI_GetVarType(InputVarNames[i]) << "\n";
			oss << "     Units: " << phreeqc_rm.BMI_GetVarUnits(InputVarNames[i]) << "\n";
			oss << "  Itemsize: " << phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
			oss << "    NBytes: " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) << "\n";
		}
	}
	// OutputVarNames
	{
		std::vector<std::string> OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
		int count = phreeqc_rm.BMI_GetOutputItemCount();
		assert(OutputVarNames.size() == (size_t)count);
		std::vector<std::string> OutputVarNames2;
		phreeqc_rm.BMI_GetValue("OutputVarNames", OutputVarNames2);
		assert(OutputVarNames == OutputVarNames2);
		oss << "BMI_GetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << OutputVarNames[i] << "\n";
			oss << "     Type:  " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
			oss << "     Units: " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
			oss << "  Itemsize: " << phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
			oss << "    NBytes: " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) << "\n";
		}
	}
	phreeqc_rm.OutputMessage(oss.str());
	// GetValue("Components")
	{
		int ncomps;
		phreeqc_rm.BMI_GetValue("ComponentCount", ncomps);
		assert(ncomps == (size_t)phreeqc_rm.GetComponentCount());
		std::vector<std::string> bmi_comps;
		phreeqc_rm.BMI_GetValue("Components", bmi_comps);
		std::vector<std::string> rm_comps;
		rm_comps = phreeqc_rm.GetComponents();
		assert(bmi_comps == rm_comps);
	}
	// GetValue("Concentrations")
	{
		int ncomps;
		phreeqc_rm.BMI_GetValue("ComponentCount", ncomps);
		int conc_nbytes = phreeqc_rm.BMI_GetVarNbytes("Concentrations");
		int conc_itemsize = phreeqc_rm.BMI_GetVarItemsize("Concentrations");
		int conc_dim = conc_nbytes / conc_itemsize;
		std::vector<double> bmi_conc;
		phreeqc_rm.BMI_GetValue("Concentrations", bmi_conc);
		assert(conc_dim == (int)bmi_conc.size());
		std::vector<double> rm_conc;
		phreeqc_rm.GetConcentrations(rm_conc);
		assert(rm_conc == bmi_conc);
	}
	// GetValue("Density")
	// GetDensity and BMI_GetValue("Density) always return 
	// the calculated solution density
	{
		std::vector<double> rm_density;
		phreeqc_rm.GetDensity(rm_density);
		std::vector<double> bmi_density;
		phreeqc_rm.BMI_GetValue("Density", bmi_density);
		assert(bmi_density == rm_density);
	}
	// GetValue("ErrorString")
	{
		std::string bmi_err;
		phreeqc_rm.BMI_GetValue("ErrorString", bmi_err);
		std::string rm_err = phreeqc_rm.GetErrorString();
		assert(bmi_err == rm_err);
	}
	//FilePrefix
	{
		std::string str = "NewPrefix";
		phreeqc_rm.BMI_SetValue("FilePrefix", str);
		std::string rm_prefix = phreeqc_rm.GetFilePrefix();
		assert(str == rm_prefix);
		std::string prefix;
		phreeqc_rm.BMI_GetValue("FilePrefix", prefix);
		assert(prefix == rm_prefix);
	}
	// GetValue("Gfw")
	{
		std::vector<double> bmi_gfw;
		phreeqc_rm.BMI_GetValue("Gfw", bmi_gfw);
		std::vector<double> rm_gfw = phreeqc_rm.GetGfw();
		assert(bmi_gfw == rm_gfw);
	}
	// GetValue("GridCellCount")
	{
		int bmi_ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", bmi_ngrid);
		int rm_ngrid = phreeqc_rm.GetGridCellCount();
		assert(bmi_ngrid == rm_ngrid);
	}
	// SetValue("NthSelectedOutput") 
	{
		// tested in "SelectedOutput"
		// tested in "SelectedOutputHeadings"
	}
	// GetValue("Porosity")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_porosity(ngrid, 0.3);
		phreeqc_rm.BMI_SetValue("Porosity", bmi_porosity);
		std::vector<double> rm_porosity;
		rm_porosity = phreeqc_rm.GetPorosity();
		assert(bmi_porosity == rm_porosity);
		phreeqc_rm.BMI_GetValue("Porosity", bmi_porosity);
		assert(bmi_porosity == rm_porosity);
	}
	// GetValue("Pressure")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_pressure(ngrid, 1.1);
		phreeqc_rm.BMI_SetValue("Pressure", bmi_pressure);
		std::vector<double> rm_pressure;
		rm_pressure = phreeqc_rm.GetPressure();
		assert(bmi_pressure == rm_pressure);
		phreeqc_rm.BMI_GetValue("Pressure", bmi_pressure);
		assert(bmi_pressure == rm_pressure);
	}
	// GetValue("Saturation")
	// Always returns solution_volume / (rv * porosity) for each cell
	{
		std::vector<double> bmi_sat;
		phreeqc_rm.BMI_GetValue("Saturation", bmi_sat);
		std::vector<double> rm_sat;
		phreeqc_rm.GetSaturation(rm_sat);
		assert(bmi_sat == rm_sat);
	}
	// GetValue("SolutionVolume")
	// Always returns solution_volume / (rv * porosity) for each cell
	{
		std::vector<double> bmi_vol;
		phreeqc_rm.BMI_GetValue("SolutionVolume", bmi_vol);
		std::vector<double> rm_vol;
		rm_vol = phreeqc_rm.GetSolutionVolume();
		assert(bmi_vol == rm_vol);
	}
	// GetValue("Temperature")
	{
		int ngrid;
		phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_temperature(ngrid, 50.);
		phreeqc_rm.BMI_SetValue("Temperature", bmi_temperature);
		std::vector<double> rm_temperature;
		rm_temperature = phreeqc_rm.GetTemperature();
		assert(bmi_temperature == rm_temperature);
		phreeqc_rm.BMI_GetValue("Temperature", bmi_temperature);
		assert(bmi_temperature == rm_temperature);
	}
	// GetValue("SelectedOutput")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			phreeqc_rm.BMI_SetValue("NthSelectedOutput", i);
			int bmi_nuser;
			phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", bmi_nuser);
			int rm_nuser = phreeqc_rm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_col_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", bmi_col_count);
			int rm_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			int bmi_row_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", bmi_row_count);
			int rm_row_count = phreeqc_rm.GetSelectedOutputRowCount();
			assert(bmi_row_count == rm_row_count);

			std::vector<double> bmi_so;
			phreeqc_rm.BMI_GetValue("SelectedOutput", bmi_so);
			std::vector<double> rm_so;
			phreeqc_rm.GetSelectedOutput(rm_so);
			assert(bmi_so == rm_so);
		}
	}
	// GetValue("SelectedOutputColumnCount")
	{
		int bmi_so_col_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", bmi_so_col_count);
		int rm_so_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
		assert(bmi_so_col_count == rm_so_col_count);
	}
	// GetValue("SelectedOutputCount")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
	}
	// GetValue("SelectedOutputHeadings")
	{
		int bmi_so_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = phreeqc_rm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			phreeqc_rm.BMI_SetValue("NthSelectedOutput", i);
			int bmi_nuser;
			phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", bmi_nuser);
			int rm_nuser = phreeqc_rm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_col_count;
			phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", bmi_col_count);
			int rm_col_count = phreeqc_rm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			std::vector< std::string > bmi_headings;
			phreeqc_rm.BMI_GetValue("SelectedOutputHeadings", bmi_headings);
			std::vector< std::string > rm_headings;
			IRM_RESULT status = phreeqc_rm.GetSelectedOutputHeadings(rm_headings);
			assert(bmi_headings == rm_headings);
		}
	}
	// GetValue("SelectedOutputRowCount")
	{
		int bmi_so_row_count;
		phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", bmi_so_row_count);
		int rm_so_row_count = phreeqc_rm.GetSelectedOutputRowCount();
		assert(bmi_so_row_count == rm_so_row_count);
	}
}
#endif // YAML