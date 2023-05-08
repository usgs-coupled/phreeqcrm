#ifdef USE_YAML
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "BMIPhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
#include <algorithm>
#include "yaml-cpp/yaml.h"
#include "YAMLPhreeqcRM.h"

class Ptrs 
{
public:
	Ptrs() 
	{
		brm = NULL;
		ComponentCount_ptr = NULL;
		Concentrations_ptr = NULL;
		Density_ptr = NULL;
		Gfw_ptr = NULL;
		GridCellCount_ptr = NULL;
		Saturation_ptr = NULL;
		SolutionVolume_ptr = NULL;
		Time_ptr = NULL;
		TimeStep_ptr = NULL;
		Porosity_ptr = NULL;
		Pressure_ptr = NULL;
		SelectedOutputOn_ptr = NULL;
		Temperature_ptr = NULL;
		Viscosity_ptr = NULL;
	};
	BMIPhreeqcRM* brm;
	int* ComponentCount_ptr;
	double* Concentrations_ptr;
	double* Density_ptr;
	double* Gfw_ptr;
	int* GridCellCount_ptr;
	double* Saturation_ptr;
	double* SolutionVolume_ptr;
	double* Time_ptr;
	double* TimeStep_ptr;
	double* Porosity_ptr;
	double* Pressure_ptr;
	bool* SelectedOutputOn_ptr;
	double* Temperature_ptr;
	double* Viscosity_ptr;
};

void advectionbmi_cpp_test(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);
void testing(BMIPhreeqcRM& brm, Ptrs ptrs);
void compare_ptrs(struct Ptrs ptrs);
void testbmi_cpp();
int AdvectBMI_cpp_test()
{

	// --------------------------------------------------------------------------
	// Write file to initialize PhreeqcRM
	// Based on PHREEQC Example 11
	// --------------------------------------------------------------------------
	try
	{
		std::string yaml_file("AdvectBMI_cpp_test.yaml");
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
		//int nxyz = BMIPhreeqcRM::GetGridCellCountYAML(yaml_file.c_str());
		//int nxyz = 40;

#ifdef USE_MPI
		// MPI
		BMIPhreeqcRM brm(nxyz, MPI_COMM_WORLD);
		some_data.brm_ptr = &brm;
		MP_TYPE comm = MPI_COMM_WORLD;
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			brm.SetMpiWorkerCallbackC(bmi_worker_tasks_cc);
			brm.SetMpiWorkerCallbackCookie(&some_data);
			brm.MpiWorker();
			return EXIT_SUCCESS;
		}
#else
		// OpenMP
		BMIPhreeqcRM brm;
#endif
		// Use YAML file to initialize
		brm.Initialize(yaml_file);
		int nxyz;
		brm.GetValue("GridCellCount", nxyz);
		// set pointers
		Ptrs ptrs;
		ptrs.brm = &brm;
		ptrs.ComponentCount_ptr = (int*)brm.GetValuePtr("ComponentCount");
		ptrs.Concentrations_ptr = (double*)brm.GetValuePtr("Concentrations");
		ptrs.Density_ptr = (double*)brm.GetValuePtr("Density");
		ptrs.Gfw_ptr = (double*)brm.GetValuePtr("Gfw");
		ptrs.GridCellCount_ptr = (int*)brm.GetValuePtr("GridCellCount");
		ptrs.Saturation_ptr = (double*)brm.GetValuePtr("Saturation");
		ptrs.SolutionVolume_ptr = (double*)brm.GetValuePtr("SolutionVolume");
		ptrs.Time_ptr = (double*)brm.GetValuePtr("Time");
		ptrs.TimeStep_ptr = (double*)brm.GetValuePtr("TimeStep");
		ptrs.Porosity_ptr = (double*)brm.GetValuePtr("Porosity");
		ptrs.Pressure_ptr = (double*)brm.GetValuePtr("Pressure");
		ptrs.SelectedOutputOn_ptr = (bool*)brm.GetValuePtr("SelectedOutputOn");
		ptrs.Temperature_ptr = (double*)brm.GetValuePtr("Temperature");
		ptrs.Viscosity_ptr = (double*)brm.GetValuePtr("Viscosity");
		{
			std::ostringstream oss;
			// InputVarNames
			std::vector<std::string> InputVarNames = brm.GetInputVarNames();
			oss << "SetValues variables and units:\n";
			for (size_t i = 0; i < InputVarNames.size(); i++)
			{
				oss << "  " << i << "  " << std::setw(50) << InputVarNames[i] << "   "
					<< brm.GetVarUnits(InputVarNames[i]) << "\n";
			}
			// OutputVarNames
			std::vector<std::string> OutputVarNames = brm.GetOutputVarNames();
			oss << "GetValues variables and units:\n";
			for (size_t i = 0; i < OutputVarNames.size(); i++)
			{
				oss << "  " << i << "  " << std::setw(50) << OutputVarNames[i] << "   "
					<< brm.GetVarUnits(OutputVarNames[i]) << "\n";
			}
			oss << std::endl;
			brm.OutputMessage(oss.str());
		}
		// Get number of components
		int ncomps;
		brm.GetValue("ComponentCount", ncomps);
		// Print some of the reaction module information
		{
			std::ostringstream oss;
			oss << "Database:                                         " << brm.GetDatabaseFileName().c_str() << "\n";
			oss << "Number of threads:                                " << brm.GetThreadCount() << "\n";
			oss << "Number of MPI processes:                          " << brm.GetMpiTasks() << "\n";
			oss << "MPI task number:                                  " << brm.GetMpiMyself() << "\n";
			std::string prefix;
			brm.GetValue("FilePrefix", prefix);
			oss << "File prefix:                                      " << prefix << "\n";
			int ngrid;
			brm.GetValue("GridCellCount", ngrid);
			oss << "Number of grid cells in the user's model:         " << ngrid << "\n";
			oss << "Number of chemistry cells in the reaction module: " << brm.GetChemistryCellCount() << "\n";
			int ncomps;
			brm.GetValue("ComponentCount", ncomps);
			oss << "Number of components for transport:               " << ncomps << "\n";
			oss << "Partioning of UZ solids:                          " << brm.GetPartitionUZSolids() << "\n";
			oss << "Error handler mode:                               " << brm.GetErrorHandlerMode() << std::endl;
			brm.OutputMessage(oss.str());
		}
		const std::vector<int>& f_map = brm.GetForwardMapping();
		// Get components
		std::vector<std::string> components;
		brm.GetValue("Components", components);
		// Get gfw
		std::vector<double> gfw(ncomps, 0);
		brm.GetValue("Gfw", gfw);
		// write component information
		for (int i = 0; i < ncomps; i++)
		{
			std::ostringstream strm;
			strm.width(10);
			strm << components[i] << "    " << gfw[i] << "\n";
			brm.OutputMessage(strm.str());
		}
		brm.OutputMessage("\n");
		// Get temperatures (unchanged by PhreeqcRM)
		std::vector<double> tempc;
		brm.GetValue("Temperature", tempc);
		// Get initial saturation
		IRM_RESULT status;
		std::vector<double> sat;
		brm.GetValue("Saturation", sat);
		//Get initial porosity
		std::vector<double> por;
		brm.GetValue("Porosity", por);
		//Get initial solution volume
		std::vector<double> volume;
		brm.GetValue("SolutionVolume", volume);
		//Get initial concentrations
		std::vector<double> c;
		brm.GetValue("Concentrations", c);

		// Set density, temperature, and pressure
		std::vector<double> density(nxyz, 1.0);
		std::vector<double> temperature(nxyz, 20.0);
		std::vector<double> pressure(nxyz, 2.0);
		brm.SetValue("Density", density);
		brm.SetValue("Temperature", temperature);
		brm.SetValue("Pressure", pressure);
		// --------------------------------------------------------------------------
		// Set boundary condition
		// --------------------------------------------------------------------------
		std::vector<double> bc_conc, bc_f1;
		std::vector<int> bc1, bc2;
		int nbound = 1;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
		bc2.resize(nbound, -1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = brm.InitialPhreeqc2Concentrations(bc_conc, bc1, bc2, bc_f1);
		// --------------------------------------------------------------------------
		// Transient loop
		// --------------------------------------------------------------------------
		int nsteps = 10;
		double time = 0.0;
		brm.SetValue("Time", time);
		double time_step = 86400;
		brm.SetValue("TimeStep", time_step);

		for (int steps = 0; steps < nsteps; steps++)
		{
			{
				std::ostringstream strm;
				strm << "Beginning transport calculation             " << brm.GetTime() * brm.GetTimeConversion() << " days\n";
				strm << "          Time step                         " << brm.GetTimeStep() * brm.GetTimeConversion() << " days\n";
				brm.LogMessage(strm.str());
				brm.SetScreenOn(true);
				brm.ScreenMessage(strm.str());
			}
			compare_ptrs(ptrs);
			// Transport calculation here
			advectionbmi_cpp_test(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			brm.SetValue("SelectedOutputOn", print_selected_output_on);
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = brm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			brm.SetValue("Concentrations", c);        // Transported concentrations
			brm.SetValue("Porosity", por);            // If pororosity changes due to compressibility
			brm.SetValue("Saturation", sat);          // If saturation changes
			brm.SetValue("Temperature", temperature); // If temperature changes
			brm.SetValue("Pressure", pressure);       // If pressure changes 
			brm.SetValue("TimeStep", time_step);      // Time step for kinetic reactions
			time += time_step;
			brm.SetValue("Time", time);
			// Run cells with transported conditions
			{
				std::ostringstream strm;
				strm << "Beginning reaction calculation              " << time * brm.GetTimeConversion() << " days\n";
				brm.LogMessage(strm.str());
				brm.ScreenMessage(strm.str());
			}
			compare_ptrs(ptrs);
			//Run chemistry
			brm.Update();
			compare_ptrs(ptrs);
			// Transfer data from PhreeqcRM for transport
			brm.GetValue("Concentrations", c);
			brm.GetValue("Density", density);
			brm.GetValue("SolutionVolume", volume);
			// Print results at last time step
			if (print_chemistry_on != 0)
			{
				{
					std::ostringstream oss;
					oss << "Current distribution of cells for workers\n";
					oss << "Worker      First cell        Last Cell\n";
					int n;
					n = brm.GetThreadCount() * brm.GetMpiTasks();
					for (int i = 0; i < n; i++)
					{
						oss << i << "           " << brm.GetStartCell()[i] << "                 "
							<< brm.GetEndCell()[i] << "\n";
					}
					brm.OutputMessage(oss.str());
				}
				int selected_output_count;
				brm.GetValue("SelectedOutputCount", selected_output_count);
				for (int isel = 0; isel < selected_output_count; isel++)
				{
					std::ostringstream oss;
					// Loop through possible multiple selected output definitions
					brm.SetValue("NthSelectedOutput", isel);
					oss << "Selected output sequence number: " << isel << "\n";
					int n_user;
					brm.GetValue("CurrentSelectedOutputUserNumber", n_user);
					oss << "Selected output user number:     " << n_user << "\n";
					// Get array of selected output values
					int col;
					brm.GetValue("SelectedOutputColumnCount", col);
					int bmi_row_count;
					brm.GetValue("SelectedOutputRowCount", bmi_row_count);

					// Get selected output
					std::vector<double> so;
					brm.GetValue("SelectedOutput", so);
					// Get headings
					std::vector<std::string> headings;
					brm.GetValue("SelectedOutputHeadings", headings);
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
					brm.OutputMessage(oss.str());

					{
						// Use GetValue to extract exchange composition and pH
						// YAMLAddOutputVars was called in YAML
						// to select additional OutputVarNames variables
						std::vector<double> CaX2, KX, NaX, pH, SAR;
						brm.GetValue("solution_ph", pH);
						brm.GetValue("exchange_X_species_log_molality_CaX2", CaX2);
						brm.GetValue("exchange_X_species_log_molality_KX", KX);
						brm.GetValue("exchange_X_species_log_molality_NaX", NaX);
						brm.GetValue("calculate_value_sar", SAR);
						std::ostringstream xoss;
						xoss << "      pH      CaX2     KX       NaX       SAR\n";
						for (size_t i = 0; i < nxyz; i++)
						{
							xoss << std::setw(10) << std::setprecision(5) << std::fixed <<
								pH[i] << "  " << std::pow(10.0, CaX2[i]) << "  " << 
								std::pow(10.0, KX[i]) << "  " << std::pow(10.0, NaX[i]) << 
								"  " << SAR[i] << "\n";
						}
						std::cerr << xoss.str();
					}
				}
			}
		}
		testing(brm, ptrs); // Tests GetValue
		// Clean up
		brm.Finalize();
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
advectionbmi_cpp_test(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim)
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
		PhreeqcRM brm(nxyz, MPI_COMM_WORLD);
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself > 0)
		{
			brm.MpiWorker();
			return EXIT_SUCCESS;
		}
#else
		// OpenMP
		int nthreads = 3;
		PhreeqcRM brm(nxyz, nthreads);
#endif
		IRM_RESULT status;
		// Set properties
		status = brm.SetErrorOn(true);
		status = brm.SetErrorHandlerMode(1);
		status = brm.SetFilePrefix("Advectcpp");
		if (brm.GetMpiMyself() == 0)
		{
			brm.OpenFiles();
		}
		// Set concentration units
		status = brm.SetUnitsSolution(1);      // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = brm.SetUnitsPPassemblage(2);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = brm.SetUnitsExchange(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = brm.SetUnitsSurface(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = brm.SetUnitsGasPhase(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = brm.SetUnitsSSassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = brm.SetUnitsKinetics(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		// Set representative volume
		std::vector<double> rv;
		rv.resize(nxyz, 1.0);
		status = brm.SetRepresentativeVolume(rv);
		// Set current porosity
		std::vector<double> por;
		por.resize(nxyz, 0.2);
		status = brm.SetPorosity(por);
		// Set saturation
		std::vector<double> sat;
		sat.resize(nxyz, 1.0);
		status = brm.SetSaturation(sat);
		// Set printing of chemistry file
		status = brm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility

		// --------------------------------------------------------------------------
		// Set initial conditions
		// --------------------------------------------------------------------------

		// Load database
		status = brm.LoadDatabase("phreeqc.dat");
		// Run file to define solutions and reactants for initial conditions, selected output
		bool workers = true;
		bool initial_phreeqc = true;
		bool utility = false;
		status = brm.RunFile(workers, initial_phreeqc, utility, "units.pqi");
		{
			std::string input = "DELETE; -all";
			status = brm.RunString(true, false, true, input.c_str());
		}
		//status = brm.SetFilePrefix("Units_InitialPhreeqc_2");
		if (brm.GetMpiMyself() == 0)
		{
			brm.OpenFiles();
		}
		// Set reference to components
		int ncomps = brm.FindComponents();
		const std::vector<std::string>& components = brm.GetComponents();
		// Set initial conditions
		std::vector < int > cell_numbers;
		cell_numbers.push_back(0);
		status = brm.InitialPhreeqcCell2Module(1, cell_numbers);
		cell_numbers[0] = 1;
		status = brm.InitialPhreeqcCell2Module(2, cell_numbers);
		cell_numbers[0] = 2;
		status = brm.InitialPhreeqcCell2Module(3, cell_numbers);
		// Retrieve concentrations
		std::vector<double> c;
		status = brm.SetFilePrefix("Advectcpp_units_worker");
		if (brm.GetMpiMyself() == 0)
		{
			brm.OpenFiles();
		}
		std::vector < int > print_mask;
		print_mask.resize(3, 1);
		brm.SetPrintChemistryMask(print_mask);
		brm.SetPrintChemistryOn(true, true, true);
		status = brm.RunCells();
		status = brm.GetConcentrations(c);
		std::vector<double> so;
		status = brm.GetSelectedOutput(so);
		std::vector<std::string> headings;
		{
			std::string heading;
			brm.GetSelectedOutputHeading(0, heading);
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
		IPhreeqc* util_ptr = brm.Concentrations2Utility(c, tc, p_atm);
		std::string input;
		input = "RUN_CELLS; -cells 0-2";
		// Output goes to new file
		int iphreeqc_result;
		util_ptr->SetOutputFileName("Advectcpp_units_utility.txt");
		util_ptr->SetOutputFileOn(true);
		iphreeqc_result = util_ptr->RunString(input.c_str());
		status = brm.MpiWorkerBreak();
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

void testing(BMIPhreeqcRM& brm, Ptrs ptrs)
{
	std::ostringstream oss;

	// ComponentName
	oss << brm.GetComponentName() << "\n";
	// Time
	{
		double rm_time = 3600.0;
		brm.SetValue("Time", rm_time);
		double bmi_time = 0.0;
		brm.GetValue("Time", bmi_time);
		assert(rm_time == bmi_time);
		rm_time = brm.GetTime();
		assert(rm_time == bmi_time);
		bmi_time = brm.GetCurrentTime();
		assert(rm_time == bmi_time);
		brm.SetTime(3601.0);
	}
	// TimeStep
	{
		double rm_timestep = 60.;
		brm.SetValue("TimeStep", rm_timestep);
		double bmi_timestep = 0.0;
		brm.GetValue("TimeStep", bmi_timestep);
		assert(rm_timestep == bmi_timestep);
		rm_timestep = brm.GetTimeStep();
		assert(rm_timestep == bmi_timestep);
		bmi_timestep = brm.GetTimeStep();
		assert(rm_timestep == bmi_timestep);
		brm.SetTime(61.);
	}
	compare_ptrs(ptrs);
	// EndTime
	oss << "GetEndTime:     " << brm.GetEndTime() << "\n";
	// TimeUnits
	oss << "GetTimeUnits:   " << brm.GetTimeUnits() << "\n";
	// TimeStep
	oss << "GetTimeStep:    " << brm.GetTimeStep() << "\n";
	// InputVarNames
	{
		std::vector<std::string> InputVarNames = brm.GetInputVarNames();
		int count = brm.GetInputItemCount();
		assert(InputVarNames.size() == (size_t)count);
		oss << "SetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << InputVarNames[i] << "\n";
			oss << "     Type:  " << brm.GetVarType(InputVarNames[i]) << "\n";
			oss << "     Units: " << brm.GetVarUnits(InputVarNames[i]) << "\n";
			oss << "  Itemsize: " << brm.GetVarItemsize(InputVarNames[i]) << "\n";
			oss << "    NBytes: " << brm.GetVarNbytes(InputVarNames[i]) << "\n";
		}
	}
	// OutputVarNames
	{
		std::vector<std::string> OutputVarNames = brm.GetOutputVarNames();
		int count = brm.GetOutputItemCount();
		assert(OutputVarNames.size() == (size_t)count);
		oss << "GetValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << OutputVarNames[i] << "\n";
			oss << "     Type:  " << brm.GetVarType(OutputVarNames[i]) << "\n";
			oss << "     Units: " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
			oss << "  Itemsize: " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
			oss << "    NBytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
		}
	}
	brm.OutputMessage(oss.str());
	// PointableVarNames
	{
		std::vector<std::string> PointableVarNames = brm.GetPointableVarNames();
		int count = brm.GetPointableItemCount();
		assert(PointableVarNames.size() == (size_t)count);
		oss << "PointableValues variables:\n";
		for (size_t i = 0; i < count; i++)
		{
			oss << "  " << i << "  " << PointableVarNames[i] << "\n";
			oss << "     Type:  " << brm.GetVarType(PointableVarNames[i]) << "\n";
			oss << "     Units: " << brm.GetVarUnits(PointableVarNames[i]) << "\n";
			oss << "  Itemsize: " << brm.GetVarItemsize(PointableVarNames[i]) << "\n";
			oss << "    NBytes: " << brm.GetVarNbytes(PointableVarNames[i]) << "\n";
		}
	}
	brm.OutputMessage(oss.str());
	// GetValue("Components")
	{
		int ncomps;
		brm.GetValue("ComponentCount", ncomps);
		std::vector<std::string> bmi_comps;
		brm.GetValue("Components", bmi_comps);
		std::vector<std::string> rm_comps;
		rm_comps = brm.GetComponents();
		assert(bmi_comps == rm_comps);
	}
	// GetValue("Concentrations")
	{
		int ncomps;
		brm.GetValue("ComponentCount", ncomps);
		double* concentrations_ptr = (double*)brm.GetValuePtr("Concentrations");
		int conc_nbytes = brm.GetVarNbytes("Concentrations");
		int conc_itemsize = brm.GetVarItemsize("Concentrations");
		int conc_dim = conc_nbytes / conc_itemsize;
		std::vector<double> bmi_conc;
		brm.GetValue("Concentrations", bmi_conc);
		assert(conc_dim == (int)bmi_conc.size());
		std::vector<double> rm_conc;
		brm.GetConcentrations(rm_conc);
		assert(rm_conc == bmi_conc);
		compare_ptrs(ptrs);
		for (int i = 0; i < conc_dim; i++)
		{
			bmi_conc[i] = bmi_conc[i] * 1.00001;
		}
		brm.SetConcentrations(bmi_conc);
	}
	compare_ptrs(ptrs);
	// GetValue("Density")
	// GetDensity and GetValue("Density) always return 
	// the calculated solution density
	{
		double* density_ptr = (double*)brm.GetValuePtr("Density");
		int dim = brm.GetVarNbytes("Density") / brm.GetVarItemsize("Density");
		std::vector<double> rm_density;
		brm.GetDensity(rm_density);
		std::vector<double> bmi_density;
		brm.GetValue("Density", bmi_density);
		assert(bmi_density == rm_density);
		for (int i = 0; i < dim; i++)
		{
			assert(density_ptr[i] == rm_density[i]);
			rm_density[i] *= 1.001;
		}
		brm.SetDensity(rm_density);
	}
	compare_ptrs(ptrs);
	// GetValue("ErrorString")
	{
		std::string bmi_err;
		brm.GetValue("ErrorString", bmi_err);
		std::string rm_err = brm.GetErrorString();
		assert(bmi_err == rm_err);
	}
	//FilePrefix
	{
		std::string str = "NewPrefix";
		brm.SetValue("FilePrefix", str);
		std::string rm_prefix = brm.GetFilePrefix();
		assert(str == rm_prefix);
		std::string prefix;
		brm.GetValue("FilePrefix", prefix);
		assert(prefix == rm_prefix);
	}
	// GetValue("Gfw")
	{
		double* gfw_ptr = (double*)brm.GetValuePtr("Gfw");
		std::vector<double> bmi_gfw;
		brm.GetValue("Gfw", bmi_gfw);
		std::vector<double> rm_gfw = brm.GetGfw();
		assert(bmi_gfw == rm_gfw);
		for (int i = 0; i < *ptrs.ComponentCount_ptr; i++)
		{
			assert(gfw_ptr[i] == rm_gfw[i]);
		}
	}
	// GetValue("GridCellCount")
	{
		int bmi_ngrid;
		brm.GetValue("GridCellCount", bmi_ngrid);
		int rm_ngrid = brm.GetGridCellCount();
		assert(bmi_ngrid == rm_ngrid);
		assert(bmi_ngrid == *ptrs.GridCellCount_ptr);
	}
	// SetValue("NthSelectedOutput") 
	{
		// tested in "SelectedOutput"
		// tested in "SelectedOutputHeadings"
	}
	// GetValue("Porosity")
	{
		int ngrid;
		brm.GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_porosity(ngrid, 0.3);
		brm.SetValue("Porosity", bmi_porosity);
		std::vector<double> rm_porosity;
		rm_porosity = brm.GetPorosity();
		assert(bmi_porosity == rm_porosity);
		brm.GetValue("Porosity", bmi_porosity);
		assert(bmi_porosity == rm_porosity);
	}
	compare_ptrs(ptrs);
	// testing of pointers and allocated values
	{
		double* por_ptr = (double*)brm.GetValuePtr("Porosity");
		double* por_ptr1 = (double*)brm.GetValuePtr("Porosity");
		assert(por_ptr == por_ptr1);
		double* my_porosity_alloc;
		my_porosity_alloc = (double*)malloc(*ptrs.GridCellCount_ptr * sizeof(double));
		double my_porosity_dim[20];
		std::vector<double> my_por;
		brm.GetValue("Porosity", my_por);
		brm.GetValue("Porosity", my_porosity_dim);
		brm.GetValue("Porosity", my_porosity_alloc);
		for (size_t i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(my_por[i] == my_porosity_dim[i]);
			assert(my_por[i] == my_porosity_alloc[i]);
			assert(my_por[i] == por_ptr[i]);
			assert(por_ptr[i] == por_ptr[i]);
		}
		my_por.resize(*ptrs.GridCellCount_ptr, 0.4);
		brm.SetValue("Porosity", my_por);
		brm.GetValue("Porosity", my_porosity_dim);
		brm.GetValue("Porosity", my_porosity_alloc);
		for (size_t i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(my_por[i] == my_porosity_dim[i]);
			assert(my_por[i] == my_porosity_alloc[i]);
			assert(my_por[i] == por_ptr[i]);
			assert(por_ptr[i] == por_ptr[i]);
		}
		my_por.resize(*ptrs.GridCellCount_ptr, 0.5);
		brm.SetPorosity(my_por);
		brm.GetValue("Porosity", my_porosity_dim);
		brm.GetValue("Porosity", my_porosity_alloc);
		for (size_t i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(my_por[i] == my_porosity_dim[i]);
			assert(my_por[i] == my_porosity_alloc[i]);
			assert(my_por[i] == por_ptr[i]);
			assert(por_ptr[i] == por_ptr[i]);
		}
		compare_ptrs(ptrs);
		free((void*)my_porosity_alloc);
	}
	// GetValue("Pressure")
	{
		int ngrid;
		brm.GetValue("GridCellCount", ngrid);
		int* nxyz_ptr = (int*)brm.GetValuePtr("GridCellCount");
		assert(*nxyz_ptr == ngrid);
		double* pressure_ptr = (double*)brm.GetValuePtr("Pressure");
		std::vector<double> bmi_pressure(ngrid, 1.1);
		brm.SetPressure(bmi_pressure);
		std::vector<double> rm_pressure;
		rm_pressure = brm.GetPressure();
		assert(bmi_pressure == rm_pressure);
		brm.GetValue("Pressure", bmi_pressure);
		assert(bmi_pressure == rm_pressure);
		for (size_t i = 0; i < *nxyz_ptr; i++)
		{
			assert(pressure_ptr[i] == rm_pressure[i]);
		}
		compare_ptrs(ptrs);
	}
	// GetValue("Saturation")
	// Always returns solution_volume / (rv * porosity) for each cell
	{
		double* saturation_ptr = (double*)brm.GetValuePtr("Saturation");
		std::vector<double> bmi_sat;
		brm.GetValue("Saturation", bmi_sat);
		std::vector<double> rm_sat;
		brm.GetSaturation(rm_sat);
		assert(bmi_sat == rm_sat);
		for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(saturation_ptr[i] == bmi_sat[i]);
		}
		rm_sat.resize(*ptrs.GridCellCount_ptr, 0.8);
		brm.SetValue("Saturation", rm_sat);
		brm.GetValue("Saturation", bmi_sat);
		assert(bmi_sat == rm_sat);
		for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(saturation_ptr[i] == bmi_sat[i]);
			bmi_sat[i] *= 1.01;
		}
		compare_ptrs(ptrs);
	}
	// GetValue("SolutionVolume")
	// Always returns solution_volume / (rv * porosity) for each cell
	{
		double* solution_volume_ptr = (double*)brm.GetValuePtr("SolutionVolume");
		std::vector<double> bmi_vol;
		brm.GetValue("SolutionVolume", bmi_vol);
		std::vector<double> rm_vol;
		rm_vol = brm.GetSolutionVolume();
		assert(bmi_vol == rm_vol);
		for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(solution_volume_ptr[i] == bmi_vol[i]);
		}
		compare_ptrs(ptrs);
	}
	// GetValue("Temperature")
	{
		double* temperature_ptr = (double*)brm.GetValuePtr("Temperature");
		int ngrid;
		brm.GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_temperature(ngrid, 50.);
		brm.SetValue("Temperature", bmi_temperature);
		std::vector<double> rm_temperature;
		rm_temperature = brm.GetTemperature();
		assert(bmi_temperature == rm_temperature);
		brm.GetValue("Temperature", bmi_temperature);
		assert(bmi_temperature == rm_temperature);
		for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(temperature_ptr[i] == bmi_temperature[i]);
		}
		rm_temperature.resize(*ptrs.GridCellCount_ptr, 22.0);
		compare_ptrs(ptrs);
	}
	// GetValue("Viscosity")
	{
		double* viscosity_ptr = (double*)brm.GetValuePtr("Viscosity");
		int ngrid;
		brm.GetValue("GridCellCount", ngrid);
		std::vector<double> bmi_viscosity;
		std::vector<double> rm_viscosity;
		brm.GetViscosity(rm_viscosity);
		brm.GetValue("Viscosity", bmi_viscosity);
		assert(bmi_viscosity == rm_viscosity);
		for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
		{
			assert(viscosity_ptr[i] == bmi_viscosity[i]);
		}
		compare_ptrs(ptrs);
	}
	// GetValue("SelectedOutput")
	{
		int bmi_so_count;
		brm.GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = brm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			brm.SetValue("NthSelectedOutput", i);
			int bmi_nuser;
			brm.GetValue("CurrentSelectedOutputUserNumber", bmi_nuser);
			int rm_nuser = brm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_col_count;
			brm.GetValue("SelectedOutputColumnCount", bmi_col_count);
			int rm_col_count = brm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			int bmi_row_count;
			brm.GetValue("SelectedOutputRowCount", bmi_row_count);
			int rm_row_count = brm.GetSelectedOutputRowCount();
			assert(bmi_row_count == rm_row_count);

			std::vector<double> bmi_so;
			brm.GetValue("SelectedOutput", bmi_so);
			std::vector<double> rm_so;
			brm.GetSelectedOutput(rm_so);
			assert(bmi_so == rm_so);
		}
	}
	// GetValue("SelectedOutputColumnCount")
	{
		int bmi_so_col_count;
		brm.GetValue("SelectedOutputColumnCount", bmi_so_col_count);
		int rm_so_col_count = brm.GetSelectedOutputColumnCount();
		assert(bmi_so_col_count == rm_so_col_count);
	}
	// GetValue("SelectedOutputCount")
	{
		int bmi_so_count;
		brm.GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = brm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
	}
	// GetValue("SelectedOutputHeadings")
	{
		int bmi_so_count;
		brm.GetValue("SelectedOutputCount", bmi_so_count);
		int rm_so_count = brm.GetSelectedOutputCount();
		assert(bmi_so_count == rm_so_count);
		for (int i = 0; i < bmi_so_count; i++)
		{
			brm.SetValue("NthSelectedOutput", i);
			int bmi_nuser;
			brm.GetValue("CurrentSelectedOutputUserNumber", bmi_nuser);
			int rm_nuser = brm.GetNthSelectedOutputUserNumber(i);
			assert(bmi_nuser == rm_nuser);

			int bmi_col_count;
			brm.GetValue("SelectedOutputColumnCount", bmi_col_count);
			int rm_col_count = brm.GetSelectedOutputColumnCount();
			assert(bmi_col_count == rm_col_count);

			std::vector< std::string > bmi_headings;
			brm.GetValue("SelectedOutputHeadings", bmi_headings);
			std::vector< std::string > rm_headings;
			IRM_RESULT status = brm.GetSelectedOutputHeadings(rm_headings);
			assert(bmi_headings == rm_headings);
		}
	}
	// GetValue("SelectedOutputOn")
	{
		bool* selectedoutput_on_ptr = (bool*)brm.GetValuePtr("SelectedOutputOn");
		bool bmi_so_on;
		brm.GetValue("SelectedOutputOn", bmi_so_on);
		int rm_so_row_count = brm.GetSelectedOutputRowCount();
		assert(bmi_so_on == *selectedoutput_on_ptr);
		bmi_so_on = false;
		brm.SetValue("SelectedOutputOn", bmi_so_on);
		assert(*selectedoutput_on_ptr == bmi_so_on);
		bmi_so_on = true;
		brm.SetValue("SelectedOutputOn", bmi_so_on);
		assert(*selectedoutput_on_ptr == bmi_so_on);
		brm.SetSelectedOutputOn(false);
		compare_ptrs(ptrs);
	}
	// GetValue("SelectedOutputRowCount")
	{
		int bmi_so_row_count;
		brm.GetValue("SelectedOutputRowCount", bmi_so_row_count);
		int rm_so_row_count = brm.GetSelectedOutputRowCount();
		assert(bmi_so_row_count == rm_so_row_count);
	}
	// Time
	{
		double time = 1.0;
		double* time_ptr = (double*)brm.GetValuePtr("Time");
		brm.SetTime(time);
		compare_ptrs(ptrs);
		double time1;
		brm.GetValue("Time", time1);
		assert(time1 == time);
		assert(*time_ptr = time);
		time = 2.0;
		brm.SetValue("Time", time);
		compare_ptrs(ptrs);

	}
	// TimeStep
	{
		double timestep = 1.1;
		double* timestep_ptr = (double*)brm.GetValuePtr("TimeStep");
		brm.SetValue("TimeStep", timestep);
		compare_ptrs(ptrs);
		assert(*timestep_ptr == timestep);
		timestep = 1.2;
		brm.SetTimeStep(timestep);
		compare_ptrs(ptrs);
		assert(*timestep_ptr == timestep);
		double timestep1;
		brm.GetValue("TimeStep", timestep1);
		assert(timestep1 == timestep);
		assert(*timestep_ptr = timestep);
	}
}
void compare_ptrs(struct Ptrs ptrs)
{
	assert(*ptrs.ComponentCount_ptr == ptrs.brm->GetComponentCount());
	std::vector<double> gfw = ptrs.brm->GetGfw();
	for (int i = 0; i < *ptrs.ComponentCount_ptr; i++)
	{
		assert(ptrs.Gfw_ptr[i] == gfw[i]);
	}
	assert(*ptrs.GridCellCount_ptr == ptrs.brm->GetGridCellCount());
	std::vector<double> density, saturation, SolutionVolume, 
		Porosity, Pressure, Temperature, Viscosity;
	ptrs.brm->GetDensity(density);
	ptrs.brm->GetSaturation(saturation);
	SolutionVolume = ptrs.brm->GetSolutionVolume();
	Porosity = ptrs.brm->GetPorosity();
	Pressure = ptrs.brm->GetPressure();
	Temperature = ptrs.brm->GetTemperature();
	ptrs.brm->GetViscosity(Viscosity);

	for (int i = 0; i < *ptrs.GridCellCount_ptr; i++)
	{
		assert(ptrs.Density_ptr[i] == density[i]);
		assert(ptrs.Saturation_ptr[i] == saturation[i]);
		assert(ptrs.SolutionVolume_ptr[i] == SolutionVolume[i]);
		assert(ptrs.Porosity_ptr[i] == Porosity[i]);
		assert(ptrs.Pressure_ptr[i] == Pressure[i]);
		assert(ptrs.Temperature_ptr[i] == Temperature[i]);
		assert(ptrs.Viscosity_ptr[i] == Viscosity[i]);
	}
	std::vector<double> Concentrations;
	ptrs.brm->GetConcentrations(Concentrations);
	int dim = *ptrs.ComponentCount_ptr * *ptrs.GridCellCount_ptr;
	for (int i = 0; i < dim; i++)
	{
		assert(ptrs.Concentrations_ptr[i] == Concentrations[i]);
	}
	assert(*ptrs.Time_ptr == ptrs.brm->GetTime());
	assert(*ptrs.TimeStep_ptr == ptrs.brm->GetTimeStep());
	assert(*ptrs.SelectedOutputOn_ptr == ptrs.brm->GetSelectedOutputOn());
}
void testbmi_cpp()
{
	YAMLPhreeqcRM yrm;
	IRM_RESULT status;
	int nxyz = 40;
	// Set GridCellCount
	yrm.YAMLSetGridCellCount(nxyz);
	yrm.YAMLThreadCount(3);
	std::string YAML_filename = "testbmi_cpp.yaml";
	yrm.WriteYAMLDoc(YAML_filename);
	//
	// Use all BMIPhreeqcRM methods roughly in order of use
	// 
	//
	BMIPhreeqcRM bmi;
	std::cerr << "BMIPhreeqcRM\n";
	//-------
	bmi.Initialize(YAML_filename);   // void function
	std::cerr << "Initialize\n";
	//-------
	bmi.GetValue("GridCellCount", nxyz);
	std::cerr << "GetValue('GridCellCount')\n";
	//-------
	nxyz = bmi.GetGridCellCount();
	std::cerr << "GetGridCellCount \n";
	//-------
	int n = bmi.GetThreadCount();
	std::cerr << "GetThreadCount " << n << "\n";
	//------- LoadDatabase removes definitions in instance	//-------
	// Inactive cells or symmetry
	std::vector<int> grid2chem(nxyz, -1);
	for (size_t i = 0; i < nxyz / 2; i++)
	{
		grid2chem[i] = i;
	}
	status = bmi.CreateMapping(grid2chem);
	std::cerr << "CreateMapping \n";
	bmi.LoadDatabase("phreeqc.dat");
	std::cerr << "LoadDatabase\n";
	//
	// Set properties
	// 
	//-------
	status = bmi.SetComponentH2O(false);
	std::cerr << "SetComponentH2O \n";
	//-------
	bmi.SetSpeciesSaveOn(true);
	std::cerr << "SetSpeciesSaveOn \n";
	//-------
	status = bmi.SetErrorOn(true);
	std::cerr << "SetErrorOn \n";
	//-------
	status = bmi.SetErrorHandlerMode(1);
	std::cerr << "SetErrorHandlerMode \n";
	//-------
	status = bmi.SetDumpFileName("AdvectBMI_test_py.dump");
	std::cerr << "SetDumpFileName \n";
	//-------
	status = bmi.SetFilePrefix("AdvectBMI_test_py");
	std::cerr << "SetFilePrefix \n";
	//-------
	status = bmi.OpenFiles();
	std::cerr << "OpenFiles \n";
	//-------
	status = bmi.SetPartitionUZSolids(false);
	std::cerr << "SetPartitionUZSolids \n";
	//-------
	status = bmi.SetRebalanceByCell(true);
	std::cerr << "SetRebalanceByCell \n";
	//-------
	status = bmi.SetRebalanceFraction(0.5);
	std::cerr << "SetRebalanceFraction \n";
	//-------
	status = bmi.SetScreenOn(true);
	std::cerr << "SetScreenOn \n";
	//-------
	status = bmi.SetSelectedOutputOn(true);
	std::cerr << "SetSelectedOutputOn \n";
	//-------
	std::cerr << "SetSelectedOutputOn \n";	status = bmi.SetUnitsExchange(1);
	std::cerr << "SetUnitsExchange \n";
	//-------
	status = bmi.SetUnitsGasPhase(1);
	std::cerr << "SetUnitsGasPhase \n";
	//-------
	status = bmi.SetUnitsKinetics(1);
	std::cerr << "SetUnitsKinetics \n";
	//-------
	status = bmi.SetUnitsPPassemblage(1);
	std::cerr << "SetUnitsPPassemblage \n";
	//-------
	status = bmi.SetUnitsSolution(2);
	std::cerr << "SetUnitsSolution \n";
	//-------
	status = bmi.SetUnitsSSassemblage(1);
	std::cerr << "SetUnitsSSassemblage \n";
	//-------
	status = bmi.SetUnitsSurface(1);
	std::cerr << "SetUnitsSurface \n";
	//-------
	bmi.UseSolutionDensityVolume(false);
	std::cerr << "UseSolutionDensityVolume \n";
	//-------
	double time_conversion = 1.0 / 86400.0;
	bmi.SetTimeConversion(time_conversion);
	std::cerr << "SetTimeConversion \n";
	//-------
	std::vector<double> v(nxyz, 1.0);
	status = bmi.SetRepresentativeVolume(v);
	std::cerr << "SetRepresentativeVolume \n";

	//-------Chemistry cells may be fewer than GridCellCount
	int nchem = bmi.GetChemistryCellCount();
	std::cerr << "GetChemistryCellCount \n";
	//
	// Begin initialization
	//
	status = bmi.SetPrintChemistryOn(false, true, false);
	std::cerr << "SetPrintChemistryOn \n";
	//-------
	status = bmi.RunFile(true, true, true, "advectBMI_test.pqi");
	std::cerr << "RunFile \n";
	//-------
	//bmi.AddOutputVars("SolutionActivities", "true");
	//bmi.AddOutputVars("SolutionMolalities", "true");
	//bmi.AddOutputVars("SaturationIndices", "true");
	bmi.AddOutputVars("SolutionActivities", "H+ Ca+2 Na+");
	bmi.AddOutputVars("SolutionMolalities", "OH- Cl-");
	bmi.AddOutputVars("SaturationIndices", "Calcite Dolomite");
	std::cerr << "AddOutputVars \n";
	//-------
	int ncomps = bmi.FindComponents();
	std::cerr << "FindComponents \n";
	// Component names
	ncomps = bmi.GetComponentCount();
	std::cerr << "GetComponentCount \n";
	//-------
	bmi.GetValue("ComponentCount", ncomps);
	std::cerr << "GetValue('ComponentCount')" << ncomps << "\n";
	//-------
	std::vector<std::string> str_vector = bmi.GetComponents();
	std::cerr << "GetComponents \n";
	// Species info
	n = bmi.GetSpeciesCount();
	std::cerr << "GetSpeciesCount \n";
	//-------
	str_vector = bmi.GetSpeciesNames();
	std::cerr << "GetSpeciesNames \n";
	//-------
	v = bmi.GetSpeciesD25();
	std::cerr << "GetSpeciesD25 \n";
	//-------
	v = bmi.GetSpeciesZ();
	std::cerr << "GetSpeciesZ \n";
	// Reactant lists
	std::vector<std::string> equiuilibrium_phases = bmi.GetEquilibriumPhases();
	std::cerr << "GetEquilibriumPhases \n";
	//-------
	n = bmi.GetEquilibriumPhasesCount();
	std::cerr << "GetEquilibriumPhasesCount \n";
	//-------
	str_vector = bmi.GetExchangeNames();
	std::cerr << "GetExchangeNames \n";
	//-------
	str_vector = bmi.GetExchangeSpecies();
	std::cerr << "GetExchangeSpecies \n";
	//-------
	n = bmi.GetExchangeSpeciesCount();
	std::cerr << "GetExchangeSpeciesCount \n";
	//-------
	str_vector = bmi.GetGasComponents();
	std::cerr << "GetGasComponents \n";
	//-------
	n = bmi.GetGasComponentsCount();
	std::cerr << "GetGasComponentsCount \n";
	//-------
	str_vector = bmi.GetKineticReactions();
	std::cerr << "GetKineticReactions \n";
	//-------
	n = bmi.GetKineticReactionsCount();
	std::cerr << "GetKineticReactionsCount \n";
	//-------
	n = bmi.GetSICount();
	std::cerr << "GetSICount \n";
	//-------
	str_vector = bmi.GetSINames();
	std::cerr << "GetSINames \n";
	//-------
	str_vector = bmi.GetSolidSolutionComponents();
	std::cerr << "GetSolidSolutionComponents \n";
	//-------
	n = bmi.GetSolidSolutionComponentsCount();
	std::cerr << "GetSolidSolutionComponentsCount \n";
	//-------
	str_vector = bmi.GetSolidSolutionNames();
	std::cerr << "GetSolidSolutionNames \n";
	//-------
	str_vector = bmi.GetSurfaceNames();
	std::cerr << "GetSurfaceNames \n";
	//-------
	str_vector = bmi.GetSurfaceSpecies();
	std::cerr << "GetSurfaceSpecies \n";
	//-------
	n = bmi.GetSurfaceSpeciesCount();
	std::cerr << "GetSurfaceSpeciesCount \n";
	//-------
	str_vector = bmi.GetSurfaceTypes();
	std::cerr << "GetSurfaceTypes \n";
	//
	// Remove any reactants in workers 
	//
	std::string input = "DELETE; -all";
	bmi.RunString(true, false, false, input);
	std::cerr << "RunString \n";
	//-------
	// 
	// Transfer initial conditions
	//
	std::vector<int> vi(nxyz, 1);
	status = bmi.InitialEquilibriumPhases2Module(vi);
	std::cerr << "InitialEquilibriumPhases2Module \n";
	//-------
	status = bmi.InitialExchanges2Module(vi);
	std::cerr << "InitialExchanges2Module \n";
	//-------
	status = bmi.InitialGasPhases2Module(vi);
	std::cerr << "InitialGasPhases2Module \n";
	//-------
	status = bmi.InitialKinetics2Module(vi);
	std::cerr << "InitialKinetics2Module \n";
	//-------
	status = bmi.InitialSolutions2Module(vi);
	std::cerr << "InitialSolutions2Module \n";
	//-------
	status = bmi.InitialSolidSolutions2Module(vi);
	std::cerr << "InitialSolidSolutions2Module \n";
	//-------
	status = bmi.InitialSurfaces2Module(vi);
	std::cerr << "InitialSurfaces2Module \n";
	//-------
	// Alternative A.to the previous seven methods
	std::vector<int> ic(nxyz * 7, 1);
	status = bmi.InitialPhreeqc2Module(ic);
	std::cerr << "InitialPhreeqc2Module \n";
	//-------
	// Alternative B.to the previous seven methods, possible mixing
	std::vector<int> v1(nxyz * 7, 1);
	std::vector<int> v2(nxyz * 7, -1);
	std::vector<double> f1(nxyz * 7, 1.0);
	status = bmi.InitialPhreeqc2Module(v1, v2, f1);
	std::cerr << "InitialPhreeqc2Modul mix \n";
	//-------
	// Alternative C.to the previous seven methods, initialize cells 18 and 19
	std::vector<int> cells(2);
	cells[0] = 18;
	cells[1] = 19;
	status = bmi.InitialPhreeqcCell2Module(1, cells);
	std::cerr << "InitialPhreeqcCell2Module \n";
	//
	// Boundary conditions
	// 
	std::vector<double> bc;
	vi.resize(1, 1);
	status = bmi.InitialPhreeqc2SpeciesConcentrations(bc, vi);
	std::cerr << "InitialPhreeqc2SpeciesConcentrations \n";
	//-------
	std::vector<double> bc_species;
	std::vector<int> vi1(1, 1);
	std::vector<int> vi2(1, -1);
	f1.resize(1, 1.0);
	status = bmi.InitialPhreeqc2SpeciesConcentrations(bc_species, vi1, vi2, f1);
	std::cerr << "InitialPhreeqc2SpeciesConcentrations mix \n";
	//
	// Get/Set methods for time steping
	//
	status = bmi.SetTime(0.0);
	std::cerr << "SetTime \n";
	//-------
	status = bmi.SetTimeStep(0.0);
	std::cerr << "SetTimeStep \n";
	//-------
	status = bmi.GetGasCompMoles(v);
	std::cerr << "GetGasCompMoles \n";
	//-------
	status = bmi.SetGasCompMoles(v);
	std::cerr << "SetGasCompMoles \n";
	//-------
	vi.resize(nxyz, 1);
	status = bmi.SetPrintChemistryMask(vi);
	std::cerr << "SetPrintChemistryMask \n";
	//
	// Get/Set methods for time stepping
	// 
	std::vector<double> c;
	status = bmi.GetConcentrations(c);
	std::cerr << "GetConcentrations \n";
	//-------
	status = bmi.SetConcentrations(c);
	std::cerr << "SetConcentrations \n";
	//-------
	v.resize(nxyz, 0);
	status = bmi.GetDensity(v);
	std::cerr << "GetDensity \n";
	//-------
	v.resize(nxyz, 1.1);
	status = bmi.SetDensity(v);
	std::cerr << "SetDensity \n";
	//-------
	status = bmi.GetGasCompMoles(v);
	std::cerr << "GetGasCompMoles \n";
	//-------
	status = bmi.SetGasCompMoles(v);
	std::cerr << "GetGasCompMoles \n";
	//-------
	status = bmi.GetGasCompPressures(v);
	std::cerr << "GetGasCompPressures \n";
	//-------
	status = bmi.GetGasPhaseVolume(v);
	std::cerr << "GetGasPhaseVolume \n";
	//-------
	v.resize(nxyz, 1.0);
	status = bmi.SetGasPhaseVolume(v);
	std::cerr << "SetGasPhaseVolume \n";
	//-------
	v.resize(nxyz, 0.21);
	status = bmi.SetPorosity(v);
	std::cerr << "SetPorosity \n";
	//-------
	v = bmi.GetPressure();
	std::cerr << "GetPressure \n";
	//-------
	v.resize(nxyz, 3.0);
	status = bmi.SetPressure(v);
	std::cerr << "SetPressure \n";
	//-------
	v.resize(nxyz, 1.0);
	status = bmi.SetSaturation(v);
	std::cerr << "SetSaturation \n";
	//-------
	v = bmi.GetSolutionVolume();
	std::cerr << "GetSolutionVolume \n";
	//-------
	status = bmi.GetSpeciesConcentrations(v);
	std::cerr << "GetSpeciesConcentrations \n";
	//-------
	status = bmi.SpeciesConcentrations2Module(v);
	std::cerr << "SpeciesConcentrations2Module \n";
	//-------
	status = bmi.GetSpeciesLog10Gammas(v);
	std::cerr << "GetSpeciesLog10Gammas \n";
	//-------
	status = bmi.GetSpeciesLog10Molalities(v);
	std::cerr << "GetSpeciesLog10Molalities \n";
	//-------
	v.resize(nxyz, 26.0);
	status = bmi.SetTemperature(v);
	std::cerr << "SetTemperature \n";
	//
	// Take a time step
	//
	bmi.Update();      // void function
	std::cerr << "Update\n";
	//-------
	status = bmi.RunCells();
	std::cerr << "RunCells\n";
	//
	// Selected output
	//
	status = bmi.SetNthSelectedOutput(0);
	std::cerr << "SetNthSelectedOutput \n";
	//-------
	int n_user = bmi.GetCurrentSelectedOutputUserNumber();
	std::cerr << "GetCurrentSelectedOutputUserNumber \n";
	//-------
	status = bmi.SetCurrentSelectedOutputUserNumber(333);
	std::cerr << "SetCurrentSelectedOutputUserNumber \n";
	//-------
	n = bmi.GetNthSelectedOutputUserNumber(0);
	std::cerr << "GetNthSelectedOutputUserNumber \n";
	//-------
	status = bmi.GetSelectedOutput(v);
	std::cerr << "GetSelectedOutput \n";
	//-------
	n = bmi.GetSelectedOutputColumnCount();
	std::cerr << "GetSelectedOutputColumnCount \n";
	//-------
	n = bmi.GetSelectedOutputRowCount();
	std::cerr << "GetSelectedOutputRowCount \n";
	//
	// Getters
	// 
	std::vector< std::vector<int> > back_map = bmi.GetBackwardMapping();
	std::cerr << "GetBackwardMapping \n";
	//-------
	std::string db_name = bmi.GetDatabaseFileName();
	std::cerr << "GetDatabaseFileName \n";
	//-------
	vi = bmi.GetEndCell();
	std::cerr << "GetEndCell\n";
	//-------
	n = bmi.GetErrorHandlerMode();
	std::cerr << "GetErrorHandlerMode \n";
	//-------
	std::string str = bmi.GetErrorString();
	std::cerr << "GetErrorString \n";
	//-------
	str = bmi.GetFilePrefix();
	std::cerr << "GetFilePrefix \n";
	//-------
	vi = bmi.GetForwardMapping();
	std::cerr << "GetForwardMapping \n";
	//-------
	status = bmi.GetGasCompPhi(v);
	std::cerr << "GetGasCompPhi \n";
	//-------
	v = bmi.GetGfw();
	std::cerr << "GetGfw ";
	//-------
	IPhreeqc* ipq = bmi.GetIPhreeqcPointer(0);
	std::cerr << "GetIPhreeqcPointer \n";
	//-------
	n = bmi.GetMpiMyself();
	std::cerr << "GetMpiMyself \n";
	//-------
	n = bmi.GetMpiTasks();
	std::cerr << "GetMpiTasks \n";
	//-------
	bool b = bmi.GetPartitionUZSolids();
	std::cerr << "GetPartitionUZSolids \n";
	//-------
	v = bmi.GetPorosity();
	std::cerr << "GetPorosity \n";
	//-------
	vi = bmi.GetPrintChemistryMask();
	std::cerr << "GetPrintChemistryMask \n";
	//-------
	std::vector<bool> vb = bmi.GetPrintChemistryOn();
	std::cerr << "GetPrintChemistryOn \n";
	//-------
	b = bmi.GetRebalanceByCell();
	std::cerr << "GetRebalanceByCell \n";
	//-------
	double d = bmi.GetRebalanceFraction();
	std::cerr << "GetRebalanceFraction \n";
	//-------
	status = bmi.GetSaturation(v);
	std::cerr << "GetSaturation \n";
	//-------
	status = bmi.GetSelectedOutputHeadings(str_vector);
	std::cerr << "GetSelectedOutputHeadings \n";
	//-------
	b = bmi.GetSelectedOutputOn();
	std::cerr << "GetSelectedOutputOn \n";
	//-------
	b = bmi.GetSpeciesSaveOn();
	std::cerr << "GetSpeciesSaveOn \n";
	//-------
	std::vector< cxxNameDouble > s = bmi.GetSpeciesStoichiometry();
	std::cerr << "GetSpeciesStoichiometry \n";
	//-------
	vi = bmi.GetStartCell();
	std::cerr << "GetStartCell \n";
	//-------
	v = bmi.GetTemperature();
	std::cerr << "GetTemperature \n";
	//-------
	d = bmi.GetTime();
	std::cerr << "GetTime \n";
	//-------
	d = bmi.GetTimeConversion();
	std::cerr << "GetTimeConversion \n";
	//-------
	d = bmi.GetTimeStep();
	std::cerr << "GetTimeStep \n";
	//-------
	n = bmi.GetUnitsExchange();
	std::cerr << "GetUnitsExchange \n";
	//-------
	n = bmi.GetUnitsGasPhase();
	std::cerr << "GetUnitsGasPhase \n";
	//-------
	n = bmi.GetUnitsKinetics();
	std::cerr << "GetUnitsKinetics \n";
	//-------
	n = bmi.GetUnitsPPassemblage();
	std::cerr << "GetUnitsPPassemblage \n";
	//-------
	n = bmi.GetUnitsSolution();
	std::cerr << "GetUnitsSolution \n";
	n = bmi.GetUnitsSSassemblage();
	std::cerr << "GetUnitsSSassemblage \n";
	//-------
	n = bmi.GetUnitsSurface();
	std::cerr << "GetUnitsSurface \n";
	//-------
	std::vector<IPhreeqcPhast *> w = bmi.GetWorkers();
	std::cerr << "GetWorkers \n";
	//
	// Utilities
	//
	vi.resize(1, 1);
	status = bmi.InitialPhreeqc2Concentrations(bc, vi);
	std::vector<double> tc(1, 30.0);
	std::vector<double> p_atm(1, 1.5);
	IPhreeqc* utility_ptr = bmi.Concentrations2Utility(bc, tc, p_atm);
	std::cerr << "Concentrations2Utility \n";
	//-------
	bmi.DecodeError(-2);	         // void function
	std::cerr << "DecodeError \n";
	//-------
	status = bmi.DumpModule(true);
	std::cerr << "DumpModule \n";
	//-------
	bmi.ErrorHandler(0, "string"); // void function
	std::cerr << "OK, just a test: ErrorHandler \n";
	//-------
	bmi.ErrorMessage("my error");  // void function
	std::cerr << "OK, just a test: ErrorMessage \n";
	//-------
	bmi.LogMessage("Log message");  // void method
	std::cerr << "LogMessage \n";
	//-------
	bmi.OutputMessage("Output message");  // void method
	std::cerr << "OutputMessage \n";
	//-------
	bmi.ScreenMessage("Screen message\n");  // void method
	std::cerr << "ScreenMessage \n";
	//-------
	status = bmi.StateSave(1);
	std::cerr << "StateSave \n";
	//-------
	status = bmi.StateApply(1);
	std::cerr << "StateApply \n";
	//-------
	status = bmi.StateDelete(1);
	std::cerr << "StateDelete \n";
	//-------
	bmi.WarningMessage("Warning message");  // void method
	std::cerr << "WarningMessage \n";
	//-------
	// BMI Methods
	str = bmi.GetComponentName();
	std::cerr << "GetComponentName \n";
	//-------
	d = bmi.GetCurrentTime();
	std::cerr << "GetCurrentTime \n";
	//-------
	d = bmi.GetEndTime();
	std::cerr << "GetEndTime \n";
	//-------
	n = bmi.GetInputItemCount();
	std::cerr << "GetInputItemCount \n";
	//-------
	str_vector = bmi.GetInputVarNames();
	std::cerr << "GetInputVarNames \n";
	//-------
	n = bmi.GetOutputItemCount();
	std::cerr << "GetOutputItemCount \n";
	//-------
	str_vector = bmi.GetOutputVarNames();
	std::cerr << "GetOutputVarNames \n";
	//-------
	d = bmi.GetTimeStep();
	std::cerr << "GetTimeStep \n";
	//-------
	str = bmi.GetTimeUnits();
	std::cerr << "GetTimeUnits \n";
	//-------
	bmi.GetValue("solution_saturation_index_Calcite", v);
	std::cerr << "GetValue \n";
	//-------
	n = bmi.GetVarItemsize("solution_saturation_index_Calcite");
	std::cerr << "GetVarItemsize \n";
	//-------
	n = bmi.GetVarNbytes("solution_saturation_index_Calcite");
	std::cerr << "GetVarNbytes \n";
	//-------
	str = bmi.GetVarType("solution_saturation_index_Calcite");
	std::cerr << "GetVarType \n";
	//-------
	str = bmi.GetVarUnits("solution_saturation_index_Calcite");
	std::cerr << "GetVarUnits \n";
	//status = bmi.Initialize("AdvectBMI_test_py.yaml");
	// See above
	bmi.SetValue("Time", 1.0);    // void method
	std::cerr << "SetValue";
	//-------
	bmi.Update();    // void method
	std::cerr << "Update";
	//-------	
	status = bmi.CloseFiles(); // not a BMI method, but needs to be last
	std::cerr << "CloseFiles \n";
	bmi.Finalize();    // void method
	std::cerr << "Finalize \n";
	//Should be private: status = bmi.ReturnHandler();
	//TODO status = bmi.MpiAbort();
	//TODO status = bmi.MpiWorker();
	//TODO status = bmi.MpiWorkerBreak();
	//TODO status = bmi.SetMpiWorkerCallbackC();
	//TODO status = bmi.SetMpiWorkerCallbackCookie();
	std::cerr << "Success.\n";
	return;
}
#endif // YAML