#ifdef USE_YAML
#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <stdlib.h>
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
		DensityCalculated_ptr = NULL;
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
	double* DensityCalculated_ptr;
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


void advectionbmi_cpp(std::vector<double>& c, std::vector<double> bc_conc, int ncomps, int nxyz, int dim);

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
		//int nxyz = 40;
		int nxyz = 0;
#ifdef USE_MPI
		// MPI
		int mpi_myself;
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
		{
			exit(4);
		}
		if (mpi_myself == 0)
		{
			nxyz = BMIPhreeqcRM::GetGridCellCountYAML(yaml_file.c_str());
		}
		MPI_Bcast(&nxyz, 1, MPI_INT, 0, MPI_COMM_WORLD);
		BMIPhreeqcRM brm(nxyz, MPI_COMM_WORLD);
		if (mpi_myself > 0)
		{
			brm.MpiWorker();
			brm.Finalize();
			return EXIT_SUCCESS;
		};
#else
		// OpenMP
		nxyz = BMIPhreeqcRM::GetGridCellCountYAML(yaml_file.c_str());
		BMIPhreeqcRM brm;
#endif
		// Use YAML file to initialize
		brm.Initialize(yaml_file);
		// InputVarNames
		{
			std::ostringstream oss;
			std::vector<std::string> InputVarNames = brm.GetInputVarNames();
			int count = brm.GetInputItemCount();
			assert(InputVarNames.size() == (size_t)count);
			oss << "SetValues variables and units:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << std::setw(50) << InputVarNames[i] << "   "
					<< brm.GetVarUnits(InputVarNames[i]) << "\n";
			}
			brm.OutputMessage(oss.str());
		}
		// OutputVarNames
		{
			std::ostringstream oss;
			std::vector<std::string> OutputVarNames = brm.GetOutputVarNames();
			int count = brm.GetOutputItemCount();
			assert(OutputVarNames.size() == (size_t)count);
			oss << "GetValues variables and units:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << std::setw(50) << OutputVarNames[i] << "   "
					<< brm.GetVarUnits(OutputVarNames[i]) << "\n";
			}
			brm.OutputMessage(oss.str());
		}
		// set pointers demonstration
		Ptrs ptrs;
		ptrs.brm = &brm;
		ptrs.ComponentCount_ptr = (int*)brm.GetValuePtr("ComponentCount");
		ptrs.Concentrations_ptr = (double*)brm.GetValuePtr("Concentrations");
		ptrs.DensityCalculated_ptr = (double*)brm.GetValuePtr("DensityCalculated");
		ptrs.Gfw_ptr = (double*)brm.GetValuePtr("Gfw");
		ptrs.GridCellCount_ptr = (int*)brm.GetValuePtr("GridCellCount");
		ptrs.Saturation_ptr = (double*)brm.GetValuePtr("SaturationCalculated");
		ptrs.SolutionVolume_ptr = (double*)brm.GetValuePtr("SolutionVolume");
		ptrs.Time_ptr = (double*)brm.GetValuePtr("Time");
		ptrs.TimeStep_ptr = (double*)brm.GetValuePtr("TimeStep");
		ptrs.Porosity_ptr = (double*)brm.GetValuePtr("Porosity");
		ptrs.Pressure_ptr = (double*)brm.GetValuePtr("Pressure");
		ptrs.SelectedOutputOn_ptr = (bool*)brm.GetValuePtr("SelectedOutputOn");
		ptrs.Temperature_ptr = (double*)brm.GetValuePtr("Temperature");
		ptrs.Viscosity_ptr = (double*)brm.GetValuePtr("Viscosity");

		// Get number of components
		int ncomps;
		brm.GetValue("ComponentCount", ncomps);
		// Print some of the reaction module information 
		{
			std::ostringstream oss;
			std::string prefix;
			brm.GetValue("FilePrefix", prefix);
			oss << "File prefix:                                      " << prefix << "\n";
			int ngrid;
			brm.GetValue("GridCellCount", ngrid);
			oss << "Number of grid cells in the user's model:         " << ngrid << "\n";
			int ncomps;
			brm.GetValue("ComponentCount", ncomps);
			oss << "Number of components for transport:               " << ncomps << "\n";
			brm.OutputMessage(oss.str());
		}

		// Get components
		std::vector<std::string> components;
		brm.GetValue("Components", components);

		// Get gfw
		std::vector<double> gfw(ncomps, 0);
		brm.GetValue("Gfw", gfw);
		// write component information using PhreeqcRM OutputMessage
		for (int i = 0; i < ncomps; i++)
		{
			std::ostringstream strm;
			strm.width(10);
			strm << components[i] << "    " << gfw[i] << "\n";
			brm.OutputMessage(strm.str());
		}
		brm.OutputMessage("\n");

		// Get temperatures (note, unchanged by BMIPhreeqcRM)
		std::vector<double> tempc;
		brm.GetValue("Temperature", tempc);

		// Get initial saturation
		IRM_RESULT status;
		std::vector<double> sat;
		brm.GetValue("SaturationCalculated", sat);

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
		brm.SetValue("DensityUser", density);
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
			// Transport calculation here
			advectionbmi_cpp(c, bc_conc, ncomps, nxyz, nbound);
			// Transfer data to PhreeqcRM for reactions
			bool print_selected_output_on = (steps == nsteps - 1) ? true : false;
			brm.SetValue("SelectedOutputOn", print_selected_output_on);
			bool print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = brm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			brm.SetValue("Concentrations", c);        // Transported concentrations
			brm.SetValue("Porosity", por);            // If pororosity changes due to compressibility
			brm.SetValue("SaturationUser", sat);          // If saturation changes
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

			//Run chemistry
			brm.Update();
			// Transfer data from PhreeqcRM for transport
			brm.GetValue("Concentrations", c);
			brm.GetValue("DensityCalculated", density);
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
					if (n_user == 333) continue; // Skip internal SELECTED_OUTPUT
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
						oss << "     Calculated Density: " << density[i] << "\n";
						oss << "     Calculated Volume:  " << volume[i] << "\n";
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
						// YAMLAddOutputVars can be used to select groups of
						// variables that are available through GetValue
						std::vector<double> CaX2, KX, NaX, pH;
						brm.GetValue("solution_ph", pH);
						brm.GetValue("exchange_X_species_log_molality_CaX2", CaX2);
						brm.GetValue("exchange_X_species_log_molality_KX", KX);
						brm.GetValue("exchange_X_species_log_molality_NaX", NaX);
						std::ostringstream xoss;
						xoss << "      pH      CaX2     KX       NaX\n";
						for (size_t i = 0; i < nxyz / 2; i++)
						{
							xoss << std::setw(10) << std::setprecision(5) << std::fixed <<
								pH[i] << "  " << std::pow(10.0, CaX2[i]) << "  " << 
								std::pow(10.0, KX[i]) << "  " << std::pow(10.0, NaX[i]) << "\n";
						}
						brm.OutputMessage(xoss.str()); 
					}
				}
			}
		}
		// Clean up
		status = brm.MpiWorkerBreak();
		brm.Finalize();
	}
	catch (PhreeqcRMStop)
	{
		std::string e_string = "AdvectBMI_cpp failed with an error in PhreeqcRM.";
		std::cerr << e_string << std::endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return IRM_FAIL;
	}
#if 0
	catch (...)
	{
		std::string e_string = "AdvectBMI_cpp failed with an unhandled exception.";
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

#endif // YAML