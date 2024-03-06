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

void TestAllMethods_cpp()
{
	YAMLPhreeqcRM yrm;
	IRM_RESULT status;
	int nxyz = 40;
	// Set GridCellCount
	yrm.YAMLSetGridCellCount(nxyz);
	yrm.YAMLThreadCount(3);
	std::string YAML_filename = "TestAllMethods_cpp.yaml";
#ifdef USE_MPI
	// MPI
	int mpi_myself;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		exit(4);
	}
	if (mpi_myself == 0)
	{
		yrm.WriteYAMLDoc(YAML_filename);
		assert(PhreeqcRM::GetGridCellCountYAML(YAML_filename.c_str()) == nxyz);
	}
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)  //  make sure yaml is finished being written
	{
		exit(4);
	}
	BMIPhreeqcRM bmi(nxyz, MPI_COMM_WORLD);
	if (mpi_myself > 0)
	{
		bmi.MpiWorker();
		bmi.Finalize();
		return;
	}
#else
	// OpenMP or serial
	yrm.WriteYAMLDoc(YAML_filename);
	BMIPhreeqcRM bmi;
#endif
	// Use YAML file to initialize
	bmi.Initialize(YAML_filename);   // void function
	bmi.InitializeYAML(YAML_filename);
	std::cerr << "Initialize\n";
	//
	// Use all BMIPhreeqcRM methods roughly in order of use
	// 

	//-------
	bmi.GetValue("GridCellCount", nxyz);
	nxyz = bmi.GetGridCellCount();
	int *i_ptr = (int*)bmi.GetValuePtr("GridCellCount");
	std::cerr << "GetValue('GridCellCount') \n";
	//-------
	int n = bmi.GetThreadCount();
	std::cerr << "GetThreadCount " << n << "\n";
	//-------
	// Inactive cells or symmetry
	std::vector<int> grid2chem(nxyz, -1);
	for (size_t i = 0; i < nxyz / 2; i++)
	{
		grid2chem[i] = i;
	}
	status = bmi.CreateMapping(grid2chem);
	std::cerr << "CreateMapping \n";
	//-------
	bmi.LoadDatabase("phreeqc.dat");
	std::cerr << "LoadDatabase\n";
	//
	// Set properties
	// 
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
	status = bmi.SetDumpFileName("TestAllMethods_cpp.dump");
	std::cerr << "SetDumpFileName \n";
	//-------
	status = bmi.SetFilePrefix("TestAllMethods_cpp");
	bmi.SetValue("FilePrefix", "TestAllMethods_cpp");
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
	bmi.SetValue("SelectedOutputOn", true);
	std::cerr << "SetSelectedOutputOn \n";
	//-------
	status = bmi.SetUnitsExchange(1);
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
	std::vector<int> vi(nxyz, 1);
	status = bmi.SetPrintChemistryMask(vi);
	std::cerr << "SetPrintChemistryMask \n";
	//-------
	status = bmi.SetPrintChemistryOn(false, true, false);
	std::cerr << "SetPrintChemistryOn \n";
	//
	// Define reactants available for initial 
	// and boundary conditions in this file
	//
	status = bmi.RunFile(true, true, true, "all_reactants.pqi");
	std::cerr << "RunFile \n";
	//-------
	bmi.AddOutputVars("AddOutputVars", "True");
	bmi.AddOutputVars("SolutionProperties", "True");
	bmi.AddOutputVars("SolutionTotalMolalities", "True");
	bmi.AddOutputVars("ExchangeMolalities", "True");
	bmi.AddOutputVars("SurfaceMolalities", "True");
	bmi.AddOutputVars("EquilibriumPhases", "True");
	bmi.AddOutputVars("Gases", "True");
	bmi.AddOutputVars("KineticReactants", "True");
	bmi.AddOutputVars("SolidSolutions", "True");
	bmi.AddOutputVars("CalculateValues", "True");
	bmi.AddOutputVars("SolutionActivities", "H+ Ca+2 Na+");
	bmi.AddOutputVars("SolutionMolalities", "OH- Cl-");
	bmi.AddOutputVars("SaturationIndices", "Calcite Dolomite");
	std::cerr << "AddOutputVars \n";
	//-------
	int ncomps = bmi.FindComponents();
	std::cerr << "FindComponents \n";
	//
	// Methods up to this point are useful 
	// in a YAML initialization file
	// 
	// Lists of reactants found by FindComponents follow
	// 
	int nchem = bmi.GetChemistryCellCount();
	std::cerr << "GetChemistryCellCount \n";
	//-------
	ncomps = bmi.GetComponentCount();
	bmi.GetValue("ComponentCount", ncomps);
	i_ptr = (int*) bmi.GetValuePtr("ComponentCount");
	std::cerr << "GetComponentCount \n";
	//-------
	std::vector<std::string> str_vector = bmi.GetComponents();
	bmi.GetValue("Components", str_vector);
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
	v = bmi.GetGfw();
	bmi.GetValue("Gfw", v);
	double *d_ptr = (double*)bmi.GetValuePtr("Gfw");
	std::cerr << "GetGfw ";
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
	// before populating cells with reactants
	//
	std::string input = "DELETE; -all";
	bmi.RunString(true, false, false, input);
	std::cerr << "RunString \n";
	//-------
	// 
	// Transfer initial conditions to model
	//
	vi.clear();
	vi.resize(nxyz, 1);
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
	vi.clear();
	vi.resize(1, 1);
	std::vector<double> bc;
	status = bmi.InitialPhreeqc2Concentrations(bc, vi);
	//-------
	std::vector<int> vi1(1, 1);
	std::vector<int> vi2(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = bmi.InitialPhreeqc2Concentrations(bc, vi1, vi2, f1);
	//-------
	vi.clear();
	vi.resize(1, 1);
	status = bmi.InitialPhreeqc2SpeciesConcentrations(bc, vi);
	std::cerr << "InitialPhreeqc2SpeciesConcentrations \n";
	//-------
	std::vector<double> bc_species;
	vi1.clear();
	vi1.resize(1, 1);
	vi2.clear();
	vi2.resize(1, -1);
	f1.clear();
	f1.resize(1, 1.0);
	status = bmi.InitialPhreeqc2SpeciesConcentrations(bc_species, vi1, vi2, f1);
	std::cerr << "InitialPhreeqc2SpeciesConcentrations mix \n";
	//
	// Get/Set methods for time steping
	//
	double d = bmi.GetTime();
	bmi.GetValue("Time", d);
	d = bmi.GetCurrentTime();
    d = bmi.GetStartTime();
	d_ptr = (double*)bmi.GetValuePtr("Time");
	std::cerr << "GetTime \n";
	//-------
	status = bmi.SetTime(0.0);
	bmi.SetValue("Time", 0.0);
	std::cerr << "SetTime \n";
	//-------
	d = bmi.GetTimeStep();
	bmi.GetValue("TimeStep", v);
	d_ptr = (double*)bmi.GetValuePtr("TimeStep");
	std::cerr << "GetTimeStep \n";
	//-------
	status = bmi.SetTimeStep(0.0);
	bmi.SetValue("TimeStep", 0.0);
	std::cerr << "SetTimeStep \n";
	//-------
	std::vector<double> c;
	status = bmi.GetConcentrations(c);
	bmi.GetValue("Concentrations", c);
	d_ptr = (double*)bmi.GetValuePtr("Concentrations");
	std::cerr << "GetConcentrations \n";
	//-------
	status = bmi.SetConcentrations(c);
	bmi.SetValue("Concentrations", c);
	std::cerr << "SetConcentrations \n";
	//-------
	status = bmi.GetDensityCalculated(v);
	bmi.GetValue("DensityCalculated", v);
	d_ptr = (double*)bmi.GetValuePtr("DensityCalculated");
	std::cerr << "GetDensityCalculated \n";
	//-------
	status = bmi.SetDensityUser(v);
	bmi.SetValue("DensityUser", v);
	std::cerr << "SetDensityUser \n";
	//-------
	status = bmi.GetGasCompMoles(v);
	std::cerr << "GetGasCompMoles \n";
	//-------
	status = bmi.SetGasCompMoles(v);
	std::cerr << "SetGasCompMoles \n";
	//-------
	status = bmi.GetGasCompPhi(v);
	std::cerr << "GetGasCompPhi \n";
	//-------
	status = bmi.GetGasCompPressures(v);
	std::cerr << "GetGasCompPressures \n";
	//-------
	status = bmi.GetGasPhaseVolume(v);
	std::cerr << "GetGasPhaseVolume \n";
	//-------
	status = bmi.SetGasPhaseVolume(v);
	std::cerr << "SetGasPhaseVolume \n";
	//-------
	for (size_t i = 0; i < bmi.GetComponentCount(); i++)
	{
		status = bmi.GetIthConcentration(i, v);
		//-------
		status = bmi.SetIthConcentration(i, v);
	}
	std::cerr << "GetIthConcentration \n";
	std::cerr << "SetIthConcentration \n";
	//-------
	for (int i = 0; i < bmi.GetSpeciesCount(); i++)
	{
		status = bmi.GetIthSpeciesConcentration(i, v);
		//-------
		status = bmi.SetIthSpeciesConcentration(i, v);
	}
	std::cerr << "GetIthSpeciesConcentration \n";
	std::cerr << "SetIthSpeciesConcentration \n";
	//-------
	v = bmi.GetPorosity();
	bmi.GetValue("Porosity", v);
	d_ptr = (double*)bmi.GetValuePtr("Porosity");
	std::cerr << "GetPorosity \n";
	//-------
	status = bmi.SetPorosity(v);
	bmi.SetValue("Porosity", v);
	std::cerr << "SetPorosity \n";
	//-------
	v = bmi.GetPressure();
	bmi.GetValue("Pressure", v);
	d_ptr = (double*)bmi.GetValuePtr("Pressure");
	std::cerr << "GetPressure \n";
	//-------
	status = bmi.SetPressure(v);
	bmi.SetValue("Pressure", v);
	std::cerr << "SetPressure \n";
	//-------
	status = bmi.GetSaturationCalculated(v);
	bmi.GetValue("SaturationCalculated", v);
	d_ptr = (double*)bmi.GetValuePtr("SaturationCalculated");
	std::cerr << "GetSaturationCalculated \n";
	//-------
	status = bmi.SetSaturationUser(v);
	bmi.SetValue("SaturationUser", v);
	std::cerr << "SetSaturationUser \n";
	//-------
	v = bmi.GetSolutionVolume();
	bmi.GetValue("SolutionVolume", v);
	d_ptr = (double*)bmi.GetValuePtr("SolutionVolume");
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
	v = bmi.GetTemperature();
	bmi.GetValue("Temperature", v);
	d_ptr = (double*)bmi.GetValuePtr("Temperature");
	std::cerr << "GetTemperature \n";
	//-------
	status = bmi.SetTemperature(v);
	bmi.SetValue("Temperature", v);
	std::cerr << "SetTemperature \n";
	//-------
	v = bmi.GetViscosity();
	bmi.GetValue("Viscosity", v);
	d_ptr = (double*)bmi.GetValuePtr("Viscosity");	
	std::cerr << "GetViscosity \n";
	//
	// Take a time step
	//
	bmi.Update();      // void function
	std::cerr << "Update\n";
	//-------
	status = bmi.RunCells();
	std::cerr << "RunCells\n";
	//-------
	bmi.UpdateUntil(86400.0);      // void function
	std::cerr << "UpdateUntil\n";
	//
	// Selected output
	//
	status = bmi.SetNthSelectedOutput(0);
	bmi.SetValue("NthSelectedOutput", 0);
	std::cerr << "SetNthSelectedOutput \n";
	//-------
	int n_user = bmi.GetCurrentSelectedOutputUserNumber();
	bmi.GetValue("CurrentSelectedOutputUserNumber", n_user);
	std::cerr << "GetCurrentSelectedOutputUserNumber \n";
	//-------
	n = bmi.GetNthSelectedOutputUserNumber(0);
	std::cerr << "GetNthSelectedOutputUserNumber \n";
	//-------
	status = bmi.GetSelectedOutput(v);
	bmi.GetValue("SelectedOutput", v);
	std::cerr << "GetSelectedOutput \n";
	//-------
	n = bmi.GetSelectedOutputColumnCount();
	bmi.GetValue("SelectedOutputColumnCount", n);
	std::cerr << "GetSelectedOutputColumnCount \n";
	//-------
	n = bmi.GetSelectedOutputCount();
	bmi.GetValue("SelectedOutputCount", n);
	std::cerr << "GetSelectedOutputCount \n";
	//-------
	status = bmi.GetSelectedOutputHeadings(str_vector);
	bmi.GetValue("SelectedOutputHeadings", str_vector);
	std::cerr << "GetSelectedOutputHeadings \n";
	//-------
	bool b = bmi.GetSelectedOutputOn();
	bmi.GetValue("SelectedOutputOn", b);
	bool *b_ptr = (bool*)bmi.GetValuePtr("SelectedOutputOn");
	std::cerr << "GetSelectedOutputOn \n";
	//-------
	n = bmi.GetSelectedOutputRowCount(); 
	bmi.GetValue("SelectedOutputRowCount", n);
	std::cerr << "GetSelectedOutputRowCount \n";
	//-------
	status = bmi.SetCurrentSelectedOutputUserNumber(333);
	std::cerr << "SetCurrentSelectedOutputUserNumber \n";
	//-------
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
	bmi.GetValue("ErrorString", str);
	std::cerr << "GetErrorString \n";
	//-------
	str = bmi.GetFilePrefix();
	bmi.GetValue("FilePrefix", str);
	std::cerr << "GetFilePrefix \n";
	//-------
	vi = bmi.GetForwardMapping();
	std::cerr << "GetForwardMapping \n";
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
	b = bmi.GetPartitionUZSolids();
	std::cerr << "GetPartitionUZSolids \n";
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
	d = bmi.GetRebalanceFraction();
	std::cerr << "GetRebalanceFraction \n";
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
	d = bmi.GetTimeConversion();
	std::cerr << "GetTimeConversion \n";
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
	//-------
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
#ifndef USE_MPI
	BMIPhreeqcRM *bmi2 = new BMIPhreeqcRM(10, 2); // Make another instance
	std::cerr << "Make a new instance with new. \n";
	//-------
	status = bmi2->CloseFiles(); 
	std::cerr << "CloseFiles \n";
	//-------
	bmi2->Finalize();
	delete bmi2; // delete new instance
#endif
	//-------
	vi.clear();
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
	//
	// BMI Methods
	//
	str = bmi.GetComponentName();
	std::cerr << "GetComponentName \n";
	//-------
	d = bmi.GetCurrentTime();
	std::cerr << "GetCurrentTime \n";
	//-------
	d = bmi.GetEndTime();
	std::cerr << "GetEndTime \n";
	//-------
    n = bmi.GetGridRank(0);
	std::cerr << "GetGridRank \n";
	//-------
    n = bmi.GetGridSize(0);
	std::cerr << "GetGridSize \n";
	//-------
    std::string gtype = bmi.GetGridType(0);
	std::cerr << "GetGridType \n";
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
	n = bmi.GetPointableItemCount();
	std::cerr << "GetPointableItemCount \n";
	//-------
	str_vector = bmi.GetPointableVarNames();
	std::cerr << "GetPointableVarNames \n";
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
	//-------
	//bmi.Initialize(YAML_filename);
	// See above
	bmi.SetValue("Time", 1.0);    // void method
	std::cerr << "SetValue\n";
	//-------
	bmi.Update();    // void method
	std::cerr << "Update\n";
	//-------	
 	bmi.UpdateUntil(864000.0);      // void function
	std::cerr << "UpdateUntil\n";

	std::cerr << "AddOutputVars\n";
	str_vector = bmi.GetOutputVarNames();
	for (size_t i = 0; i < str_vector.size(); i++)
	{
		int	itemsize = bmi.GetVarItemsize(str_vector[i]);
		int	nbytes = bmi.GetVarNbytes(str_vector[i]);
		std::string	vtype = bmi.GetVarType(str_vector[i]);
		if (itemsize == 0) itemsize++;
		if (nbytes == 0) nbytes++;
		int dim = nbytes / itemsize;
		if (vtype == "double")
		{
			if (dim == 1)
			{
				double dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<double> dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
		else if (vtype == "int")
		{
			if (dim == 1)
			{
				int dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<int> dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		} else if (vtype == "bool")
		{
			if (dim == 1)
			{
				bool dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest << std::endl;
			}
		}
		else if (vtype == "std::string")
		{
			if (dim == 1)
			{
				std::string dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest << std::endl;
			}
			else
			{
				std::vector<std::string> dest;
				bmi.GetValue(str_vector[i], dest);
				std::cerr << "     " << str_vector[i] << "  " << dest[0] << std::endl;
			}
		}
	}


	//-------	
	bmi.MpiWorkerBreak();
	bmi.Finalize();    // void method
	std::cerr << "Finalize \n";
	//Should be private: status = bmi.ReturnHandler();
	//TODO status = bmi.MpiAbort();
	//TODO status = bmi.SetMpiWorkerCallbackC();
	//TODO status = bmi.SetMpiWorkerCallbackCookie();
	std::cerr << "Success.\n";
	return;
}
#endif // YAML
