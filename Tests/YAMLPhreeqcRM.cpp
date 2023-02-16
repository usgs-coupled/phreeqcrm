#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "YAMLPhreeqcRM.h"
/* ---------------------------------------------------------------------- */
void
YAMLPhreeqcRM::MakeMethodMap()
/* ---------------------------------------------------------------------- */
{
	method_map["PhreeqcRM"] = METHODS::PHREEQCRM;
	method_map["CloseFiles"] = METHODS::CLOSEFILES;
	method_map["Concentrations2Utility"] = METHODS::CONCENTRATIONS2UTILITY;
	method_map["CreateMapping"] = METHODS::CREATEMAPPING;
	method_map["DecodeError"] = METHODS::DECODEERROR;
	method_map["DumpModule"] = METHODS::DUMPMODULE;
	method_map["ErrorHandler"] = METHODS::ERRORHANDLER;
	method_map["ErrorMessage"] = METHODS::ERRORMESSAGE;
	method_map["FindComponents"] = METHODS::FINDCOMPONENTS;
	method_map["GetBackwardMapping"] = METHODS::GETBACKWARDMAPPING;
	method_map["GetChemistryCellCount"] = METHODS::GETCHEMISTRYCELLCOUNT;
	method_map["GetComponentCount"] = METHODS::GETCOMPONENTCOUNT;
	method_map["GetComponents"] = METHODS::GETCOMPONENTS;
	method_map["GetConcentrations"] = METHODS::GETCONCENTRATIONS;
	method_map["GetDatabaseFileName"] = METHODS::GETDATABASEFILENAME;
	method_map["GetDensity"] = METHODS::GETDENSITY;
	method_map["GetEndCell"] = METHODS::GETENDCELL;
	method_map["GetEquilibriumPhases"] = METHODS::GETEQUILIBRIUMPHASES;
	method_map["GetEquilibriumPhasesCount"] = METHODS::GETEQUILIBRIUMPHASESCOUNT;
	method_map["GetErrorHandlerMode"] = METHODS::GETERRORHANDLERMODE;
	method_map["GetErrorString"] = METHODS::GETERRORSTRING;
	method_map["GetExchangeNames"] = METHODS::GETEXCHANGENAMES;
	method_map["GetExchangeSpecies"] = METHODS::GETEXCHANGESPECIES;
	method_map["GetExchangeSpeciesCount"] = METHODS::GETEXCHANGESPECIESCOUNT;
	method_map["GetFilePrefix"] = METHODS::GETFILEPREFIX;
	method_map["GetForwardMapping"] = METHODS::GETFORWARDMAPPING;
	method_map["GetGasComponents"] = METHODS::GETGASCOMPONENTS;
	method_map["GetGasComponentsCount"] = METHODS::GETGASCOMPONENTSCOUNT;
	method_map["GetGasCompMoles"] = METHODS::GETGASCOMPMOLES;
	method_map["GetGasCompPressures"] = METHODS::GETGASCOMPPRESSURES;
	method_map["GetGasCompPhi"] = METHODS::GETGASCOMPPHI;
	method_map["GetGasPhaseVolume"] = METHODS::GETGASPHASEVOLUME;
	method_map["GetGfw"] = METHODS::GETGFW;
	method_map["GetGridCellCount"] = METHODS::GETGRIDCELLCOUNT;
	method_map["GetIPhreeqcPointer"] = METHODS::GETIPHREEQCPOINTER;
	method_map["GetKineticReactions"] = METHODS::GETKINETICREACTIONS;
	method_map["GetKineticReactionsCount"] = METHODS::GETKINETICREACTIONSCOUNT;
	method_map["GetMpiMyself"] = METHODS::GETMPIMYSELF;
	method_map["GetMpiTasks"] = METHODS::GETMPITASKS;
	method_map["GetNthSelectedOutputUserNumber"] = METHODS::GETNTHSELECTEDOUTPUTUSERNUMBER;
	method_map["GetPartitionUZSolids"] = METHODS::GETPARTITIONUZSOLIDS;
	method_map["GetPressure"] = METHODS::GETPRESSURE;
	method_map["GetPrintChemistryMask"] = METHODS::GETPRINTCHEMISTRYMASK;
	method_map["GetPrintChemistryOn"] = METHODS::GETPRINTCHEMISTRYON;
	method_map["GetRebalanceByCell"] = METHODS::GETREBALANCEBYCELL;
	method_map["GetRebalanceFraction"] = METHODS::GETREBALANCEFRACTION;
	method_map["GetSaturation"] = METHODS::GETSATURATION;
	method_map["GetSelectedOutput"] = METHODS::GETSELECTEDOUTPUT;
	method_map["GetSelectedOutputColumnCount"] = METHODS::GETSELECTEDOUTPUTCOLUMNCOUNT;
	method_map["GetSelectedOutputCount"] = METHODS::GETSELECTEDOUTPUTCOUNT;
	method_map["GetSelectedOutputHeading"] = METHODS::GETSELECTEDOUTPUTHEADING;
	method_map["GetSelectedOutputOn"] = METHODS::GETSELECTEDOUTPUTON;
	method_map["GetSelectedOutputRowCount"] = METHODS::GETSELECTEDOUTPUTROWCOUNT;
	method_map["GetSICount"] = METHODS::GETSICOUNT;
	method_map["GetSINames"] = METHODS::GETSINAMES;
	method_map["GetSolidSolutionComponents"] = METHODS::GETSOLIDSOLUTIONCOMPONENTS;
	method_map["GetSolidSolutionComponentsCount"] = METHODS::GETSOLIDSOLUTIONCOMPONENTSCOUNT;
	method_map["GetSolidSolutionNames"] = METHODS::GETSOLIDSOLUTIONNAMES;
	method_map["GetSolutionVolume"] = METHODS::GETSOLUTIONVOLUME;
	method_map["GetSpeciesConcentrations"] = METHODS::GETSPECIESCONCENTRATIONS;
	method_map["GetSpeciesCount"] = METHODS::GETSPECIESCOUNT;
	method_map["GetSpeciesD25"] = METHODS::GETSPECIESD25;
	method_map["GetSpeciesLog10Gammas"] = METHODS::GETSPECIESLOG10GAMMAS;
	method_map["GetSpeciesLog10Molalities"] = METHODS::GETSPECIESLOG10MOLALITIES;
	method_map["GetSpeciesNames"] = METHODS::GETSPECIESNAMES;
	method_map["GetSpeciesSaveOn"] = METHODS::GETSPECIESSAVEON;
	method_map["GetSpeciesStoichiometry"] = METHODS::GETSPECIESSTOICHIOMETRY;
	method_map["GetSpeciesZ"] = METHODS::GETSPECIESZ;
	method_map["GetStartCell"] = METHODS::GETSTARTCELL;
	method_map["GetSurfaceNames"] = METHODS::GETSURFACENAMES;
	method_map["GetSurfaceSpecies"] = METHODS::GETSURFACESPECIES;
	method_map["GetSurfaceSpeciesCount"] = METHODS::GETSURFACESPECIESCOUNT;
	method_map["GetSurfaceTypes"] = METHODS::GETSURFACETYPES;
	method_map["GetTemperature"] = METHODS::GETTEMPERATURE;
	method_map["GetThreadCount"] = METHODS::GETTHREADCOUNT;
	method_map["GetTime"] = METHODS::GETTIME;
	method_map["GetTimeConversion"] = METHODS::GETTIMECONVERSION;
	method_map["GetTimeStep"] = METHODS::GETTIMESTEP;
	method_map["GetUnitsExchange"] = METHODS::GETUNITSEXCHANGE;
	method_map["GetUnitsGasPhase"] = METHODS::GETUNITSGASPHASE;
	method_map["GetUnitsKinetics"] = METHODS::GETUNITSKINETICS;
	method_map["GetUnitsPPassemblage"] = METHODS::GETUNITSPPASSEMBLAGE;
	method_map["GetUnitsSolution"] = METHODS::GETUNITSSOLUTION;
	method_map["GetUnitsSSassemblage"] = METHODS::GETUNITSSSASSEMBLAGE;
	method_map["GetUnitsSurface"] = METHODS::GETUNITSSURFACE;
	method_map["GetWorkers"] = METHODS::GETWORKERS;
	method_map["InitialPhreeqc2Concentrations"] = METHODS::INITIALPHREEQC2CONCENTRATIONS;
	method_map["InitialPhreeqc2Concentrations"] = METHODS::INITIALPHREEQC2CONCENTRATIONS;
	method_map["InitialPhreeqc2Module"] = METHODS::INITIALPHREEQC2MODULE;
	method_map["InitialPhreeqc2Module"] = METHODS::INITIALPHREEQC2MODULE;
	method_map["InitialPhreeqc2SpeciesConcentrations"] = METHODS::INITIALPHREEQC2SPECIESCONCENTRATIONS;
	method_map["InitialPhreeqc2SpeciesConcentrations"] = METHODS::INITIALPHREEQC2SPECIESCONCENTRATIONS;
	method_map["InitialPhreeqcCell2Module"] = METHODS::INITIALPHREEQCCELL2MODULE;
	method_map["LoadDatabase"] = METHODS::LOADDATABASE;
	method_map["LogMessage"] = METHODS::LOGMESSAGE;
	method_map["MpiAbort"] = METHODS::MPIABORT;
	method_map["MpiWorker"] = METHODS::MPIWORKER;
	method_map["MpiWorkerBreak"] = METHODS::MPIWORKERBREAK;
	method_map["OpenFiles"] = METHODS::OPENFILES;
	method_map["OutputMessage"] = METHODS::OUTPUTMESSAGE;
	method_map["RunCells"] = METHODS::RUNCELLS;
	method_map["ReturnHandler"] = METHODS::RETURNHANDLER;
	method_map["RunFile"] = METHODS::RUNFILE;
	method_map["RunString"] = METHODS::RUNSTRING;
	method_map["ScreenMessage"] = METHODS::SCREENMESSAGE;
	method_map["SetComponentH2O"] = METHODS::SETCOMPONENTH2O;
	method_map["SetConcentrations"] = METHODS::SETCONCENTRATIONS;
	method_map["SetCurrentSelectedOutputUserNumber"] = METHODS::SETCURRENTSELECTEDOUTPUTUSERNUMBER;
	method_map["SetDensity"] = METHODS::SETDENSITY;
	method_map["SetDumpFileName"] = METHODS::SETDUMPFILENAME;
	method_map["SetErrorHandlerMode"] = METHODS::SETERRORHANDLERMODE;
	method_map["SetErrorOn"] = METHODS::SETERRORON;
	method_map["SetFilePrefix"] = METHODS::SETFILEPREFIX;
	method_map["SetGasCompMoles"] = METHODS::SETGASCOMPMOLES;
	method_map["SetGasPhaseVolume"] = METHODS::SETGASPHASEVOLUME;
	method_map["SetMpiWorkerCallbackC"] = METHODS::SETMPIWORKERCALLBACKC;
	method_map["SetMpiWorkerCallbackCookie"] = METHODS::SETMPIWORKERCALLBACKCOOKIE;
	method_map["SetMpiWorkerCallbackFortran"] = METHODS::SETMPIWORKERCALLBACKFORTRAN;
	method_map["SetPartitionUZSolids"] = METHODS::SETPARTITIONUZSOLIDS;
	method_map["SetPorosity"] = METHODS::SETPOROSITY;
	method_map["SetPressure"] = METHODS::SETPRESSURE;
	method_map["SetPrintChemistryMask"] = METHODS::SETPRINTCHEMISTRYMASK;
	method_map["SetPrintChemistryOn"] = METHODS::SETPRINTCHEMISTRYON;
	method_map["SetRebalanceByCell"] = METHODS::SETREBALANCEBYCELL;
	method_map["SetRebalanceFraction"] = METHODS::SETREBALANCEFRACTION;
	method_map["SetRepresentativeVolume"] = METHODS::SETREPRESENTATIVEVOLUME;
	method_map["SetSaturation"] = METHODS::SETSATURATION;
	method_map["SetScreenOn"] = METHODS::SETSCREENON;
	method_map["SetSelectedOutputOn"] = METHODS::SETSELECTEDOUTPUTON;
	method_map["SetSpeciesSaveOn"] = METHODS::SETSPECIESSAVEON;
	method_map["SetTemperature"] = METHODS::SETTEMPERATURE;
	method_map["SetTime"] = METHODS::SETTIME;
	method_map["SetTimeConversion"] = METHODS::SETTIMECONVERSION;
	method_map["SetTimeStep"] = METHODS::SETTIMESTEP;
	method_map["SetUnitsExchange"] = METHODS::SETUNITSEXCHANGE;
	method_map["SetUnitsGasPhase"] = METHODS::SETUNITSGASPHASE;
	method_map["SetUnitsKinetics"] = METHODS::SETUNITSKINETICS;
	method_map["SetUnitsPPassemblage"] = METHODS::SETUNITSPPASSEMBLAGE;
	method_map["SetUnitsSolution"] = METHODS::SETUNITSSOLUTION;
	method_map["SetUnitsSSassemblage"] = METHODS::SETUNITSSSASSEMBLAGE;
	method_map["SetUnitsSurface"] = METHODS::SETUNITSSURFACE;
	method_map["SpeciesConcentrations2Module"] = METHODS::SPECIESCONCENTRATIONS2MODULE;
	method_map["StateSave"] = METHODS::STATESAVE;
	method_map["StateApply"] = METHODS::STATEAPPLY;
	method_map["StateDelete"] = METHODS::STATEDELETE;
	method_map["UseSolutionDensityVolume"] = METHODS::USESOLUTIONDENSITYVOLUME;
	method_map["WarningMessage"] = METHODS::WARNINGMESSAGE;
}
YAMLPhreeqcRM::YAMLPhreeqcRM()
{
}
void YAMLPhreeqcRM::clear()
{
	YAML::Node empty;
	YAML_doc = empty;
}
void YAMLPhreeqcRM::YAMLCloseFiles (void)
{
	YAML_doc["CloseFiles"] = true;
	return;
};
void YAMLPhreeqcRM::YAMLCreateMapping(std::vector< int >& grid2chem)
{
	YAML_doc["CreateMapping"] = grid2chem;
	return;
};
void YAMLPhreeqcRM::YAMLFindComponents()
{
	YAML_doc["FindComponents"] = true;
	return;
}
void YAMLPhreeqcRM::YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1) 
{
	YAML_doc["InitialPhreeqc2Module"] = initial_conditions1;
	return;
};
void YAMLPhreeqcRM::YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1) 
{
	YAML::Node node;
	node[initial_conditions1] = initial_conditions1;
	node[initial_conditions2] = initial_conditions2;
	node[fraction1] = fraction1;
	YAML_doc["InitialPhreeqc2Module"] = node;
};
//void YAMLPhreeqcRM::YAMLInitialPhreeqc2SpeciesConcentrations(std::vector< double > destination_c, std::vector< int > boundary_solution1) {};
//void YAMLPhreeqcRM::YAMLInitialPhreeqc2SpeciesConcentrations(std::vector< double > destination_c, std::vector< int > boundary_solution1, std::vector< int > boundary_solution2, std::vector< double > fraction1) {};
//void YAMLPhreeqcRM::YAMLInitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers) {};
//void YAMLPhreeqcRM::YAMLLoadDatabase(std::string database) {};
void YAMLPhreeqcRM::YAMLOpenFiles(void)
{
	YAML_doc["OpenFiles"] = true;
	return;
};
//void YAMLPhreeqcRM::YAMLOutputMessage(std::string str) {};
//void YAMLPhreeqcRM::YAMLRunCells(void) {};
//void YAMLPhreeqcRM::YAMLRunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name) {};
//void YAMLPhreeqcRM::YAMLRunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string) {};
//void YAMLPhreeqcRM::YAMLScreenMessage(std::string str) {};
void YAMLPhreeqcRM::YAMLSetComponentH2O(bool tf)
{
	YAML_doc["SetComponentH2O"] = tf;
	return;
};
//void YAMLPhreeqcRM::YAMLSetConcentrations(std::vector< double > c) {};
//void YAMLPhreeqcRM::YAMLSetCurrentSelectedOutputUserNumber(int n_user) {};
//void YAMLPhreeqcRM::YAMLSetDensity(std::vector< double > density) {};
//void YAMLPhreeqcRM::YAMLSetDumpFileName(std::string dump_name) {};
void YAMLPhreeqcRM::YAMLSetErrorHandlerMode(int mode) 
{
	YAML_doc["SetErrorHandlerMode"] = mode;
	return;
};
//void YAMLPhreeqcRM::YAMLSetErrorOn(bool tf) {};
void YAMLPhreeqcRM::YAMLSetFilePrefix(std::string prefix) 
{
	YAML_doc["SetFilePrefix"] = prefix;
	return;
};
//void YAMLPhreeqcRM::YAMLSetGasCompMoles(std::vector< double > gas_moles) {};
//void YAMLPhreeqcRM::YAMLSetGasPhaseVolume(std::vector< double > gas_volume) {};
void YAMLPhreeqcRM::YAMLSetPartitionUZSolids(bool tf)  
{
	YAML_doc["SetPartitionUZSolids"] = tf;
	return;
};
//void YAMLPhreeqcRM::YAMLSetPorosity(std::vector< double > por) {};
//void YAMLPhreeqcRM::YAMLSetPressure(std::vector< double > p) {};
//void YAMLPhreeqcRM::YAMLSetPrintChemistryMask(std::vector< int > cell_mask) {};
//void YAMLPhreeqcRM::YAMLSetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility) {};
void YAMLPhreeqcRM::YAMLSetRebalanceByCell(bool tf) 
{
	YAML_doc["SetRebalanceByCell"] = tf;
	return;
};
void YAMLPhreeqcRM::YAMLSetRebalanceFraction(double f) 
{
	YAML_doc["SetRebalanceFraction"] = f;
	return;
};
//void YAMLPhreeqcRM::YAMLSetRepresentativeVolume(std::vector< double > rv) {};
//void YAMLPhreeqcRM::YAMLSetSaturation(std::vector< double > sat) {};
//void YAMLPhreeqcRM::YAMLSetScreenOn(bool tf) {};
//void YAMLPhreeqcRM::YAMLSetSelectedOutputOn(bool tf) {};
//void YAMLPhreeqcRM::YAMLSetSpeciesSaveOn(bool save_on) {};
//void YAMLPhreeqcRM::YAMLSetTemperature(std::vector< double > t) {};
//void YAMLPhreeqcRM::YAMLSetTime(double time) {};
//void YAMLPhreeqcRM::YAMLSetTimeConversion(double conv_factor) {};
//void YAMLPhreeqcRM::YAMLSetTimeStep(double time_step) {};
void YAMLPhreeqcRM::YAMLSetUnitsExchange(int option)
{
	YAML_doc["SetUnitsExchange"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsGasPhase(int option)
{
	YAML_doc["SetUnitsGasPhase"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsKinetics(int option)
{
	YAML_doc["SetUnitsKinetics"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsPPassemblage(int option)
{
	YAML_doc["SetUnitsPPassemblage"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSolution(int option)
{
	YAML_doc["SetUnitsSolution"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSSassemblage(int option)
{
	YAML_doc["SetUnitsSSassemblage"] = option;
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSurface(int option)
{
	YAML_doc["SetUnitsSurface"] = option;
	return;
};
//void YAMLPhreeqcRM::YAMLSpeciesConcentrations2Module(std::vector< double > species_conc) {};
//void YAMLPhreeqcRM::YAMLStateSave(int istate) {};
//void YAMLPhreeqcRM::YAMLStateApply(int istate) {};
//void YAMLPhreeqcRM::YAMLStateDelete(int istate) {};
void YAMLPhreeqcRM::YAMLUseSolutionDensityVolume(bool tf)  
{
	YAML_doc["UseSolutionDensityVolume"] = tf;
	return;
};
//void YAMLPhreeqcRM::YAMLWarningMessage(std::string warnstr) {};
//

