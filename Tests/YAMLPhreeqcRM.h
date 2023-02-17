#include <map>
#include <string>
#include "yaml-cpp/yaml.h"
#pragma once
class YAMLPhreeqcRM
{

public:
	YAMLPhreeqcRM();
	void clear();
	//void MakeMethodMap();
	const YAML::Node &GetYAML_doc() { return this->YAML_doc; };
	YAML::Node YAML_doc;
	// methods
	void YAMLCloseFiles(void);
	void YAMLCreateMapping(std::vector< int >& grid2chem);
	void YAMLFindComponents();
	void YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1);
	void YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
	void YAMLInitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
	void YAMLLoadDatabase(std::string database);
	void YAMLOpenFiles(void);
	void YAMLOutputMessage(std::string str);
	void YAMLRunCells(void);
	void YAMLRunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
	void YAMLRunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
	void YAMLScreenMessage(std::string str);
	void YAMLSetComponentH2O(bool tf);
	void YAMLSetConcentrations(std::vector< double > c);
	void YAMLSetCurrentSelectedOutputUserNumber(int n_user);
	void YAMLSetDensity(std::vector< double > density);
	void YAMLSetDumpFileName(std::string dump_name);
	void YAMLSetErrorHandlerMode(int mode);
	void YAMLSetErrorOn(bool tf);
	void YAMLSetFilePrefix(std::string prefix);
	void YAMLSetGasCompMoles(std::vector< double > gas_moles);
	void YAMLSetGasPhaseVolume(std::vector< double > gas_volume);
	void YAMLSetPartitionUZSolids(bool tf);
	void YAMLSetPorosity(std::vector< double > por);
	const YAML::Node& GetYAMLDoc() { return this->YAML_doc; }
	void YAMLSetPressure(std::vector< double > p);
	void YAMLSetPrintChemistryMask(std::vector< int > cell_mask);
	void YAMLSetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
	void YAMLSetRebalanceByCell(bool tf);
	void YAMLSetRebalanceFraction(double f);
	void YAMLSetRepresentativeVolume(std::vector< double > rv);
	void YAMLSetSaturation(std::vector< double > sat);
	void YAMLSetScreenOn(bool tf);
	void YAMLSetSelectedOutputOn(bool tf);
	void YAMLSetSpeciesSaveOn(bool save_on);
	void YAMLSetTemperature(std::vector< double > t);
	void YAMLSetTime(double time);
	void YAMLSetTimeConversion(double conv_factor);
	void YAMLSetTimeStep(double time_step);
	void YAMLSetUnitsExchange(int option);
	void YAMLSetUnitsGasPhase(int option);
	void YAMLSetUnitsKinetics(int option);
	void YAMLSetUnitsPPassemblage(int option);
	void YAMLSetUnitsSolution(int option);
	void YAMLSetUnitsSSassemblage(int option);
	void YAMLSetUnitsSurface(int option);
	void YAMLSpeciesConcentrations2Module(std::vector< double > species_conc);
	void YAMLStateSave(int istate);
	void YAMLStateApply(int istate);
	void YAMLStateDelete(int istate);
	void YAMLUseSolutionDensityVolume(bool tf);
	void YAMLWarningMessage(std::string warnstr);
	// data
	std::map<std::string, int> method_map;
};

