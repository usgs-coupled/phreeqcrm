#if !defined(VARMANAGER_H_INCLUDED)
#define VARMANAGER_H_INCLUDED
#include <vector>
#include <list>
#include <set>
#include <map>
#include <string>
#include <assert.h>
#include "PhreeqcRM.h"
#include "BMI_var.h"
class IRM_DLL_EXPORT VarManager
{
public:
	std::map<std::string, Variant> VariantMap;
	std::set<std::string> UpdateSet;
	VarManager(PhreeqcRM* rm_ptr);
	enum class VAR_TASKS {
		Init,
		Update,
		GetPtr,
		Info,
		GetVar,
		SetVar,
		no_op
	};
	VAR_TASKS task;
	enum class VARIABLES {
		ComponentCount,
		Components,
		Concentrations,
		Density,
		ErrorString,
		FilePrefix,
		Gfw,
		GridCellCount,
		InputVarNames,
		NthSelectedOutput,
		OutputVarNames,
		Saturation,
		SelectedOutput,
		SelectedOutputColumnCount,
		SelectedOutputCount,
		SelectedOutputHeadings,
		SelectedOutputRowCount,
		SolutionVolume,
		Time,
		TimeStep,
		CurrentSelectedOutputUserNumber,
		Porosity,
		Pressure,
		SelectedOutputOn,
		Temperature
	};
	void ComponentCount_var(PhreeqcRM* rm_ref);
#ifdef SKIP
	void Components_var(PhreeqcRM* rm_ref);
	void Concentrations_var(PhreeqcRM* rm_ref);
	void Density_var(PhreeqcRM* rm_ref);
	void ErrorString_var(PhreeqcRM* rm_ref);
	void FilePrefix_var(PhreeqcRM* rm_ref);
	void Gfw_var(PhreeqcRM* rm_ref);
	void GridCellCount_var(PhreeqcRM* rm_ref);
	void InputVarNames_var(PhreeqcRM* rm_ref);
	void NthSelectedOutput_var(PhreeqcRM* rm_ref);
	void OutputVarNames_var(PhreeqcRM* rm_ref);
	void Saturation_var(PhreeqcRM* rm_ref);
	void SelectedOutput_var(PhreeqcRM* rm_ref);
	void SelectedOutputColumnCount_var(PhreeqcRM* rm_ref);
	void SelectedOutputCount_var(PhreeqcRM* rm_ref);
	void SelectedOutputHeadings_var(PhreeqcRM* rm_ref);
	void SelectedOutputRowCount_var(PhreeqcRM* rm_ref);
	void SolutionVolume_var(PhreeqcRM* rm_ref);
	void Time_var(PhreeqcRM* rm_ref);
	void TimeStep_var(PhreeqcRM* rm_ref);
	void CurrentSelectedOutputUserNumber_var(PhreeqcRM* rm_ref);
	void Porosity_var(PhreeqcRM* rm_ref);
	void Pressure_var(PhreeqcRM* rm_ref);
	void SelectedOutputOn_var(PhreeqcRM* rm_ref);
	void Temperature_var(PhreeqcRM* rm_ref);
#endif
	typedef void (VarManager::* VarFunction)(PhreeqcRM* rm_ptr);
	//typedef VarManager* (*NewDogFunction)(void);
	//typedef void (VarManager::* VarFunction)(PhreeqcRM* rm_ptr); // function pointer type
	//typedef std::map<std::string, VarFunction> VarFunction_map;
	//VarFunction_map FnMap;
	//void test() { VarFunction x = VarManager::Concentrations_var; 
	//VarVariant vv;
	//vv.SetFn(VarFux);
};
#endif

