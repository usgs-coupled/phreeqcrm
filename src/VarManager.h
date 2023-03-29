#if !defined(VARMANAGER_H_INCLUDED)
#define VARMANAGER_H_INCLUDED
#include <vector>
#include <list>
#include <set>
#include <map>
#include <string>
#include <assert.h>
#include "PhreeqcRM.h"
//class PhreeqcRM;
#include "BMI_var.h"
class VARS;

class IRM_DLL_EXPORT VarManager
{
public:
	enum class VAR_TASKS {
		RMUpdate,
		UpdateState,
		GetPtr,  
		GetVar,
		SetVar,
		no_op
	};
	enum class VARS {
		NotFound,
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
	// Constructor
	VarManager(PhreeqcRM* rm_ptr);
	// Data
	PhreeqcRM* rm_ptr;
public:
	BMIVariant VarExchange;
	std::set<VARS> PointerSet;
	VARS CurrentVar;
	std::map < std::string, VARS> EnumMap;
	VAR_TASKS task;
	std::map<VARS, BMIVariant> VariantMap;
	// Methods
	VARS GetEnum(std::string name);
	void RM2BMIUpdate(std::string);

	//std::map<std::string, BMIVariant> VariantMap;
	VARS GetCurrentVar() { return this->CurrentVar; }
	void SetCurrentVar(VarManager::VARS v) { this->CurrentVar = v; }
	// Function pointer definition
	typedef void (VarManager::* VarFunction)(void);
	//!typedef void (VarManager::* VarFunction)(PhreeqcRM* rm_ptr);
	//typedef VarManager* (*NewDogFunction)(void);
	//typedef void (VarManager::* VarFunction)(PhreeqcRM* rm_ptr); // function pointer type
	//typedef std::map<std::string, VarFunction> VarFunction_map;
	//VarFunction_map FnMap;
	//void test() { VarFunction x = VarManager::Concentrations_var; 
	//VarVariant vv;
	//vv.SetFn(VarFux);

	// Var functions
	void ComponentCount_var();
	//void Porosity_var();
};
#endif

