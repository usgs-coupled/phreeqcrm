#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "YAMLPhreeqcRM.h"

std::map<size_t, YAMLPhreeqcRM*> YAMLPhreeqcRM::Instances;
std::mutex YAMLPhreeqcRM::InstancesLock;
size_t YAMLPhreeqcRM::InstancesIndex = 0;

YAMLPhreeqcRM::YAMLPhreeqcRM()
{
	InstancesLock.lock();
	this->Index = YAMLPhreeqcRM::InstancesIndex++;
	std::map<size_t, YAMLPhreeqcRM*>::value_type instance(this->Index, this);
	YAMLPhreeqcRM::Instances.insert(instance);
	InstancesLock.unlock(); 
	this->style = YAML::EmitterStyle::value::Default; // Default, Flow, Block
}
YAMLPhreeqcRM::~YAMLPhreeqcRM()
{
	InstancesLock.lock();
	std::map<size_t, YAMLPhreeqcRM*>::iterator it = YAMLPhreeqcRM::Instances.find(this->Index);
	if (it != YAMLPhreeqcRM::Instances.end())
	{
		YAMLPhreeqcRM::Instances.erase(it);
	}
	InstancesLock.unlock();
}
void YAMLPhreeqcRM::Clear()
{
	YAML::Node empty;
	YAML_doc = empty;
}
int YAMLPhreeqcRM::GetId(void)const
{
	return (int)this->Index;
}
void YAMLPhreeqcRM::WriteYAMLDoc(std::string file_name)
{
	std::ofstream ofs = std::ofstream(file_name.c_str(), std::ofstream::out);
	ofs << this->GetYAMLDoc();
	ofs.close();
}
void YAMLPhreeqcRM::YAMLAddOutputVars(std::string option, std::string definition)
{
	YAML::Node node;
	node["key"] = "AddOutputVars";
	node["option"] = option;
	node["definition"] = definition;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLCloseFiles (void)
{
	YAML::Node node;
	node["key"] = "CloseFiles";
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLCreateMapping(std::vector< int >& grid2chem)
{
	YAML::Node node;
	node["key"] = "CreateMapping";
	node["grid2chem"] = grid2chem;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLDumpModule(bool dump_on, bool append)
{
	YAML::Node node;
	node["key"] = "DumpModule";
	node["dump_on"] = dump_on;
	node["append"] = append;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLFindComponents()
{
	YAML::Node node;
	node["key"] = "FindComponents";
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLInitialSolutions2Module(std::vector< int > solutions)
{
	YAML::Node node;
	node["key"] = "InitialSolutions2Module";
	node["solutions"] = solutions;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};	
void YAMLPhreeqcRM::YAMLInitialEquilibriumPhases2Module(std::vector< int > equilibrium_phases)
{
	YAML::Node node;
	node["key"] = "InitialEquilibriumPhases2Module";
	node["equilibrium_phases"] = equilibrium_phases;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLInitialExchanges2Module(std::vector< int > exchanges)
{
	YAML::Node node;
	node["key"] = "InitialExchanges2Module";
	node["exchanges"] = exchanges;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLInitialSurfaces2Module(std::vector< int > surfaces)
{
	YAML::Node node;
	node["key"] = "InitialSurfaces2Module";
	node["surfaces"] = surfaces;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLInitialGasPhases2Module(std::vector< int > gas_phases)
{
	YAML::Node node;
	node["key"] = "InitialGasPhases2Module";
	node["gas_phases"] = gas_phases;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLInitialSolidSolutions2Module(std::vector< int > solid_solutions)
{
	YAML::Node node;
	node["key"] = "InitialSolidSolutions2Module";
	node["solid_solutions"] = solid_solutions;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLInitialKinetics2Module(std::vector< int > kinetics)
{
	YAML::Node node;
	node["key"] = "InitialKinetics2Module";
	node["kinetics"] = kinetics;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};

void YAMLPhreeqcRM::YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1) 
{
	YAML::Node node;
	node["key"] = "InitialPhreeqc2Module";
	node["initial_conditions1"] = initial_conditions1;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};

void YAMLPhreeqcRM::YAMLInitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1) 
{
	YAML::Node node;
	node["key"] = "InitialPhreeqc2Module";
	node["initial_conditions1"] = initial_conditions1;
	node["initial_conditions2"] = initial_conditions2;
	node["fraction1"] = fraction1;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
};

void YAMLPhreeqcRM::YAMLInitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers) 
{
	YAML::Node node;
	node["key"] = "InitialPhreeqcCell2Module";
	node["n"] = n;
	node["cell_numbers"] = cell_numbers;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}

void YAMLPhreeqcRM::YAMLLoadDatabase(std::string database) 
{
	YAML::Node node;
	node["key"] = "LoadDatabase";
	node["database"] = database;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLLogMessage(std::string str)
{
	YAML::Node node;
	node["key"] = "LogMessage";
	node["str"] = str;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLOpenFiles(void)
{
	YAML::Node node;
	node["key"] = "OpenFiles";
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLOutputMessage(std::string str) 
{
	YAML::Node node;
	node["key"] = "OutputMessage";
	node["str"] = str;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLRunCells(void) 
{
	YAML::Node node;
	node["key"] = "RunCells";
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLRunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name) 
{
	YAML::Node node;
	node["key"] = "RunFile";
	node["workers"] = workers;
	node["initial_phreeqc"] = initial_phreeqc;
	node["utility"] = utility;
	node["chemistry_name"] = chemistry_name;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLRunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string) 
{
	YAML::Node node;
	node["key"] = "RunString";
	node["workers"] = workers;
	node["initial_phreeqc"] = initial_phreeqc;
	node["utility"] = utility;
	node["input_string"] = input_string;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLScreenMessage(std::string str) 
{
	YAML::Node node;
	node["key"] = "ScreenMessage";
	node["str"] = str;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetComponentH2O(bool tf)
{
	YAML::Node node;
	node["key"] = "SetComponentH2O";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetConcentrations(std::vector< double >& c) 
{
	YAML::Node node;
	node["key"] = "SetConcentrations";
	node["c"] = c;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetCurrentSelectedOutputUserNumber(int n_user) 
{
	YAML::Node node;
	node["key"] = "SetCurrentSelectedOutputUserNumber";
	node["n_user"] = n_user;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetDensityUser(std::vector< double > density)
{
	YAML::Node node;
	node["key"] = "SetDensityUser";
	node["density"] = density;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetDumpFileName(std::string dump_name) 
{
	YAML::Node node;
	node["key"] = "SetDumpFileName";
	node["dump_name"] = dump_name;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetErrorHandlerMode(int mode) 
{
	YAML::Node node;
	node["key"] = "SetErrorHandlerMode";
	node["mode"] = mode;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetErrorOn(bool tf) 
{
	YAML::Node node;
	node["key"] = "SetErrorOn";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetFilePrefix(std::string prefix) 
{
	YAML::Node node;
	node["key"] = "SetFilePrefix";
	node["prefix"] = prefix;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetGasCompMoles(std::vector< double > gas_moles) 
{
	YAML::Node node;
	node["key"] = "SetGasCompMoles";
	node["gas_moles"] = gas_moles;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetGasPhaseVolume(std::vector< double > gas_volume) 
{
	YAML::Node node;
	node["key"] = "SetGasPhaseVolume";
	node["gas_volume"] = gas_volume;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetGridCellCount(int count)
{
	YAML::Node node;
	node["key"] = "SetGridCellCount";
	node["count"] = count;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetNthSelectedOutput(int n)
{
	YAML::Node node;
	node["key"] = "SetNthSelectedOutput";
	node["n"] = n;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetPartitionUZSolids(bool tf)  
{
	YAML::Node node;
	node["key"] = "SetPartitionUZSolids";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetPorosity(std::vector< double > por) 
{
	YAML::Node node;
	node["key"] = "SetPorosity";
	node["por"] = por;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetPressure(std::vector< double > p) 
{
	YAML::Node node;
	node["key"] = "SetPressure";
	node["p"] = p;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetPrintChemistryMask(std::vector< int > cell_mask) 
{
	YAML::Node node;
	node["key"] = "SetPrintChemistryMask";
	node["cell_mask"] = cell_mask;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility) 
{
	YAML::Node node;
	node["key"] = "SetPrintChemistryOn";
	node["workers"] = workers;
	node["initial_phreeqc"] = initial_phreeqc;
	node["utility"] = utility;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetRebalanceByCell(bool tf) 
{
	YAML::Node node;
	node["key"] = "SetRebalanceByCell";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetRebalanceFraction(double f) 
{
	YAML::Node node;
	node["key"] = "SetRebalanceFraction";
	node["f"] = f;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetRepresentativeVolume(std::vector< double > rv) 
{
	YAML::Node node;
	node["key"] = "SetRepresentativeVolume";
	node["rv"] = rv;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetSaturationUser(std::vector< double > sat)
{
	YAML::Node node;
	node["key"] = "SetSaturationUser";
	node["sat"] = sat;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetScreenOn(bool tf) 
{
	YAML::Node node;
	node["key"] = "SetScreenOn";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetSelectedOutputOn(bool tf) 
{
	YAML::Node node;
	node["key"] = "SetSelectedOutputOn";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetSpeciesSaveOn(bool save_on) 
{
	YAML::Node node;
	node["key"] = "SetSpeciesSaveOn";
	node["save_on"] = save_on;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetTemperature(std::vector< double > t)
{
	YAML::Node node;
	node["key"] = "SetTemperature";
	node["t"] = t;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetTime(double time)
{
	YAML::Node node;
	node["key"] = "SetTime";
	node["time"] = time;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetTimeConversion(double conv_factor)
{
	YAML::Node node;
	node["key"] = "SetTimeConversion";
	node["conv_factor"] = conv_factor;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetTimeStep(double time_step)
{
	YAML::Node node;
	node["key"] = "SetTimeStep";
	node["time_step"] = time_step;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSetUnitsExchange(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsExchange";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsGasPhase(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsGasPhase";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsKinetics(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsKinetics";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsPPassemblage(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsPPassemblage";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSolution(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsSolution";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSSassemblage(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsSSassemblage";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLSetUnitsSurface(int option)
{
	YAML::Node node;
	node["key"] = "SetUnitsSurface";
	node["option"] = option;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLSpeciesConcentrations2Module(std::vector< double > species_conc) 
{
	YAML::Node node;
	node["key"] = "SpeciesConcentrations2Module";
	node["species_conc"] = species_conc;
	node.SetStyle(this->style);
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLStateSave(int istate)
{
	YAML::Node node;
	node["key"] = "StateSave";
	node["istate"] = istate;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLStateApply(int istate)
{
	YAML::Node node;
	node["key"] = "StateApply";
	node["istate"] = istate;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLStateDelete(int istate)
{
	YAML::Node node;
	node["key"] = "StateDelete";
	node["istate"] = istate;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLThreadCount(int nthreads)
{
	YAML::Node node;
	node["key"] = "ThreadCount";
	node["nthreads"] = nthreads;
	YAML_doc.push_back(node);
	return;
}
void YAMLPhreeqcRM::YAMLUseSolutionDensityVolume(bool tf)  
{
	YAML::Node node;
	node["key"] = "UseSolutionDensityVolume";
	node["tf"] = tf;
	YAML_doc.push_back(node);
	return;
};
void YAMLPhreeqcRM::YAMLWarningMessage(std::string warnstr) 
{
	YAML::Node node;
	node["key"] = "WarningMessage";
	node["warnstr"] = warnstr;
	YAML_doc.push_back(node);
	return;
}
//
// helper functions
//
int
YAMLPhreeqcRMLib::CreateYAMLPhreeqcRM(void)
{
	int n = IRM_OUTOFMEMORY;
	YAMLPhreeqcRM* YAMLPhreeqcRMPtr;
	try
	{
		YAMLPhreeqcRMPtr = new YAMLPhreeqcRM;
		n = (int)YAMLPhreeqcRMPtr->Index;
	}
	catch (const std::bad_alloc&)
	{
		return IRM_OUTOFMEMORY;
	}
	return n;
}

IRM_RESULT
YAMLPhreeqcRMLib::DestroyYAMLPhreeqcRM(int id)
{
	IRM_RESULT retval = IRM_BADINSTANCE;
	if (id >= 0)
	{
		if (YAMLPhreeqcRM* ptr = YAMLPhreeqcRMLib::GetInstance(id))
		{
			delete ptr;
			retval = IRM_OK;
		}
	}
	return retval;
}

YAMLPhreeqcRM*
YAMLPhreeqcRMLib::GetInstance(int id)
{
	YAMLPhreeqcRM* instance = 0;
	YAMLPhreeqcRM::InstancesLock.lock();
	std::map<size_t, YAMLPhreeqcRM*>::iterator it = YAMLPhreeqcRM::Instances.find(size_t(id));
	if (it != YAMLPhreeqcRM::Instances.end())
	{
		instance = (*it).second;
	}
	YAMLPhreeqcRM::InstancesLock.unlock();
	return instance;
}
#endif