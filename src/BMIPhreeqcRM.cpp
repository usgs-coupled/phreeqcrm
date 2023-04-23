#include "BMIPhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <queue>

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif

#include "Phreeqc.h"
#include "IPhreeqcPhast.h"
#include "PhreeqcRM.h"
#include "VarManager.h"

std::map<size_t, BMIPhreeqcRM*> BMIPhreeqcRM::Instances;
size_t BMIPhreeqcRM::InstancesIndex = 0;

//// static BMIPhreeqcRM methods
/* ---------------------------------------------------------------------- */
void
BMIPhreeqcRM::CleanupBMIModuleInstances(void)
/* ---------------------------------------------------------------------- */
{
	std::map<size_t, BMIPhreeqcRM*>::iterator it = BMIPhreeqcRM::Instances.begin();
	std::vector<BMIPhreeqcRM*> bmirm_list;
	for (; it != BMIPhreeqcRM::Instances.end(); it++)
	{
		bmirm_list.push_back(it->second);
	}
	for (size_t i = 0; i < bmirm_list.size(); i++)
	{
		delete bmirm_list[i];
	}
}
/* ---------------------------------------------------------------------- */
int
BMIPhreeqcRM::CreateBMIModule(int nxyz, MP_TYPE nthreads)
/* ---------------------------------------------------------------------- */
{
	//_CrtSetDbgbool ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_crtBreakAlloc = 5144;
	int n = IRM_OUTOFMEMORY;
	try
	{
		BMIPhreeqcRM* bmirm_ptr = new BMIPhreeqcRM(nxyz, nthreads);
		if (bmirm_ptr)
		{
			n = (int)bmirm_ptr->Index;
			BMIPhreeqcRM::Instances[n] = bmirm_ptr;
			bmirm_ptr->language = "F90";
			return n;
		}
	}
	catch (...)
	{
		return IRM_OUTOFMEMORY;
	}
	return IRM_OUTOFMEMORY;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMIPhreeqcRM::DestroyBMIModule(int id)
/* ---------------------------------------------------------------------- */
{
	IRM_RESULT retval = IRM_BADINSTANCE;
	std::map<size_t, BMIPhreeqcRM*>::iterator it = BMIPhreeqcRM::Instances.find(size_t(id));
	if (it != BMIPhreeqcRM::Instances.end())
	{
		delete (*it).second;
		retval = IRM_OK;
	}
	return retval;
}
/* ---------------------------------------------------------------------- */
BMIPhreeqcRM*
BMIPhreeqcRM::GetInstance(int id)
/* ---------------------------------------------------------------------- */
{
	std::map<size_t, BMIPhreeqcRM*>::iterator it = BMIPhreeqcRM::Instances.find(size_t(id));
	if (it != BMIPhreeqcRM::Instances.end())
	{
		return (*it).second;
	}
	return 0;
}
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM() :
PhreeqcRM(PhreeqcRM::default_nxyz, PhreeqcRM::default_data_for_parallel_processing, nullptr, true)
{
	this->language = "cpp";
}
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads, nullptr, true) 
{
	this->language = "cpp";
}
// Destructor
BMIPhreeqcRM::~BMIPhreeqcRM()
{
}
void BMIPhreeqcRM::Construct(PhreeqcRM::Initializer i)
{
	this->PhreeqcRM::Construct(i);
	//std::map<size_t, BMIPhreeqcRM*>::value_type instance(this->GetWorkers()[0]->Get_Index(), this);
	std::map<size_t, BMIPhreeqcRM*>::value_type instance(this->Index, this);
	BMIPhreeqcRM::Instances.insert(instance);
	this->var_man = new VarManager((PhreeqcRM*)this);
	//this->language = "cpp";
}

// Model control functions.
void BMIPhreeqcRM::Initialize(std::string config_file)
{
#ifdef USE_YAML
	YAML::Node yaml = YAML::LoadFile(config_file);
	std::string keyword;
	YAML::Node node;
	if (yaml["SetGridCellCount"].IsDefined())
	{
		this->initializer.nxyz_arg = yaml["SetGridCellCount"].as<int>();
	}
	if (yaml["ThreadCount"].IsDefined())
	{
		this->initializer.data_for_parallel_processing = yaml["ThreadCount"].as<int>();
	}
#endif

	this->Construct(this->initializer);

#ifdef USE_YAML
	this->InitializeYAML(config_file);
#endif
}
void BMIPhreeqcRM::Update()
{
	this->RunCells();
	this->SetTime(this->GetTime() + this->GetTimeStep());
	this->UpdateVariables();
}
void BMIPhreeqcRM::UpdateVariables()
{
	this->var_man->task = VarManager::VAR_TASKS::Update;
	std::set<RMVARS>& UpdateSet = this->var_man->UpdateSet;
	for (auto it = this->var_man->UpdateSet.begin(); it != UpdateSet.end(); it++)
	{
		VarManager::VarFunction fn = this->var_man->GetFn(*it);
		//		((*this).*f)(); 
		// ((*this->var_man).*fn)();
		((*this->var_man).*fn)();
	}
}
void BMIPhreeqcRM::UpdateUntil(double time)
{
	double time_step = time - this->GetTime();
	if (time_step >= 0)
	{
		this->SetTimeStep(time_step);
		this->RunCells();
		this->SetTime(time);
		this->UpdateVariables();
	}
}
void BMIPhreeqcRM::Finalize()
{
	this->CloseFiles();
}
int BMIPhreeqcRM::GetInputItemCount()
{
	int count = 0;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			//((*this->var_man).*bv.GetFn())();
			((*this->var_man).*bv.GetFn())();
		}
		if (bv.GetHasSetter())
		{
			count++;
		}
	}
	return count;
}
int BMIPhreeqcRM::GetOutputItemCount()
{
	int count = 0;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}	
		if (bv.GetHasGetter())
		{
			count++;
		}
	}
	return count;
}
int BMIPhreeqcRM::GetPointableItemCount()
{
	int count = 0;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		if (bv.GetHasPtr())
		{
			count++;
		}
	}
	return count;
}
std::vector<std::string> BMIPhreeqcRM::GetInputVarNames()
{
	std::vector <std::string> names;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		if (bv.GetHasSetter()) 
		{
			names.push_back(bv.GetName());
		}
	}
	return names;
}
std::vector<std::string>  BMIPhreeqcRM::GetOutputVarNames()
{ 
	std::vector <std::string> names;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		if (bv.GetHasGetter())
		{
			names.push_back(bv.GetName());
		}
	}
	return names;
}
std::vector<std::string> BMIPhreeqcRM::GetPointableVarNames()
{
	std::vector <std::string> names;
	for (auto it = this->var_man->VariantMap.begin();
		it != this->var_man->VariantMap.end(); it++)
	{
		BMIVariant& bv = it->second;
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		if (it->first == RMVARS::InputVarNames)
		{
			continue;
		}
		if (it->first == RMVARS::OutputVarNames)
		{
			continue;
		}
		if (bv.GetHasPtr())
		{
			names.push_back(bv.GetName());
		}
	}
	return names;
}
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		if (this->language == "cpp")
		{
			return bv.GetCType();
		}
		else if (this->language == "F90")
		{
			return bv.GetFType();
		}
		else if (this->language == "Py")
		{
			return bv.GetPType();
		}
	}
	assert(false);
	return "Unknown language.";
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		return bv.GetUnits();
	}
	assert(false);
	return "";
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		return bv.GetItemsize();
	}
	assert(false);
	return 0;
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		return bv.GetNbytes();
	}
	assert(false);
	return 0;
}
double BMIPhreeqcRM::GetCurrentTime()
{
	return this->GetTime();
}
double BMIPhreeqcRM::GetStartTime()
{
	return this->GetTime();
}
double BMIPhreeqcRM::GetEndTime()
{
	return this->GetTime() + this->GetTimeStep();
}
//double BMIPhreeqcRM::GetTimeStep()
//{
//	return this->GetTimeStep();
//}
void BMIPhreeqcRM::GetValue(const std::string name, void* dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		int Nbytes = this->var_man->VarExchange.GetNbytes();
		int dim = this->var_man->VarExchange.GetDim();
		if (this->var_man->VarExchange.GetCType() == "bool" && dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetBVarPtr(), Nbytes);
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "int" && dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetIVarPtr(), Nbytes);
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "double" && dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetDVarPtr(), Nbytes);
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "std::vector<std::string>")
		{
			int itemsize = this->GetVarItemsize(name);
			//int nbytes = this->GetVarNbytes(name);
			std::stringstream all;
			for (size_t i = 0; i < this->var_man->VarExchange.GetStringVectorRef().size(); i++)
			{
				all << std::left << std::setfill(' ') << std::setw(itemsize) << this->var_man->VarExchange.GetStringVectorRef()[i];
			}
			memcpy( dest, all.str().data(), all.str().size());
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "std::string" && dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetStringRef().data(), Nbytes);
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "double" && dim > 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetDoubleVectorPtr(), Nbytes);
			return;
		}
		if (this->var_man->VarExchange.GetCType() == "int" && dim > 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetIntVectorPtr(), Nbytes);
			return;
		}
	}
	std::ostringstream oss;
	oss << "BMI GetValue void* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "bool");
		dest = this->var_man->VarExchange.GetBVar();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool* dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "bool");
		int dim = this->var_man->VarExchange.GetDim();
		int nbytes = this->var_man->VarExchange.GetNbytes();
		if (dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetBVarPtr(), nbytes);
			return;
		}
	}
	std::ostringstream oss;
	oss << "BMI GetValue bool* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "double");
		dest = this->var_man->VarExchange.GetDVar();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double* dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "double");
		int dim = this->var_man->VarExchange.GetDim();
		int nbytes = this->var_man->VarExchange.GetNbytes();
		if (dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetDVarPtr(), nbytes);
			return;
		}
		else if (dim > 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetDoubleVectorPtr(), nbytes);
			return;
		}
	}
	std::ostringstream oss;
	oss << "BMI GetValue double* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		dest = this->var_man->VarExchange.GetIVar();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int* dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "int");
		int dim = this->var_man->VarExchange.GetDim();
		int nbytes = this->var_man->VarExchange.GetNbytes();
		if (dim == 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetIVarPtr(), nbytes);
			return;
		}
		else if (dim > 1)
		{
			memcpy(dest, this->var_man->VarExchange.GetIntVectorPtr(), nbytes);
			return;
		}
		std::ostringstream oss;
		oss << "BMI GetValue int* failed for variable " << name << std::endl;
		this->ErrorMessage(oss.str(), true);
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::string& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "std::string");
		dest = this->var_man->VarExchange.GetStringVar();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<double>& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		assert(this->var_man->VarExchange.GetCType() == "double");
		dest = this->var_man->VarExchange.GetDoubleVectorRef();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<int>& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		dest = this->var_man->VarExchange.GetIntVectorRef();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<std::string>& dest)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->task = VarManager::VAR_TASKS::GetVar;
		((*this->var_man).*bv.GetFn())();
		dest = this->var_man->VarExchange.GetStringVectorRef();
		return;
	}
	assert(false);
	return;
}
void* BMIPhreeqcRM::GetValuePtr(const std::string name)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (bv.GetVoidPtr() == NULL)
		{
			this->var_man->task = VarManager::VAR_TASKS::GetPtr;
			((*this->var_man).*bv.GetFn())();
		}
		return bv.GetVoidPtr();
	}
	assert(false);
	return NULL;
}
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store the variable in var_man->VarExchange
		int Nbytes = bv.GetNbytes();
		int itemsize = bv.GetItemsize();
		int dim = Nbytes / itemsize;
		if (bv.GetCType() == "bool" && dim == 1)
		{
			memcpy(this->var_man->VarExchange.GetBVarPtr(), src, Nbytes);
		}
		else if (bv.GetCType() == "int" && dim == 1)
		{
			memcpy(this->var_man->VarExchange.GetIVarPtr(), src, Nbytes);
		}
		else if (bv.GetCType() == "double" && dim == 1)
		{
			memcpy(this->var_man->VarExchange.GetDVarPtr(), src, Nbytes);
		}
		else if (bv.GetCType() == "std::string")
		{
			this->var_man->VarExchange.GetStringRef() = (char*)src;
		}
		else if (bv.GetCType() == "double" && dim > 1)
		{
			this->var_man->VarExchange.GetDoubleVectorRef().resize(dim);
			memcpy(this->var_man->VarExchange.GetDoubleVectorPtr(), src, Nbytes);
		}
		else if (bv.GetCType() == "int" && dim > 1)
		{
			this->var_man->VarExchange.GetIntVectorRef().resize(dim);
			memcpy(this->var_man->VarExchange.GetIntVectorPtr(), src, Nbytes);
		}
		//if (this->var_man->VarExchange.GetType() == "StringVector")
		//{
		//	// Don't think this is possible
		//	//int itemsize = this->GetVarItemsize(name);
		//	//int nbytes = this->GetVarNbytes(name);
		//	//std::stringstream all;
		//	//for (size_t i = 0; i < this->var_man->VarExchange.GetStringVectorRef().size(); i++)
		//	//{
		//	//	all << std::left << std::setfill(' ') << std::setw(itemsize) << this->var_man->VarExchange.GetStringVectorRef()[i];
		//	//}
		//	//memcpy( src, all.str().size());
		//}
		else
		{
			std::ostringstream oss;
			oss << "BMI failed in SetValue void* for variable " << name << std::endl;
			this->ErrorMessage(oss.str(), true);
			throw PhreeqcRMStop();
		}
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, bool src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.SetCType("bool");
		this->var_man->VarExchange.SetBVar(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, char* src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.SetStringVar(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, double src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.SetDVar(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, int src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.SetIVar(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, const std::string src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.SetStringVar(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<double> src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}	
		// Check dimension
		int dim = bv.GetDim(); 
		if (dim != src.size())
		{
			std::ostringstream oss;
			oss << "Dimension error in SetValue: " << name;
			this->ErrorMessage(oss.str());
			return;
		}
		// Store in var_man->VarExchange
		this->var_man->VarExchange.GetDoubleVectorRef().resize(bv.GetDim());
		this->var_man->VarExchange.SetDoubleVector(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<int> src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}	// Store in var_man->VarExchange
		this->var_man->VarExchange.GetIntVectorRef().resize(bv.GetDim());
		this->var_man->VarExchange.SetIntVector(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<std::string> src)
{
	RMVARS v_enum = this->var_man->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (!bv.GetInitialized())
		{
			this->var_man->task = VarManager::VAR_TASKS::Info;
			((*this->var_man).*bv.GetFn())();
		}
		this->var_man->VarExchange.SetCType("std::vector<std::string>");
		this->var_man->VarExchange.SetStringVector(src);
		// Set the variable
		this->var_man->task = VarManager::VAR_TASKS::SetVar;
		((*this->var_man).*bv.GetFn())();
		return;
	}
	assert(false);
	return;
}
//BMIPhreeqcRM::VarFunction BMIPhreeqcRM::GetFn(const std::string name)
//{
//	this->var_man->VarExchange.Clear();
//	std::string name_lc = name;
//	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
//	auto it = varfn_map.find(name_lc);
//	if (it == varfn_map.end())
//	{
//		std::ostringstream oss;
//		oss << "Unknown variable: " << name;
//		this->ErrorMessage(oss.str());
//		return NULL;
//	}
//
//	VarManager::VAR_TASKS task_save = this->var_man->task;
//	this->var_man->task = VarManager::VAR_TASKS::Info;
//	it->second(this);
//	this->var_man->task = task_save;
//	if (this->var_man->VarExchange.GetNotImplementedRef())
//	{
//		std::ostringstream oss;
//		oss << "Not implemented for variable: " << name;
//		this->ErrorMessage(oss.str());
//		return NULL;
//	}
//	if (this->var_man->task == VarManager::VAR_TASKS::GetVar)
//	{
//		if (!this->var_man->VarExchange.GetHasGetter() || this->var_man->VarExchange.GetNotImplementedRef())
//		{
//			std::ostringstream oss;
//			oss << "Cannot get variable: " << name;
//			this->ErrorMessage(oss.str());
//			return NULL;
//		}
//	}
//	if (this->var_man->task == VarManager::VAR_TASKS::SetVar)
//	{
//		if (!this->var_man->VarExchange.GetHasSetter() || this->var_man->VarExchange.GetNotImplementedRef())
//		{
//			std::ostringstream oss;
//			oss << "Cannot set variable: " << name;
//			this->ErrorMessage(oss.str());
//			return NULL;
//		}
//	}
//	if (this->var_man->task == VarManager::VAR_TASKS::GetPtr)
//	{
//		if (this->var_man->VarExchange.GetNotImplementedRef())
//		{
//			std::ostringstream oss;
//			oss << "Cannot get a pointer to variable: " << name;
//			this->ErrorMessage(oss.str());
//			return NULL;
//		}
//	}
//	return it->second;
//}

int BMIPhreeqcRM::GetGridRank(const int grid)
{
	if (grid == 0)
	{
		return 1;
	}
	return 0;
};

int BMIPhreeqcRM::GetGridSize(const int grid)
{
	if (grid == 0)
	{
		return this->GetGridCellCount();
	}
	return 0;
};

std::string BMIPhreeqcRM::GetGridType(const int grid)
{
	if (grid == 0)
	{
		return "points";
	}
	return "Undefined grid identifier";
};


//////////////////

IRM_RESULT BMIPhreeqcRM::BMIGenerateSelectedOutput()
{
	IRM_RESULT return_value = IRM_OK;
	int line_no = 10;
	int Itemsize = (int)sizeof(double);
	int Nbytes = this->GetGridCellCount() * Itemsize;
	std::ostringstream headings;
	std::ostringstream code;
	BMISelectedOutputVars.clear();
	auto it = BMISelecteOutputDefs.begin();
	for (; it != BMISelecteOutputDefs.end(); it++)
	{
		if (it->first == "output_aqueous_ph")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_ph", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_ph";
				BMIVariant bv(name, "-", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH -LA('H+')" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_pe")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_pe", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_pe";
				BMIVariant bv(name, "-", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH -LA('e-')" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_alkalinity")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_alkalinity", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_alkalinity";
				BMIVariant bv(name, "eq kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH ALK" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_ionic_strength")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_ionic_strength", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_ionic_strength";
				BMIVariant bv(name, "mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH MU" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_water_mass")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_water_mass", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_water_mass";
				BMIVariant bv(name, "kg", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH TOT('water')" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_charge_balance")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_charge_balance", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_charge_balance";
				BMIVariant bv(name, "eq kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH CHARGE_BALANCE / TOT('water')" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_percent_error")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				ErrorMessage("Unknown input for output_aqueous_percent_error", true);
				return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "aqueous_percent_error";
				BMIVariant bv(name, "-", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH PERCENT_ERROR" << std::endl;
				line_no += 10;
				break;
			}
			}
		}
		else if (it->first == "output_aqueous_total_molalities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
				case 0:
					continue;
				case 1:
				{
					item_set = ElementRedoxSet;
					break;
				}
				case 2:
				{
					item_set = tokenize(it->second);
					break;
				}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "aqueous_total_molality_" + *item_it;
				BMIVariant bv(name, "mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH TOT('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_aqueous_molalities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < species_names.size(); i++)
				{
					item_set.insert(species_names[i]);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "aqueous_species_log_molality_" + *item_it;
				BMIVariant bv(name, "log mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH LM('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_aqueous_activities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < species_names.size(); i++)
				{
					item_set.insert(species_names[i]);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "aqueous_species_log_activity_" + *item_it;
				BMIVariant bv(name, "log -", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH LA('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_exchange_molalities")
		{
			std::set<std::string> item_set;
			std::map<std::string, std::string> item_map;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < ExchangeSpeciesNamesList.size(); i++)
				{
					item_set.insert(ExchangeSpeciesNamesList[i]);
					item_map[ExchangeSpeciesNamesList[i]] = ExchangeNamesList[i];
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "exchange_species_log_molality_" + *item_it;
				if (item_map.size() > 0)
				{
					std::string xname = item_map[*item_it];
					name = "exchange_" + xname + "_species_log_molality_" + *item_it;
				}
				BMIVariant bv(name, "log mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH LM('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_surface_molalities")
		{
			std::set<std::string> item_set;
			std::map<std::string, std::string> item_map;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < SurfaceNamesList.size(); i++)
				{
					item_set.insert(SurfaceNamesList[i]);
					item_map[SurfaceNamesList[i]] = SurfaceTypesList[i];
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "surface_species_log_molality_" + *item_it;
				if (item_map.size() > 0)
				{
					std::string type = item_map[*item_it];
					name = "surface_" + type + "_species_log_molality_" + *item_it;
				}
				BMIVariant bv(name, "log mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH LM('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_equilibrium_phases")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < EquilibriumPhasesList.size(); i++)
				{
					item_set.insert(EquilibriumPhasesList[i]);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				{
					std::string name = "equilibrium_phases_moles_" + *item_it;
					BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH EQUI('" << name << "')\n";
					line_no += 10;
				}
				{
					std::string name = "equilibrium_phases_delta_moles_" + *item_it;
					BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH EQUI_DELTA('" << name << "')\n";
					line_no += 10;
				}
			}
		}
		else if (it->first == "output_saturation_indices")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < SINamesList.size(); i++)
				{
					item_set.insert(SINamesList[i]);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "aqueous_saturation_index_" + *item_it;
				BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH SI('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_gases")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
				case 0:
					continue;
				case 1:
				{
					for (size_t i = 0; i < GasComponentsList.size(); i++)
					{
						item_set.insert(GasComponentsList[i]);
					}
					break;
				}
				case 2:
				{
					item_set = tokenize(it->second);
					break;
				}
			}
			{
				std::string name = "gas_phase_volume";
				BMIVariant bv(name, "L", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH SYS('gas') * GAS_VM\n";
				line_no += 10;
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				{
					std::string name = "gas_phase_moles_" + *item_it;
					BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH GAS('" << name << "')\n";
					line_no += 10;
				}
				{
					std::string name = "gas_phase_pressure_" + *item_it;
					BMIVariant bv(name, "atm", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH PR_P('" << name << "')\n";
					line_no += 10;
				}
				{
					std::string name = "gas_phase_phi_" + *item_it;
					BMIVariant bv(name, "atm-1", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH PR_PHI('" << name << "')\n";
					line_no += 10;
				}
			}
		}
		else if (it->first == "output_kinetic_reactants")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < KineticReactionsList.size(); i++)
				{
					item_set.insert(KineticReactionsList[i]);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				{
					std::string name = "kinetic_reaction_moles_" + *item_it;
					BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH KIN('" << name << "')\n";
					line_no += 10;
				}
				{
					std::string name = "kinetic_reaction_delta_moles_" + *item_it;
					BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
					bv.SetColumn((int)BMISelectedOutputVars.size());
					BMISelectedOutputVars[name] = bv;
					headings << name << "\t";
					code << line_no << " PUNCH KIN_DELTA('" << name << "')\n";
					line_no += 10;
				}
			}
		}
		else if (it->first == "output_solid_solutions")
		{
			std::set<std::string> item_set;
			std::map<std::string, std::string> item_map;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < SolidSolutionComponentsList.size(); i++)
				{
					item_set.insert(SolidSolutionComponentsList[i]);
					item_map[SolidSolutionComponentsList[i]] = SolidSolutionNamesList[i];
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "solid_solution_moles_" + *item_it;
				if (item_map.size() > 0)
				{
					std::string xname = item_map[*item_it];
					name = "solid_solution_" + xname + "_moles_" + *item_it;
				}
				BMIVariant bv(name, "mol", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH S_S('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_calculate_values")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				auto it = this->workers[0]->Get_PhreeqcPtr()->GetCalculateValueMap().begin();
				for (; it != this->workers[0]->Get_PhreeqcPtr()->GetCalculateValueMap().end(); it++)
				{
					item_set.insert(it->first);
				}
				break;
			}
			case 2:
			{
				item_set = tokenize(it->second);
				break;
			}
			}
			auto item_it = item_set.begin();
			for (; item_it != item_set.end(); item_it++)
			{
				std::string name = "calculate_value_" + *item_it;
				BMIVariant bv(name, "unknown", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH CALC_VALUE('" << name << "')\n";
				line_no += 10;
			}
		}
		else
		{
			std::ostringstream oss;
			oss << "Unknown output request " << it->first;
			this->ErrorMessage(oss.str(), true);
			throw PhreeqcStop();
		}
	}
	std::ostringstream data_block;
	data_block << "SELECTED_OUTPUT 777777777; USER_PUNCH 777777777;" << std::endl;
	data_block << headings.str() << std::endl;
	data_block << code.str() << std::endl;
	this->RunString(true, false, false, data_block.str());
	BMISelecteOutputDefs.clear();
	return return_value;
}
int BMIPhreeqcRM::BMICheckSelectedOutputDef(bool tf_only, std::string& def)
{
	std::string def_lc = def;
	std::transform(def_lc.begin(), def_lc.end(), def_lc.begin(),
		tolower);
	if (def_lc == "false")
	{
		return 0;
	}
	if (def_lc == "true")
	{
		return 1;
	}
	if (tf_only)
	{
		return -1;
	}
	return 2;
}
std::set<std::string> BMIPhreeqcRM::tokenize(const std::string& def_in)
{
	std::set<std::string> item_set;
	std::string def, a_token;
	def = def_in;
	for (size_t i = 0; i < def.size(); i++)
	{
		// check for c1 and replace
		if (def[i] == '\t')
			def[i] = ' ';
		if (def[i] == ',')
			def[i] = ' ';
	}
	std::stringstream ss(def);
	while (std::getline(ss, a_token, ' '))
	{
		a_token = trim(a_token);
		if (a_token.size() > 0)
		{
			item_set.insert(a_token);
		}
	}
	return item_set;
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int
PhreeqcRM::FindComponents(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *   Counts components in any defined solution, gas_phase, exchanger,
	 *   surface, or pure_phase_assemblage
	 *
	 *   Returns
	 *		n_comp, which is total, including H, O, elements, and Charge
	 *      names, which contains character strings with names of components
	 */
	bool clear = false;
	this->phreeqcrm_error_string.clear();
	try
	{
#ifdef USE_MPI
		if (this->mpi_myself == 0)
		{
			int method = METHOD_FINDCOMPONENTS;
			MPI_Bcast(&method, 1, MPI_INT, 0, phreeqcrm_comm);
		}
#endif
		// Always include H, O, Charge

		std::set<std::string> component_set;

		size_t fixed_components = 3;
		if (this->component_h2o)
			fixed_components = 4;

		// save old components
		for (size_t i = fixed_components; i < this->components.size(); i++)
		{
			component_set.insert(this->components[i]);
		}

		// Get other components
		IPhreeqcPhast* phast_iphreeqc_worker = this->GetWorkers()[this->nthreads];
		size_t count_components = phast_iphreeqc_worker->GetComponentCount();

		size_t i;
		for (i = 0; i < count_components; i++)
		{
			std::string comp(phast_iphreeqc_worker->GetComponent((int)i));
			assert(comp != "H");
			assert(comp != "O");
			assert(comp != "Charge");
			assert(comp != "charge");

			component_set.insert(comp);
		}
		// clear and refill components in vector
		this->components.clear();

		// Always include H, O, Charge
		if (this->component_h2o)
			this->components.push_back("H2O");
		this->components.push_back("H");
		this->components.push_back("O");
		this->components.push_back("Charge");
		for (std::set<std::string>::iterator it = component_set.begin(); it != component_set.end(); it++)
		{
			this->components.push_back(*it);
		}
		// Calculate gfw for components
		this->gfw.clear();
		for (i = 0; i < components.size(); i++)
		{
			if (components[i] == "Charge")
			{
				this->gfw.push_back(1.0);
			}
			else
			{
				this->gfw.push_back(phast_iphreeqc_worker->Get_gfw(components[i].c_str()));
			}
		}
		// Get list of species
		if (this->species_save_on)
		{
			phast_iphreeqc_worker->PhreeqcPtr->save_species = true;
		}
		// Make lists regardless of species_save_on
		{
			int next = phast_iphreeqc_worker->PhreeqcPtr->next_user_number(Keywords::KEY_SOLUTION);
			{
				std::ostringstream in;
				in << "SOLUTION " << next << "\n";
				for (i = 0; i < components.size(); i++)
				{
					if (components[i] == "H") continue;
					if (components[i] == "O") continue;
					if (components[i] == "H2O") continue;
					if (components[i] == "Charge") continue;
					in << components[i] << " 1e-6\n";
				}
				int status = phast_iphreeqc_worker->RunString(in.str().c_str());
				if (status != 0)
				{
					this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
					throw PhreeqcRMStop();
				}
			}
			species_names.clear();
			species_z.clear();
			s_num2rm_species_num.clear();
			species_stoichiometry.clear();
			for (int i = 0; i < (int)phast_iphreeqc_worker->PhreeqcPtr->s_x.size(); i++)
			{
				species_names.push_back(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->name);
				species_z.push_back(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->z);
				species_d_25.push_back(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->dw);
				s_num2rm_species_num[phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->number] = i;
				cxxNameDouble nd(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->next_elt);
				nd.add("Charge", phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->z);
				species_stoichiometry.push_back(nd);
			}
			ElementRedoxSet.clear();
			for (size_t i = 0; i < phast_iphreeqc_worker->PhreeqcPtr->master.size(); i++)
			{
				if (phast_iphreeqc_worker->PhreeqcPtr->master[i]->in != FALSE)
				{
					std::string e = phast_iphreeqc_worker->PhreeqcPtr->master[i]->elt->name;
					if (e != "E")
					{
						ElementRedoxSet.insert(e);
					}
					e = phast_iphreeqc_worker->PhreeqcPtr->master[i]->elt->primary->elt->name;
					if (e != "E")
					{
						ElementRedoxSet.insert(e);
					}
				}
			}
			for (int i = 0; i < (int)phast_iphreeqc_worker->PhreeqcPtr->phases.size(); i++)
			{
				if (phast_iphreeqc_worker->PhreeqcPtr->phases[i]->in == TRUE)
				{
					SINamesList.push_back(phast_iphreeqc_worker->PhreeqcPtr->phases[i]->name);
				}
			}
			{
				std::ostringstream in;
				in << "DELETE; -solution " << next << "\n";
				int status = phast_iphreeqc_worker->RunString(in.str().c_str());
				if (status != 0)
				{
					this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
					throw PhreeqcRMStop();
				}
			}
		}
		// Make all lists
		{
			// Make set for surfaces
			std::set<std::string> surface_types_set;
			std::map<std::string, std::string> surface_names_map;
			if (!clear)
			{
				for (size_t ii = 0; ii < this->SurfaceSpeciesNamesList.size(); ii++)
				{
					surface_types_set.insert(this->SurfaceTypesList[ii]);
					surface_names_map[this->SurfaceTypesList[ii]] = this->SurfaceNamesList[ii];
				}
			}
			// add new surface types 
			{
				const std::list<std::string>& surftype = phast_iphreeqc_worker->GetSurfaceTypeList();
				const std::list<std::string>& surfnames = phast_iphreeqc_worker->GetSurfaceNamesList();
				{
					std::list<std::string>::const_iterator surftype_it = surftype.begin();
					std::list<std::string>::const_iterator surfnames_it = surfnames.begin();
					for (; surftype_it != surftype.end(); surftype_it++)
					{
						surface_types_set.insert(*surftype_it);
						surface_names_map[*surftype_it] = *surfnames_it++;
					}
				}
			}
			// make set for exchange
			std::set<std::string> ex_set;
			if (!clear)
			{
				for (size_t ii = 0; ii < this->ExchangeNamesList.size(); ii++)
				{
					ex_set.insert(this->ExchangeNamesList[ii]);
				}
			}
			// add new exchange sites
			{
				const std::list<std::string>& ex = phast_iphreeqc_worker->GetExchangeNamesList();
				{
					std::list<std::string>::const_iterator ex_it = ex.begin();
					for (; ex_it != ex.end(); ex_it++)
					{
						ex_set.insert(*ex_it);
					}
				}
			}
			// write solution
			int next = phast_iphreeqc_worker->PhreeqcPtr->next_user_number(Keywords::KEY_SOLUTION);
			if (ex_set.size() > 0 || surface_types_set.size() > 0)
			{
				std::ostringstream in;
				in << "SOLUTION " << next << "\n";
				for (i = 0; i < components.size(); i++)
				{
					if (components[i] == "H") continue;
					if (components[i] == "O") continue;
					if (components[i] == "H2O") continue;
					if (components[i] == "Charge") continue;
					in << components[i] << " 1e-6\n";
				}
				in << "END\n";
				int status = phast_iphreeqc_worker->RunString(in.str().c_str());
				if (status != 0)
				{
					this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
					throw PhreeqcRMStop();
				}
			}
			// write surface and save vectors
			int next_surf = phast_iphreeqc_worker->PhreeqcPtr->next_user_number(Keywords::KEY_SURFACE);
			this->SurfaceSpeciesNamesList.clear();
			this->SurfaceTypesList.clear();
			this->SurfaceNamesList.clear();
			if (surface_types_set.size() > 0)
			{
				std::set<std::string>::iterator cit = surface_types_set.begin();
				for (; cit != surface_types_set.end(); cit++)
				{
					std::ostringstream in;
					in << "SURFACE " << next_surf << "\n";
					in << "  -eq " << next << "\n";
					in << "  " << *cit << "  0.001  1   1\n";
					int status = phast_iphreeqc_worker->RunString(in.str().c_str());
					if (status != 0)
					{
						this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
						throw PhreeqcRMStop();
					}
					// fill surface vectors
					for (int i = 0; i < (int)phast_iphreeqc_worker->PhreeqcPtr->s_x.size(); i++)
					{
						if (phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->type == SURF)
						{
							this->SurfaceSpeciesNamesList.push_back(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->name);
							this->SurfaceTypesList.push_back(*cit);
							this->SurfaceNamesList.push_back(surface_names_map[*cit]);
						}
					}
					{
						std::ostringstream in1;
						in1 << "DELETE -surface " << next_surf << "\n";
						status = phast_iphreeqc_worker->RunString(in1.str().c_str());
						if (status != 0)
						{
							this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
							throw PhreeqcRMStop();
						}
					}
				}
			}
			// Exchange species
			int next_ex = phast_iphreeqc_worker->PhreeqcPtr->next_user_number(Keywords::KEY_EXCHANGE);
			this->ExchangeSpeciesNamesList.clear();
			this->ExchangeNamesList.clear();
			if (ex_set.size() > 0)
			{
				std::set<std::string>::iterator cit = ex_set.begin();
				for (; cit != ex_set.end(); cit++)
				{
					std::ostringstream in;
					in << "EXCHANGE " << next_ex << "\n";
					in << "  -eq " << next << "\n";
					in << "  " << *cit << "  0.001\n";
					int status = phast_iphreeqc_worker->RunString(in.str().c_str());
					if (status != 0)
					{
						this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
						throw PhreeqcRMStop();
					}
					for (int i = 0; i < (int)phast_iphreeqc_worker->PhreeqcPtr->s_x.size(); i++)
					{
						if (phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->type == EX)
						{
							this->ExchangeSpeciesNamesList.push_back(phast_iphreeqc_worker->PhreeqcPtr->s_x[i]->name);
							this->ExchangeNamesList.push_back(*cit);
						}
					}
					{
						std::ostringstream in1;
						in1 << "DELETE -exchange " << next_ex << "\n";
						status = phast_iphreeqc_worker->RunString(in1.str().c_str());
						if (status != 0)
						{
							this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
							throw PhreeqcRMStop();
						}
					}
				}
			}
			if (ex_set.size() > 0 || surface_types_set.size() > 0)
			{
				std::ostringstream in;
				in << "DELETE; -solution " << next << "\n";
				int status = phast_iphreeqc_worker->RunString(in.str().c_str());
				if (status != 0)
				{
					this->ErrorMessage(phast_iphreeqc_worker->GetErrorString());
					throw PhreeqcRMStop();
				}
			}
			// equilibrium_phases
			{
				// Move to set
				std::set<std::string> eq_set;
				if (!clear)
				{
					for (size_t i = 0; i < this->EquilibriumPhasesList.size(); i++)
					{
						eq_set.insert(this->EquilibriumPhasesList[i]);
					}
				}
				// add new equilibrium phases to set
				std::list<std::string>::const_iterator cit = phast_iphreeqc_worker->GetEquilibriumPhasesList().begin();
				for (; cit != phast_iphreeqc_worker->GetEquilibriumPhasesList().end(); cit++)
				{
					eq_set.insert(*cit);
				}
				// move set to vector
				this->EquilibriumPhasesList.clear();
				std::set<std::string>::const_iterator eqit = eq_set.begin();
				for (; eqit != eq_set.end(); eqit++)
				{
					this->EquilibriumPhasesList.push_back(*eqit);
				}
			}
			// gas phase components
			{
				// Move to set
				std::set<std::string> g_set;
				if (!clear)
				{
					for (size_t i = 0; i < this->GasComponentsList.size(); i++)
					{
						g_set.insert(this->GasComponentsList[i]);
					}
				}
				// add new gas components to set
				std::list<std::string>::const_iterator cit = phast_iphreeqc_worker->GetGasComponentsList().begin();
				for (; cit != phast_iphreeqc_worker->GetGasComponentsList().end(); cit++)
				{
					g_set.insert(*cit);
				}
				// move set to vector
				this->GasComponentsList.clear();
				std::set<std::string>::const_iterator git = g_set.begin();
				for (; git != g_set.end(); git++)
				{
					this->GasComponentsList.push_back(*git);
				}
			}
			// Kinetics
			{
				// Move to set
				std::set<std::string> k_set;
				if (!clear)
				{
					for (size_t i = 0; i < this->KineticReactionsList.size(); i++)
					{
						k_set.insert(this->KineticReactionsList[i]);
					}
				}
				// add new kinetic reactions to set
				std::list<std::string>::const_iterator cit = phast_iphreeqc_worker->GetKineticReactionsList().begin();
				for (; cit != phast_iphreeqc_worker->GetKineticReactionsList().end(); cit++)
				{
					k_set.insert(*cit);
				}
				// move set to vector
				this->KineticReactionsList.clear();
				std::set<std::string>::const_iterator kit = k_set.begin();
				for (; kit != k_set.end(); kit++)
				{
					this->KineticReactionsList.push_back(*kit);
				}
			}
			// Solid solutions
			{
				// move existing component names to set and solid solution names to map
				std::set<std::string> sscomp_set;
				std::map<std::string, std::string> ssnames_map;
				if (!clear)
				{
					for (size_t i = 0; i < this->SolidSolutionComponentsList.size(); i++)
					{
						sscomp_set.insert(this->SolidSolutionComponentsList[i]);
						ssnames_map[this->SolidSolutionComponentsList[i]] = this->SolidSolutionNamesList[i];
					}
				}
				// add new component names set and solid solution names to map
				std::list<std::string>::const_iterator cit = phast_iphreeqc_worker->GetSolidSolutionComponentsList().begin();
				std::list<std::string>::const_iterator nit = phast_iphreeqc_worker->GetSolidSolutionNamesList().begin();
				for (; cit != phast_iphreeqc_worker->GetSolidSolutionComponentsList().end(); cit++)
				{
					if (sscomp_set.find(*cit) == sscomp_set.end())
					{
						sscomp_set.insert(*cit);
						ssnames_map[*cit] = *nit++;
					}
				}
				// Move from set and map to vectors
				this->SolidSolutionComponentsList.clear();
				this->SolidSolutionNamesList.clear();
				std::set<std::string>::const_iterator set_it = sscomp_set.begin();
				for (; set_it != sscomp_set.end(); set_it++)
				{
					this->SolidSolutionComponentsList.push_back(*set_it);
					this->SolidSolutionNamesList.push_back(ssnames_map[*set_it]);
				}
			}
		}
	}
	catch (...)
	{
		return this->ReturnHandler(IRM_FAIL, "PhreeqcRM::FindComponents");
	}
	return (int)this->components.size();
}
#endif