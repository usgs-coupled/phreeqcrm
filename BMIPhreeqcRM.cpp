#include "BMIPhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>

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
			n = (int)bmirm_ptr->GetWorkers()[0]->Get_Index();
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
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	std::map<size_t, BMIPhreeqcRM*>::value_type instance(this->GetWorkers()[0]->Get_Index(), this);
	BMIPhreeqcRM::Instances.insert(instance);

	this->var_man = new VarManager((PhreeqcRM*)this);
	this->language = "cpp";
}
// Model control functions.
void BMIPhreeqcRM::Initialize(std::string config_file)
{
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


