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
	for (auto it = UpdateMap.begin(); it != UpdateMap.end(); it++)
	{
		std::string name_lc = *it;
		std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
		auto m_it = varfn_map.find(name_lc);
		if (m_it != varfn_map.end())
		{
			m_it->second(this);
		}
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
	this->var_man->task = VarManager::VAR_TASKS::no_op;
	int count = 0;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		if (it->first == "inputvarnames")
		{
			continue;
		}
		if (it->first == "outputvarnames")
		{
			continue;
		}
		//this->var_man->VarExchange.SetSet(false);
		it->second(this);
		//if (this->var_man->VarExchange.GetSet()) count++;
		if (this->var_man->VarExchange.GetHasSetter()) count++;
	}
	return count;
}
int BMIPhreeqcRM::GetOutputItemCount()
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	int count = 0;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		if (it->first == "inputvarnames")
		{
			count++;
			continue;
		}
		if (it->first == "outputvarnames")
		{
			count++;
			continue;
		}
		//this->var_man->VarExchange.SetGet(false);
		it->second(this);
		if (this->var_man->VarExchange.GetHasGetter()) count++;
	}
	return count;
}
int BMIPhreeqcRM::GetPointableItemCount()
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	int count = 0;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		it->second(this);
		if (this->var_man->VarExchange.GetHasPtr()) count++;
	}
	return count;
}
std::vector<std::string> BMIPhreeqcRM::GetInputVarNames()
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	std::vector<std::string> names;
	VarFunction_map::iterator it;
	for (it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		if (it->first == "inputvarnames")
		{
			continue;
		}
		if (it->first == "outputvarnames")
		{
			continue;
		}
		it->second(this);
		if (this->var_man->VarExchange.GetHasSetter()) 
		{
			names.push_back(this->var_man->VarExchange.GetName());
		}
	}
	return names;
}
std::vector<std::string>  BMIPhreeqcRM::GetOutputVarNames()
{ 
	this->var_man->task = VarManager::VAR_TASKS::Info;
	std::vector<std::string> names;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		if (it->first == "inputvarnames")
		{
			names.push_back("InputVarNames");
			continue;
		}
		if (it->first == "outputvarnames")
		{
			names.push_back("OutputVarNames");
			continue;
		}
		if (this->var_man->VarExchange.GetHasGetter())
		it->second(this);
		if (this->var_man->VarExchange.GetHasGetter())
		{
			names.push_back(this->var_man->VarExchange.GetName());
		}
	}
	return names;
}
std::vector<std::string> BMIPhreeqcRM::GetPointableVarNames()
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	std::vector<std::string> names;
	VarFunction_map::iterator it;
	for (it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		it->second(this);
		if (this->var_man->VarExchange.GetHasPtr())
		{
			names.push_back(this->var_man->VarExchange.GetName());
		}
	}
	return names;
}
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn(this);
	if (this->language == "cpp")
	{
		return this->var_man->VarExchange.GetCType();
	} 
	else if (this->language == "F90")
	{
		return this->var_man->VarExchange.GetFType();
	}
	else if (this->language == "Py")
	{
		return this->var_man->VarExchange.GetPType();
	}
	return "Unknown language.";
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn(this);
	return this->var_man->VarExchange.GetUnits();
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn(this);
	return this->var_man->VarExchange.GetItemsize();
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	this->var_man->task = VarManager::VAR_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn(this);
	return this->var_man->VarExchange.GetNbytes();
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
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	int Nbytes = this->var_man->VarExchange.GetNbytes();
	//int dim = this->var_man->VarExchange.GetDim();
	int dim = this->var_man->VarExchange.GetDim();
	if (this->var_man->VarExchange.GetCType() == "bool" && dim == 1)
	{
		//memcpy(dest, this->var_man->VarExchange.GetBVarPtr(), Nbytes);
		memcpy(dest, this->var_man->VarExchange.GetBVarPtr(), Nbytes);
		return;
	}
	if (this->var_man->VarExchange.GetCType() == "int" && dim == 1)
	{
		memcpy(dest, &this->var_man->VarExchange.GetIVarRef(), Nbytes);
		return;
	}
	if (this->var_man->VarExchange.GetCType() == "double" && dim == 1)
	{
		memcpy(dest, &this->var_man->VarExchange.GetDVarRef(), Nbytes);
		return;
	}
	//if (this->var_man->VarExchange.GetCType() == "std::string" && dim == 1)
	//{
	//	memcpy(dest, this->var_man->VarExchange.GetStringRef().data(), Nbytes);
	//}
	if (this->var_man->VarExchange.GetCType() == "double" && dim > 1)
	{
		//memcpy(dest, this->var_man->VarExchange.GetDoubleVectorRef().data(), Nbytes);
		memcpy(dest, this->var_man->VarExchange.GetDoubleVectorRef().data(), Nbytes);
		return;
	}
	if (this->var_man->VarExchange.GetCType() == "int" && dim > 1)
	{
		memcpy(dest, this->var_man->VarExchange.GetIntVectorRef().data(), Nbytes);
		return;
	}
	//if (this->var_man->VarExchange.GetCType() == "std::vector<std::string>")
	//{
	//	int itemsize = this->var_man->VarExchange.GetItemsize();
	//	std::stringstream all;
	//	for (size_t i = 0; i < this->var_man->VarExchange.GetStringVectorRef().size(); i++)
	//	{
	//		all << std::left << std::setfill(' ') << std::setw(itemsize) << this->var_man->VarExchange.GetStringVectorRef()[i];
	//	}
	//	memcpy(dest, all.str().c_str(), all.str().size());
	//}
	std::ostringstream oss;
	oss << "BMI GetValue void* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->var_man->VarExchange.GetCType() == "bool");
	dest = this->var_man->VarExchange.GetBVarRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool* dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->var_man->VarExchange.GetCType() == "bool");
	int dim = this->var_man->VarExchange.GetDim();
	int nbytes = this->var_man->VarExchange.GetNbytes();
	if (dim == 1)
	{
		memcpy(dest, this->var_man->VarExchange.GetBVarPtr(), nbytes);
		return;
	}
	//else if (dim > 1)
	//{
	//	memcpy(dest, this->var_man->VarExchange.BoolVector.data(), nbytes);
	//	return;
	//}
	std::ostringstream oss;
	oss << "BMI GetValue bool* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->var_man->VarExchange.GetCType() == "double");
	dest = this->var_man->VarExchange.GetDVarRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double* dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->var_man->VarExchange.GetCType() == "double");
	int dim = this->var_man->VarExchange.GetDim();
	int nbytes = this->var_man->VarExchange.GetNbytes();
	if (dim == 1)
	{
		//memcpy(dest, &this->var_man->VarExchange.GetDVarRef(), nbytes);
		memcpy(dest, &this->var_man->VarExchange.GetDVarRef(), nbytes);
		return;
	}
	else if (dim > 1)
	{
		memcpy(dest, this->var_man->VarExchange.GetDoubleVectorRef().data(), nbytes);
		return;
	}
	std::ostringstream oss;
	oss << "BMI GetValue double* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->var_man->VarExchange.GetIVarRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int* dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->var_man->VarExchange.GetCType() == "int");
	int dim = this->var_man->VarExchange.GetDim();
	int nbytes = this->var_man->VarExchange.GetNbytes();
	if (dim == 1)
	{
		memcpy(dest, &this->var_man->VarExchange.GetIVarRef(), nbytes);
		return;
	}
	else if (dim > 1)
	{
		memcpy(dest, this->var_man->VarExchange.GetIntVectorRef().data(), nbytes);
		return;
	}
	std::ostringstream oss;
	oss << "BMI GetValue int* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::string& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	//dest = this->var_man->VarExchange.GetStringRef();
	dest = this->var_man->VarExchange.GetStringRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<double>& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->var_man->VarExchange.GetDoubleVectorRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<int>& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->var_man->VarExchange.GetIntVectorRef();
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<std::string>& dest)
{
	this->var_man->task = VarManager::VAR_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	//dest = this->var_man->VarExchange.GetStringVectorRef();
	dest = this->var_man->VarExchange.GetStringVectorRef();
	return;
}
void* BMIPhreeqcRM::GetValuePtr(const std::string name)
{
	this->var_man->task = VarManager::VAR_TASKS::GetPtr;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return NULL;
	fn(this);
	return this->var_man->VarExchange.GetVoidPtr();
}
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store the variable in var_man->VarExchange
	int Nbytes = this->var_man->VarExchange.GetNbytes();
	int itemsize = this->var_man->VarExchange.GetItemsize();
	int dim = Nbytes / itemsize;
	if (this->var_man->VarExchange.GetCType() == "bool" && dim == 1)
	{
		memcpy(this->var_man->VarExchange.GetBVarPtr(), src, Nbytes);
	} 
	else if (this->var_man->VarExchange.GetCType() == "int" && dim == 1)
	{
		memcpy(&this->var_man->VarExchange.GetIVarRef(), src, Nbytes);
	} 
	else if (this->var_man->VarExchange.GetCType() == "double" && dim == 1)
	{
		memcpy(&this->var_man->VarExchange.GetDVarRef(), src, Nbytes);
	} 
	else if (this->var_man->VarExchange.GetCType() == "std::string")
	{ 
		this->var_man->VarExchange.GetStringRef() = (char*)src;
	}
	else if (this->var_man->VarExchange.GetCType() == "double" && dim > 1)
	{
		this->var_man->VarExchange.GetDoubleVectorRef().resize(dim);
		memcpy(this->var_man->VarExchange.GetDoubleVectorRef().data(), src, Nbytes);
	}
	else if (this->var_man->VarExchange.GetCType() == "int" && dim > 1)
	{
		//this->var_man->VarExchange.GetIntVectorRef().resize(dim);
		this->var_man->VarExchange.GetIntVectorRef().resize(dim);
		memcpy(&this->var_man->VarExchange.GetIntVectorRef(), src, Nbytes);
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
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, bool src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("bool");
	this->var_man->VarExchange.GetBVarRef() = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, char* src)
{
	std::string str = src;
	SetValue(name, str);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, double src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("double");
	this->var_man->VarExchange.GetDVarRef() = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, int src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("int");
	this->var_man->VarExchange.GetIVarRef() = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, const std::string src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("std::string");
	this->var_man->VarExchange.GetStringRef() = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<double> src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Check dimension
	int dim = this->GetVarNbytes(name) / this->GetVarItemsize(name);
	if (dim != src.size())
	{
		std::ostringstream oss;
		oss << "Dimension error in SetValue: " << name;
		this->ErrorMessage(oss.str());
		return;
	}
	// Store in var_man->VarExchange
	this->var_man->VarExchange.GetDoubleVectorRef() = src;
	// Set the variable
	this->var_man->VarExchange.SetCType("std::vector<double>");
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<int> src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("std::vector<int>");
	this->var_man->VarExchange.GetIntVectorRef() = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<std::string> src)
{
	this->var_man->task = VarManager::VAR_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	//// Check dimension
	//int dim = this->GetVarNbytes(name) / this->GetVarItemsize(name);
	//if (dim != src.size())
	//{
	//	std::ostringstream oss;
	//	oss << "Dimension error in SetValue: " << name;
	//	this->ErrorMessage(oss.str());
	//	return;
	//}
	// Store in var_man->VarExchange
	this->var_man->VarExchange.SetCType("std::vector<std::string>");
	this->var_man->VarExchange.GetStringVectorRef() = src;
	// Set the variable
	fn(this);
	return;
}
BMIPhreeqcRM::VarFunction BMIPhreeqcRM::GetFn(const std::string name)
{
	this->var_man->VarExchange.Clear();
	std::string name_lc = name;
	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
	auto it = varfn_map.find(name_lc);
	if (it == varfn_map.end())
	{
		std::ostringstream oss;
		oss << "Unknown variable: " << name;
		this->ErrorMessage(oss.str());
		return NULL;
	}

	VarManager::VAR_TASKS task_save = this->var_man->task;
	this->var_man->task = VarManager::VAR_TASKS::Info;
	it->second(this);
	this->var_man->task = task_save;
	if (this->var_man->VarExchange.GetNotImplementedRef())
	{
		std::ostringstream oss;
		oss << "Not implemented for variable: " << name;
		this->ErrorMessage(oss.str());
		return NULL;
	}
	if (this->var_man->task == VarManager::VAR_TASKS::GetVar)
	{
		if (!this->var_man->VarExchange.GetHasGetter() || this->var_man->VarExchange.GetNotImplementedRef())
		{
			std::ostringstream oss;
			oss << "Cannot get variable: " << name;
			this->ErrorMessage(oss.str());
			return NULL;
		}
	}
	if (this->var_man->task == VarManager::VAR_TASKS::SetVar)
	{
		if (!this->var_man->VarExchange.GetHasSetter() || this->var_man->VarExchange.GetNotImplementedRef())
		{
			std::ostringstream oss;
			oss << "Cannot set variable: " << name;
			this->ErrorMessage(oss.str());
			return NULL;
		}
	}
	if (this->var_man->task == VarManager::VAR_TASKS::GetPtr)
	{
		if (this->var_man->VarExchange.GetNotImplementedRef())
		{
			std::ostringstream oss;
			oss << "Cannot get a pointer to variable: " << name;
			this->ErrorMessage(oss.str());
			return NULL;
		}
	}
	return it->second;
}

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


