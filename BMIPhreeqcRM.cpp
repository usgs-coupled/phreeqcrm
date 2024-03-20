#include "BMIPhreeqcRM.h"
#include "BMIVariant.h"
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

//// static BMIPhreeqcRM methods
/* ---------------------------------------------------------------------- */
void
BMIPhreeqcRM::CleanupBMIModuleInstances(void)
/* ---------------------------------------------------------------------- */
{
	// This is necessary until template methods are fixed in StaticIndexer<>
	// see TestBMIdtor built as a dll on windows
	// ie PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx) causes missing symbol
	const std::lock_guard<std::mutex> lock(PhreeqcRM::_InstancesLock);
	std::list<BMIPhreeqcRM*> derived_items;
	for (auto pair : PhreeqcRM::StaticIndexer::_Instances)
	{
		if (BMIPhreeqcRM* derived = dynamic_cast<BMIPhreeqcRM*>(pair.second))
		{
			derived_items.push_back(derived);
		}
	}
	for (auto item : derived_items)
	{
		delete item;
	}
}
/* ---------------------------------------------------------------------- */
int
BMIPhreeqcRM::CreateBMIModule()
/* ---------------------------------------------------------------------- */
{
	//_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_crtBreakAlloc = 5144;
	try
	{
		BMIPhreeqcRM* bmirm_ptr = new BMIPhreeqcRM();
		if (bmirm_ptr)
		{
			bmirm_ptr->language = "F90";
			return bmirm_ptr->GetIndex();
		}
	}
	catch (...)
	{
		return IRM_OUTOFMEMORY;
	}
	return IRM_OUTOFMEMORY;
}
/* ---------------------------------------------------------------------- */
int
BMIPhreeqcRM::CreateBMIModule(int nxyz, MP_TYPE nthreads)
/* ---------------------------------------------------------------------- */
{
	//_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_crtBreakAlloc = 5144;
	try
	{
		BMIPhreeqcRM* bmirm_ptr = new BMIPhreeqcRM(nxyz, nthreads);
		if (bmirm_ptr)
		{
			bmirm_ptr->language = "F90";
			return bmirm_ptr->GetIndex();
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
	//return PhreeqcRM::Destroy<BMIPhreeqcRM>(id);
	if (BMIPhreeqcRM::GetInstance(id))
	{
		return PhreeqcRM::Destroy(id);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
BMIPhreeqcRM*
BMIPhreeqcRM::GetInstance(int id)
/* ---------------------------------------------------------------------- */
{
	//return PhreeqcRM::GetInstance<BMIPhreeqcRM>(id);
	return dynamic_cast<BMIPhreeqcRM*>(PhreeqcRM::GetInstance(id));
}
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM()
: PhreeqcRM(PhreeqcRM::default_nxyz, PhreeqcRM::default_data_for_parallel_processing, nullptr, true)
, var_man( nullptr )
, constructed( false )
{
	this->language = "cpp";
#if defined(WITH_PYBIND11)
	this->_initialized = false;
	this->language = "Py";
#endif
#if defined(swig_python_EXPORTS)
	this->language = "Py";
#endif
}
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, MP_TYPE nthreads)
: PhreeqcRM(nxyz, nthreads, nullptr, true) 
, var_man( nullptr )
, constructed( false )
{
	this->language = "cpp";
#if defined(WITH_PYBIND11)
	this->_initialized = false;
	this->language = "Py";
#endif
#if defined(swig_python_EXPORTS)
	this->language = "Py";
#endif
}
// Destructor
BMIPhreeqcRM::~BMIPhreeqcRM()
{
	delete(this->var_man);
}
void BMIPhreeqcRM::AddOutputVars(std::string option, std::string def)
{
	assert(this->var_man);
	this->var_man->AddOutputVars(option, def);
}
void BMIPhreeqcRM::ClearBMISelectedOutput(void)
{
	assert(this->var_man);
	this->var_man->BMISelectedOutput.clear();
}
void BMIPhreeqcRM::Construct(PhreeqcRM::Initializer i)
{
	if (constructed) return;
	this->PhreeqcRM::Construct(i);
	this->var_man = new VarManager((PhreeqcRM*)this);
#if defined(WITH_PYBIND11)
	this->_initialized = true;
#endif
#if defined(swig_python_EXPORTS) || defined(WITH_PYBIND11)
	phreeqcrm_io->Set_screen_on(false);
#endif
}

// Model control functions.
void BMIPhreeqcRM::Initialize(std::string config_file)
{
#ifdef USE_YAML
	if (config_file.size() != 0)
	{
		YAML::Node yaml = YAML::LoadFile(config_file);

		bool found_nxyz = false;
		bool found_threads = false;

		for (auto it = yaml.begin(); it != yaml.end(); it++)
		{
			YAML::Node node1 = *it;
			auto it1 = node1.begin();
			std::string keyword = it1++->second.as<std::string>();
			if (keyword == "SetGridCellCount")
			{
				this->initializer.nxyz_arg = it1++->second.as<int>();
				found_nxyz = true;
			}
			if (keyword == "ThreadCount")
			{
#if !defined(USE_MPI)
				this->initializer.data_for_parallel_processing = it1++->second.as<int>();
#endif
				found_threads = true;
			}
			if (found_threads && found_nxyz) break;
		}
	}
#endif

	this->Construct(this->initializer);
	constructed = true;
#ifdef USE_YAML
	if (config_file.size() != 0)
	{
		this->InitializeYAML(config_file);
	}
#endif
}
void BMIPhreeqcRM::Update()
{
	this->RunCells();
	this->SetTime(this->GetTime() + this->GetTimeStep());
	this->UpdateVariables();
}
void BMIPhreeqcRM::UpdateBMI(RMVARS v_enum)
{
	assert(this->var_man);
	this->var_man->RM2BMIUpdate(v_enum);
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
#ifdef USE_MPI
//	this->MpiWorkerBreak();
#endif
	this->CloseFiles();
#if defined(WITH_PYBIND11)
	delete this->var_man;
	this->var_man = nullptr;
	this->_initialized = false;
#endif
}
void BMIPhreeqcRM::GenerateAutoOutputVars()
{
#ifdef USE_MPI
	if (this->mpi_myself == 0)
	{
		if (var_man != nullptr)
		{
			var_man->GenerateAutoOutputVars();
			this->SetCurrentSelectedOutputUserNumber(var_man->BMISelectedOutputUserNumber);
			if (var_man->NeedInitialRun)
			{
				bool current = this->phreeqcrm_io->Get_screen_on();
				this->SetScreenOn(false);
				this->RunCells();
				this->SetScreenOn(current);
			}
			// Initialize BMI variables
			var_man->task = VarManager::VAR_TASKS::Info;
			for (auto it = this->var_man->VariantMap.begin();
				it != this->var_man->VariantMap.end(); it++)
			{
				BMIVariant& bv = it->second;
				bv.SetInitialized(false);
				((*this->var_man).*bv.GetFn())();
			}
		}
	}
#else
	if (var_man != nullptr)
	{ 
		var_man->GenerateAutoOutputVars();
		this->SetCurrentSelectedOutputUserNumber(var_man->BMISelectedOutputUserNumber);
		//if (var_man->NeedInitialRun)
		//{
		//	bool current = this->phreeqcrm_io->Get_screen_on();
		//	this->SetScreenOn(false);
		//	this->RunCells();
		//	this->SetScreenOn(current);
		//}
		// Initialize BMI variables
		var_man->task = VarManager::VAR_TASKS::Info;
		for (auto it = this->var_man->VariantMap.begin();
			it != this->var_man->VariantMap.end(); it++)
		{
			BMIVariant& bv = it->second;
			bv.SetInitialized(false);
			((*this->var_man).*bv.GetFn())();
		}
	}
#endif
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
	count += (int)this->var_man->AutoOutputVars.size();
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
	for(auto it = var_man->AutoOutputVars.begin();
		it != var_man->AutoOutputVars.end(); it++)
	{ 
		names.push_back(it->second.GetName());
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

std::vector<std::string> BMIPhreeqcRM::GetReadOnlyVarNames()
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
		if (!bv.GetHasSetter())
		{
			names.push_back(bv.GetName());
		}
	}
	return names;
}

std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	RMVARS v_enum = this->GetEnum(name);
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
		else if (this->language == "C")
		{
			return bv.GetClangType();
		}
	}
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			if (this->language == "cpp")
			{
				return it->second.GetCType();
			}
			else if (this->language == "F90")
			{
				return it->second.GetFType();
			}
			else if (this->language == "Py")
			{
				return it->second.GetPType();
			}
			else if (this->language == "C")
			{
				return it->second.GetClangType();
			}
		}
	}
	throw std::runtime_error("Failed in GetVarType.");
	return "Failed in GetVarType.";
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	RMVARS v_enum = this->GetEnum(name);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			return it->second.GetUnits();
		}
	}
	throw std::runtime_error("Failed in GetVarUnits.");
	return "Failed in GetVarUnits.";
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	RMVARS v_enum = this->GetEnum(name);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			return it->second.GetItemsize();
		}
	}
	throw std::runtime_error("Failed in GetVarItemsize.");
	return 0;
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	RMVARS v_enum = this->GetEnum(name);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			return it->second.GetNbytes();
		}
	}
	throw std::runtime_error("Failed in GetVarNbytes.");
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
	RMVARS v_enum = this->GetEnum(name);
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
			bool tf = *this->var_man->VarExchange.GetBVarPtr();
			int tf_int = 1;
			if (!tf) tf_int = 0;
			assert(tf == this->var_man->VarExchange.GetBVar());
			memcpy(dest, &tf_int, sizeof(int));
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
		if (this->var_man->VarExchange.GetCType() == "std::string" && dim == 0)
		{
			memcpy(dest, this->var_man->VarExchange.GetStringRef().data(), Nbytes);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			if (var_man->BMISelectedOutput.size() == 0)
			{
				int n_user = GetCurrentSelectedOutputUserNumber();
				SetCurrentSelectedOutputUserNumber(var_man->BMISelectedOutputUserNumber);
				this->GetSelectedOutput(var_man->BMISelectedOutput);
				SetCurrentSelectedOutputUserNumber(n_user);
			}
			int column = it->second.GetColumn();
			int nxyz = GetGridCellCount();
			void* ptr = (void*)&(var_man->BMISelectedOutput[column * nxyz]);
			memcpy(dest, ptr, it->second.GetNbytes());
			return;
		}
	}	
	std::ostringstream oss;
	oss << "BMI GetValue void* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool* dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double* dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			if (var_man->BMISelectedOutput.size() == 0)
			{
				int n_user = GetCurrentSelectedOutputUserNumber();
				SetCurrentSelectedOutputUserNumber(var_man->BMISelectedOutputUserNumber);
				this->GetSelectedOutput(var_man->BMISelectedOutput);
				SetCurrentSelectedOutputUserNumber(n_user);
			}
			int column = it->second.GetColumn();
			int nxyz = GetGridCellCount();
			void* ptr = (void*)&(var_man->BMISelectedOutput[column * nxyz]);
			memcpy(dest, ptr, it->second.GetNbytes());
			return;
		}
	}
	std::ostringstream oss;
	oss << "BMI GetValue double* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int* dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
		if (this->language == "C")
		{
			assert(this->var_man->VarExchange.GetClangType() == "int");
		}
		else
		{
			assert(this->var_man->VarExchange.GetCType() == "int");
		}
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::string& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<double>& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	{
		std::string name_lc = name;
		std::transform(name_lc.begin(), name_lc.end(),
			name_lc.begin(), tolower);
		auto it = var_man->AutoOutputVars.find(name_lc);
		if (it != var_man->AutoOutputVars.end())
		{
			if (var_man->BMISelectedOutput.size() == 0)
			{
				int n_user = GetCurrentSelectedOutputUserNumber();
				SetCurrentSelectedOutputUserNumber(var_man->BMISelectedOutputUserNumber);
				this->GetSelectedOutput(var_man->BMISelectedOutput);
				SetCurrentSelectedOutputUserNumber(n_user);
			}
			int column = it->second.GetColumn();
			int nxyz = GetGridCellCount();
			void* ptr = (void*)&(var_man->BMISelectedOutput[column * nxyz]);
			dest.resize(nxyz);
			memcpy(dest.data(), ptr, it->second.GetNbytes());
			return;
		}
	}
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<int>& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<std::string>& dest)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in GetValue.");
	return;
}
void* BMIPhreeqcRM::GetValuePtr(const std::string name)
{
	RMVARS v_enum = this->GetEnum(name);
	if (v_enum != RMVARS::NotFound)
	{
		this->var_man->SetLanguage(this->language);
		BMIVariant& bv = this->var_man->VariantMap[v_enum];
		//VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
		if (bv.GetVoidPtr() == NULL)
		{
			this->var_man->task = VarManager::VAR_TASKS::GetPtr;
			((*this->var_man).*bv.GetFn())();
		}
		return bv.GetVoidPtr();
	}
	throw std::runtime_error("Failed in GetValuePtr.");
	return NULL;
}
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	RMVARS v_enum = this->GetEnum(name);
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
			int src_int = *(int*)src;
			bool tf = true;
			if (src_int == 0) tf = false;
			memcpy(this->var_man->VarExchange.GetBVarPtr(), &tf, Nbytes);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, bool src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, const char* src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, double src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, int src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, const std::string src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<double> src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<int> src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<std::string> src)
{
	RMVARS v_enum = this->GetEnum(name);
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
	throw std::runtime_error("Failed in SetValue.");
	return;
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

RMVARS BMIPhreeqcRM::GetEnum(const std::string name)
{
	if (this->var_man != nullptr)
	{
		return this->var_man->GetEnum(name);
	}
	std::cerr << "BMIPhreeqcRM has not been initialized." << std::endl;
	return RMVARS::NotFound;
};

//////////////////


