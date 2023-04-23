#include "VarManager.h"
#include "IPhreeqcPhast.h"
#include "Phreeqc.h"
#include "RMVARS.h"
#include <assert.h>
#include <algorithm>
VarManager::VarManager(PhreeqcRM* rm_ptr_in) 
{
	this->VariantMap[RMVARS::ComponentCount] = 
		BMIVariant(&VarManager::ComponentCount_Var, "ComponentCount");
	this->VariantMap[RMVARS::Components] =
		BMIVariant(&VarManager::Components_Var, "Components");
	this->VariantMap[RMVARS::Concentrations] =
		BMIVariant(&VarManager::Concentrations_Var, "Concentrations");
	this->VariantMap[RMVARS::CurrentSelectedOutputUserNumber] =
		BMIVariant(&VarManager::CurrentSelectedOutputUserNumber_Var, "CurrentSelectedOutputUserNumber");
	this->VariantMap[RMVARS::Density] =
		BMIVariant(&VarManager::Density_Var, "Density");
	this->VariantMap[RMVARS::ErrorString] =
		BMIVariant(&VarManager::ErrorString_Var, "ErrorString");
	this->VariantMap[RMVARS::FilePrefix] =
		BMIVariant(&VarManager::FilePrefix_Var, "FilePrefix");
	this->VariantMap[RMVARS::Gfw] =
		BMIVariant(&VarManager::Gfw_Var, "Gfw");
	this->VariantMap[RMVARS::GridCellCount] =
		BMIVariant(&VarManager::GridCellCount_Var, "GridCellCount");
	this->VariantMap[RMVARS::NthSelectedOutput] =
		BMIVariant(&VarManager::NthSelectedOutput_Var, "NthSelectedOutput");
	this->VariantMap[RMVARS::Porosity] =
		BMIVariant(&VarManager::Porosity_Var, "Porosity");
	this->VariantMap[RMVARS::Pressure] =
		BMIVariant(&VarManager::Pressure_Var, "Pressure");
	this->VariantMap[RMVARS::Saturation] =
		BMIVariant(&VarManager::Saturation_Var, "Saturation");
	this->VariantMap[RMVARS::SelectedOutput] =
		BMIVariant(&VarManager::SelectedOutput_Var, "SelectedOutput");
	this->VariantMap[RMVARS::SelectedOutputColumnCount] =
		BMIVariant(&VarManager::SelectedOutputColumnCount_Var, "SelectedOutputColumnCount");
	this->VariantMap[RMVARS::SelectedOutputCount] =
		BMIVariant(&VarManager::SelectedOutputCount_Var, "SelectedOutputCount");
	this->VariantMap[RMVARS::SelectedOutputHeadings] =
		BMIVariant(&VarManager::SelectedOutputHeadings_Var, "SelectedOutputHeadings");
	this->VariantMap[RMVARS::SelectedOutputOn] =
		BMIVariant(&VarManager::SelectedOutputOn_Var, "SelectedOutputOn");
	this->VariantMap[RMVARS::SelectedOutputRowCount] =
		BMIVariant(&VarManager::SelectedOutputRowCount_Var, "SelectedOutputRowCount");
	this->VariantMap[RMVARS::SolutionVolume] =
		BMIVariant(&VarManager::SolutionVolume_Var, "SolutionVolume");
	this->VariantMap[RMVARS::Temperature] =
		BMIVariant(&VarManager::Temperature_Var, "Temperature");
	this->VariantMap[RMVARS::Time] =
		BMIVariant(&VarManager::Time_Var, "Time");
	this->VariantMap[RMVARS::TimeStep] =
		BMIVariant(&VarManager::TimeStep_Var, "TimeStep");
	this->VariantMap[RMVARS::Viscosity] =
		BMIVariant(&VarManager::Viscosity_Var, "Viscosity");
	///!!!VarFunction x = &VarManager::ComponentCount_var;
	///!!! (this->*x)(rm_ptr); // Remember this !!!///
	//auto it = VariantMap.begin();
	//VarFunction x = it->second.GetFn();
	//(this->*x)();
	rm_ptr = rm_ptr_in;
	this->task = VarManager::VAR_TASKS::no_op;
	this->CurrentVar = RMVARS::NotFound;
	//auto it = VarVariantMap.begin();
	//VarFunction x = it->second.GetFn(rm_ptr);
	//x(rm_ptr);
	// Make EnumMap
	for (auto it = VariantMap.begin(); it != VariantMap.end(); it++)
	{
		std::string name_lc = it->second.GetName();
		std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
		EnumMap[name_lc] = it->first;
	}
}
void VarManager::RM2BMIUpdate(RMVARS v_enum)
{
	if (this->PointerSet.size() == 0) return;
	//RMVARS v_enum = this->GetEnum(name);
	if (this->GetCurrentVar() == v_enum) return;
	auto it = this->VariantMap.find(v_enum);
	if (it != VariantMap.end())
	{
		this->task = VarManager::VAR_TASKS::RMUpdate;
		VarFunction f =  it->second.GetFn(); 
		(this->*f)();
	}
	return;
}
RMVARS VarManager::GetEnum(const std::string name)
{
	std::string name_lc = name;
	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
	auto m_it = EnumMap.find(name_lc);
	if (m_it != EnumMap.end())
	{
		return m_it->second;
	}
	return RMVARS::NotFound;
}
VarManager::VarFunction VarManager::GetFn(RMVARS v_enum)
{
	auto it = this->GetVariantMap().find(v_enum);
	if (it != this->GetVariantMap().end())
	{
		//VarManager::VarFunction f = it->second.GetFn();
		//VarFunction*(f)();
		return it->second.GetFn();
	}
	return NULL;
}

//// Start_var
void VarManager::ComponentCount_Var()
{
	RMVARS VARS_myself = RMVARS::ComponentCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("count", false, true, true, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetComponentCount());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		int v = rm_ptr->GetComponentCount();
		bv.SetIVar(v);
		bv.SetVoidPtr((void*)(bv.GetIVarPtr()));
		this->PointerSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetComponentCount();
		bv.SetIVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}	
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Components_Var()
{
	RMVARS VARS_myself = RMVARS::Components;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		std::vector<std::string> comps = rm_ptr->GetComponents();
		size_t size = 0;
		for (size_t i = 0; i < comps.size(); i++)
		{
			if (comps[i].size() > size) size = comps[i].size();
		}
		int Itemsize = (int)size;
		int Nbytes = (int)(size * comps.size());
		//name, std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("names", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
		bv.SetStringVector(rm_ptr->GetComponents());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		std::vector<const char*> CharVector;
		std::vector<std::string>&  Components = bv.GetStringVectorRef();
		for (size_t i = 0; i < Components.size(); i++)
		{
			CharVector.push_back(Components[i].c_str());
		}
		bv.SetCharVector(CharVector);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		const std::vector<std::string>& v = rm_ptr->GetComponents();
		this->VarExchange.SetStringVector(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Concentrations_Var()
{
	RMVARS VARS_myself = RMVARS::Concentrations;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(double);
		int Nbytes = (int)sizeof(double) *
			rm_ptr->GetGridCellCount() * rm_ptr->GetComponentCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize
		int option = rm_ptr->GetUnitsSolution();
		std::string units;
		switch (option)
		{
		case 1:
			units = "mg/L";
			break;
		case 2:
			units = "mol/L";
			break;
		case 3:
			units = "kg/kgs";
			break;
		}
		bv.SetBasic(units, true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		rm_ptr->GetConcentrations(bv.GetDoubleVectorRef());
		rm_ptr->GetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		rm_ptr->GetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		rm_ptr->GetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetConcentrations(this->VarExchange.GetDoubleVectorRef());
		rm_ptr->GetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Density_Var()
{
	this->SetCurrentVar(RMVARS::Density);
	BMIVariant& bv = this->VariantMap[RMVARS::Density];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("kg L-1", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		rm_ptr->GetDensity(this->VarExchange.GetDoubleVectorRef());
		rm_ptr->GetDensity(bv.GetDoubleVectorRef());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		rm_ptr->GetDensity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(RMVARS::Density);
		this->UpdateSet.insert(RMVARS::Density);
		break;
	}
	case VarManager::VAR_TASKS::RMUpdate:
	{
		// Concentrations may change when SetDensity is called
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	{
		rm_ptr->GetDensity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		rm_ptr->SetDensity(this->VarExchange.GetDoubleVectorRef());
		// don't update density vector, only get solution densities
		//bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		// Concentrations may change when SetDensity is called
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::ErrorString_Var()
{
	RMVARS VARS_myself = RMVARS::ErrorString;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)rm_ptr->GetErrorString().size();
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("error", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("std::string", "character", "");
		this->VarExchange.GetStringRef() = rm_ptr->GetErrorString(); 
		bv.GetStringRef() = rm_ptr->GetErrorString();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetStringRef() = rm_ptr->GetErrorString();
		bv.GetStringRef() = rm_ptr->GetErrorString();
		bv.SetItemsize(bv.GetStringRef().size());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::FilePrefix_Var()
{
	RMVARS VARS_myself = RMVARS::FilePrefix;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	//if (!bv.GetInitialized())
	//{
		int Itemsize = (int)rm_ptr->GetFilePrefix().size();
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("prefix", true, true, false, Nbytes, Itemsize);
		bv.SetTypes("std::string", "character", "");
		//this->VarExchange.GetStringRef() = rm_ptr->GetFilePrefix();
		bv.GetStringRef() = rm_ptr->GetFilePrefix();
		//bv.SetInitialized(true);
	//}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		//int Itemsize = (int)rm_ptr->GetFilePrefix().size();
		//int Nbytes = Itemsize;
		////name, std::string units, set, get, ptr, Nbytes, Itemsize  
		//bv.SetBasic("prefix", true, true, false, Nbytes, Itemsize);
		//bv.SetTypes("std::string", "character", "");
		this->VarExchange.GetStringRef() = rm_ptr->GetFilePrefix();
		bv.GetStringRef() = rm_ptr->GetFilePrefix();
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		int Itemsize = (int)this->VarExchange.GetStringRef().size();
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize 
		bv.SetBasic("prefix", true, true, false, Nbytes, Itemsize);
		//bv.SetTypes("std::string", "character", "");
		rm_ptr->SetFilePrefix(this->VarExchange.GetStringRef());
		bv.GetStringRef() = this->VarExchange.GetStringRef();
		break;
	}
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Gfw_Var()
{
	RMVARS VARS_myself = RMVARS::Gfw;
	this->SetCurrentVar(RMVARS::Gfw);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetComponentCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("g mol-1", false, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetGfw();
		bv.GetDoubleVectorRef() = rm_ptr->GetGfw();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetGfw();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetGfw();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::GridCellCount_Var()
{
	RMVARS VARS_myself = RMVARS::GridCellCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("count", false, true, true, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetGridCellCount());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		int v = rm_ptr->GetGridCellCount();
		bv.SetIVar(v);
		bv.SetVoidPtr((void*)(bv.GetIVarPtr()));
		this->PointerSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetGridCellCount();
		bv.SetIVar(v);
		this->PointerSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
//
//void VarManager::InputVarNames_Var() // Implement in BRMPhreeqcRM TODO???
//
void VarManager::NthSelectedOutput_Var()
{
	RMVARS VARS_myself = RMVARS::NthSelectedOutput;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("id", true, false, false, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(-1);
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	case VarManager::VAR_TASKS::GetVar:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		int v = *this->VarExchange.GetIVarPtr();
		rm_ptr->SetNthSelectedOutput(v);
		break;
	}
	case VarManager::VAR_TASKS::RMUpdate:
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
//
//void VarManager::OutputVarNames_Var() // Implement in BRMPhreeqcRM TODO???
//
void VarManager::Saturation_Var()
{
	RMVARS VARS_myself = RMVARS::Saturation;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("unitless", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		rm_ptr->GetSaturation(this->VarExchange.GetDoubleVectorRef());
		rm_ptr->GetSaturation(bv.GetDoubleVectorRef());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		rm_ptr->GetSaturation(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::RMUpdate:
	{
		// GetSaturation does not change with SetSaturation
		// But Concentrations may change
		//rm_ptr->GetSaturation(this->VarExchange.GetDoubleVectorRef());
		//bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		// Concentrations change when Saturation changes
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	{
		rm_ptr->GetSaturation(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		rm_ptr->SetSaturation(this->VarExchange.GetDoubleVectorRef());
		// Saturation is solution volume / (porosity*rv)
		// GetSaturation does not change with SetSaturation
		// bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		// Concentrations may change when Saturation changes
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutput_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutput;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetSelectedOutputRowCount() *
			rm_ptr->GetSelectedOutputColumnCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("user specified", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		rm_ptr->GetSelectedOutput(this->VarExchange.GetDoubleVectorRef());
		rm_ptr->GetSelectedOutput(bv.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
		assert(false);
		break;
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutputColumnCount_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutputColumnCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("count", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetSelectedOutputColumnCount());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetSelectedOutputColumnCount();
		bv.SetIVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutputCount_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutputCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("count", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetSelectedOutputCount());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetSelectedOutputCount();
		bv.SetIVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutputHeadings_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutputHeadings;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		rm_ptr->GetSelectedOutputHeadings(bv.GetStringVectorRef());
		this->VarExchange.GetStringVectorRef() = bv.GetStringVectorRef();
		std::vector<std::string> headings = bv.GetStringVectorRef();
		size_t size = 0;
		for (size_t i = 0; i < headings.size(); i++)
		{
			if (headings[i].size() > size) size = headings[i].size();
		}
		int Itemsize = (int)size;
		int Nbytes = (int)(size * headings.size());
		//name, std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("names", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		//std::vector<const char*> CharVector;
		//std::vector<std::string>& Components = bv.GetStringVectorRef();
		//for (size_t i = 0; i < Components.size(); i++)
		//{
		//	CharVector.push_back(Components[i].c_str());
		//}
		//bv.SetCharVector(CharVector);
		//this->PointerSet.insert(VARS_myself);
		//this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		rm_ptr->GetSelectedOutputHeadings(bv.GetStringVectorRef());
		this->VarExchange.SetStringVector(bv.GetStringVectorRef());
		std::vector<std::string>& headings = bv.GetStringVectorRef();
		size_t size = 0;
		for (size_t i = 0; i < headings.size(); i++)
		{
			if (headings[i].size() > size) size = headings[i].size();
		}
		bv.SetItemsize(size);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutputRowCount_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutputRowCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("count", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetSelectedOutputRowCount());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetSelectedOutputRowCount();
		bv.SetIVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SolutionVolume_Var()
{
	RMVARS VARS_myself = RMVARS::SolutionVolume;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("L", false, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetSolutionVolume();
		bv.GetDoubleVectorRef() = rm_ptr->GetSolutionVolume();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetSolutionVolume();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetSolutionVolume();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Time_Var()
{
	RMVARS VARS_myself = RMVARS::Time;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("s", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.SetDVar(rm_ptr->GetTime());
		bv.SetDVar(rm_ptr->GetTime());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.SetDVar(rm_ptr->GetTime());
		bv.SetDVar(rm_ptr->GetTime());
		bv.SetVoidPtr((void*)(bv.GetDVarPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.SetDVar(rm_ptr->GetTime());
		bv.SetDVar(rm_ptr->GetTime());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetTime(*this->VarExchange.GetDVarPtr());
		bv.SetDVar(*this->VarExchange.GetDVarPtr());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::TimeStep_Var()
{
	RMVARS VARS_myself = RMVARS::TimeStep;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("s", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.SetDVar(rm_ptr->GetTimeStep());
		bv.SetDVar(rm_ptr->GetTimeStep());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.SetDVar(rm_ptr->GetTimeStep());
		bv.SetDVar(rm_ptr->GetTimeStep());
		bv.SetVoidPtr((void*)(bv.GetDVarPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.SetDVar(rm_ptr->GetTimeStep());
		bv.SetDVar(rm_ptr->GetTimeStep());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetTimeStep(*this->VarExchange.GetDVarPtr());
		bv.SetDVar(*this->VarExchange.GetDVarPtr());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::CurrentSelectedOutputUserNumber_Var()
{
	RMVARS VARS_myself = RMVARS::CurrentSelectedOutputUserNumber;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(int);
		int Nbytes = (int)sizeof(int);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("id", false, true, false, Nbytes, Itemsize);
		bv.SetTypes("int", "integer", "int");
		bv.SetIVar(rm_ptr->GetCurrentSelectedOutputUserNumber());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetCurrentSelectedOutputUserNumber();
		bv.SetIVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::RMUpdate:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Porosity_Var()
{
	RMVARS VARS_myself = RMVARS::Porosity;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("unitless", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}	
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		// Concentrations change when Saturation changes
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
		// Saturation changes with change in porosity
		std::vector<double> sat;
		rm_ptr->GetSaturation(sat);
		RMVARS VARS_sat = RMVARS::Saturation;
		BMIVariant& bv_sat = this->VariantMap[VARS_sat];
		bv_sat.SetDoubleVector(sat);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		rm_ptr->SetPorosity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		// Concentrations change when Saturation changes
		std::vector<double> c;
		rm_ptr->GetConcentrations(c);
		RMVARS VARS_c = RMVARS::Concentrations;
		BMIVariant& bv_c = this->VariantMap[VARS_c];
		bv_c.SetDoubleVector(c);
		// Saturation changes with change in porosity
		std::vector<double> sat;
		rm_ptr->GetSaturation(sat);
		RMVARS VARS_sat = RMVARS::Saturation;
		BMIVariant& bv_sat = this->VariantMap[VARS_sat];
		bv_sat.SetDoubleVector(sat);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Pressure_Var()
{
	RMVARS VARS_myself = RMVARS::Pressure;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("atm", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPressure();
		bv.GetDoubleVectorRef() = rm_ptr->GetPressure();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPressure();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPressure();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetPressure(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::SelectedOutputOn_Var()
{
	RMVARS VARS_myself = RMVARS::SelectedOutputOn;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = (int)sizeof(bool);
		int Nbytes = (int)sizeof(bool);
		//std::string units, set, get, ptr, Nbytes, Itemsize
		bv.SetBasic("bool", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("bool", "logical", "");
		bv.SetBVar(rm_ptr->GetSelectedOutputOn());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		int v = rm_ptr->GetSelectedOutputOn();
		bv.SetBVar(v);
		bv.SetVoidPtr((void*)(bv.GetBVarPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	case VarManager::VAR_TASKS::GetVar:
	{
		bool v = rm_ptr->GetSelectedOutputOn();
		bv.SetBVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		bool v = this->VarExchange.GetBVar();
		bv.SetBVar(v);
		rm_ptr->SetSelectedOutputOn(v);
		break;
	}
	case VarManager::VAR_TASKS::Info:
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Temperature_Var()
{
	RMVARS VARS_myself = RMVARS::Temperature;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("C", true, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetTemperature();
		bv.GetDoubleVectorRef() = rm_ptr->GetTemperature();
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetTemperature();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetTemperature();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetTemperature(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
void VarManager::Viscosity_Var()
{
	RMVARS VARS_myself = RMVARS::Viscosity;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		int Itemsize = sizeof(double);
		int Nbytes = Itemsize * rm_ptr->GetGridCellCount();
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("mPa s", false, true, true, Nbytes, Itemsize);
		bv.SetTypes("double", "real(kind=8)", "");
		rm_ptr->GetViscosity(this->VarExchange.GetDoubleVectorRef());
		rm_ptr->GetViscosity(bv.GetDoubleVectorRef());
		bv.SetInitialized(true);
	}
	switch (this->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		rm_ptr->GetViscosity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		bv.SetVoidPtr((void*)(bv.GetDoubleVectorPtr()));
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::Update:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		rm_ptr->GetViscosity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		assert(false);
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(RMVARS::NotFound);
}
/// end_
////////////////////////////////

void VarManager::BMIGenerateSelectedOutput()
{
	int line_no = 10;
	int Itemsize = (int)sizeof(double);
	int Nbytes = rm_ptr->GetGridCellCount() * Itemsize;
	std::ostringstream headings;
	std::ostringstream code;
	BMISelectedOutputVars.clear();
	auto it = BMISelecteOutputDefs.begin();
	for (; it != BMISelecteOutputDefs.end(); it++)
	{
		if (it->first == "output_solution_properties")
		{
			switch (BMICheckSelectedOutputDef(true, it->second))
			{
			case -1:
				rm_ptr->ErrorMessage("Unknown input for output_solution_properties", true);
				//return_value = IRM_INVALIDARG;
				continue;
			case 0:
				continue;
			case 1:
			{
				std::string name = "solution_ph";
				BMIVariant bv(name, "-", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH -LA('H+')" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_pe";
				BMIVariant bv(name, "-", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH -LA('e-')" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_alkalinity";
				BMIVariant bv(name, "eq kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH ALK" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_ionic_strength";
				BMIVariant bv(name, "mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH MU" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_water_mass";
				BMIVariant bv(name, "kg", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH TOT('water')" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_charge_balance";
				BMIVariant bv(name, "eq kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH CHARGE_BALANCE / TOT('water')" << std::endl;
				line_no += 10;
				break;
			}
			{
				std::string name = "solution_percent_error";
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
		else if (it->first == "output_solution_total_molalities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				item_set = rm_ptr->ElementRedoxSet;
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
				std::string name = "solution_total_molality_" + *item_it;
				BMIVariant bv(name, "mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH TOT('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_solution_molalities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < rm_ptr->species_names.size(); i++)
				{
					item_set.insert(rm_ptr->species_names[i]);
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
				std::string name = "solution_species_log_molality_" + *item_it;
				BMIVariant bv(name, "log mol kgw-1", false, true, false, Nbytes, Itemsize);
				bv.SetColumn((int)BMISelectedOutputVars.size());
				BMISelectedOutputVars[name] = bv;
				headings << name << "\t";
				code << line_no << " PUNCH LM('" << name << "')\n";
				line_no += 10;
			}
		}
		else if (it->first == "output_solution_activities")
		{
			std::set<std::string> item_set;
			switch (BMICheckSelectedOutputDef(false, it->second))
			{
			case 0:
				continue;
			case 1:
			{
				for (size_t i = 0; i < rm_ptr->species_names.size(); i++)
				{
					item_set.insert(rm_ptr->species_names[i]);
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
				std::string name = "solution_species_log_activity_" + *item_it;
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
				for (size_t i = 0; i < rm_ptr->ExchangeSpeciesNamesList.size(); i++)
				{
					item_set.insert(rm_ptr->ExchangeSpeciesNamesList[i]);
					item_map[rm_ptr->ExchangeSpeciesNamesList[i]] = rm_ptr->ExchangeNamesList[i];
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
				for (size_t i = 0; i < rm_ptr->SurfaceNamesList.size(); i++)
				{
					item_set.insert(rm_ptr->SurfaceNamesList[i]);
					item_map[rm_ptr->SurfaceNamesList[i]] = rm_ptr->SurfaceTypesList[i];
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
				for (size_t i = 0; i < rm_ptr->EquilibriumPhasesList.size(); i++)
				{
					item_set.insert(rm_ptr->EquilibriumPhasesList[i]);
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
				for (size_t i = 0; i < rm_ptr->SINamesList.size(); i++)
				{
					item_set.insert(rm_ptr->SINamesList[i]);
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
				for (size_t i = 0; i < rm_ptr->GasComponentsList.size(); i++)
				{
					item_set.insert(rm_ptr->GasComponentsList[i]);
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
				for (size_t i = 0; i < rm_ptr->KineticReactionsList.size(); i++)
				{
					item_set.insert(rm_ptr->KineticReactionsList[i]);
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
				for (size_t i = 0; i < rm_ptr->SolidSolutionComponentsList.size(); i++)
				{
					item_set.insert(rm_ptr->SolidSolutionComponentsList[i]);
					item_map[rm_ptr->SolidSolutionComponentsList[i]] = rm_ptr->SolidSolutionNamesList[i];
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
				auto it = rm_ptr->GetWorkers()[0]->Get_PhreeqcPtr()->GetCalculateValueMap().begin();
				for (; it != rm_ptr->workers[0]->Get_PhreeqcPtr()->GetCalculateValueMap().end(); it++)
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
			rm_ptr->ErrorMessage(oss.str(), true);
			throw PhreeqcStop();
		}
	}
	std::ostringstream data_block;
	data_block << "SELECTED_OUTPUT 777777777; USER_PUNCH 777777777;" << std::endl;
	data_block << headings.str() << std::endl;
	data_block << code.str() << std::endl;
	rm_ptr->RunString(true, false, false, data_block.str());
	BMISelecteOutputDefs.clear();
	return;
}
int VarManager::BMICheckSelectedOutputDef(bool tf_only, std::string& def)
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
std::set<std::string> VarManager::tokenize(const std::string& def_in)
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


//////////////////
