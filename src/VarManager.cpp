#include "VarManager.h"
#include <assert.h>
#include <algorithm>
VarManager::VarManager(PhreeqcRM* rm_ptr_in) 
{
	this->VariantMap[VARS::ComponentCount] = 
		BMIVariant(&VarManager::ComponentCount_Var, "ComponentCount");
	this->VariantMap[VARS::Components] =
		BMIVariant(&VarManager::Components_Var, "Components");
	this->VariantMap[VARS::Concentrations] =
		BMIVariant(&VarManager::Concentrations_Var, "Concentrations");
	this->VariantMap[VARS::Density] =
		BMIVariant(&VarManager::Density_Var, "Density");
	this->VariantMap[VARS::ErrorString] =
		BMIVariant(&VarManager::ErrorString_Var, "ErrorString");
	this->VariantMap[VARS::FilePrefix] =
		BMIVariant(&VarManager::FilePrefix_Var, "FilePrefix");
	this->VariantMap[VARS::Gfw] =
		BMIVariant(&VarManager::Gfw_Var, "Gfw");
	this->VariantMap[VARS::GridCellCount] =
		BMIVariant(&VarManager::GridCellCount_Var, "GridCellCount");
	//this->VariantMap[VARS::InputVarNames] =
	//	BMIVariant(&VarManager::InputVarNames_Var, "InputVarNames");
	this->VariantMap[VARS::NthSelectedOutput] =
		BMIVariant(&VarManager::NthSelectedOutput_Var, "NthSelectedOutput");
	//this->VariantMap[VARS::OutputVarNames] =
	//	BMIVariant(&VarManager::OutputVarNames_Var, "OutputVarNames");
	this->VariantMap[VARS::Saturation] =
		BMIVariant(&VarManager::Saturation_Var, "Saturation");
	this->VariantMap[VARS::SelectedOutput] =
		BMIVariant(&VarManager::SelectedOutput_Var, "SelectedOutput");
	this->VariantMap[VARS::SelectedOutputColumnCount] =
		BMIVariant(&VarManager::SelectedOutputColumnCount_Var, "SelectedOutputColumnCount");
	this->VariantMap[VARS::SelectedOutputCount] =
		BMIVariant(&VarManager::SelectedOutputCount_Var, "SelectedOutputCount");
	this->VariantMap[VARS::SelectedOutputHeadings] =
		BMIVariant(&VarManager::SelectedOutputHeadings_Var, "SelectedOutputHeadings");
	this->VariantMap[VARS::SelectedOutputRowCount] =
		BMIVariant(&VarManager::SelectedOutputRowCount_Var, "SelectedOutputRowCount");
	this->VariantMap[VARS::SolutionVolume] =
		BMIVariant(&VarManager::SolutionVolume_Var, "SolutionVolume");
	this->VariantMap[VARS::Time] =
		BMIVariant(&VarManager::Time_Var, "Time");
	this->VariantMap[VARS::TimeStep] =
		BMIVariant(&VarManager::TimeStep_Var, "TimeStep");
	this->VariantMap[VARS::CurrentSelectedOutputUserNumber] =
		BMIVariant(&VarManager::CurrentSelectedOutputUserNumber_Var, "CurrentSelectedOutputUserNumber");
	this->VariantMap[VARS::Porosity] =
		BMIVariant(&VarManager::Porosity_Var, "Porosity");
	this->VariantMap[VARS::Pressure] =
		BMIVariant(&VarManager::Pressure_Var, "Pressure");
	this->VariantMap[VARS::SelectedOutputOn] =
		BMIVariant(&VarManager::SelectedOutputOn_Var, "SelectedOutputOn");
	this->VariantMap[VARS::Temperature] =
		BMIVariant(&VarManager::Temperature_Var, "Temperature");
	///!!!VarFunction x = &VarManager::ComponentCount_var;
	///!!! (this->*x)(rm_ptr); // Remember this !!!///
	//auto it = VariantMap.begin();
	//VarFunction x = it->second.GetFn();
	//(this->*x)();
	rm_ptr = rm_ptr_in;
	this->task = VarManager::VAR_TASKS::no_op;
	this->CurrentVar = VarManager::VARS::NotFound;
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
void VarManager::RM2BMIUpdate(VarManager::VARS v_enum)
{
	if (this->PointerSet.size() == 0) return;
	//VarManager::VARS v_enum = this->GetEnum(name);
	if (this->GetCurrentVar() != v_enum) return;
	auto it = this->VariantMap.find(v_enum);
	if (it != VariantMap.end())
	{
		this->task = VarManager::VAR_TASKS::RMUpdate;
		VarFunction f =  it->second.GetFn(); 
		(this->*f)();
	}
	return;
}
VarManager::VARS VarManager::GetEnum(const std::string name)
{
	std::string name_lc = name;
	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
	auto m_it = EnumMap.find(name_lc);
	if (m_it != EnumMap.end())
	{
		return m_it->second;
	}
	return VarManager::VARS::NotFound;
}
//// Start_var
void VarManager::ComponentCount_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::ComponentCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Components_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Components;
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
		const std::vector<std::string>& v = rm_ptr->GetComponents();
		bv.SetStringVector(v);
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Concentrations_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Concentrations;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{

		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		rm_ptr->GetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetConcentrations(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Density_Var()
{
	this->SetCurrentVar(VarManager::VARS::Density);
	BMIVariant& bv = this->VariantMap[VarManager::VARS::Density];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VarManager::VARS::Density];
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
		this->PointerSet.insert(VarManager::VARS::Density);
		this->UpdateSet.insert(VarManager::VARS::Density);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		rm_ptr->GetDensity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetDensity(this->VarExchange.GetDoubleVectorRef());
		// don't update density vector, only get solution densities
		//bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::ErrorString_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::ErrorString;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetStringRef() = rm_ptr->GetErrorString();
		bv.GetStringRef() = rm_ptr->GetErrorString();
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::FilePrefix_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::FilePrefix;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
		int Itemsize = (int)rm_ptr->GetFilePrefix().size();
		int Nbytes = Itemsize;
		//name, std::string units, set, get, ptr, Nbytes, Itemsize  
		bv.SetBasic("error", true, true, false, Nbytes, Itemsize);
		bv.SetTypes("std::string", "character", "");
		this->VarExchange.GetStringRef() = rm_ptr->GetFilePrefix();
		bv.GetStringRef() = rm_ptr->GetFilePrefix();
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
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetStringRef() = rm_ptr->GetFilePrefix();
		bv.GetStringRef() = rm_ptr->GetFilePrefix();
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetFilePrefix(this->VarExchange.GetStringRef());
		bv.GetStringRef() = this->VarExchange.GetStringRef();
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Gfw_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Gfw;
	this->SetCurrentVar(VarManager::VARS::Gfw);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		//
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::GridCellCount_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::GridCellCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
//
//void VarManager::InputVarNames_Var() // Implement in BRMPhreeqcRM TODO???
//
void VarManager::NthSelectedOutput_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::NthSelectedOutput;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
//
//void VarManager::OutputVarNames_Var() // Implement in BRMPhreeqcRM TODO???
//
void VarManager::Saturation_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Saturation;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		rm_ptr->GetSaturation(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetSaturation(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutput_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutput;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutputColumnCount_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutputColumnCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutputCount_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutputCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{

		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutputHeadings_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutputCount;
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
		std::vector<const char*> CharVector;
		std::vector<std::string>& Components = bv.GetStringVectorRef();
		for (size_t i = 0; i < Components.size(); i++)
		{
			CharVector.push_back(Components[i].c_str());
		}
		bv.SetCharVector(CharVector);
		this->PointerSet.insert(VARS_myself);
		this->UpdateSet.insert(VARS_myself);
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutputRowCount_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutputRowCount;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SolutionVolume_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SolutionVolume;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Time_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Time;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::TimeStep_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::TimeStep;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::CurrentSelectedOutputUserNumber_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::CurrentSelectedOutputUserNumber;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	{
		assert(false);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Porosity_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Porosity;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::GetVar:
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	{
		this->VarExchange.GetDoubleVectorRef() = rm_ptr->GetPorosity();
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->SetPorosity(this->VarExchange.GetDoubleVectorRef());
		bv.SetDoubleVector(this->VarExchange.GetDoubleVectorRef());
		break;
	case VarManager::VAR_TASKS::no_op:
	case VarManager::VAR_TASKS::Info:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Pressure_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Pressure;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::SelectedOutputOn_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::SelectedOutputOn;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
	case VarManager::VAR_TASKS::RMUpdate:
	case VarManager::VAR_TASKS::GetVar:
	{
		int v = rm_ptr->GetSelectedOutputOn();
		bv.SetBVar(v);
		break;
	}
	case VarManager::VAR_TASKS::SetVar:
	{
		bool v = this->VarExchange.GetBVarPtr();
		bv.SetBVar(v);
		rm_ptr->SetSelectedOutputOn(v);
		break;
	}
	case VarManager::VAR_TASKS::no_op:
		break;
	}
	this->VarExchange.CopyScalars(bv);
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
void VarManager::Temperature_Var()
{
	VarManager::VARS VARS_myself = VarManager::VARS::Temperature;
	this->SetCurrentVar(VARS_myself);
	BMIVariant& bv = this->VariantMap[VARS_myself];
	if (!bv.GetInitialized())
	{
		BMIVariant& bv = this->VariantMap[VARS_myself];
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
	case VarManager::VAR_TASKS::UpdateState:
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
	this->SetCurrentVar(VarManager::VARS::NotFound);
}
/// end_
////////////////////////////////

#ifdef SKIP


void VarManager::InputVarNames_var(PhreeqcRM* rm_ptr)
{
	static bool initialized = false;
	static std::vector<std::string> InputVarNames;
	static std::vector<const char*> CharVector;
	//
	PhreeqcRM::BMI_TASKS task_save = rm_ptr->task;
	std::vector<std::string> names = rm_ptr->GetInputVarNames();
	rm_ptr->task = task_save;
	size_t size = 0;
	for (size_t i = 0; i < names.size(); i++)
	{
		if (names[i].size() > size) size = names[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * names.size());
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("InputVarNames", "names", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	rm_ptr->bmi_variant.bmi_var = bv;
	switch (rm_ptr->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		if (!initialized)
		{
			InputVarNames = rm_ptr->GetInputVarNames();
			initialized = true;
			for (size_t i = 0; i < InputVarNames.size(); i++)
			{
				CharVector.push_back(InputVarNames[i].c_str());
			}
		}
		rm_ptr->bmi_variant.SetCharVector(CharVector);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
		rm_ptr->bmi_variant.StringVector = rm_ptr->GetInputVarNames();
		rm_ptr->bmi_variant.bmi_var = bv;
		break;
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->bmi_variant.NotImplemented = true;
		break;
	case VarManager::VAR_TASKS::Info:
		break;
	case VarManager::VAR_TASKS::no_op:
		break;
	}
}
void VarManager::OutputVarNames_var(PhreeqcRM* rm_ptr)
{
	static bool initialized = false;
	static std::vector<std::string> OutputVarNames;
	static std::vector<const char*> CharVector;
	//
	PhreeqcRM::BMI_TASKS task_save = rm_ptr->task;
	std::vector<std::string> names = rm_ptr->GetOutputVarNames();
	rm_ptr->task = task_save;
	size_t size = 0;
	for (size_t i = 0; i < names.size(); i++)
	{
		if (names[i].size() > size) size = names[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * names.size());
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("OutputVarNames", "names", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	rm_ptr->bmi_variant.bmi_var = bv;
	switch (rm_ptr->task)
	{
	case VarManager::VAR_TASKS::GetPtr:
	{
		if (!initialized)
		{
			OutputVarNames = rm_ptr->GetOutputVarNames();
			initialized = true;
			for (size_t i = 0; i < OutputVarNames.size(); i++)
			{
				CharVector.push_back(OutputVarNames[i].c_str());
			}
		}
		rm_ptr->bmi_variant.SetCharVector(CharVector);
		break;
	}
	case VarManager::VAR_TASKS::GetVar:
		rm_ptr->bmi_variant.StringVector = rm_ptr->GetOutputVarNames();
		rm_ptr->bmi_variant.bmi_var = bv;
		break;
	case VarManager::VAR_TASKS::SetVar:
		rm_ptr->bmi_variant.NotImplemented = true;
		break;
	case VarManager::VAR_TASKS::Info:
		break;
	case VarManager::VAR_TASKS::no_op:
		break;
	}
}
#endif


//////////////////
