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

void ComponentCount_var(BMIPhreeqcRM* brm_ref);
void Components_var(BMIPhreeqcRM* brm_ref);
void Concentrations_var(BMIPhreeqcRM* brm_ref);
void Density_var(BMIPhreeqcRM* brm_ref);
void ErrorString_var(BMIPhreeqcRM* brm_ref);
void FilePrefix_var(BMIPhreeqcRM* brm_ref);
void Gfw_var(BMIPhreeqcRM* brm_ref);
void GridCellCount_var(BMIPhreeqcRM* brm_ref);
void InputVarNames_var(BMIPhreeqcRM* brm_ref);
void NthSelectedOutput_var(BMIPhreeqcRM* brm_ref);
void OutputVarNames_var(BMIPhreeqcRM* brm_ref);
void Saturation_var(BMIPhreeqcRM* brm_ref);
void SelectedOutput_var(BMIPhreeqcRM* brm_ref);
void SelectedOutputColumnCount_var(BMIPhreeqcRM* brm_ref);
void SelectedOutputCount_var(BMIPhreeqcRM* brm_ref);
void SelectedOutputHeadings_var(BMIPhreeqcRM* brm_ref);
void SelectedOutputRowCount_var(BMIPhreeqcRM* brm_ref);
void SolutionVolume_var(BMIPhreeqcRM* brm_ref);
void Time_var(BMIPhreeqcRM* brm_ref);
void TimeStep_var(BMIPhreeqcRM* brm_ref);
void CurrentSelectedOutputUserNumber_var(BMIPhreeqcRM* brm_ref);
void Porosity_var(BMIPhreeqcRM* brm_ref);
void Pressure_var(BMIPhreeqcRM* brm_ref);
void SelectedOutputOn_var(BMIPhreeqcRM* brm_ref);
void Temperature_var(BMIPhreeqcRM* brm_ref);
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
void BMI_Variant::Clear()
{
	bmi_var = BMI_Var("", "", false, false, false, 0, 0);
	b_var = false;
	i_var = -1;
	d_var = -1;
	string_var = "";
	IntVector.clear();
	DoubleVector.clear();
	StringVector.clear();
	NotImplemented = false;
	double* double_ptr = NULL;
	int* int_ptr = NULL;
	bool* bool_ptr = NULL;
	char* char_ptr = NULL;
}
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	std::map<size_t, BMIPhreeqcRM*>::value_type instance(this->GetWorkers()[0]->Get_Index(), this);
	BMIPhreeqcRM::Instances.insert(instance);

	this->var_man = new VarManager((PhreeqcRM*)this);
	
	varfn_map["components"] = &Components_var;
	varfn_map["componentcount"] = &ComponentCount_var;
	varfn_map["concentrations"] = &Concentrations_var;
	varfn_map["density"] = &Density_var;
	varfn_map["errorstring"] = &ErrorString_var;
	varfn_map["fileprefix"] = &FilePrefix_var;
	varfn_map["gfw"] = &Gfw_var;
	varfn_map["gridcellcount"] = &GridCellCount_var;
	varfn_map["inputvarnames"] = &InputVarNames_var;
	varfn_map["nthselectedoutput"] = &NthSelectedOutput_var;
	varfn_map["outputvarnames"] = &OutputVarNames_var;
	varfn_map["saturation"] = &Saturation_var;
	varfn_map["selectedoutput"] = &SelectedOutput_var;
	varfn_map["selectedoutputcolumncount"] = &SelectedOutputColumnCount_var;
	varfn_map["selectedoutputcount"] = &SelectedOutputCount_var;
	varfn_map["selectedoutputheadings"] = &SelectedOutputHeadings_var;
	varfn_map["selectedoutputrowcount"] = &SelectedOutputRowCount_var;
	varfn_map["solutionvolume"] = &SolutionVolume_var;
	varfn_map["time"] = &Time_var;
	varfn_map["timestep"] = &TimeStep_var;
	varfn_map["currentselectedoutputusernumber"] = &CurrentSelectedOutputUserNumber_var;
	varfn_map["porosity"] = &Porosity_var;
	varfn_map["pressure"] = &Pressure_var;
	varfn_map["selectedoutputon"] = &SelectedOutputOn_var;
	varfn_map["temperature"] = &Temperature_var;
	this->task = BMIPhreeqcRM::BMI_TASKS::no_op;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::Update;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
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
		this->bmi_variant.SetSet(false);
		it->second(this);
		if (this->bmi_variant.GetSet()) count++;
	}
	return count;
}
int BMIPhreeqcRM::GetOutputItemCount()
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
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
		this->bmi_variant.SetGet(false);
		it->second(this);
		if (this->bmi_variant.GetGet()) count++;
	}
	return count;
}
int BMIPhreeqcRM::GetPointableItemCount()
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	int count = 0;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		this->bmi_variant.SetHasPtr(false);
		it->second(this);
		if (this->bmi_variant.GetHasPtr()) count++;
	}
	return count;
}
std::vector<std::string> BMIPhreeqcRM::GetInputVarNames()
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
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
		this->bmi_variant.SetSet(false);
		it->second(this);
		if (this->bmi_variant.GetSet())
		{
			names.push_back(this->bmi_variant.GetName());
		}
	}
	return names;
}
std::vector<std::string>  BMIPhreeqcRM::GetOutputVarNames()
{ 
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
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
		this->bmi_variant.SetGet(false);
		it->second(this);
		if (this->bmi_variant.GetGet())
		{
			names.push_back(this->bmi_variant.GetName());
		}
	}
	return names;
}
std::vector<std::string> BMIPhreeqcRM::GetPointableVarNames()
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	std::vector<std::string> names;
	VarFunction_map::iterator it;
	for (it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		this->bmi_variant.SetHasPtr(false);
		it->second(this);
		if (this->bmi_variant.GetHasPtr())
		{
			names.push_back(this->bmi_variant.GetName());
		}
	}
	return names;
}
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn(this);
	if (this->language == "cpp")
	{
		return this->bmi_variant.GetCType();
	} 
	else if (this->language == "F90")
	{
		return this->bmi_variant.GetFType();
	}
	else if (this->language == "Py")
	{
		return this->bmi_variant.GetPType();
	}
	return "Unknown language.";
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn(this);
	return this->bmi_variant.GetUnits();
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn(this);
	return this->bmi_variant.GetItemsize();
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn(this);
	return this->bmi_variant.GetNbytes();
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
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	int Nbytes = this->bmi_variant.GetNbytes();
	int dim = this->bmi_variant.bmi_var.GetDim();
	if (this->bmi_variant.GetCType() == "bool" && dim == 1)
	{
		memcpy(dest, &this->bmi_variant.b_var, Nbytes);
		return;
	}
	if (this->bmi_variant.GetCType() == "int" && dim == 1)
	{
		memcpy(dest, &this->bmi_variant.i_var, Nbytes);
		return;
	}
	if (this->bmi_variant.GetCType() == "double" && dim == 1)
	{
		memcpy(dest, &this->bmi_variant.d_var, Nbytes);
		return;
	}
	//if (this->bmi_variant.GetCType() == "std::string" && dim == 1)
	//{
	//	memcpy(dest, this->bmi_variant.string_var.data(), Nbytes);
	//}
	if (this->bmi_variant.GetCType() == "double" && dim > 1)
	{
		memcpy(dest, this->bmi_variant.DoubleVector.data(), Nbytes);
		return;
	}
	if (this->bmi_variant.GetCType() == "int" && dim > 1)
	{
		memcpy(dest, this->bmi_variant.IntVector.data(), Nbytes);
		return;
	}
	//if (this->bmi_variant.GetCType() == "std::vector<std::string>")
	//{
	//	int itemsize = this->bmi_variant.GetItemsize();
	//	std::stringstream all;
	//	for (size_t i = 0; i < this->bmi_variant.StringVector.size(); i++)
	//	{
	//		all << std::left << std::setfill(' ') << std::setw(itemsize) << this->bmi_variant.StringVector[i];
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
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->bmi_variant.GetCType() == "bool");
	dest = this->bmi_variant.b_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool* dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->bmi_variant.GetCType() == "bool");
	int dim = this->bmi_variant.bmi_var.GetDim();
	int nbytes = this->bmi_variant.bmi_var.GetNbytes();
	if (dim == 1)
	{
		memcpy(dest, &this->bmi_variant.b_var, nbytes);
		return;
	}
	//else if (dim > 1)
	//{
	//	memcpy(dest, this->bmi_variant.BoolVector.data(), nbytes);
	//	return;
	//}
	std::ostringstream oss;
	oss << "BMI GetValue bool* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->bmi_variant.GetCType() == "double");
	dest = this->bmi_variant.d_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double* dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->bmi_variant.GetCType() == "double");
	int dim = this->bmi_variant.bmi_var.GetDim();
	int nbytes = this->bmi_variant.bmi_var.GetNbytes();
	if (dim == 1)
	{
		memcpy(dest, &this->bmi_variant.d_var, nbytes);
		return;
	}
	else if (dim > 1)
	{
		memcpy(dest, this->bmi_variant.DoubleVector.data(), nbytes);
		return;
	}
	std::ostringstream oss;
	oss << "BMI GetValue double* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.i_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int* dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	assert(this->bmi_variant.GetCType() == "int");
	int dim = this->bmi_variant.bmi_var.GetDim();
	int nbytes = this->bmi_variant.bmi_var.GetNbytes();
	if (dim == 1)
	{
		memcpy(dest, &this->bmi_variant.i_var, nbytes);
		return;
	}
	else if (dim > 1)
	{
		memcpy(dest, this->bmi_variant.IntVector.data(), nbytes);
		return;
	}
	std::ostringstream oss;
	oss << "BMI GetValue int* failed for variable " << name << std::endl;
	this->ErrorMessage(oss.str(), true);
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::string& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.string_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<double>& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.DoubleVector;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<int>& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.IntVector;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<std::string>& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.StringVector;
	return;
}
void* BMIPhreeqcRM::GetValuePtr(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetPtr;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return NULL;
	fn(this);
	return this->bmi_variant.GetVoidPtr();
}
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store the variable in bmi_variant
	int Nbytes = this->bmi_variant.GetNbytes();
	int itemsize = this->bmi_variant.GetItemsize();
	int dim = Nbytes / itemsize;
	if (this->bmi_variant.GetCType() == "bool" && dim == 1)
	{
		memcpy(&this->bmi_variant.b_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "int" && dim == 1)
	{
		memcpy(&this->bmi_variant.i_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "double" && dim == 1)
	{
		memcpy(&this->bmi_variant.d_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "std::string")
	{ 
		this->bmi_variant.string_var = (char*)src;
	}
	else if (this->bmi_variant.GetCType() == "double" && dim > 1)
	{
		this->bmi_variant.DoubleVector.resize(dim);
		memcpy(this->bmi_variant.DoubleVector.data(), src, Nbytes);
	}
	else if (this->bmi_variant.GetCType() == "int" && dim > 1)
	{
		this->bmi_variant.IntVector.resize(dim);
		memcpy(&this->bmi_variant.IntVector, src, Nbytes);
	}
	//if (this->bmi_variant.GetType() == "StringVector")
	//{
	//	// Don't think this is possible
	//	//int itemsize = this->GetVarItemsize(name);
	//	//int nbytes = this->GetVarNbytes(name);
	//	//std::stringstream all;
	//	//for (size_t i = 0; i < this->bmi_variant.StringVector.size(); i++)
	//	//{
	//	//	all << std::left << std::setfill(' ') << std::setw(itemsize) << this->bmi_variant.StringVector[i];
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
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetCType("bool");
	this->bmi_variant.b_var = src;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetCType("double");
	this->bmi_variant.d_var = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, int src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetCType("int");
	this->bmi_variant.i_var = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, const std::string src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetCType("std::string");
	this->bmi_variant.string_var = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<double> src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
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
	// Store in bmi_variant
	this->bmi_variant.DoubleVector = src;
	// Set the variable
	this->bmi_variant.SetCType("std::vector<double>");
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<int> src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetCType("std::vector<int>");
	this->bmi_variant.IntVector = src;
	// Set the variable
	fn(this);
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<std::string> src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
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
	// Store in bmi_variant
	this->bmi_variant.SetCType("std::vector<std::string>");
	this->bmi_variant.StringVector = src;
	// Set the variable
	fn(this);
	return;
}
BMIPhreeqcRM::VarFunction BMIPhreeqcRM::GetFn(const std::string name)
{
	this->bmi_variant.Clear();
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

	BMIPhreeqcRM::BMI_TASKS task_save = this->task;
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	it->second(this);
	this->task = task_save;
	if (this->bmi_variant.NotImplemented)
	{
		std::ostringstream oss;
		oss << "Not implemented for variable: " << name;
		this->ErrorMessage(oss.str());
		return NULL;
	}
	if (task == BMIPhreeqcRM::BMI_TASKS::GetVar)
	{
		if (!this->bmi_variant.GetGet() || this->bmi_variant.NotImplemented)
		{
			std::ostringstream oss;
			oss << "Cannot get variable: " << name;
			this->ErrorMessage(oss.str());
			return NULL;
		}
	}
	if (task == BMIPhreeqcRM::BMI_TASKS::SetVar)
	{
		if (!this->bmi_variant.GetSet() || this->bmi_variant.NotImplemented)
		{
			std::ostringstream oss;
			oss << "Cannot set variable: " << name;
			this->ErrorMessage(oss.str());
			return NULL;
		}
	}
	if (task == BMIPhreeqcRM::BMI_TASKS::GetPtr)
	{
		if (this->bmi_variant.NotImplemented)
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

//// Start_var
////////////////////////////////
void Components_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<std::string> Components;
	static std::vector<const char*> CharVector;
	//
	std::vector<std::string> comps = brm_ptr->GetComponents();
	size_t size = 0;
	for (size_t i = 0; i < comps.size(); i++)
	{
		if (comps[i].size() > size) size = comps[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * comps.size());
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Components", "names", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			Components = brm_ptr->GetComponents();
			initialized = true;
			for (size_t i = 0; i < Components.size(); i++)
			{
				CharVector.push_back(Components[i].c_str());
			}
		}
		brm_ptr->bmi_variant.SetCharVector(CharVector);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetComponents();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void ComponentCount_var(BMIPhreeqcRM* brm_ptr)
{
	static int ComponentCount;
	//
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("ComponentCount", "names", false, true, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		ComponentCount = brm_ptr->GetComponentCount();
		brm_ptr->bmi_variant.SetVoidPtr((void*) & ComponentCount);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetComponentCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Concentrations_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Concentrations;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->GetConcentrations(brm_ptr->bmi_variant.DoubleVector);
		assert(Concentrations.size() ==  brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Concentrations.data(), 
			brm_ptr->bmi_variant.DoubleVector.data(), 
			Concentrations.size()*sizeof(double));
		return;
	}
	int Itemsize = (int)sizeof(double);
	int Nbytes = (int)sizeof(double) *
		brm_ptr->GetGridCellCount() * brm_ptr->GetComponentCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Concentrations", "mol L-1", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			brm_ptr->GetConcentrations(Concentrations);
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Concentrations");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Concentrations.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetConcentrations(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		if (initialized) 
		{
			assert(Concentrations.size() == brm_ptr->bmi_variant.DoubleVector.size());
			memcpy(Concentrations.data(),
				brm_ptr->bmi_variant.DoubleVector.data(),
				Concentrations.size() * sizeof(double));
		}
		brm_ptr->SetConcentrations(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Density_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Density;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->GetDensity(brm_ptr->bmi_variant.DoubleVector);
		assert(Density.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Density.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			Density.size() * sizeof(double));
		return;
	}
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Density", "kg L-1", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			brm_ptr->GetDensity(Density);
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Density");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Density.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetDensity(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		// SetDensity does not affect GetDensity
		brm_ptr->SetDensity(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void ErrorString_var(BMIPhreeqcRM* brm_ptr)
{

	int Itemsize = (int)brm_ptr->GetErrorString().size();
	int Nbytes = Itemsize;
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("ErrorString", "error", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::string", "character", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.string_var = brm_ptr->GetErrorString();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void FilePrefix_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = brm_ptr->GetFilePrefix().size();
	int Nbytes = Itemsize;
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("FilePrefix", "name", true, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::string", "character", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.string_var = brm_ptr->GetFilePrefix();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetFilePrefix(brm_ptr->bmi_variant.string_var);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Gfw_var(BMIPhreeqcRM* brm_ptr)
{
	static std::vector<double> Gfw;
	bool initialized = false;
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetComponentCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Gfw", "g mol-1", false, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			Gfw = brm_ptr->GetGfw();
			initialized = true;
		}
		brm_ptr->bmi_variant.SetVoidPtr(Gfw.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetGfw();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void GridCellCount_var(BMIPhreeqcRM* brm_ptr)
{
	static int GridCellCount;
	//
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("GridCellCount", "count", false, true, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		GridCellCount = brm_ptr->GetGridCellCount();
		brm_ptr->bmi_variant.SetVoidPtr((void*) & GridCellCount);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetGridCellCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void InputVarNames_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<std::string> InputVarNames;
	static std::vector<const char*> CharVector;
	//
	BMIPhreeqcRM::BMI_TASKS task_save = brm_ptr->task;
	std::vector<std::string> names = brm_ptr->GetInputVarNames();
	brm_ptr->task = task_save;
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
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			InputVarNames = brm_ptr->GetInputVarNames();
			initialized = true;
			for (size_t i = 0; i < InputVarNames.size(); i++)
			{
				CharVector.push_back(InputVarNames[i].c_str());
			}
		}
		brm_ptr->bmi_variant.SetCharVector(CharVector);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetInputVarNames();
		brm_ptr->bmi_variant.bmi_var = bv;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void NthSelectedOutput_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("NthSelectedOutput", "id", true, false, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetNthSelectedOutput(brm_ptr->bmi_variant.i_var);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void OutputVarNames_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<std::string> OutputVarNames;
	static std::vector<const char*> CharVector;
	//
	BMIPhreeqcRM::BMI_TASKS task_save = brm_ptr->task;
	std::vector<std::string> names = brm_ptr->GetOutputVarNames();
	brm_ptr->task = task_save;
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
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			OutputVarNames = brm_ptr->GetOutputVarNames();
			initialized = true;
			for (size_t i = 0; i < OutputVarNames.size(); i++)
			{
				CharVector.push_back(OutputVarNames[i].c_str());
			}
		}
		brm_ptr->bmi_variant.SetCharVector(CharVector);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetOutputVarNames();
		brm_ptr->bmi_variant.bmi_var = bv;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Saturation_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Saturation;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->GetSaturation(brm_ptr->bmi_variant.DoubleVector);
		assert(Saturation.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Saturation.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			Saturation.size() * sizeof(double));
		return;
	}
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Saturation", "unitless", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			brm_ptr->GetSaturation(Saturation);
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Saturation");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Saturation.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSaturation(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		// SetSaturation does not affect GetSaturation
		brm_ptr->SetSaturation(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutput_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetSelectedOutputRowCount() *
		brm_ptr->GetSelectedOutputColumnCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutput", "user", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSelectedOutput(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutputColumnCount_var(BMIPhreeqcRM* brm_ptr)
{
	static int SelectedOutputColumnCount;
	//
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputColumnCount", "count", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputColumnCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutputCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputCount", "count", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutputHeadings_var(BMIPhreeqcRM* brm_ptr)
{
	std::vector<std::string> headings;
	brm_ptr->GetSelectedOutputHeadings(headings);
	size_t size = 0;
	for (size_t i = 0; i < headings.size(); i++)
	{
		if (headings[i].size() > size) size = headings[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * headings.size());
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputHeadings", "names", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSelectedOutputHeadings(brm_ptr->bmi_variant.StringVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutputRowCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputRowCount", "count", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputRowCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SolutionVolume_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> SolutionVolume;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetSolutionVolume();
		assert(SolutionVolume.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(SolutionVolume.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			SolutionVolume.size() * sizeof(double));
		return;
	}
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("SolutionVolume", "L", false, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			SolutionVolume = brm_ptr->GetSolutionVolume();
			initialized = true;
			brm_ptr->GetUpdateMap().insert("SolutionVolume");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)SolutionVolume.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetSolutionVolume();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Time_var(BMIPhreeqcRM* brm_ptr)
{
	static double Time;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		Time = brm_ptr->GetTime();
		return;
	}
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Time", "s", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "float");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		Time = brm_ptr->GetTime();
		brm_ptr->bmi_variant.SetVoidPtr(&Time);
		brm_ptr->GetUpdateMap().insert("Time");
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.d_var = brm_ptr->GetTime();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		Time = brm_ptr->bmi_variant.d_var;
		brm_ptr->SetTime(brm_ptr->bmi_variant.d_var);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void TimeStep_var(BMIPhreeqcRM* brm_ptr)
{
	static double TimeStep;
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		TimeStep = brm_ptr->GetTimeStep();
		return;
	}
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("TimeStep", "s", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "float");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		TimeStep = brm_ptr->GetTimeStep();
		brm_ptr->bmi_variant.SetVoidPtr((void*) & TimeStep);
		brm_ptr->GetUpdateMap().insert("TimeStep");
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.d_var = brm_ptr->GetTimeStep();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		TimeStep = brm_ptr->bmi_variant.d_var;
		brm_ptr->SetTimeStep(brm_ptr->bmi_variant.d_var);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void CurrentSelectedOutputUserNumber_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("CurrentSelectedOutputUserNumber", "id", false, true, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetCurrentSelectedOutputUserNumber();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Porosity_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Porosity;
	// 
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPorosity();
		assert(Porosity.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Porosity.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			Porosity.size() * sizeof(double));
		return;
	}
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Porosity", "unitless", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			Porosity = brm_ptr->GetPorosity();
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Porosity");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Porosity.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPorosity();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		if (initialized)
		{
			assert(Porosity.size() == brm_ptr->bmi_variant.DoubleVector.size());
			memcpy(Porosity.data(),
				brm_ptr->bmi_variant.DoubleVector.data(),
				Nbytes);
		}
		brm_ptr->SetPorosity(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Pressure_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Pressure;
	// 
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPressure();
		assert(Pressure.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Pressure.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			Pressure.size() * sizeof(double));
		return;
	}
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Pressure", "atm", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			Pressure = brm_ptr->GetPressure();
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Pressure");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Pressure.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPressure();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		if (initialized)
		{
			assert(Pressure.size() == brm_ptr->bmi_variant.DoubleVector.size());
			memcpy(Pressure.data(),
				brm_ptr->bmi_variant.DoubleVector.data(),
				Pressure.size() * sizeof(double));
		}
		brm_ptr->SetPressure(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void SelectedOutputOn_var(BMIPhreeqcRM* brm_ptr)
{
	static bool SelectedOutputOn;
	//
	//
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		SelectedOutputOn = brm_ptr->GetSelectedOutputOn();
		return;
	}
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, std::string units, set, get, ptr, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputOn", "bool", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("bool", "logical", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		SelectedOutputOn = brm_ptr->GetSelectedOutputOn();
		brm_ptr->bmi_variant.SetVoidPtr((void*) & SelectedOutputOn);
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.b_var = brm_ptr->GetSelectedOutputOn();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		SelectedOutputOn = brm_ptr->bmi_variant.b_var;
		brm_ptr->SetSelectedOutputOn(brm_ptr->bmi_variant.b_var);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}
void Temperature_var(BMIPhreeqcRM* brm_ptr)
{
	static bool initialized = false;
	static std::vector<double> Temperature;
	// 
	if (brm_ptr->task == BMIPhreeqcRM::BMI_TASKS::Update)
	{
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetTemperature();
		assert(Temperature.size() == brm_ptr->bmi_variant.DoubleVector.size());
		memcpy(Temperature.data(),
			brm_ptr->bmi_variant.DoubleVector.data(),
			Temperature.size() * sizeof(double));
		return;
	}
	//
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, std::string units, set, get, ptr, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Temperature", "C", true, true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetPtr:
	{
		if (!initialized)
		{
			Temperature = brm_ptr->GetTemperature();
			initialized = true;
			brm_ptr->GetUpdateMap().insert("Temperature");
		}
		brm_ptr->bmi_variant.SetVoidPtr((void*)Temperature.data());
		break;
	}
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetTemperature();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		if (initialized)
		{
			assert(Temperature.size() == brm_ptr->bmi_variant.DoubleVector.size());
			memcpy(Temperature.data(),
				brm_ptr->bmi_variant.DoubleVector.data(),
				Temperature.size() * sizeof(double));
		}
		brm_ptr->SetTemperature(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::Info:
		break;
	case BMIPhreeqcRM::BMI_TASKS::no_op:
		break;
	}
}


//////////////////


