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
BMIPhreeqcRM::CleanupBmiModuleInstances(void)
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
BMIPhreeqcRM::CreateBmiModule(int nxyz, MP_TYPE nthreads)
/* ---------------------------------------------------------------------- */
{
	//_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_crtBreakAlloc = 5144;
	int n = IRM_OUTOFMEMORY;
	try
	{
		BMIPhreeqcRM* bmirm_ptr = new BMIPhreeqcRM(nxyz, nthreads);
		if (bmirm_ptr)
		{
			n = (int)bmirm_ptr->GetWorkers()[0]->Get_Index();
			BMIPhreeqcRM::Instances[n] = bmirm_ptr;
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
BMIPhreeqcRM::DestroyBmiModule(int id)
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
	bmi_var = BMI_Var("", "", false, false, 0, 0);
	//Nbytes = 0;
	//Itemsize = 0;
	b_var = false;
	i_var = -1;
	d_var = -1;
	string_var = "";
	IntVector.clear();
	DoubleVector.clear();
	StringVector.clear();
	NotImplemented = false;
}
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	std::map<size_t, BMIPhreeqcRM*>::value_type instance(this->GetWorkers()[0]->Get_Index(), this);
	BMIPhreeqcRM::Instances.insert(instance);
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
}
// Model control functions.
void BMIPhreeqcRM::Initialize(std::string config_file)
{
	this->InitializeYAML(config_file);
}
void BMIPhreeqcRM::Update()
{
	this->RunCells();
	this->SetTime(this->GetTime() + this->GetTimeStep());
}
void BMIPhreeqcRM::UpdateUntil(double time)
{
	double time_step = time - this->GetTime();
	if (time_step >= 0)
	{
		this->SetTimeStep(time_step);
		this->RunCells();
		this->SetTime(time);
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
		this->bmi_variant.SetGet(false);
		it->second(this);
		if (this->bmi_variant.GetGet()) count++;
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
			names.push_back("InputVarNames");
			continue;
		}
		if (it->first == "outputvarnames")
		{
			names.push_back("OutputVarNames");
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
			continue;
		}
		if (it->first == "outputvarnames")
		{
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
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn(this);
	return this->bmi_variant.GetCType();
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
	if (this->bmi_variant.GetCType() == "bool")
	{
		memcpy(dest, &this->bmi_variant.b_var, Nbytes);
	}
	if (this->bmi_variant.GetCType() == "int")
	{
		memcpy(dest, &this->bmi_variant.i_var, Nbytes);
	}
	if (this->bmi_variant.GetCType() == "double")
	{
		memcpy(dest, &this->bmi_variant.d_var, Nbytes);
	}
	if (this->bmi_variant.GetCType() == "std::string")
	{
		memcpy(dest, this->bmi_variant.string_var.data(), Nbytes);
	}
	if (this->bmi_variant.GetCType() == "std::vector<double>")
	{
		memcpy(dest, this->bmi_variant.DoubleVector.data(), Nbytes);
	}
	if (this->bmi_variant.GetCType() == "std::vector<int>")
	{
		memcpy(dest, this->bmi_variant.IntVector.data(), Nbytes);
	}
	if (this->bmi_variant.GetCType() == "std::vector<std::string>")
	{
		int itemsize = this->bmi_variant.GetItemsize();
		std::stringstream all;
		for (size_t i = 0; i < this->bmi_variant.StringVector.size(); i++)
		{
			all << std::left << std::setfill(' ') << std::setw(itemsize) << this->bmi_variant.StringVector[i];
		}
		memcpy(dest, all.str().c_str(), all.str().size());
	}
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, bool& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.b_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::GetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn(this);
	dest = this->bmi_variant.d_var;
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
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store the variable in bmi_variant
	int Nbytes = this->bmi_variant.GetNbytes();
	int itemsize = this->bmi_variant.GetItemsize();
	int dim = Nbytes / itemsize;
	if (this->bmi_variant.GetCType() == "bool")
	{
		memcpy(&this->bmi_variant.b_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "int")
	{
		memcpy(&this->bmi_variant.i_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "double")
	{
		memcpy(&this->bmi_variant.d_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetCType() == "std::string")
	{ 
		this->bmi_variant.string_var = (char*)src;
	}
	else if (this->bmi_variant.GetCType() == "std::vector<double>")
	{
		this->bmi_variant.DoubleVector.resize(dim);
		memcpy(this->bmi_variant.DoubleVector.data(), src, Nbytes);
	}
	else if (this->bmi_variant.GetCType() == "std::vector<int>")
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
		this->ErrorMessage("Unknown input to SetValue", true);
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
	return it->second;
}
////////////////////////////////
void Components_var(BMIPhreeqcRM* brm_ptr)
{
	std::vector<std::string> comps = brm_ptr->GetComponents();
	size_t size = 0;
	for (size_t i = 0; i < comps.size(); i++)
	{
		if (comps[i].size() > size) size = comps[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * comps.size());
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Components", "names", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetComponents();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void ComponentCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("ComponentCount", "names", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetComponentCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void Concentrations_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(double);
	int Nbytes = (int)sizeof(double) *
		brm_ptr->GetGridCellCount() * brm_ptr->GetComponentCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Concentrations", "mol L-1", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:,:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		 brm_ptr->GetConcentrations(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		 brm_ptr->SetConcentrations(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}
void Density_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Density", "kg L-1", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetDensity(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetDensity(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}
void ErrorString_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)brm_ptr->GetErrorString().size();
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("ErrorString", "error", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::string", "character", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.string_var = brm_ptr->GetErrorString();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void FilePrefix_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = brm_ptr->GetFilePrefix().size();
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("FilePrefix", "name", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::string", "character", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.string_var = brm_ptr->GetFilePrefix();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetFilePrefix(brm_ptr->bmi_variant.string_var);
		break;
	}
}
void Gfw_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetComponentCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Gfw", "g mol-1", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetGfw();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void GridCellCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("GridCellCount", "count", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetGridCellCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void InputVarNames_var(BMIPhreeqcRM* brm_ptr)
{
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
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("InputVarNames", "names", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetInputVarNames();
		brm_ptr->bmi_variant.bmi_var = bv;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void NthSelectedOutput_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("NthSelectedOutput", "id", true, false, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetNthSelectedOutput(brm_ptr->bmi_variant.i_var);
		break;
	}
}
void OutputVarNames_var(BMIPhreeqcRM* brm_ptr)
{
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
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("OutputVarNames", "names", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.StringVector = brm_ptr->GetOutputVarNames();
		brm_ptr->bmi_variant.bmi_var = bv;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void Saturation_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Saturation", "unitless", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSaturation(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetSaturation(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}
void SelectedOutput_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetSelectedOutputRowCount() * 
		brm_ptr->GetSelectedOutputColumnCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutput", "user", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:,:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSelectedOutput(brm_ptr->bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputColumnCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputColumnCount", "count", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputColumnCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputCount", "count", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
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
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputHeadings", "names", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<std::string>", "character(len=:),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->GetSelectedOutputHeadings(brm_ptr->bmi_variant.StringVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputRowCount_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputRowCount", "count", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetSelectedOutputRowCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void SolutionVolume_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("SolutionVolume", "L", false, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetSolutionVolume();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void Time_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Time", "s", true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "float");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.d_var = brm_ptr->GetTime();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetTime(brm_ptr->bmi_variant.d_var);
		break;
	}
}
void TimeStep_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("TimeStep", "s", true, true, Nbytes, Itemsize);
	bv.SetTypes("double", "real(kind=8)", "float");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.d_var = brm_ptr->GetTimeStep();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetTimeStep(brm_ptr->bmi_variant.d_var);
		break;
	}
}
void CurrentSelectedOutputUserNumber_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("CurrentSelectedOutputUserNumber", "id", false, true, Nbytes, Itemsize);
	bv.SetTypes("int", "integer", "int");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.i_var = brm_ptr->GetCurrentSelectedOutputUserNumber();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->bmi_variant.NotImplemented = true;
		break;
	}
}
void Porosity_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Porosity", "unitless", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPorosity();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetPorosity(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}
void Pressure_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Pressure", "atm", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetPressure();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetPressure(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}
void SelectedOutputOn_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputOn", "flag", true, true, Nbytes, Itemsize);
	bv.SetTypes("bool", "logical", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.b_var = brm_ptr->GetSelectedOutputOn();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetSelectedOutputOn(brm_ptr->bmi_variant.b_var);
		break;
	}
}
void Temperature_var(BMIPhreeqcRM* brm_ptr)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ptr->GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Temperature", "C", true, true, Nbytes, Itemsize);
	bv.SetTypes("std::vector<double>", "real(kind=8),allocatable,dimension(:)", "");
	brm_ptr->bmi_variant.bmi_var = bv;
	switch (brm_ptr->task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ptr->bmi_variant.DoubleVector = brm_ptr->GetTemperature();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ptr->SetTemperature(brm_ptr->bmi_variant.DoubleVector);
		break;
	}
}

//////////////////


