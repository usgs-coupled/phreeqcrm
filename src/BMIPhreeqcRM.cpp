#include "BMIPhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>
void Components_var(BMIPhreeqcRM& brm_ref);
void ComponentCount_var(BMIPhreeqcRM& brm_ref);
void Concentrations_var(BMIPhreeqcRM& brm_ref);
void Density_var(BMIPhreeqcRM& brm_ref);
void ErrorString_var(BMIPhreeqcRM& brm_ref);
void FilePrefix_var(BMIPhreeqcRM& brm_ref);
void Gfw_var(BMIPhreeqcRM& brm_ref);
void GridCellCount_var(BMIPhreeqcRM& brm_ref);
void NthSelectedOutput_var(BMIPhreeqcRM& brm_ref);
void Saturation_var(BMIPhreeqcRM& brm_ref);
void SelectedOutput_var(BMIPhreeqcRM& brm_ref);
void SelectedOutputColumnCount_var(BMIPhreeqcRM& brm_ref);
void SelectedOutputCount_var(BMIPhreeqcRM& brm_ref);
void SelectedOutputHeadings_var(BMIPhreeqcRM& brm_ref);
void SelectedOutputRowCount_var(BMIPhreeqcRM& brm_ref);
void SolutionVolume_var(BMIPhreeqcRM& brm_ref);
void Time_var(BMIPhreeqcRM& brm_ref);
void TimeStep_var(BMIPhreeqcRM& brm_ref);
void CurrentSelectedOutputUserNumber_var(BMIPhreeqcRM& brm_ref);
void Porosity_var(BMIPhreeqcRM& brm_ref);
void Pressure_var(BMIPhreeqcRM& brm_ref);
void SelectedOutputOn_var(BMIPhreeqcRM& brm_ref);
void Temperature_var(BMIPhreeqcRM& brm_ref);
void InputVarNames_var(BMIPhreeqcRM& brm_ref);
void OutputVarNames_var(BMIPhreeqcRM& brm_ref);


void BMI_Variant::Clear()
{
	bmi_var = BMI_Var("", "", "", false, false, 0, 0);
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
void FilePrefix_var(BMIPhreeqcRM& brm_ref);
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	varfn_map["Components"] = &Components_var;
	varfn_map["ComponentCount"] = &ComponentCount_var;
	varfn_map["Concentrations"] = &Concentrations_var;
	varfn_map["Density"] = &Density_var;
	varfn_map["ErrorString"] = &ErrorString_var;
	varfn_map["FilePrefix"] = &FilePrefix_var;
	varfn_map["Gfw"] = &Gfw_var;
	varfn_map["GridCellCount"] = &GridCellCount_var;
	varfn_map["NthSelectedOutput"] = &NthSelectedOutput_var;
	varfn_map["Saturation"] = &Saturation_var;
	varfn_map["SelectedOutput"] = &SelectedOutput_var;
	varfn_map["SelectedOutputColumnCount"] = &SelectedOutputColumnCount_var;
	varfn_map["SelectedOutputCount"] = &SelectedOutputCount_var;
	varfn_map["SelectedOutputHeadings"] = &SelectedOutputHeadings_var;
	varfn_map["SelectedOutputRowCount"] = &SelectedOutputRowCount_var;
	varfn_map["SolutionVolume"] = &SolutionVolume_var;
	varfn_map["Time"] = &Time_var;
	varfn_map["TimeStep"] = &TimeStep_var;
	varfn_map["CurrentSelectedOutputUserNumber"] = &CurrentSelectedOutputUserNumber_var;
	varfn_map["Porosity"] = &Porosity_var;
	varfn_map["Pressure"] = &Pressure_var;
	varfn_map["SelectedOutputOn"] = &SelectedOutputOn_var;
	varfn_map["Temperature"] = &Temperature_var;
	varfn_map["InputVarNames"] = &InputVarNames_var;
	varfn_map["OutputVarNames"] = &OutputVarNames_var;
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
		it->second;
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
		it->second;
		if (this->bmi_variant.GetGet()) count++;
	}
	return count;
}
std::vector<std::string> BMIPhreeqcRM::GetInputVarNames()
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	std::vector<std::string> names;
	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		this->bmi_variant.SetSet(false);
		it->second;
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
		this->bmi_variant.SetGet(false);
		it->second;
		if (this->bmi_variant.GetGet())
		{
			names.push_back(this->bmi_variant.GetName());
		}
	}
	return names;
}
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{

	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn;
	return this->bmi_variant.GetType();
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn;
	return this->bmi_variant.GetUnits();
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn;
	return this->bmi_variant.GetItemsize();
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Info;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn;
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
double BMIPhreeqcRM::GetTimeStep()
{
	return this->GetTimeStep();
}
void BMIPhreeqcRM::GetValue(const std::string name, void* dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	int Nbytes = this->GetVarNbytes(name);
	if (this->bmi_variant.GetType() == "bool")
	{
		memcpy(dest, &this->bmi_variant.b_var, Nbytes);
	}
	if (this->bmi_variant.GetType() == "int")
	{
		memcpy(dest, &this->bmi_variant.i_var, Nbytes);
	}
	if (this->bmi_variant.GetType() == "double")
	{
		memcpy(dest, &this->bmi_variant.d_var, Nbytes);
	}
	if (this->bmi_variant.GetType() == "string")
	{
		memcpy(dest, this->bmi_variant.string_var.data(), Nbytes);
	}
	if (this->bmi_variant.GetType() == "double,1d")
	{
		memcpy(dest, this->bmi_variant.DoubleVector.data(), Nbytes);
	}
	if (this->bmi_variant.GetType() == "int,1d")
	{
		memcpy(dest, &this->bmi_variant.IntVector, Nbytes);
	}
	if (this->bmi_variant.GetType() == "StringVector")
	{
		int itemsize = this->GetVarItemsize(name);
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
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.b_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, double& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.d_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, int& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.i_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::string& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.string_var;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<double>& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.DoubleVector;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<int>& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.IntVector;
	return;
}
void BMIPhreeqcRM::GetValue(const std::string name, std::vector<std::string>& dest)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	fn;
	dest = this->bmi_variant.StringVector;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, void* src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store the variable in bmi_variant
	int Nbytes = this->bmi_variant.GetNbytes();
	int itemsize = this->bmi_variant.GetNbytes();
	int dim = Nbytes / itemsize;
	if (this->bmi_variant.GetType() == "bool")
	{
		memcpy(&this->bmi_variant.b_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetType() == "int")
	{
		memcpy(&this->bmi_variant.i_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetType() == "double")
	{
		memcpy(&this->bmi_variant.d_var, src, Nbytes);
	} 
	else if (this->bmi_variant.GetType() == "string")
	{ 
		this->bmi_variant.string_var = (char*)src;
	}
	else if (this->bmi_variant.GetType() == "double,1d")
	{
		this->bmi_variant.DoubleVector.resize(dim);
		memcpy(this->bmi_variant.DoubleVector.data(), src, Nbytes);
	}
	else if (this->bmi_variant.GetType() == "int,1d")
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
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, bool src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetType("bool");
	this->bmi_variant.b_var = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, double src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetType("double");
	this->bmi_variant.d_var = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, int src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetType("int");
	this->bmi_variant.i_var = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}void BMIPhreeqcRM::SetValue(const std::string name, std::string src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetType("string");
	this->bmi_variant.string_var = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<double> src)
{
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
	this->bmi_variant.SetType("double,1d");
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<int> src)
{
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return;
	// Store in bmi_variant
	this->bmi_variant.SetType("int,1d");
	this->bmi_variant.IntVector = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
	return;
}
void BMIPhreeqcRM::SetValue(const std::string name, std::vector<std::string> src)
{
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
	this->bmi_variant.SetType("string,1d");
	this->bmi_variant.StringVector = src;
	// Set the variable
	this->task = BMIPhreeqcRM::BMI_TASKS::SetVar;
	fn;
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
	it->second;
	this->task = task_save;
	if (this->bmi_variant.NotImplemented)
	{
		std::ostringstream oss;
		oss << "Not implemented for variable: " << name;
		this->ErrorMessage(oss.str());
		return NULL;
	}
	if (task_save == BMIPhreeqcRM::BMI_TASKS::GetVar)
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
		if (!this->bmi_variant.GetGet() || this->bmi_variant.NotImplemented)
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
void Components_var(BMIPhreeqcRM& brm_ref)
{
	std::vector<std::string> comps = brm_ref.GetComponents();
	size_t size = 0;
	for (size_t i = 0; i < comps.size(); i++)
	{
		if (comps[i].size() > size) size = comps[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * comps.size());
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Components", "character,1d", "names", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.StringVector = brm_ref.GetComponents();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void ComponentCount_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("ComponentCount", "integer", "names", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetComponentCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void Concentrations_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(double);
	int Nbytes = (int)sizeof(double) *
		brm_ref.GetGridCellCount() * brm_ref.GetComponentCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("Concentrations", "double,2d", "mol L-1", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		 brm_ref.GetConcentrations(brm_ref.bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		 brm_ref.SetConcentrations(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void Density_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Density", "double,1d", "kg L-1", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.GetDensity(brm_ref.bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetDensity(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void ErrorString_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)brm_ref.GetErrorString().size();
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("ErrorString", "character", "error", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.string_var = brm_ref.GetErrorString();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void FilePrefix_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = brm_ref.GetFilePrefix().size();
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("FilePrefix", "character", "name", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.string_var = brm_ref.GetFilePrefix();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetFilePrefix(brm_ref.bmi_variant.string_var);
		break;
	}
}
void Gfw_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetComponentCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Gfw", "double,1d", "g mol-1", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.DoubleVector = brm_ref.GetGfw();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void GridCellCount_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("GridCellCount", "integer", "count", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetGridCellCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void NthSelectedOutput_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("NthSelectedOutput", "integer", "id", true, false, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetNthSelectedOutput(brm_ref.bmi_variant.i_var);
		break;
	}
}
void Saturation_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Saturation", "double,1d", "unitless", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.GetSaturation(brm_ref.bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetSaturation(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void SelectedOutput_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetSelectedOutputRowCount() * 
		brm_ref.GetSelectedOutputColumnCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutput", "double,2d", "user", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.GetSelectedOutput(brm_ref.bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputColumnCount_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputColumnCount", "integer", "count", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetSelectedOutputColumnCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputCount_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputCount", "integer", "count", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetSelectedOutputCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputHeadings_var(BMIPhreeqcRM& brm_ref)
{
	std::vector<std::string> headings;
	brm_ref.GetSelectedOutputHeadings(headings);
	size_t size = 0;
	for (size_t i = 0; i < headings.size(); i++)
	{
		if (headings[i].size() > size) size = headings[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * headings.size());
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputHeadings", "character,1d", "names", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.GetSelectedOutputHeadings(brm_ref.bmi_variant.StringVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void SelectedOutputRowCount_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputRowCount", "integer", "count", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetSelectedOutputRowCount();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void SolutionVolume_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("SolutionVolume", "double,1d", "L", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.DoubleVector = brm_ref.GetSolutionVolume();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void Time_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Time", "double", "s", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.d_var = brm_ref.GetTime();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetTime(brm_ref.bmi_variant.d_var);
		break;
	}
}
void TimeStep_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize;
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("TimeStep", "double", "s", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.d_var = brm_ref.GetTimeStep();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetTimeStep(brm_ref.bmi_variant.d_var);
		break;
	}
}
void CurrentSelectedOutputUserNumber_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("CurrentSelectedOutputUserNumber", "integer", "id", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.i_var = brm_ref.GetCurrentSelectedOutputUserNumber();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
void Porosity_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Porosity", "double,1d", "unitless", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.DoubleVector = brm_ref.GetPorosity();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetPorosity(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void Pressure_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Pressure", "double,1d", "atm", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.DoubleVector = brm_ref.GetPressure();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetPressure(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void SelectedOutputOn_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = (int)sizeof(int);
	int Nbytes = (int)sizeof(int);
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("SelectedOutputOn", "logical", "flag", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.b_var = brm_ref.GetSelectedOutputOn();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetSelectedOutputOn(brm_ref.bmi_variant.b_var);
		break;
	}
}
void Temperature_var(BMIPhreeqcRM& brm_ref)
{
	int Itemsize = sizeof(double);
	int Nbytes = Itemsize * brm_ref.GetGridCellCount();
	//name, type, std::string units, set, get, Nbytes, Itemsize  
	BMI_Var bv = BMI_Var("Temperature", "double,1d", "C", true, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.DoubleVector = brm_ref.GetTemperature();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetTemperature(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void InputVarNames_var(BMIPhreeqcRM& brm_ref)
{
	std::vector<std::string> headings = brm_ref.GetInputVarNames();
	size_t size = 0;
	for (size_t i = 0; i < headings.size(); i++)
	{
		if (headings[i].size() > size) size = headings[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * headings.size());
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("InputVarNames", "character,1d", "string", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.StringVector = brm_ref.GetInputVarNames();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}

void OutputVarNames_var(BMIPhreeqcRM& brm_ref)
{
	std::vector<std::string> headings = brm_ref.GetInputVarNames();
	size_t size = 0;
	for (size_t i = 0; i < headings.size(); i++)
	{
		if (headings[i].size() > size) size = headings[i].size();
	}
	int Itemsize = (int)size;
	int Nbytes = (int)(size * headings.size());
	//name, type, std::string units, set, get, Nbytes, Itemsize
	BMI_Var bv = BMI_Var("OutputVarNames", "character,1d", "string", false, true, Nbytes, Itemsize);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.StringVector = brm_ref.GetOutputVarNames();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.bmi_variant.NotImplemented = true;
		break;
	}
}
//////////////////


