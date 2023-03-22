#include "BMIPhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>
void BMI_Variant::Clear()
{
	bmi_var = BMI_Var("", "", "", false, false);
	Nbytes = 0;
	Itemsize = 0;
	b_var = false;
	i_var = -1;
	d_var = -1;
	string_var = "";
	IntVector.clear();
	DoubleVector.clear();
	StringVector.clear();
	NotImplemented = true;
}
void FilePrefix_var(BMIPhreeqcRM& brm_ref);
// Constructor
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	varfn_map["fileprefix"] = &FilePrefix_var;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::CountInput;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::CountOutput;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::CountInput;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::CountOutput;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::Units;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return "";
	fn;
	return this->bmi_variant.GetUnits();
}

int BMIPhreeqcRM::GetVarItemsize(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Itemsize;
	BMIPhreeqcRM::VarFunction fn = GetFn(name);
	if (fn == NULL) return 0;
	fn;
	return this->bmi_variant.GetItemsize();
}

int BMIPhreeqcRM::GetVarNbytes(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Nbytes;
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
	this->task = BMIPhreeqcRM::BMI_TASKS::Type;
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
	//name, type, std::string units, set, get
	BMI_Var bv = BMI_Var("Components", "character,1d", "names", false, true);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::CountInput:
	case BMIPhreeqcRM::BMI_TASKS::CountOutput:
	case BMIPhreeqcRM::BMI_TASKS::Type:
	case BMIPhreeqcRM::BMI_TASKS::Units:
		break;
	case BMIPhreeqcRM::BMI_TASKS::Itemsize:
	case BMIPhreeqcRM::BMI_TASKS::Nbytes:
	{
		std::vector<std::string> comps = brm_ref.GetComponents();
		size_t size = 0;
		for (size_t i = 0; i < comps.size(); i++)
		{
			if (comps[i].size() > size) size = comps[i].size();
		}
		brm_ref.bmi_variant.Itemsize = (int)size;
		brm_ref.bmi_variant.Nbytes = (int)(size * comps.size());
		break;
	}
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
	//name, type, std::string units, set, get
	BMI_Var bv = BMI_Var("ComponentCount", "integer", "names", false, true);
	brm_ref.bmi_variant.bmi_var = bv;
	brm_ref.bmi_variant.Itemsize = (int)sizeof(int);
	brm_ref.bmi_variant.Nbytes = (int)sizeof(int);
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::CountInput:
	case BMIPhreeqcRM::BMI_TASKS::CountOutput:
	case BMIPhreeqcRM::BMI_TASKS::Type:
	case BMIPhreeqcRM::BMI_TASKS::Units:
	case BMIPhreeqcRM::BMI_TASKS::Itemsize:
	case BMIPhreeqcRM::BMI_TASKS::Nbytes:
		break;
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
	//name, type, std::string units, set, get
	BMI_Var bv = BMI_Var("Concentrations", "double,2d", "mol L-1", true, true);
	brm_ref.bmi_variant.bmi_var = bv;
	brm_ref.bmi_variant.Itemsize = (int)sizeof(double);
	brm_ref.bmi_variant.Nbytes = (int)sizeof(double) * 
		brm_ref.GetGridCellCount() * brm_ref.GetComponentCount();
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::CountInput:
	case BMIPhreeqcRM::BMI_TASKS::CountOutput:
	case BMIPhreeqcRM::BMI_TASKS::Type:
	case BMIPhreeqcRM::BMI_TASKS::Units:
	case BMIPhreeqcRM::BMI_TASKS::Itemsize:
	case BMIPhreeqcRM::BMI_TASKS::Nbytes:
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		 brm_ref.GetConcentrations(brm_ref.bmi_variant.DoubleVector);
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		 brm_ref.SetConcentrations(brm_ref.bmi_variant.DoubleVector);
		break;
	}
}
void FilePrefix_var(BMIPhreeqcRM& brm_ref)
{
	//name, type, std::string units, set, get
	BMI_Var bv = BMI_Var("FilePrefix", "character", "name", true, true);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::CountInput:
	case BMIPhreeqcRM::BMI_TASKS::CountOutput:
	case BMIPhreeqcRM::BMI_TASKS::Type:
	case BMIPhreeqcRM::BMI_TASKS::Units:
		break;
	case BMIPhreeqcRM::BMI_TASKS::Itemsize:
		brm_ref.bmi_variant.Itemsize = brm_ref.GetFilePrefix().size();
	case BMIPhreeqcRM::BMI_TASKS::Nbytes:
		brm_ref.bmi_variant.Nbytes = brm_ref.GetFilePrefix().size();
		break;
	case BMIPhreeqcRM::BMI_TASKS::GetVar:
		brm_ref.bmi_variant.string_var = brm_ref.GetFilePrefix();
		break;
	case BMIPhreeqcRM::BMI_TASKS::SetVar:
		brm_ref.SetFilePrefix(brm_ref.bmi_variant.string_var);
		break;
	}
}

