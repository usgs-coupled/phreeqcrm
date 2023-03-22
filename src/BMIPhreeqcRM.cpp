#include "BMIPhreeqcRM.h"
#include "BMI_Var.h"
void FilePrefix_var(BMIPhreeqcRM& brm_ref);
BMIPhreeqcRM::BMIPhreeqcRM(int nxyz, int nthreads) :
PhreeqcRM(nxyz, nthreads) 
{
	varfn_map["fileprefix"] = &FilePrefix_var;
	this->task = BMIPhreeqcRM::BMI_TASKS::no_op;
}
int BMIPhreeqcRM::GetInputItemCount()
{
	this->bmi_variant.Clear();
	int count = 0;

	for (auto it = varfn_map.begin(); it != varfn_map.end(); it++)
	{
		it->second;

	}
	return count;
}
std::string BMIPhreeqcRM::GetVarType(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Type;
	std::string name_lc = name;
	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
	auto it = varfn_map.find(name_lc);
	if (it == varfn_map.end())
	{
		return "";
	}
	it->second;
	return this->bmi_variant.GetType();
}
std::string BMIPhreeqcRM::GetVarUnits(const std::string name)
{
	this->task = BMIPhreeqcRM::BMI_TASKS::Units;
	std::string name_lc = name;
	std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
	auto it = varfn_map.find(name_lc);
	if (it == varfn_map.end())
	{
		return "";
	}
	return this->bmi_variant.GetUnits();
}

void FilePrefix_var(BMIPhreeqcRM& brm_ref)
{
	//name, type, std::string units, set, get
	BMI_Var bv = BMI_Var("FilePrefix", "character", "name", true, true);
	brm_ref.bmi_variant.bmi_var = bv;
	switch (brm_ref.task)
	{
	case BMIPhreeqcRM::BMI_TASKS::count_input:
	case BMIPhreeqcRM::BMI_TASKS::count_output:
	case BMIPhreeqcRM::BMI_TASKS::Type:
	case BMIPhreeqcRM::BMI_TASKS::Units:
		break;
	case BMIPhreeqcRM::BMI_TASKS::Itemsize:
		brm_ref.bmi_variant.Itemsize = sizeof(char);
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
void BMI_Variant::Clear()
{
	bmi_var = BMI_Var("","","",false,false);
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
