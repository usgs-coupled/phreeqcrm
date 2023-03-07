#ifdef USE_MPI
#include "mpi.h"
#endif
#include "PhreeqcRM.h"
#include "RM_interface_F.h"
#include "IPhreeqcPhastLib.h"
#include "Phreeqc.h"
#include "PHRQ_io.h"
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
static void
rmpadfstring(char* dest, const char* src, unsigned int len)
{
	size_t sofar;

	for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
		*dest++ = *src++;

	while (sofar++ < len)
		*dest++ = ' ';
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetComponentName(int* id, char* chem_name, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns "PhreeqcRM"
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (chem_name != NULL)
		{
			if (*l1 > 0)
			{
				rmpadfstring(chem_name, Reaction_module_ptr->BMI_GetComponentName().c_str(), (unsigned int)*l1);
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetCurrentTime(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current simulation time, in seconds
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTime();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetEndTime(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current simulation time, in seconds
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->BMI_GetEndTime();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int        
RMF_BMI_GetInputItemCount(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be set
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->BMI_GetInputItemCount();
	}
	return IRM_BADINSTANCE;
}
IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetInputVarNames(int* id, char* names, int* l1)
{	
	// Retrieves names of variables that can be set
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector< std::string > VarNames = Reaction_module_ptr->BMI_GetInputVarNames();
		size_t len = 0;
		for (size_t i = 0; i < VarNames.size(); i++)
		{
			if (VarNames[i].size() > len) len = VarNames[i].size();
		}
		size_t total_len = VarNames.size() * len;
		if (*l1 < (int)total_len)
		{
			std::stringstream all;
			for (size_t i = 0; i < VarNames.size(); i++)
			{
				all << std::left << std::setfill(' ') << std::setw(len) << VarNames[i];
				memcpy(names, all.str().c_str(), all.str().size());
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetOutputItemCount(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be retrieved
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->BMI_GetOutputItemCount();
	}
	return IRM_BADINSTANCE;
}
IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetOutputVarNames(int* id, char* names, int* l1)
{
	// Retrieves names of variables that can be retrieved
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector< std::string > VarNames = Reaction_module_ptr->BMI_GetOutputVarNames();
		size_t len = 0;
		for (size_t i = 0; i < VarNames.size(); i++)
		{
			if (VarNames[i].size() > len) len = VarNames[i].size();
		}
		size_t total_len = VarNames.size() * len;
		if (*l1 < (int)total_len)
		{
			std::stringstream all;
			for (size_t i = 0; i < VarNames.size(); i++)
			{
				all << std::left << std::setfill(' ') << std::setw(len) << VarNames[i];
				memcpy(names, all.str().c_str(), all.str().size());
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetTimeStep(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current time step, in seconds
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTimeStep();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetTimeUnits(int* id, char* units, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns time units
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (units != NULL)
		{
			if (*l1 > 0)
			{
				rmpadfstring(units, Reaction_module_ptr->BMI_GetTimeUnits().c_str(), (unsigned int)*l1);
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetValue(int* id, char* var, void* dest)
/* ---------------------------------------------------------------------- */
{
	// Returns value(s) for var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (var != NULL)
		{
			std::string str_var = var;
			size_t end = str_var.find_last_not_of(' ');
			str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);				
			Reaction_module_ptr->BMI_GetValue(str_var, dest);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetVarItemsize(int* id, char* var)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of bytes needed for one item
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		return Reaction_module_ptr->BMI_GetVarItemsize(str_var);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetVarNbytes(int* id, char* var)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number total number of bytes needed for the buffer
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		return Reaction_module_ptr->BMI_GetVarNbytes(str_var);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetVarType(int* id, char* var, char* vtype, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns type of variable var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		std::string type_cpp = Reaction_module_ptr->BMI_GetVarType(str_var);
		if (*l1 > 0)
		{
			rmpadfstring(vtype, type_cpp.c_str(), (unsigned int)*l1);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetVarUnits(int* id, char* var, char* units, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		std::string units_cpp = Reaction_module_ptr->BMI_GetVarUnits(str_var);
		if (*l1 > 0)
		{
			rmpadfstring(units, units_cpp.c_str(), (unsigned int)*l1);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
//IRM_DLL_EXPORT IRM_RESULT RMF_BMI_Initialize(int* id, char* config_file);
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_Initialize(int* id, char* config_file)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->BMI_Initialize(config_file);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_SetValue(int* id, char* var, void* src)
/* ---------------------------------------------------------------------- */
{
	// Returns value(s) for var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (var != NULL)
		{
			std::string str_var = var;
			size_t end = str_var.find_last_not_of(' ');
			str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
			Reaction_module_ptr->BMI_SetValue(str_var, src);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
// IRM_RESULT RMF_BMI_Update(int* id);
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_Update(int* id)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	PhreeqcRM* Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->BMI_Update();
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetGridCellCountYAML(const char* config)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	return GetGridCellCountYAML(config);
}