#ifdef USE_MPI
#include "mpi.h"
#endif
#include "BMIVariant.h"
#include "BMIPhreeqcRM.h"
//#include "RM_interface_F.h"
#include "BMI_interface_F.h"
#include "IPhreeqcPhastLib.h"
#include "Phreeqc.h"
#include "PHRQ_io.h"
#include "BMIVariant.h"
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
#ifdef USE_MPI
/* ---------------------------------------------------------------------- */
int
RM_BMI_Create(int* nxyz, int* nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::CreateBMIModule(*nxyz, MPI_Comm_f2c(*nthreads));
}
#else
/* ---------------------------------------------------------------------- */
int
RM_BMI_Create(int* nxyz, int* nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::CreateBMIModule(*nxyz, *nthreads);
}
#endif
/* ---------------------------------------------------------------------- */
int
RM_BMI_Create_default()
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::CreateBMIModule();
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_Destroy(int* id)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::DestroyBMIModule(*id);
}
IRM_RESULT        RMF_BMI_AddOutputVars(int* id, char* option_in, char* def_in)
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::string option = option_in;
		std::string def = def_in;
		bmirm_ptr->AddOutputVars(option, def);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetComponentName(int* id, char* chem_name, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns "PhreeqcRM"
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		if (chem_name != NULL)
		{
			if (*l1 > 0)
			{
				rmpadfstring(chem_name, bmirm_ptr->GetComponentName().c_str(), (unsigned int)*l1);
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
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetTime();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetEndTime(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current simulation time, in seconds
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetEndTime();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int        
RMF_BMI_GetInputItemCount(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be set
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetInputItemCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetOutputItemCount(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be retrieved
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetOutputItemCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetPointableItemCount(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be retrieved
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetPointableItemCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetNamesSize(int* id, const char* type, int* dest)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be retrieved
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::vector<std::string> v;
		std::string stype = type;
		if (stype == "inputvarnames")
		{
			v = bmirm_ptr->GetInputVarNames();
		}
		if (stype == "outputvarnames")
		{
			v = bmirm_ptr->GetOutputVarNames();
		}
		if (stype == "pointablevarnames")
		{
			v = bmirm_ptr->GetPointableVarNames();
		}
		int size = 0;
		for (size_t i = 0; i < v.size(); i++)
		{
			if (v[i].size() > size) size = v[i].size();
		}
		memcpy(dest, &size, sizeof(int));
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetNames(int* id, const char* type, char* dest)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number of variables that can be retrieved
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	std::string stype = type;
	if (bmirm_ptr)
	{
		std::vector<std::string> v;
		if (stype == "inputvarnames")
		{
			v = bmirm_ptr->GetInputVarNames();
		}
		if (stype == "outputvarnames")
		{
			v = bmirm_ptr->GetOutputVarNames();
		}
		if (stype == "pointablevarnames")
		{
			v = bmirm_ptr->GetPointableVarNames();
		}
		int size = 0;
		for (size_t i = 0; i < v.size(); i++)
		{
			if (v[i].size() > size) size = v[i].size();
		}
		int itemsize = size;
		std::stringstream all;
		for (size_t i = 0; i < v.size(); i++)
		{
			all << std::left << std::setfill(' ') << std::setw(itemsize) << v[i];
		}
		memcpy( dest, all.str().data(), all.str().size());
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetStartTime(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current time step, in seconds
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetStartTime();
	}
	return -1.0;
}
/* ---------------------------------------------------------------------- */
double
RMF_BMI_GetTimeStep(int* id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current time step, in seconds
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		return bmirm_ptr->GetTimeStep();
	}
	return -1.0;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetTimeUnits(int* id, char* units, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns time units
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		if (units != NULL)
		{
			if (*l1 > 0)
			{
				rmpadfstring(units, bmirm_ptr->GetTimeUnits().c_str(), (unsigned int)*l1);
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
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		if (var != NULL)
		{
			std::string str_var = var;
			size_t end = str_var.find_last_not_of(' ');
			str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
			std::string type = bmirm_ptr->GetVarType(var);
			bmirm_ptr->GetValue(str_var, dest);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetValuePtr(int* id, char* var, void*& dest)
/* ---------------------------------------------------------------------- */
{
	// Returns value(s) for var
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		if (var != NULL)
		{
			std::string str_var = var;
			size_t end = str_var.find_last_not_of(' ');
			str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
			std::string type = bmirm_ptr->GetVarType(var);
			dest = bmirm_ptr->GetValuePtr(str_var);
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
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		return bmirm_ptr->GetVarItemsize(str_var);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_BMI_GetVarNbytes(int* id, char* var)
/* ---------------------------------------------------------------------- */
{
	// Retrieves number total number of bytes needed for the buffer
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		return bmirm_ptr->GetVarNbytes(str_var);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_GetVarType(int* id, char* var, char* vtype, int* l1)
/* ---------------------------------------------------------------------- */
{
	// Returns type of variable var
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		std::string type_cpp = bmirm_ptr->GetVarType(str_var);
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
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		std::string str_var = var;
		size_t end = str_var.find_last_not_of(' ');
		str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
		std::string units_cpp = bmirm_ptr->GetVarUnits(str_var);
		if (*l1 > 0)
		{
			rmpadfstring(units, units_cpp.c_str(), (unsigned int)*l1);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
#ifdef USE_YAML
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_Initialize(int* id, char* config_file)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		bmirm_ptr->Initialize(config_file);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
#endif
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_SetValue(int* id, char* var, void* src)
/* ---------------------------------------------------------------------- */
{
	// Returns value(s) for var
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		if (var != NULL)
		{
			std::string str_var = var;
			size_t end = str_var.find_last_not_of(' ');
			str_var = (end == std::string::npos) ? "" : str_var.substr(0, end + 1);
			bmirm_ptr->SetValue(str_var, src);
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
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		bmirm_ptr->Update();
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_BMI_UpdateUntil(int* id, double* time)
/* ---------------------------------------------------------------------- */
{
	// Returns units of variable var
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(*id);
	if (bmirm_ptr)
	{
		bmirm_ptr->UpdateUntil(*time);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}