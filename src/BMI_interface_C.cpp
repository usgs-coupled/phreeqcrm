#ifdef USE_MPI
#include "mpi.h"
#endif
#include "BMIVariant.h"
#include "BMIPhreeqcRM.h"
#include "BMI_interface_F.h"
#include "IPhreeqcPhastLib.h"
#include "Phreeqc.h"
#include "PHRQ_io.h"
#include "BMIVariant.h"
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
//static void
//rmpadfstring(char* dest, const char* src, unsigned int len)
//{
//	size_t sofar;
//
//	for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
//		*dest++ = *src++;
//
//	while (sofar++ < len)
//		*dest++ = ' ';
//}
//#ifdef USE_MPI
///* ---------------------------------------------------------------------- */
//int
//RM_BMI_Create(int* nxyz, int* nthreads)
///* ---------------------------------------------------------------------- */
//{
//	//
//	// Creates reaction module, called by root and MPI workers
//	//
//	return BMIPhreeqcRM::CreateBMIModule(*nxyz, MPI_Comm_f2c(*nthreads));
//}
//#else
///* ---------------------------------------------------------------------- */
//int
//RM_BMI_Create_default()
///* ---------------------------------------------------------------------- */
//{
//	//
//	// Creates reaction module, called by root and MPI workers
//	//
//	return BMIPhreeqcRM::CreateBMIModule();
//}
///* ---------------------------------------------------------------------- */
//int
//RM_BMI_Create(int* nxyz, int* nthreads)
///* ---------------------------------------------------------------------- */
//{
//	//
//	// Creates reaction module, called by root and MPI workers
//	//
//	return BMIPhreeqcRM::CreateBMIModule(*nxyz, *nthreads);
//}
//#endif
int
RM_BMI_Destroy(int id)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::DestroyBMIModule(id);
}
IRM_RESULT        
RM_BMI_AddOutputVars(int id, char* option_in, char* def_in)
{
	return RMF_BMI_AddOutputVars(&id, option_in, def_in);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetComponentName(int id, char* chem_name, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetComponentName(&id, chem_name, &l1);
}
/* ---------------------------------------------------------------------- */
double
RM_BMI_GetCurrentTime(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetCurrentTime(&id);
}
/* ---------------------------------------------------------------------- */
double
RM_BMI_GetEndTime(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetCurrentTime(&id);
}

/* ---------------------------------------------------------------------- */
int        
RM_BMI_GetInputItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetInputItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetInputVarName(int id, char* name, int i)
/* ---------------------------------------------------------------------- */
{
	// Returns ith output variable name
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::vector<std::string> names = bmirm_ptr->GetInputVarNames();
		if ((size_t)i < names.size())
		{
			memcpy(name, names[i].c_str(), names[i].size());
			name[names[i].size()] = '\0';
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RM_BMI_GetInputVarNamesSize(int id)
/* ---------------------------------------------------------------------- */
{
	IRM_RESULT status;
	int l;
	status = (IRM_RESULT)RMF_BMI_GetNamesSize(&id, "inputvarnames", &l);
	l = l + 1;
	if (status != IRM_OK) l = -1;
	return l;
}


/* ---------------------------------------------------------------------- */
int
RM_BMI_GetOutputItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetOutputItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetOutputVarName(int id, char* name, int i)
/* ---------------------------------------------------------------------- */
{
	// Returns ith output variable name
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::vector<std::string> names = bmirm_ptr->GetOutputVarNames();
		if ((size_t)i < names.size())
		{
			memcpy(name, names[i].c_str(), names[i].size());
			name[names[i].size()] = '\0';
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int        
RM_BMI_GetOutputVarNamesSize(int id)
/* ---------------------------------------------------------------------- */
{
	IRM_RESULT status;
	int l;
	status = (IRM_RESULT) RMF_BMI_GetNamesSize(&id, "outputvarnames", &l);
	l = l + 1;
	if (status != IRM_OK) l = -1;
	return l;
}


/* ---------------------------------------------------------------------- */
int
RM_BMI_GetPointableItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetPointableItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetPointableVarName(int id, char* name, int i)
/* ---------------------------------------------------------------------- */
{
	// Returns ith output variable name
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::vector<std::string> names = bmirm_ptr->GetPointableVarNames();
		if ((size_t)i < names.size())
		{
			memcpy(name, names[i].c_str(), names[i].size());
			name[names[i].size()] = '\0';
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RM_BMI_GetPointableVarNamesSize(int id)
/* ---------------------------------------------------------------------- */
{
	IRM_RESULT status;
	int l;
	status = (IRM_RESULT)RMF_BMI_GetNamesSize(&id, "pointablevarnames", &l);
	l = l + 1;
	if (status != IRM_OK) l = -1;
	return l;
}



/* ---------------------------------------------------------------------- */
double
RM_BMI_GetTimeStep(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetTimeStep(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetTimeUnits(int id, char* units, int* l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetTimeUnits(&id, units, l1);
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValue_char(int id, char* var, char* dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValue_double(int id, char* var, double* dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValue_int(int id, char* var, int* dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}


/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValuePtr_char(int id, char* var, char** dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
//}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValuePtr_double(int id, char* var, double** dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetValuePtr_int(int id, char* var, int** dest)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->GetValue(name, dest);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RM_BMI_GetVarItemsize(int id, char* var)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarItemsize(&id, var);
}
/* ---------------------------------------------------------------------- */
int
RM_BMI_GetVarNbytes(int id, char* var)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarNbytes(&id, var);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetVarType(int id, char* var, char* vtype, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarType(&id, var, vtype, &l1);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_GetVarUnits(int id, char* var, char* units, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarUnits(&id, var, units, &l1);
}
#ifdef USE_YAML
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_Initialize(int id, char* config_file)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_Initialize(&id, config_file);
}
#endif


/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_SetValue_char(int id, char* var, char* src)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->SetValue(name, src);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
///* ---------------------------------------------------------------------- */
//IRM_RESULT
//RM_BMI_SetValue(int id, char* var, double src)
///* ---------------------------------------------------------------------- */
//{
//	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
//	if (bmirm_ptr)
//	{
//		std::string name = var;
//		bmirm_ptr->SetValue(name, src);
//		return IRM_OK;
//	}
//	return IRM_BADINSTANCE;
//}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_SetValue_double(int id, char* var, double* src)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->SetValue(name, src);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
///* ---------------------------------------------------------------------- */
//IRM_RESULT
//RM_BMI_SetValue(int id, char* var, int src)
///* ---------------------------------------------------------------------- */
//{
//	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
//	if (bmirm_ptr)
//	{
//		std::string name = var;
//		bmirm_ptr->SetValue(name, src);
//		return IRM_OK;
//	}
//	return IRM_BADINSTANCE;
//}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_SetValue_int(int id, char* var, int* src)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		bmirm_ptr->SetValue(name, src);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_Update(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_Update(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_BMI_UpdateUntil(int id, double time)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_UpdateUntil(&id, &time);
}