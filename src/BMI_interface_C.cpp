#ifdef USE_MPI
#include "mpi.h"
#endif
#include "BMIVariant.h"
#include "BMIPhreeqcRM.h"
#include "BMI_interface_C.h"
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
#ifdef USE_MPI
/* ---------------------------------------------------------------------- */
int
BMI_Create(int* nxyz, int* nthreads)
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
BMI_Create_default()
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::CreateBMIModule();
}
/* ---------------------------------------------------------------------- */
int
BMI_Create(int nxyz, int nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::CreateBMIModule(nxyz, nthreads);
}
#endif
int
BMI_Destroy(int id)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return BMIPhreeqcRM::DestroyBMIModule(id);
}
IRM_RESULT        
BMI_AddOutputVars(int id, char* option_in, char* def_in)
{
	return RMF_BMI_AddOutputVars(&id, option_in, def_in);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetComponentName(int id, char* chem_name, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetComponentName(&id, chem_name, &l1);
}
/* ---------------------------------------------------------------------- */
double
BMI_GetCurrentTime(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetCurrentTime(&id);
}
/* ---------------------------------------------------------------------- */
double
BMI_GetEndTime(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetCurrentTime(&id);
}

/* ---------------------------------------------------------------------- */
int        
BMI_GetInputItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetInputItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetInputVarName(int id, char* name, int i)
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
BMI_GetInputVarNamesSize(int id)
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
BMI_GetOutputItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetOutputItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetOutputVarName(int id, char* name, int i)
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
BMI_GetOutputVarNamesSize(int id)
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
BMI_GetPointableItemCount(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetPointableItemCount(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetPointableVarName(int id, char* name, int i)
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
BMI_GetPointableVarNamesSize(int id)
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
BMI_GetTimeStep(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetTimeStep(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetTimeUnits(int id, char* units, int* l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetTimeUnits(&id, units, l1);
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetValueChar(int id, char* var, char* dest)
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
BMI_GetValueDouble(int id, char* var, double* dest)
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
BMI_GetValueInt(int id, char* var, int* dest)
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
void*
BMI_GetValuePtr(int id, char* var)
/* ---------------------------------------------------------------------- */
{
	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
	if (bmirm_ptr)
	{
		std::string name = var;
		return bmirm_ptr->GetValuePtr(name);
	}
	return NULL;
}

///* ---------------------------------------------------------------------- */
//IRM_RESULT
//BMI_GetValuePtrChar(int id, char* var, char** dest)
///* ---------------------------------------------------------------------- */
//{
//	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
//	if (bmirm_ptr)
//	{
//		std::string name = var;
//		*dest = (char*)bmirm_ptr->GetValuePtr(name);
//		return IRM_OK;
//	}
//	return IRM_BADINSTANCE;
//}
///* ---------------------------------------------------------------------- */
//double*
//BMI_GetValuePtrDouble(int id, char* var)
///* ---------------------------------------------------------------------- */
//{
//	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
//	if (bmirm_ptr)
//	{
//		std::string name = var;
//		return (double*)bmirm_ptr->GetValuePtr(name);
//	}
//	return NULL;
//}
///* ---------------------------------------------------------------------- */
//int*
//BMI_GetValuePtrInt(int id, char* var)
///* ---------------------------------------------------------------------- */
//{
//	BMIPhreeqcRM* bmirm_ptr = BMIPhreeqcRM::GetInstance(id);
//	if (bmirm_ptr)
//	{
//		std::string name = var;
//		return (int*)bmirm_ptr->GetValuePtr(name);
//	}
//	return NULL;
//}

/* ---------------------------------------------------------------------- */
int
BMI_GetVarItemsize(int id, char* var)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarItemsize(&id, var);
}
/* ---------------------------------------------------------------------- */
int
BMI_GetVarNbytes(int id, char* var)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarNbytes(&id, var);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetVarType(int id, char* var, char* vtype, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarType(&id, var, vtype, &l1);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_GetVarUnits(int id, char* var, char* units, int l1)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_GetVarUnits(&id, var, units, &l1);
}
#ifdef USE_YAML
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_Initialize(int id, char* config_file)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_Initialize(&id, config_file);
}
#endif


/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_SetValueChar(int id, char* var, char* src)
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
//BMI_SetValue(int id, char* var, double src)
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
BMI_SetValueDouble(int id, char* var, double* src)
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
//BMI_SetValue(int id, char* var, int src)
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
BMI_SetValueInt(int id, char* var, int* src)
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
BMI_Update(int id)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_Update(&id);
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
BMI_UpdateUntil(int id, double time)
/* ---------------------------------------------------------------------- */
{
	return RMF_BMI_UpdateUntil(&id, &time);
}