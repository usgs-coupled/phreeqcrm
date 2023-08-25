///*! @file BMI_interface_F.h
//	@brief C/Fortran Documentation
//*/
#if !defined(BMI_INTERFACE_F_H_INCLUDED)
#define BMI_INTERFACE_F_H_INCLUDED
#include "IrmResult.h"
#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_AddOutputVars(int* id, char* option, char* def);
	IRM_DLL_EXPORT int        RM_BMI_Create(int* nxyz, int* nthreads = nullptr);
	IRM_DLL_EXPORT int        RM_BMI_Create_default();
	IRM_DLL_EXPORT int		  RMF_BMI_Destroy(int* id);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetComponentName(int* id, char* chem_name, int* l1);
	IRM_DLL_EXPORT double     RMF_BMI_GetCurrentTime(int* id);
	IRM_DLL_EXPORT double     RMF_BMI_GetEndTime(int* id);
	IRM_DLL_EXPORT int        RMF_BMI_GetInputItemCount(int* id);
	IRM_DLL_EXPORT int        RMF_BMI_GetNames(int* id, const char* type, char* dest);
	IRM_DLL_EXPORT int        RMF_BMI_GetNamesSize(int* id, const char* type, int* dest);
	IRM_DLL_EXPORT int        RMF_BMI_GetOutputItemCount(int* id);
	IRM_DLL_EXPORT int		  RMF_BMI_GetPointableItemCount(int* id);
	IRM_DLL_EXPORT double     RMF_BMI_GetTimeStep(int* id);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetTimeUnits(int* id, char* units, int* l1);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetValue(int* id, char* name, void* dest);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetValuePtr(int* id, char* var, void*& dest);
	IRM_DLL_EXPORT int        RMF_BMI_GetVarItemsize(int* id, char* name);
	IRM_DLL_EXPORT int        RMF_BMI_GetVarNbytes(int* id, char* name);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetVarType(int* id, char* name, char* vtype, int* l1);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_GetVarUnits(int* id, char* name, char* units, int* l1);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_Initialize(int* id, char* config_file);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_SetValue(int* id, char* name, void* src);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_Update(int* id);
	IRM_DLL_EXPORT IRM_RESULT RMF_BMI_UpdateUntil(int* id, double* time);

#if defined(__cplusplus)
}
#endif

#endif // BMI_INTERFACE_F_H_INCLUDED
