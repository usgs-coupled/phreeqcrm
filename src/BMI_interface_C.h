///*! @file BMI_interface_C.h
//	@brief C/Fortran Documentation
//*/
#if !defined(BMI_INTERFACE_C_H_INCLUDED)
#define BMI_INTERFACE_C_H_INCLUDED
#include "IrmResult.h"
#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif
	IRM_DLL_EXPORT IRM_RESULT BMI_AddOutputVars(int id, char* option, char* def);
	IRM_DLL_EXPORT int        BMI_Create(int nxyz, int nthreads);
	//IRM_DLL_EXPORT int        BMI_Create_default();
	IRM_DLL_EXPORT int		  BMI_Destroy(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetComponentName(int id, char* chem_name, int l1);
	IRM_DLL_EXPORT double     BMI_GetCurrentTime(int id);
	IRM_DLL_EXPORT double     BMI_GetEndTime(int id);
	IRM_DLL_EXPORT int        BMI_GetInputItemCount(int id);
	IRM_DLL_EXPORT int        BMI_GetInputVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetInputVarName(int id, char* name, int i);
	//IRM_DLL_EXPORT int        BMI_GetNames(int id, const char* type, char* dest);
	//IRM_DLL_EXPORT int        BMI_GetNamesSize(int id, const char* type, int* dest);
	IRM_DLL_EXPORT int        BMI_GetOutputItemCount(int id);
	IRM_DLL_EXPORT int        BMI_GetOutputVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetOutputVarName(int id, char* name, int i);
	IRM_DLL_EXPORT int		  BMI_GetPointableItemCount(int id);
	IRM_DLL_EXPORT int        BMI_GetPointableVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetPointableVarName(int id, char* name, int i);
	IRM_DLL_EXPORT double     BMI_GetTimeStep(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetTimeUnits(int id, char* units, int* l1);
	// GetValue
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueInt(int id, char* var, int* dest);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueDouble(int id, char* var, double* dest);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetValueChar(int id, char* var, char* dest);
	// GetValuePtr
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrInt(int id, char* var, int** dest);
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrDouble(int id, char* var, double** dest);
	//IRM_DLL_EXPORT IRM_RESULT BMI_GetValuePtrChar(int id, char* var, char** dest);
	IRM_DLL_EXPORT void* BMI_GetValuePtr(int id, char* var);

	
	IRM_DLL_EXPORT int        BMI_GetVarItemsize(int id, char* name);
	IRM_DLL_EXPORT int        BMI_GetVarNbytes(int id, char* name);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetVarType(int id, char* name, char* vtype, int l1);
	IRM_DLL_EXPORT IRM_RESULT BMI_GetVarUnits(int id, char* name, char* units, int l1);
	IRM_DLL_EXPORT IRM_RESULT BMI_Initialize(int id, char* config_file);
	// SetValue
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueChar(int id, char* name, const char* src);
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueDouble(int id, char* name, double src);
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueDoubleArray(int id, char* name, double* src);
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueInt(int id, char* name, int src);
	IRM_DLL_EXPORT IRM_RESULT BMI_SetValueIntArray(int id, char* name, int* src);

	IRM_DLL_EXPORT IRM_RESULT BMI_Update(int id);
	IRM_DLL_EXPORT IRM_RESULT BMI_UpdateUntil(int id, double time);

#if defined(__cplusplus)
}
#endif

#endif // BMI_INTERFACE_C_H_INCLUDED
