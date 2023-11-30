///*! @file BMI_interface_F.h
//	@brief C/Fortran Documentation
//*/
#if !defined(BMI_INTERFACE_C_H_INCLUDED)
#define BMI_INTERFACE_C_H_INCLUDED
#include "IrmResult.h"
#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_AddOutputVars(int id, char* option, char* def);
	//IRM_DLL_EXPORT int        RM_BMI_Create(int* nxyz, int* nthreads = nullptr);
	//IRM_DLL_EXPORT int        RM_BMI_Create_default();
	IRM_DLL_EXPORT int		  RM_BMI_Destroy(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetComponentName(int id, char* chem_name, int l1);
	IRM_DLL_EXPORT double     RM_BMI_GetCurrentTime(int id);
	IRM_DLL_EXPORT double     RM_BMI_GetEndTime(int id);
	IRM_DLL_EXPORT int        RM_BMI_GetInputItemCount(int id);
	IRM_DLL_EXPORT int        RM_BMI_GetInputVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetInputVarName(int id, char* name, int i);
	//IRM_DLL_EXPORT int        RM_BMI_GetNames(int id, const char* type, char* dest);
	//IRM_DLL_EXPORT int        RM_BMI_GetNamesSize(int id, const char* type, int* dest);
	IRM_DLL_EXPORT int        RM_BMI_GetOutputItemCount(int id);
	IRM_DLL_EXPORT int        RM_BMI_GetOutputVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetOutputVarName(int id, char* name, int i);
	IRM_DLL_EXPORT int		  RM_BMI_GetPointableItemCount(int id);
	IRM_DLL_EXPORT int        RM_BMI_GetPointableVarNamesSize(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetPointableVarName(int id, char* name, int i);
	IRM_DLL_EXPORT double     RM_BMI_GetTimeStep(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetTimeUnits(int id, char* units, int* l1);
	// GetValue
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValue_int(int id, char* var, int* dest);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValue_double(int id, char* var, double* dest);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValue_char(int id, char* var, char* dest);
	// GetValuePtr
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValuePtr_int(int id, char* var, int** dest);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValuePtr_double(int id, char* var, double** dest);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetValuePtr_char(int id, char* var, char** dest);
	
	IRM_DLL_EXPORT int        RM_BMI_GetVarItemsize(int id, char* name);
	IRM_DLL_EXPORT int        RM_BMI_GetVarNbytes(int id, char* name);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetVarType(int id, char* name, char* vtype, int l1);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_GetVarUnits(int id, char* name, char* units, int l1);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_Initialize(int id, char* config_file);
	// SetValue
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_SetValue_char(int id, char* name, char* src);
	//IRM_DLL_EXPORT IRM_RESULT RM_BMI_SetValue(int id, char* name, double src);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_SetValue_double(int id, char* name, double* src);
	//IRM_DLL_EXPORT IRM_RESULT RM_BMI_SetValue(int id, char* name, int src);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_SetValue_int(int id, char* name, int* src);

	IRM_DLL_EXPORT IRM_RESULT RM_BMI_Update(int id);
	IRM_DLL_EXPORT IRM_RESULT RM_BMI_UpdateUntil(int id, double time);

#if defined(__cplusplus)
}
#endif

#endif // BMI_INTERFACE_C_H_INCLUDED
