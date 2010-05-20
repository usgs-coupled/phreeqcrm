#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)


#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
IPQ_DLL_EXPORT int  ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
IPQ_DLL_EXPORT int  ADDERROR(int *id, char *error_msg, unsigned int len)
{
	return AddErrorF(id, error_msg, len);
}
IPQ_DLL_EXPORT int  ADDWARNING(int *id, char *warn_msg, unsigned int len)
{
	return AddWarningF(id, warn_msg, len);
}
IPQ_DLL_EXPORT int  CLEARACCUMULATEDLINES(int *id)
{
	return ClearAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT int  CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
IPQ_DLL_EXPORT int  DESTROYIPHREEQC(int *id)
{
	return DestroyIPhreeqcF(id);
}
IPQ_DLL_EXPORT void GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
IPQ_DLL_EXPORT int  GETDUMPFILEON(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
IPQ_DLL_EXPORT void GETDUMPSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGLINECOUNT(int *id)
{
	return GetDumpStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETDUMPSTRINGON(int *id)
{
	return GetDumpStringOnF(id);
}
IPQ_DLL_EXPORT int  GETERRORFILEON(int *id)
{
	return GetErrorFileOnF(id);
}
IPQ_DLL_EXPORT void GETERRORSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETERRORSTRINGLINECOUNT(int *id)
{
	return GetErrorStringLineCountF(id);
}
IPQ_DLL_EXPORT int  GETLOGFILEON(int *id)
{
	return GetLogFileOnF(id);
}
IPQ_DLL_EXPORT int  GETOUTPUTFILEON(int *id)
{
	return GetOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTFILEON(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
IPQ_DLL_EXPORT int  GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
IPQ_DLL_EXPORT void GETWARNINGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
IPQ_DLL_EXPORT int  GETWARNINGSTRINGLINECOUNT(int *id)
{
	return GetWarningStringLineCountF(id);
}
IPQ_DLL_EXPORT int  LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
IPQ_DLL_EXPORT int  LOADDATABASESTRING(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
IPQ_DLL_EXPORT void OUTPUTACCUMULATEDLINES(int *id)
{
	OutputAccumulatedLinesF(id);
}
IPQ_DLL_EXPORT void OUTPUTERRORSTRING(int *id)
{
	OutputErrorStringF(id);
}
IPQ_DLL_EXPORT void OUTPUTWARNINGSTRING(int *id)
{
	OutputWarningStringF(id);
}
IPQ_DLL_EXPORT int  RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
IPQ_DLL_EXPORT int  RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
IPQ_DLL_EXPORT int  RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
IPQ_DLL_EXPORT void SETDUMPFILEON(int *id, int *dump_on)
{
	SetDumpFileOnF(id, dump_on);
}
IPQ_DLL_EXPORT void SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
IPQ_DLL_EXPORT void SETERRORFILEON(int *id, int *error_on)
{
	SetErrorFileOnF(id, error_on);
}
IPQ_DLL_EXPORT void SETLOGFILEON(int *id, int *log_on)
{
	SetLogFileOnF(id, log_on);
}
IPQ_DLL_EXPORT void SETOUTPUTFILEON(int *id, int *output_on)
{
	SetOutputFileOnF(id, output_on);
}
IPQ_DLL_EXPORT void SETSELECTEDOUTPUTFILEON(int *id, int *selected_on)
{
	SetSelOutFileOnF(id, selected_on);
}
IPQ_DLL_EXPORT int  UNLOADDATABASE(int *id)
{
	return UnLoadDatabaseF(id);
}

#if defined(__cplusplus)
}
#endif

#endif