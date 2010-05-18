#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)


#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
DLL_EXPORT int  ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
// AddError
DLL_EXPORT int  CLEARACCUMULATEDLINES(int *id)
{
	return ClearAccumulatedLinesF(id);
}
DLL_EXPORT int  CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
DLL_EXPORT int  DESTROYIPHREEQC(int *id)
{
	return DestroyIPhreeqcF(id);
}
DLL_EXPORT void GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
DLL_EXPORT int  GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
DLL_EXPORT int  GETDUMPFILEON(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
DLL_EXPORT void GETDUMPSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  GETDUMPSTRINGLINECOUNT(int *id)
{
	return GetDumpStringLineCountF(id);
}
DLL_EXPORT int  GETDUMPSTRINGON(int *id)
{
	return GetDumpStringOnF(id);
}
DLL_EXPORT int  GETERRORFILEON(int *id)
{
	return GetErrorFileOnF(id);
}
DLL_EXPORT void GETERRORSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  GETERRORSTRINGLINECOUNT(int *id)
{
	return GetErrorStringLineCountF(id);
}
DLL_EXPORT int  GETLOGFILEON(int *id)
{
	return GetLogFileOnF(id);
}
DLL_EXPORT int  GETOUTPUTFILEON(int *id)
{
	return GetOutputFileOnF(id);
}
DLL_EXPORT int  GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
DLL_EXPORT int  GETSELECTEDOUTPUTFILEON(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
DLL_EXPORT int  GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
DLL_EXPORT int  GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
DLL_EXPORT void GETWARNINGSTRINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  GETWARNINGSTRINGLINECOUNT(int *id)
{
	return GetWarningStringLineCountF(id);
}
DLL_EXPORT int  LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
DLL_EXPORT int  LOADDATABASESTRING(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
DLL_EXPORT void OUTPUTACCUMULATEDLINES(int *id)
{
	OutputAccumulatedLinesF(id);
}
DLL_EXPORT void OUTPUTERRORSTRING(int *id)
{
	OutputErrorStringF(id);
}
DLL_EXPORT void OUTPUTWARNINGSTRING(int *id)
{
	OutputWarningStringF(id);
}
DLL_EXPORT int  RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
DLL_EXPORT int  RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
DLL_EXPORT int  RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
DLL_EXPORT void SETDUMPFILEON(int *id, int *dump_on)
{
	SetDumpFileOnF(id, dump_on);
}
DLL_EXPORT void SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
DLL_EXPORT void SETERRORFILEON(int *id, int *error_on)
{
	SetErrorFileOnF(id, error_on);
}
DLL_EXPORT void SETLOGFILEON(int *id, int *log_on)
{
	SetLogFileOnF(id, log_on);
}
DLL_EXPORT void SETOUTPUTFILEON(int *id, int *output_on)
{
	SetOutputFileOnF(id, output_on);
}
DLL_EXPORT void SETSELECTEDOUTPUTFILEON(int *id, int *selected_on)
{
	SetSelOutFileOnF(id, selected_on);
}
DLL_EXPORT int  UNLOADDATABASE(int *id)
{
	return UnLoadDatabaseF(id);
}

#if defined(__cplusplus)
}
#endif

#endif