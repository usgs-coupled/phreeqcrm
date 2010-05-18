#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cref /assume:underscore
//
DLL_EXPORT int  accumulateline_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
// AddError
DLL_EXPORT int  clearaccumulatedlines_(int *id)
{
	return ClearAccumulatedLinesF(id);
}
DLL_EXPORT int  createiphreeqc_(void)
{
	return CreateIPhreeqcF();
}
DLL_EXPORT int  destroyiphreeqc_(int *id)
{
	return DestroyIPhreeqcF(id);
}
DLL_EXPORT void getcomponent_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
DLL_EXPORT int  getcomponentcount_(int *id)
{
	return GetComponentCountF(id);
}
DLL_EXPORT int  getdumpfileon_(int *id)
{
	return GetDumpFileOnF(id);
}
// GetDumpString
DLL_EXPORT void getdumpstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  getdumpstringlinecount_(int *id)
{
	return GetDumpStringLineCountF(id);
}
DLL_EXPORT int  getdumpstringon_(int *id)
{
	return GetDumpStringOnF(id);
}
DLL_EXPORT int  geterrorfileon_(int *id)
{
	return GetErrorFileOnF(id);
}
// GetErrorString
DLL_EXPORT void geterrorstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  geterrorstringlinecount_(int *id)
{
	return GetErrorStringLineCountF(id);
}
DLL_EXPORT int  getlogfileon_(int *id)
{
	return GetLogFileOnF(id);
}
DLL_EXPORT int  getoutputfileon_(int *id)
{
	return GetOutputFileOnF(id);
}
DLL_EXPORT int  getselectedoutputcolumncount_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
DLL_EXPORT int  getselectedoutputfileon_(int *id)
{
	return GetSelectedOutputFileOnF(id);
}
DLL_EXPORT int  getselectedoutputrowcount_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
DLL_EXPORT int  getselectedoutputvalue_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
// GetWarningString
DLL_EXPORT void getwarningstringline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningStringLineF(id, n, line, line_length);
}
DLL_EXPORT int  getwarningstringlinecount_(int *id)
{
	return GetWarningStringLineCountF(id);
}
DLL_EXPORT int  loaddatabase_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
DLL_EXPORT int  loaddatabasestring_(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
DLL_EXPORT void outputaccumulatedlines_(int *id)
{
	OutputAccumulatedLinesF(id);
}
DLL_EXPORT void outputerrorstring_(int *id)
{
	OutputErrorStringF(id);
}
DLL_EXPORT void outputwarningstring_(int *id)
{
	OutputWarningStringF(id);
}
DLL_EXPORT int  runaccumulated_(int *id)
{
	return RunAccumulatedF(id);
}
DLL_EXPORT int  runfile_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
DLL_EXPORT int  runstring_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
DLL_EXPORT void setdumpfileon_(int *id, int *dump_on)
{
	SetDumpFileOnF(id, dump_on);
}
DLL_EXPORT void setdumpstringon_(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
DLL_EXPORT void seterrorfileon_(int *id, int *error_on)
{
	SetErrorFileOnF(id, error_on);
}
DLL_EXPORT void setlogfileon_(int *id, int *log_on)
{
	SetLogFileOnF(id, log_on);
}
DLL_EXPORT void setoutputfileon_(int *id, int *output_on)
{
	SetOutputFileOnF(id, output_on);
}
DLL_EXPORT void setselectedoutputfileon_(int *id, int *selected_on)
{
	SetSelOutFileOnF(id, selected_on);
}
DLL_EXPORT int  unloaddatabase_(int *id)
{
	return UnLoadDatabaseF(id);
}

#if defined(__cplusplus)
}
#endif

#endif
