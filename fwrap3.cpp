#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cref /assume:underscore
//
DLL_EXPORT int createiphreeqc_(void)
{
	return CreateIPhreeqcF();
}
DLL_EXPORT int destroyiphreeqc_(int *id)
{
	return DestroyIPhreeqcF(id);
}
DLL_EXPORT int loaddatabase_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
DLL_EXPORT int loaddatabasestring_(int *id, char *input, unsigned int len)
{
	return LoadDatabaseStringF(id, input, len);
}
DLL_EXPORT int unloaddatabase_(int *id)
{
	return UnLoadDatabaseF(id);
}
DLL_EXPORT void outputlasterror_(int *id)
{
	OutputLastErrorF(id);
}
DLL_EXPORT void outputlastwarning_(int *id)
{
	OutputLastWarningF(id);
}
DLL_EXPORT int accumulateline_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
DLL_EXPORT void setselectedoutputon_(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
DLL_EXPORT void setoutputon_(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
DLL_EXPORT void seterroron_(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
DLL_EXPORT void setlogon_(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
DLL_EXPORT void setdumpon_(int *id, int *dump_on)
{
	SetLogOnF(id, dump_on);
}
DLL_EXPORT void setdumpstringon_(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
DLL_EXPORT int getdumplinecount_(int *id)
{
	return GetDumpLineCountF(id);
}
DLL_EXPORT void getdumpline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
DLL_EXPORT int geterrorlinecount_(int *id)
{
	return GetErrorLineCountF(id);
}
DLL_EXPORT void geterrorline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
DLL_EXPORT int getwarninglinecount_(int *id)
{
	return GetWarningLineCountF(id);
}
DLL_EXPORT void getwarningline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningLineF(id, n, line, line_length);
}
DLL_EXPORT int getcomponentcount_(int *id)
{
	return GetComponentCountF(id);
}
DLL_EXPORT void getcomponent_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
DLL_EXPORT int runaccumulated_(int *id)
{
	return RunAccumulatedF(id);
}
DLL_EXPORT int runfile_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
DLL_EXPORT int runstring_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
DLL_EXPORT void outputlines_(int *id)
{
	OutputLinesF(id);
}
DLL_EXPORT int getselectedoutputrowcount_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
DLL_EXPORT int getselectedoutputcolumncount_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
DLL_EXPORT int getselectedoutputvalue_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif
