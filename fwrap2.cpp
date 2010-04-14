#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
DLL_EXPORT int CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
DLL_EXPORT int LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
DLL_EXPORT void OUTPUTLASTERROR(int *id)
{
	OutputLastErrorF(id);
}
DLL_EXPORT int ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
DLL_EXPORT void SETSELECTEDOUTPUTON(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
DLL_EXPORT void SETOUTPUTON(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
DLL_EXPORT void SETERRORON(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
DLL_EXPORT void SETLOGON(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
DLL_EXPORT void SETDUMPON(int *id, int *dump_on)
{
	SetDumpOnF(id, dump_on);
}
DLL_EXPORT void SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
DLL_EXPORT int GETDUMPLINECOUNT(int *id)
{
	return GetDumpLineCountF(id);
}
DLL_EXPORT void GETDUMPLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
DLL_EXPORT int GETERRORLINECOUNT(int *id)
{
	return GetErrorLineCountF(id);
}
DLL_EXPORT void GETERRORLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
DLL_EXPORT int GETWARNINGLINECOUNT(int *id)
{
	return GetWarningLineCountF(id);
}
DLL_EXPORT void GETWARNINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningLineF(id, n, line, line_length);
}
DLL_EXPORT int GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
DLL_EXPORT void GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
DLL_EXPORT int RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
DLL_EXPORT int RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
DLL_EXPORT int RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
DLL_EXPORT void OUTPUTLINES(int *id)
{
	OutputLinesF(id);
}
DLL_EXPORT int GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
DLL_EXPORT int GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
DLL_EXPORT int GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif