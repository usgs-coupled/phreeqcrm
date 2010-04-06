#include "IPhreeqcLib.h"  // TODO DELETE AFTER RENAMING TO IPhreeqc.h
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
int CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
int LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
void OUTPUTLASTERROR(int *id)
{
	OutputLastErrorF(id);
}
int ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
void SETSELECTEDOUTPUTON(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
void SETOUTPUTON(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
void SETERRORON(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
void SETLOGON(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
void SETDUMPON(int *id, int *dump_on)
{
	SetDumpOnF(id, dump_on);
}
void SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
int GETDUMPLINECOUNT(int *id)
{
	return GetDumpLineCountF(id);
}
void GETDUMPLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
int GETERRORLINECOUNT(int *id)
{
	return GetErrorLineCountF(id);
}
void GETERRORLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
int GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
void GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
int RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
int RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
int RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
void OUTPUTLINES(int *id)
{
	OutputLinesF(id);
}
int GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
int GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
int GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif