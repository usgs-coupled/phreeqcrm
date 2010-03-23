#include "Var.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
int LOADDATABASE(char *filename, unsigned int len)
{
	return LoadDatabaseF(filename, len);
}
void OUTPUTLASTERROR(void)
{
	OutputLastErrorF();
}
int ACCUMULATELINE(char *line, unsigned int len)
{
	return AccumulateLineF(line, len);
}
void SETSELECTEDOUTPUTON(int *selected_on)
{
	SetSelectedOutputOnF(selected_on);
}
void SETOUTPUTON(int *output_on)
{
	SetOutputOnF(output_on);
}
void SETERRORON(int *error_on)
{
	SetErrorOnF(error_on);
}
void SETLOGON(int *log_on)
{
	SetLogOnF(log_on);
}
void SETDUMPON(int *dump_on)
{
	SetDumpOnF(dump_on);
}
void SETDUMPSTRINGON(int *dump_string_on)
{
	SetDumpStringOnF(dump_string_on);
}
int GETDUMPLINECOUNT(void)
{
	return GetDumpLineCountF();
}
void GETDUMPLINE(int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(n, line, line_length);
}
int GETERRORLINECOUNT(void)
{
	return GetErrorLineCountF();
}
void GETERRORLINE(int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(n, line, line_length);
}
int GETCOMPONENTCOUNT(void)
{
	return GetComponentCountF();
}
void GETCOMPONENT(int *n, char* line, unsigned int line_length)
{
	GetComponentF(n, line, line_length);
}
int RUNACCUMULATED(void)
{
	return RunAccumulatedF();
}
int RUNFILE(char *filename, unsigned int len)
{
	return RunFileF(filename, len);
}
int RUNSTRING(char *input, unsigned int len)
{
	return RunStringF(input, len);
}
void OUTPUTLINES(void)
{
	OutputLinesF();
}
int GETSELECTEDOUTPUTROWCOUNT(void)
{
	return GetSelectedOutputRowCountF();
}
int GETSELECTEDOUTPUTCOLUMNCOUNT(void)
{
	return GetSelectedOutputColumnCountF();
}
int GETSELECTEDOUTPUTVALUE(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(row, col, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif