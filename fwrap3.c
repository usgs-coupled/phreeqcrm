#include "Var.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cref /assume:underscore
//
int loaddatabase_(char *filename, unsigned int len)
{
	return LoadDatabaseF(filename, len);
}
void outputlasterror_(void)
{
	OutputLastErrorF();
}
int accumulateline_(char *line, unsigned int len)
{
	return AccumulateLineF(line, len);
}
void setselectedoutputon_(int *selected_on)
{
	SetSelectedOutputOnF(selected_on);
}
void setoutputon_(int *output_on)
{
	SetOutputOnF(output_on);
}
void seterroron_(int *error_on)
{
	SetErrorOnF(error_on);
}
void setlogon_(int *log_on)
{
	SetLogOnF(log_on);
}
void setdumpon_(int *dump_on)
{
	SetLogOnF(dump_on);
}
void setdumpstringon_(int *dump_string_on)
{
	SetDumpStringOnF(dump_string_on);
}
int getdumplinecount_(void)
{
	return GetDumpLineCountF();
}
void getdumpline_(int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(n, line, line_length);
}
int geterrorlinecount_(void)
{
	return GetErrorLineCountF();
}
void geterrorline_(int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(n, line, line_length);
}
int getcomponentcount_(void)
{
	return GetComponentCountF();
}
void getcomponent_(int *n, char* line, unsigned int line_length)
{
	GetComponentF(n, line, line_length);
}
int run_(void)
{
	return RunF();
}
int runfile_(char *filename, unsigned int len)
{
	return RunFileF(filename, len);
}
int runstring_(char *input, unsigned int len)
{
	return RunFileF(input, len);
}
void outputlines_(void)
{
	OutputLinesF();
}
int getselectedoutputrowcount_(void)
{
	return GetSelectedOutputRowCountF() - 1;
}
int getselectedoutputcolumncount_(void)
{
	return GetSelectedOutputColumnCountF();
}
int getselectedoutputvalue_(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	int adjcol = *col - 1;
	return GetSelectedOutputValueF(row, &adjcol, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif