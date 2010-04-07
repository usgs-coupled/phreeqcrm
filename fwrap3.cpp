#include "IPhreeqc.h"
#include "fwrap.h"

#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cref /assume:underscore
//
int createiphreeqc_(void)
{
	return CreateIPhreeqcF();
}
int loaddatabase_(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
void outputlasterror_(int *id)
{
	OutputLastErrorF(id);
}
int accumulateline_(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
void setselectedoutputon_(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
void setoutputon_(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
void seterroron_(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
void setlogon_(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
void setdumpon_(int *id, int *dump_on)
{
	SetLogOnF(id, dump_on);
}
void setdumpstringon_(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
int getdumplinecount_(int *id)
{
	return GetDumpLineCountF(id);
}
void getdumpline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
int geterrorlinecount_(int *id)
{
	return GetErrorLineCountF(id);
}
void geterrorline_(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
int getcomponentcount_(int *id)
{
	return GetComponentCountF(id);
}
void getcomponent_(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
int runaccumulated_(int *id)
{
	return RunAccumulatedF(id);
}
int runfile_(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
int runstring_(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
void outputlines_(int *id)
{
	OutputLinesF(id);
}
int getselectedoutputrowcount_(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
int getselectedoutputcolumncount_(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
int getselectedoutputvalue_(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif