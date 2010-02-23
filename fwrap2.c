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
int RUN(int *output_on, int *error_on, int *log_on, int *selected_on)
{
	return RunF(output_on, error_on, log_on, selected_on);
}
int RUNFILE(char *filename, unsigned int len, int *output_on, int *error_on, int *log_on, int *selected_on)
{
	return RunFileF(output_on, error_on, log_on, selected_on, filename, len);
}
void OUTPUTLINES(void)
{
	OutputLinesF();
}
int GETSELECTEDOUTPUTROWCOUNT(void)
{
	return GetSelectedOutputRowCountF() - 1;
}
int GETSELECTEDOUTPUTCOLUMNCOUNT(void)
{
	return GetSelectedOutputColumnCountF();
}
int GETSELECTEDOUTPUTVALUE(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	int adjcol = *col - 1;
	return GetSelectedOutputValueF(row, &adjcol, vtype, dvalue, svalue, svalue_length);
}

#if defined(__cplusplus)
}
#endif

#endif