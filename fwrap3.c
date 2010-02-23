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
int run_(int *output_on, int *error_on, int *log_on, int *selected_on)
{
	return RunF(output_on, error_on, log_on, selected_on);
}
int runfile_(char *filename, unsigned int len, int *output_on, int *error_on, int *log_on, int *selected_on)
{
	return RunFileF(output_on, error_on, log_on, selected_on, filename, len);
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