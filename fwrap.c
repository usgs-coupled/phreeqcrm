#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include "phrqtype.h"
#include "IPhreeqc.h"

/*******************************
When using GNU gcc/g++/g77
compile fortran with:
	g77 -fno-second-underscore
********************************/
#if defined(__GNUC__)
#define LoadDatabaseF                  loaddatabasef_
#define AccumulateLineF                accumulatelinef_
#define RunF                           runf_
#define RunFileF                       runfilef_
#define GetSelectedOutputValueF        getselectedoutputvaluef_
#define OutputLastErrorF               outputlasterrorf_
#define OutputLinesF                   outputlinesf_
#define GetSelectedOutputRowCountF     getselectedoutputrowcountf_
#define GetSelectedOutputColumnCountF  getselectedoutputcolumncountf_
#define GetSelectedOutputValueF        getselectedoutputvaluef_
#define SystemF                        systemf_
#endif

#include "fwrap.h"

char *
f2cstring(char* fstring, int len)
{
    char *cstr, *str;
    int  i;

    str = fstring;
    for (i = len - 1; i >= 0 && !isgraph((int)str[i]); i--);
    cstr = (char *) malloc((size_t) (i + 2));
    if (!cstr) return 0;

    cstr[i + 1] = '\0';
    memcpy(cstr,str,i+1);
    return cstr;
}

void
padfstring(char *dest, const char *src, unsigned int len)
{
    unsigned int sofar;

    for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
        *dest++ = *src++;

    while (sofar++ < len)
        *dest++ = ' ';
}

int
LoadDatabaseF(char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		AddError("LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabase(cfilename);
	free(cfilename);
	return n;
}

int
LoadDatabaseStringF(char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		AddError("LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseString(cinput);
	free(cinput);
	return n;
}


VRESULT
AccumulateLineF(char *line, unsigned int line_length)
{
	VRESULT n;
	char* cline;

	cline = f2cstring(line, line_length);
	if (!cline)
	{
		AddError("AccumulateLine: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	n = AccumulateLine(cline);
	free(cline);
	return n;
}

void
SetSelectedOutputOnF(int* sel_on)
{
	::SetSelectedOutputOn(*sel_on);
}

void
SetOutputOnF(int* output_on)
{
	::SetOutputOn(*output_on);
}

void
SetErrorOnF(int* error_on)
{
	::SetErrorOn(*error_on);
}

void
SetLogOnF(int* log_on)
{
	::SetLogOn(*log_on);
}

void
SetDumpOnF(int* dump_on)
{
	::SetDumpOn(*dump_on);
}

void
SetDumpStringOnF(int* dump_string_on)
{
	::SetDumpStringOn(*dump_string_on);
}

int
GetDumpLineCountF(void)
{
	return ::GetDumpLineCount();
}

void
GetDumpLineF(int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetDumpLine((*n) - 1), line_length);
}

int
GetErrorLineCountF(void)
{
	return ::GetErrorLineCount();
}

void
GetErrorLineF(int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetErrorLine((*n) - 1), line_length);
}

int
RunF(void)
{
	return ::Run();
}

int
RunFileF(char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		AddError("RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunFile(cfilename);
	free(cfilename);
	return n;
}

int
RunStringF(char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		AddError("RunString: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunString(cinput);
	free(cinput);
	return n;
}

int
GetSelectedOutputRowCountF(void)
{
	return ::GetSelectedOutputRowCount();
}

int
GetSelectedOutputColumnCountF(void)
{
	return ::GetSelectedOutputColumnCount();
}

VRESULT
GetSelectedOutputValueF(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	VRESULT result;
	VAR v;
	VarInit(&v);
	result = ::GetSelectedOutputValue(*row, *col, &v);

	switch (v.type) {
	case TT_EMPTY:
		*vtype = v.type;
		break;
	case TT_ERROR:
		*vtype = v.type;
		break;
	case TT_LONG:
		*vtype = TT_DOUBLE;
		*dvalue = (double)v.lVal;
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		break;
	case TT_STRING:
		*vtype = v.type;
		padfstring(svalue, v.sVal, svalue_length);
		break;
	default:
		assert(0);
	}
	::VarClear(&v);
	return result;
}

void
OutputLastErrorF(void)
{
	::OutputLastError();
}

void
OutputLinesF(void)
{
	::OutputLines();
}

#if defined(__cplusplus)
extern "C" {
#endif

int
SystemF(char* command, unsigned int command_length)
{
	char* ccommand;

	ccommand = f2cstring(command, command_length);
	if (!ccommand)
	{
		AddError("System: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = system(ccommand);
	free(ccommand);
	return n;
}

#if defined(__cplusplus)
}
#endif



#if defined(_WIN32)

#if defined(__cplusplus)
extern "C" {
#endif

//
// Intel Fortran compiler 9.1 /iface:cvf
//
int __stdcall LOADDATABASE(char *filename, unsigned int len)
{
	return LoadDatabaseF(filename, len);
}
void __stdcall OUTPUTLASTERROR(void)
{
	OutputLastErrorF();
}
int __stdcall ACCUMULATELINE(char *line, unsigned int len)
{
	return AccumulateLineF(line, len);
}
void __stdcall SETSELECTEDOUTPUTON(int *selected_on)
{
	SetSelectedOutputOnF(selected_on);
}
void __stdcall SETOUTPUTON(int *output_on)
{
	SetOutputOnF(output_on);
}
void __stdcall SETERRORON(int *error_on)
{
	SetErrorOnF(error_on);
}
void __stdcall SETLOGON(int *log_on)
{
	SetLogOnF(log_on);
}
void __stdcall SETDUMPON(int *dump_on)
{
	SetDumpOnF(dump_on);
}
void __stdcall SETDUMPSTRINGON(int *dump_string_on)
{
	SetDumpStringOnF(dump_string_on);
}
int __stdcall GETDUMPLINECOUNT(void)
{
	return GetDumpLineCountF();
}
void __stdcall GETDUMPLINE(int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(n, line, line_length);
}
int __stdcall GETERRORLINECOUNT(void)
{
	return GetErrorLineCountF();
}
void __stdcall GETERRORLINE(int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(n, line, line_length);
}
int __stdcall RUN(void)
{
	return RunF();
}
int __stdcall RUNFILE(char *filename, unsigned int len)
{
	return RunFileF(filename, len);
}
int __stdcall RUNSTRING(char *input, unsigned int len)
{
	return RunStringF(input, len);
}
void __stdcall OUTPUTLINES(void)
{
	OutputLinesF();
}
int __stdcall GETSELECTEDOUTPUTROWCOUNT(void)
{
	return GetSelectedOutputRowCountF() - 1;
}
int __stdcall GETSELECTEDOUTPUTCOLUMNCOUNT(void)
{
	return GetSelectedOutputColumnCountF();
}
int __stdcall GETSELECTEDOUTPUTVALUE(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	int adjcol = *col - 1;
	return GetSelectedOutputValueF(row, &adjcol, vtype, dvalue, svalue, svalue_length);
}
#if defined(__cplusplus)
}
#endif

#endif

