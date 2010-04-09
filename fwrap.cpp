#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* sprintf */
#include "phrqtype.h"
#include "IPhreeqc.h"

/*******************************
When using GNU gcc/g++/g77
compile fortran with:
	g77 -fno-second-underscore
********************************/
#if defined(__GNUC__)
#define AccumulateLineF                accumulatelinef_
#define CreateIPhreeqcF                createiphreeqcf_
#define GetComponentCountF             getcomponentcountf_
#define GetComponentF                  getcomponentf_
#define GetDumpLineCountF              getdumplinecountf_
#define GetDumpLineF                   getdumplinef_
#define GetErrorLineCountF             geterrorlinecountf_
#define GetErrorLineF                  geterrorlinef_
#define GetSelectedOutputColumnCountF  getselectedoutputcolumncountf_
#define GetSelectedOutputRowCountF     getselectedoutputrowcountf_
#define GetSelectedOutputValueF        getselectedoutputvaluef_
#define GetWarningLineCountF           getwarninglinecountf_
#define GetWarningLineF                getwarninglinef_
#define LoadDatabaseF                  loaddatabasef_
#define LoadDatabaseStringF            loaddatabasestringf_
#define OutputLastErrorF               outputlasterrorf_
#define OutputLinesF                   outputlinesf_
#define RunAccumulatedF                runaccumulatedf_
#define RunFileF                       runfilef_
#define RunStringF                     runstringf_
#define SetDumpOnF                     setdumponf_
#define SetDumpStringOnF               setdumpstringonf_
#define SetErrorOnF                    seterroronf_
#define SetLogOnF                      setlogonf_
#define SetOutputOnF                   setoutputonf_
#define SetSelectedOutputOnF           setselectedoutputonf_
//#define SystemF                        systemf_
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
CreateIPhreeqcF(void)
{
	return ::CreateIPhreeqc();
}

int
LoadDatabaseF(int *id, char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabase(*id, cfilename);
	free(cfilename);
	return n;
}

int
LoadDatabaseStringF(int *id, char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "LoadDatabaseString: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseString(*id, cinput);
	free(cinput);
	return n;
}

IPQ_RESULT
AccumulateLineF(int *id, char *line, unsigned int line_length)
{
	IPQ_RESULT n;
	char* cline;

	cline = f2cstring(line, line_length);
	if (!cline)
	{
		::AddError(*id, "AccumulateLine: Out of memory.\n");
		return IPQ_OUTOFMEMORY;
	}

	n = ::AccumulateLine(*id, cline);
	free(cline);
	return n;
}

IPQ_RESULT
SetSelectedOutputOnF(int *id, int* sel_on)
{
	return ::SetSelectedOutputOn(*id, *sel_on);
}

IPQ_RESULT
SetOutputOnF(int *id, int* output_on)
{
	return ::SetOutputOn(*id, *output_on);
}

IPQ_RESULT
SetErrorOnF(int *id, int* error_on)
{
	return ::SetErrorOn(*id, *error_on);
}

IPQ_RESULT
SetLogOnF(int *id, int* log_on)
{
	return ::SetLogOn(*id, *log_on);
}

IPQ_RESULT
SetDumpOnF(int *id, int* dump_on)
{
	return ::SetDumpOn(*id, *dump_on);
}

IPQ_RESULT
SetDumpStringOnF(int *id, int* dump_string_on)
{
	return ::SetDumpStringOn(*id, *dump_string_on);
}

int
GetDumpLineCountF(int *id)
{
	return ::GetDumpLineCount(*id);
}

void
GetDumpLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetDumpLine(*id, (*n) - 1), line_length);
}

int
GetErrorLineCountF(int *id)
{
	return ::GetErrorLineCount(*id);
}

void
GetErrorLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetErrorLine(*id, (*n) - 1), line_length);
}

int
GetWarningLineCountF(int *id)
{
	return ::GetWarningLineCount(*id);
}

void
GetWarningLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetWarningLine(*id, (*n) - 1), line_length);
}

int
GetComponentCountF(int *id)
{
	return ::GetComponentCount(*id);
}

void
GetComponentF(int *id, int *n, char* comp, unsigned int line_length)
{
	padfstring(comp, ::GetComponent(*id, (*n) - 1), line_length);
}

int
RunAccumulatedF(int *id)
{
	return ::RunAccumulated(*id);
}

int
RunFileF(int *id, char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddError(*id, "RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunFile(*id, cfilename);
	free(cfilename);
	return n;
}

int
RunStringF(int *id, char* input, unsigned int input_length)
{
	char* cinput;

	cinput = f2cstring(input, input_length);
	if (!cinput)
	{
		::AddError(*id, "RunString: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunString(*id, cinput);
	free(cinput);
	return n;
}

int
GetSelectedOutputRowCountF(int *id)
{
	int rows = ::GetSelectedOutputRowCount(*id);
	if (rows > 0)
	{
		rows -= 1;
	}
	return rows;
}

int
GetSelectedOutputColumnCountF(int *id)
{
	return ::GetSelectedOutputColumnCount(*id);
}

IPQ_RESULT
GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	IPQ_RESULT result;
	VAR v;
	VarInit(&v);
	char buffer[100];

	int adjcol = *col - 1;
	result = ::GetSelectedOutputValue(*id, *row, adjcol, &v);

	switch (v.type)
	{
	case TT_EMPTY:
		*vtype = v.type;
		break;
	case TT_ERROR:
		*vtype = v.type;
		break;
	case TT_LONG:
		*vtype = TT_DOUBLE;
		*dvalue = (double)v.lVal;
		::sprintf(buffer, "%ld", v.lVal);
		padfstring(svalue, buffer, svalue_length);
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		::sprintf(buffer, "%23.15e", v.dVal);
		padfstring(svalue, buffer, svalue_length);
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
OutputLastErrorF(int *id)
{
	::OutputLastError(*id);
}

void
OutputLinesF(int *id)
{
	::OutputLines(*id);
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
		//AddError("System: Out of memory.\n");
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
DLL_EXPORT int __stdcall CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
DLL_EXPORT int __stdcall LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
DLL_EXPORT void __stdcall OUTPUTLASTERROR(int *id)
{
	OutputLastErrorF(id);
}
DLL_EXPORT int __stdcall ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
DLL_EXPORT void __stdcall SETSELECTEDOUTPUTON(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
DLL_EXPORT void __stdcall SETOUTPUTON(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
DLL_EXPORT void __stdcall SETERRORON(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
DLL_EXPORT void __stdcall SETLOGON(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
DLL_EXPORT void __stdcall SETDUMPON(int *id, int *dump_on)
{
	SetDumpOnF(id, dump_on);
}
DLL_EXPORT void __stdcall SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
DLL_EXPORT int __stdcall GETDUMPLINECOUNT(int *id)
{
	return GetDumpLineCountF(id);
}
DLL_EXPORT void __stdcall GETDUMPLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
DLL_EXPORT int __stdcall GETERRORLINECOUNT(int *id)
{
	return GetErrorLineCountF(id);
}
DLL_EXPORT void __stdcall GETERRORLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
DLL_EXPORT int __stdcall GETWARNINGLINECOUNT(int *id)
{
	return GetWarningLineCountF(id);
}
DLL_EXPORT void __stdcall GETWARNINGLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetWarningLineF(id, n, line, line_length);
}
DLL_EXPORT int __stdcall GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
DLL_EXPORT void __stdcall GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
DLL_EXPORT int __stdcall RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
DLL_EXPORT int __stdcall RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
DLL_EXPORT int __stdcall RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
DLL_EXPORT void __stdcall OUTPUTLINES(int *id)
{
	OutputLinesF(id);
}
DLL_EXPORT int __stdcall GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
DLL_EXPORT int __stdcall GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
DLL_EXPORT int __stdcall GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
#if defined(__cplusplus)
}
#endif

#endif

