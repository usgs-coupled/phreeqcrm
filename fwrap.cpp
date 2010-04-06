#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* sprintf */
#include "phrqtype.h"
#include "IPhreeqcLib.h"  // TODO DELETE AFTER RENAMING TO IPhreeqc.h

/*******************************
When using GNU gcc/g++/g77
compile fortran with:
	g77 -fno-second-underscore
********************************/
#if defined(__GNUC__)
#define LoadDatabaseF                  loaddatabasef_
#define AccumulateLineF                accumulatelinef_
#define RunAccumulatedF                runaccumulatedf_
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
		::AddErrorM(*id, "LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseM(*id, cfilename);
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
		::AddErrorM(*id, "LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	int n = ::LoadDatabaseStringM(*id, cinput);
	free(cinput);
	return n;
}


IPL_RESULT
AccumulateLineF(int *id, char *line, unsigned int line_length)
{
	IPL_RESULT n;
	char* cline;

	cline = f2cstring(line, line_length);
	if (!cline)
	{
		::AddErrorM(*id, "AccumulateLine: Out of memory.\n");
		return IPL_OUTOFMEMORY;
	}

	n = ::AccumulateLineM(*id, cline);
	free(cline);
	return n;
}

IPL_RESULT
SetSelectedOutputOnF(int *id, int* sel_on)
{
	return ::SetSelectedOutputOnM(*id, *sel_on);
}

IPL_RESULT
SetOutputOnF(int *id, int* output_on)
{
	return ::SetOutputOnM(*id, *output_on);
}

IPL_RESULT
SetErrorOnF(int *id, int* error_on)
{
	return ::SetErrorOnM(*id, *error_on);
}

IPL_RESULT
SetLogOnF(int *id, int* log_on)
{
	return ::SetLogOnM(*id, *log_on);
}

IPL_RESULT
SetDumpOnF(int *id, int* dump_on)
{
	return ::SetDumpOnM(*id, *dump_on);
}

IPL_RESULT
SetDumpStringOnF(int *id, int* dump_string_on)
{
	return ::SetDumpStringOnM(*id, *dump_string_on);
}

int
GetDumpLineCountF(int *id)
{
	return ::GetDumpLineCountM(*id);
}

void
GetDumpLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetDumpLineM(*id, (*n) - 1), line_length);
}

int
GetErrorLineCountF(int *id)
{
	return ::GetErrorLineCountM(*id);
}

void
GetErrorLineF(int *id, int* n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetErrorLineM(*id, (*n) - 1), line_length);
}

int
GetComponentCountF(int *id)
{
	return ::GetComponentCountM(*id);
}

void
GetComponentF(int *id, int *n, char* line, unsigned int line_length)
{
	padfstring(line, ::GetComponentM(*id, (*n) - 1), line_length);
}

int
RunAccumulatedF(int *id)
{
	return ::RunAccumulatedM(*id);
}

int
RunFileF(int *id, char* filename, unsigned int filename_length)
{
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename)
	{
		::AddErrorM(*id, "RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunFileM(*id, cfilename);
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
		::AddErrorM(*id, "RunString: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	int n = ::RunStringM(*id, cinput);
	free(cinput);
	return n;
}

int
GetSelectedOutputRowCountF(int *id)
{
	int rows = ::GetSelectedOutputRowCountM(*id);
	if (rows > 0)
	{
		rows -= 1;
	}
	return rows;
}

int
GetSelectedOutputColumnCountF(int *id)
{
	return ::GetSelectedOutputColumnCountM(*id);
}

IPL_RESULT
GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	IPL_RESULT result;
	VAR v;
	VarInit(&v);
	char buffer[100];

	int adjcol = *col - 1;
	result = ::GetSelectedOutputValueM(*id, *row, adjcol, &v);

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
	::OutputLastErrorM(*id);
}

void
OutputLinesF(int *id)
{
	::OutputLinesM(*id);
}

#if defined(__cplusplus)
extern "C" {
#endif

// COMMENT: {4/5/2010 7:12:55 PM}int
// COMMENT: {4/5/2010 7:12:55 PM}SystemF(char* command, unsigned int command_length)
// COMMENT: {4/5/2010 7:12:55 PM}{
// COMMENT: {4/5/2010 7:12:55 PM}	char* ccommand;
// COMMENT: {4/5/2010 7:12:55 PM}
// COMMENT: {4/5/2010 7:12:55 PM}	ccommand = f2cstring(command, command_length);
// COMMENT: {4/5/2010 7:12:55 PM}	if (!ccommand)
// COMMENT: {4/5/2010 7:12:55 PM}	{
// COMMENT: {4/5/2010 7:12:55 PM}		AddError("System: Out of memory.\n");
// COMMENT: {4/5/2010 7:12:55 PM}		return (int)VR_OUTOFMEMORY;
// COMMENT: {4/5/2010 7:12:55 PM}	}
// COMMENT: {4/5/2010 7:12:55 PM}
// COMMENT: {4/5/2010 7:12:55 PM}	int n = system(ccommand);
// COMMENT: {4/5/2010 7:12:55 PM}	free(ccommand);
// COMMENT: {4/5/2010 7:12:55 PM}	return n;
// COMMENT: {4/5/2010 7:12:55 PM}}

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
int __stdcall CREATEIPHREEQC(void)
{
	return CreateIPhreeqcF();
}
int __stdcall LOADDATABASE(int *id, char *filename, unsigned int len)
{
	return LoadDatabaseF(id, filename, len);
}
void __stdcall OUTPUTLASTERROR(int *id)
{
	OutputLastErrorF(id);
}
int __stdcall ACCUMULATELINE(int *id, char *line, unsigned int len)
{
	return AccumulateLineF(id, line, len);
}
void __stdcall SETSELECTEDOUTPUTON(int *id, int *selected_on)
{
	SetSelectedOutputOnF(id, selected_on);
}
void __stdcall SETOUTPUTON(int *id, int *output_on)
{
	SetOutputOnF(id, output_on);
}
void __stdcall SETERRORON(int *id, int *error_on)
{
	SetErrorOnF(id, error_on);
}
void __stdcall SETLOGON(int *id, int *log_on)
{
	SetLogOnF(id, log_on);
}
void __stdcall SETDUMPON(int *id, int *dump_on)
{
	SetDumpOnF(id, dump_on);
}
void __stdcall SETDUMPSTRINGON(int *id, int *dump_string_on)
{
	SetDumpStringOnF(id, dump_string_on);
}
int __stdcall GETDUMPLINECOUNT(int *id)
{
	return GetDumpLineCountF(id);
}
void __stdcall GETDUMPLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetDumpLineF(id, n, line, line_length);
}
int __stdcall GETERRORLINECOUNT(int *id)
{
	return GetErrorLineCountF(id);
}
void __stdcall GETERRORLINE(int *id, int *n, char* line, unsigned int line_length)
{
	GetErrorLineF(id, n, line, line_length);
}
int __stdcall GETCOMPONENTCOUNT(int *id)
{
	return GetComponentCountF(id);
}
void __stdcall GETCOMPONENT(int *id, int *n, char* line, unsigned int line_length)
{
	GetComponentF(id, n, line, line_length);
}
int __stdcall RUNACCUMULATED(int *id)
{
	return RunAccumulatedF(id);
}
int __stdcall RUNFILE(int *id, char *filename, unsigned int len)
{
	return RunFileF(id, filename, len);
}
int __stdcall RUNSTRING(int *id, char *input, unsigned int len)
{
	return RunStringF(id, input, len);
}
void __stdcall OUTPUTLINES(int *id)
{
	OutputLinesF(id);
}
int __stdcall GETSELECTEDOUTPUTROWCOUNT(int *id)
{
	return GetSelectedOutputRowCountF(id);
}
int __stdcall GETSELECTEDOUTPUTCOLUMNCOUNT(int *id)
{
	return GetSelectedOutputColumnCountF(id);
}
int __stdcall GETSELECTEDOUTPUTVALUE(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	return GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, svalue_length);
}
#if defined(__cplusplus)
}
#endif

#endif

