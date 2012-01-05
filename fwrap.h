#ifndef __FWRAP__H
#define __FWRAP__H

#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

#if defined(FC_FUNC)
#define AccumulateLineF               FC_FUNC (accumulatelinef,               ACCUMULATELINEF)
#define AddErrorF                     FC_FUNC (adderrorf,                     ADDERRORF)
#define AddWarningF                   FC_FUNC (addwarningf,                   ADDWARNINGF)
#define ClearAccumulatedLinesF        FC_FUNC (clearaccumulatedlinesf,        CLEARACCUMULATEDLINESF)
#define CreateIPhreeqcF               FC_FUNC (createiphreeqcf,               CREATEIPHREEQCF)
#define DestroyIPhreeqcF              FC_FUNC (destroyiphreeqcf,              DESTROYIPHREEQCF)
#define GetComponentCountF            FC_FUNC (getcomponentcountf,            GETCOMPONENTCOUNTF)
#define GetComponentF                 FC_FUNC (getcomponentf,                 GETCOMPONENTF)
#define GetDumpStringLineCountF       FC_FUNC (getdumpstringlinecountf,       GETDUMPSTRINGLINECOUNTF)
#define GetDumpStringLineF            FC_FUNC (getdumpstringlinef,            GETDUMPSTRINGLINEF)
#define GetDumpFileNameF              FC_FUNC (getdumpfilenamef,              GETDUMPFILENAMEF)
#define GetDumpFileOnF                FC_FUNC (getdumpfileonf,                GETDUMPFILEONF)
#define GetDumpStringOnF              FC_FUNC (getdumpstringonf,              GETDUMPSTRINGONF)
#define GetErrorStringLineCountF      FC_FUNC (geterrorstringlinecountf,      GETERRORSTRINGLINECOUNTF)
#define GetErrorStringLineF           FC_FUNC (geterrorstringlinef,           GETERRORSTRINGLINEF)
#define GetErrorFileOnF               FC_FUNC (geterrorfileonf,               GETERRORFILEONF)
#define GetLogFileNameF               FC_FUNC (getlogfilenamef,               GETLOGFILENAMEF)
#define GetLogFileOnF                 FC_FUNC (getlogfileonf,                 GETLOGFILEONF)
#define GetLogStringLineCountF        FC_FUNC (getlogstringlinecountf,        GETLOGSTRINGLINECOUNTF)
#define GetLogStringLineF             FC_FUNC (getlogstringlinef,             GETLOGSTRINGLINEF)
#define GetLogStringOnF               FC_FUNC (getlogstringonf,               GETLOGSTRINGONF)
#define GetOutputFileNameF            FC_FUNC (getoutputfilenamef,            GETOUTPUTFILENAMEF)
#define GetOutputFileOnF              FC_FUNC (getoutputfileonf,              GETOUTPUTFILEONF)
#define GetOutputStringLineF          FC_FUNC (getoutputstringlinef,          GETOUTPUTSTRINGLINEF)
#define GetOutputStringLineCountF     FC_FUNC (getoutputstringlinecountf,     GETOUTPUTSTRINGLINECOUNTF)
#define GetOutputStringOnF            FC_FUNC (getoutputstringonf,            GETOUTPUTSTRINGONF)
#define GetSelectedOutputColumnCountF FC_FUNC (getselectedoutputcolumncountf, GETSELECTEDOUTPUTCOLUMNCOUNTF)
#define GetSelectedOutputFileOnF      FC_FUNC (getselectedoutputfileonf,      GETSELECTEDOUTPUTFILEONF)
#define GetSelectedOutputRowCountF    FC_FUNC (getselectedoutputrowcountf,    GETSELECTEDOUTPUTROWCOUNTF)
#define GetSelectedOutputValueF       FC_FUNC (getselectedoutputvaluef,       GETSELECTEDOUTPUTVALUEF)
#define GetWarningStringLineCountF    FC_FUNC (getwarningstringlinecountf,    GETWARNINGSTRINGLINECOUNTF)
#define GetWarningStringLineF         FC_FUNC (getwarningstringlinef,         GETWARNINGSTRINGLINEF)
#define LoadDatabaseF                 FC_FUNC (loaddatabasef,                 LOADDATABASEF)
#define LoadDatabaseStringF           FC_FUNC (loaddatabasestringf,           LOADDATABASESTRINGF)
#define OutputErrorStringF            FC_FUNC (outputerrorstringf,            OUTPUTERRORSTRINGF)
#define OutputAccumulatedLinesF       FC_FUNC (outputaccumulatedlinesf,       OUTPUTACCUMULATEDLINESF)
#define OutputWarningStringF          FC_FUNC (outputwarningstringf,          OUTPUTWARNINGSTRINGF)
#define RunAccumulatedF               FC_FUNC (runaccumulatedf,               RUNACCUMULATEDF)
#define RunFileF                      FC_FUNC (runfilef,                      RUNFILEF)
#define RunStringF                    FC_FUNC (runstringf,                    RUNSTRINGF)
#define SetDumpFileNameF              FC_FUNC (setdumpfilenamef,              SETDUMPFILENAMEF)
#define SetDumpFileOnF                FC_FUNC (setdumpfileonf,                SETDUMPFILEONF)
#define SetDumpStringOnF              FC_FUNC (setdumpstringonf,              SETDUMPSTRINGONF)
#define SetErrorFileOnF               FC_FUNC (seterrorfileonf,               SETERRORFILEONF)
#define SetLogFileNameF               FC_FUNC (setlogfilenamef,               SETLOGFILENAMEF)
#define SetLogFileOnF                 FC_FUNC (setlogfileonf,                 SETLOGFILEONF)
#define SetLogStringOnF               FC_FUNC (setlogstringonf,               SETLOGSTRINGONF)
#define SetOutputFileNameF            FC_FUNC (setoutputfilenamef,            SETOUTPUTFILENAMEF)
#define SetOutputFileOnF              FC_FUNC (setoutputfileonf,              SETOUTPUTFILEONF)
#define SetOutputStringOnF            FC_FUNC (setoutputstringonf,            SETOUTPUTSTRINGONF)
#define SetSelOutFileOnF              FC_FUNC (setseloutfileonf,              SETSELOUTFILEONF)
#endif /* FC_FUNC */

#if defined(__cplusplus)
extern "C" {
#endif

  IPQ_RESULT AccumulateLineF(int *id, char *line, unsigned int line_length);
  int        AddErrorF(int *id, char *error_msg, unsigned int len);
  int        AddWarningF(int *id, char *warn_msg, unsigned int len);
  IPQ_RESULT ClearAccumulatedLinesF(int *id);
  int        CreateIPhreeqcF(void);
  int        DestroyIPhreeqcF(int *id);
  int        GetComponentCountF(int *id);
  void       GetComponentF(int *id, int* n, char* line, unsigned int line_length);
  int        GetDumpStringLineCountF(int *id);
  void       GetDumpStringLineF(int *id, int* n, char* line, unsigned int line_length);
  void       GetDumpFileNameF(int *id, char* filename, unsigned int filename_length);
  int        GetDumpFileOnF(int *id);
  int        GetDumpStringOnF(int *id);
  int        GetErrorStringLineCountF(int *id);
  void       GetErrorStringLineF(int *id, int* n, char* line, unsigned int line_length);
  int        GetErrorFileOnF(int *id);
  void       GetLogFileNameF(int *id, char* filename, unsigned int filename_length);
  int        GetLogFileOnF(int *id);
  int        GetLogStringLineCountF(int *id);
  void       GetLogStringLineF(int *id, int* n, char* line, unsigned int line_length);
  int        GetLogStringOnF(int *id);
  void       GetOutputFileNameF(int *id, char* filename, unsigned int filename_length);
  int        GetOutputFileOnF(int *id);
  int        GetOutputStringLineCountF(int *id);
  void       GetOutputStringLineF(int *id, int* n, char* line, unsigned int line_length);
  int        GetOutputStringOnF(int *id);
  int        GetSelectedOutputColumnCountF(int *id);
  int        GetSelectedOutputFileOnF(int *id);
  int        GetSelectedOutputRowCountF(int *id);
  IPQ_RESULT GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);
  int        GetWarningStringLineCountF(int *id);
  void       GetWarningStringLineF(int *id, int* n, char* line, unsigned int line_length);
  int        LoadDatabaseF(int *id, char* filename, unsigned int filename_length);
  int        LoadDatabaseStringF(int *id, char* input, unsigned int input_length);
  void       OutputErrorStringF(int *id);
  void       OutputAccumulatedLinesF(int *id);
  void       OutputWarningStringF(int *id);
  int        RunAccumulatedF(int *id);
  int        RunFileF(int *id, char* filename, unsigned int filename_length);
  int        RunStringF(int *id, char* input, unsigned int input_length);
  IPQ_RESULT SetDumpFileNameF(int *id, char* fname, unsigned int fname_length);
  IPQ_RESULT SetDumpFileOnF(int *id, int* dump_on);
  IPQ_RESULT SetDumpStringOnF(int *id, int* dump_string_on);
  IPQ_RESULT SetErrorFileOnF(int *id, int* error_on);
  IPQ_RESULT SetLogFileNameF(int *id, char* fname, unsigned int fname_length);
  IPQ_RESULT SetLogFileOnF(int *id, int* log_on);
  IPQ_RESULT SetLogStringOnF(int *id, int* log_string_on);
  IPQ_RESULT SetOutputFileNameF(int *id, char* fname, unsigned int fname_length);
  IPQ_RESULT SetOutputFileOnF(int *id, int* output_on);
  IPQ_RESULT SetOutputStringOnF(int *id, int* output_string_on);
  IPQ_RESULT SetSelOutFileOnF(int *id, int* selected_output_on);

#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
