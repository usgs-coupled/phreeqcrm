#ifndef __FWRAP__H
#define __FWRAP__H

#if defined(_WINDLL)
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#if defined(FC_FUNC)
#define AccumulateLineF               FC_FUNC (accumulatelinef,               ACCUMULATELINEF)
#define CreateIPhreeqcF               FC_FUNC (createiphreeqcf,               CREATEIPHREEQCF)
#define DestroyIPhreeqcF              FC_FUNC (destroyiphreeqcf,              DESTROYIPHREEQCF)
#define GetComponentCountF            FC_FUNC (getcomponentcountf,            GETCOMPONENTCOUNTF)
#define GetComponentF                 FC_FUNC (getcomponentf,                 GETCOMPONENTF)
#define GetDumpLineCountF             FC_FUNC (getdumplinecountf,             GETDUMPLINECOUNTF)
#define GetDumpLineF                  FC_FUNC (getdumplinef,                  GETDUMPLINEF)
#define GetDumpOnF                    FC_FUNC (getdumponf,                    GETDUMPONF)
#define GetDumpStringOnF              FC_FUNC (getdumpstringonf,              GETDUMPSTRINGONF)
#define GetErrorLineCountF            FC_FUNC (geterrorlinecountf,            GETERRORLINECOUNTF)
#define GetErrorLineF                 FC_FUNC (geterrorlinef,                 GETERRORLINEF)
#define GetErrorOnF                   FC_FUNC (geterroronf,                   GETERRORONF)
#define GetLogOnF                     FC_FUNC (getlogonf,                     GETLOGONF)
#define GetOutputOnF                  FC_FUNC (getoutputonf,                  GETOUTPUTONF)
#define GetSelectedOutputColumnCountF FC_FUNC (getselectedoutputcolumncountf, GETSELECTEDOUTPUTCOLUMNCOUNTF)
#define GetSelectedOutputOnF          FC_FUNC (getselectedoutputonf,          GETSELECTEDOUTPUTONF)
#define GetSelectedOutputRowCountF    FC_FUNC (getselectedoutputrowcountf,    GETSELECTEDOUTPUTROWCOUNTF)
#define GetSelectedOutputValueF       FC_FUNC (getselectedoutputvaluef,       GETSELECTEDOUTPUTVALUEF)
#define GetWarningLineCountF          FC_FUNC (getwarninglinecountf,          GETWARNINGLINECOUNTF)
#define GetWarningLineF               FC_FUNC (getwarninglinef,               GETWARNINGLINEF)
#define LoadDatabaseF                 FC_FUNC (loaddatabasef,                 LOADDATABASEF)
#define LoadDatabaseStringF           FC_FUNC (loaddatabasestringf,           LOADDATABASESTRINGF)
#define OutputErrorF                  FC_FUNC (outputerrorf,                  OUTPUTERRORF)
#define OutputLinesF                  FC_FUNC (outputlinesf,                  OUTPUTLINESF)
#define OutputWarningF                FC_FUNC (outputwarningf,                OUTPUTWARNINGF)
#define RunAccumulatedF               FC_FUNC (runaccumulatedf,               RUNACCUMULATEDF)
#define RunFileF                      FC_FUNC (runfilef,                      RUNFILEF)
#define RunStringF                    FC_FUNC (runstringf,                    RUNSTRINGF)
#define SetDumpOnF                    FC_FUNC (setdumponf,                    SETDUMPONF)
#define SetDumpStringOnF              FC_FUNC (setdumpstringonf,              SETDUMPSTRINGONF)
#define SetErrorOnF                   FC_FUNC (seterroronf,                   SETERRORONF)
#define SetLogOnF                     FC_FUNC (setlogonf,                     SETLOGONF)
#define SetOutputOnF                  FC_FUNC (setoutputonf,                  SETOUTPUTONF)
#define SetSelectedOutputOnF          FC_FUNC (setselectedoutputonf,          SETSELECTEDOUTPUTONF)
#define UnLoadDatabaseF               FC_FUNC (unloaddatabasef,               UNLOADDATABASEF)
#endif /* FC_FUNC */

#if defined(__cplusplus)
extern "C" {
#endif

  IPQ_RESULT AccumulateLineF(int *id, char *line, unsigned int line_length);
  int        CreateIPhreeqcF(void);
  int        DestroyIPhreeqcF(int *id);
  int        GetComponentCountF(int *id);
  void       GetComponentF(int *id, int* n, char* line, unsigned int line_length);
  int        GetDumpLineCountF(int *id);
  void       GetDumpLineF(int *id, int* n, char* line, unsigned int line_length);
  int        GetErrorLineCountF(int *id);
  void       GetErrorLineF(int *id, int* n, char* line, unsigned int line_length);
  int        GetErrorOnF(int *id);
  int        GetDumpOnF(int *id);
  int        GetDumpStringOnF(int *id);
  int        GetLogOnF(int *id);
  int        GetOutputOnF(int *id);
  int        GetSelectedOutputColumnCountF(int *id);
  int        GetSelectedOutputOnF(int *id);
  int        GetSelectedOutputRowCountF(int *id);
  IPQ_RESULT GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);
  int        GetWarningLineCountF(int *id);
  void       GetWarningLineF(int *id, int* n, char* line, unsigned int line_length);
  int        LoadDatabaseF(int *id, char* filename, unsigned int filename_length);
  int        LoadDatabaseStringF(int *id, char* input, unsigned int input_length);
  void       OutputErrorF(int *id);
  void       OutputLinesF(int *id);
  void       OutputWarningF(int *id);
  int        RunAccumulatedF(int *id);
  int        RunFileF(int *id, char* filename, unsigned int filename_length);
  int        RunStringF(int *id, char* input, unsigned int input_length);
  IPQ_RESULT SetDumpOnF(int *id, int* dump_on);
  IPQ_RESULT SetDumpStringOnF(int *id, int* dump_string_on);
  IPQ_RESULT SetErrorOnF(int *id, int* error_on);
  IPQ_RESULT SetLogOnF(int *id, int* log_on);
  IPQ_RESULT SetOutputOnF(int *id, int* output_on);
  IPQ_RESULT SetSelectedOutputOnF(int *id, int* selected_output_on);
  int        UnLoadDatabaseF(int *id);

#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
