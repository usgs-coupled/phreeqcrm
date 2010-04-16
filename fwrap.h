#ifndef __FWRAP__H
#define __FWRAP__H

#if defined(_WINDLL)
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#if defined(F77_FUNC)
#define CreateIPhreeqcF               F77_FUNC (createiphreeqcf,               CREATEIPHREEQCF)
#define DestroyIPhreeqcF              F77_FUNC (destroyiphreeqcf,              DESTROYIPHREEQCF)
#define LoadDatabaseF                 F77_FUNC (loaddatabasef,                 LOADDATABASEF)
#define LoadDatabaseStringF           F77_FUNC (loaddatabasestringf,           LOADDATABASESTRINGF)
#define UnLoadDatabaseF               F77_FUNC (unloaddatabasef,               UNLOADDATABASEF)
#define AccumulateLineF               F77_FUNC (accumulatelinef,               ACCUMULATELINEF)
#define RunAccumulatedF               F77_FUNC (runaccumulatedf,               RUNACCUMULATEDF)
#define RunFileF                      F77_FUNC (runfilef,                      RUNFILEF)
#define RunStringF                    F77_FUNC (runstringf,                    RUNSTRINGF)
#define GetSelectedOutputRowCountF    F77_FUNC (getselectedoutputrowcountf,    GETSELECTEDOUTPUTROWCOUNTF)
#define GetSelectedOutputColumnCountF F77_FUNC (getselectedoutputcolumncountf, GETSELECTEDOUTPUTCOLUMNCOUNTF)
#define GetSelectedOutputValueF       F77_FUNC (getselectedoutputvaluef,       GETSELECTEDOUTPUTVALUEF)
#define SetSelectedOutputOnF          F77_FUNC (setselectedoutputonf,          SETSELECTEDOUTPUTONF)
#define SetOutputOnF                  F77_FUNC (setoutputonf,                  SETOUTPUTONF)
#define SetErrorOnF                   F77_FUNC (seterroronf,                   SETERRORONF)
#define SetLogOnF                     F77_FUNC (setlogonf,                     SETLOGONF)
#define SetDumpOnF                    F77_FUNC (setdumponf,                    SETDUMPONF)
#define SetDumpStringOnF              F77_FUNC (setdumpstringonf,              SETDUMPSTRINGONF)
#define GetDumpLineCountF             F77_FUNC (getdumplinecountf,             GETDUMPLINECOUNTF)
#define GetDumpLineF                  F77_FUNC (getdumplinef,                  GETDUMPLINEF)
#define GetErrorLineCountF            F77_FUNC (geterrorlinecountf,            GETERRORLINECOUNTF)
#define GetErrorLineF                 F77_FUNC (geterrorlinef,                 GETERRORLINEF)
#define GetWarningLineCountF          F77_FUNC (getwarninglinecountf,          GETWARNINGLINECOUNTF)
#define GetWarningLineF               F77_FUNC (getwarninglinef,               GETWARNINGLINEF)
#define GetComponentCountF            F77_FUNC (getcomponentcountf,            GETCOMPONENTCOUNTF)
#define GetComponentF                 F77_FUNC (getcomponentf,                 GETCOMPONENTF)
#define OutputLastErrorF              F77_FUNC (outputlasterrorf,              OUTPUTLASTERRORF)
#define OutputLastWarningF            F77_FUNC (outputlastwarningf,            OUTPUTLASTWARNINGF)
#define OutputLinesF                  F77_FUNC (outputlinesf,                  OUTPUTLINESF)
#endif /* F77_FUNC */

#if defined(__cplusplus)
extern "C" {
#endif

  int CreateIPhreeqcF(void);

  int DestroyIPhreeqcF(int *id);

  int LoadDatabaseF(int *id, char* filename, unsigned int filename_length);

  int LoadDatabaseStringF(int *id, char* input, unsigned int input_length);

  int UnLoadDatabaseF(int *id);

  IPQ_RESULT AccumulateLineF(int *id, char *line, unsigned int line_length);

  int RunAccumulatedF(int *id);

  int RunFileF(int *id, char* filename, unsigned int filename_length);

  int RunStringF(int *id, char* input, unsigned int input_length);

  int GetSelectedOutputRowCountF(int *id);

  int GetSelectedOutputColumnCountF(int *id);

  IPQ_RESULT GetSelectedOutputValueF(int *id, int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);

  IPQ_RESULT SetSelectedOutputOnF(int *id, int* selected_output_on);

  IPQ_RESULT SetOutputOnF(int *id, int* output_on);

  IPQ_RESULT SetErrorOnF(int *id, int* error_on);

  IPQ_RESULT SetLogOnF(int *id, int* log_on);

  IPQ_RESULT SetDumpOnF(int *id, int* dump_on);

  IPQ_RESULT SetDumpStringOnF(int *id, int* dump_string_on);

  int GetDumpLineCountF(int *id);

  void GetDumpLineF(int *id, int* n, char* line, unsigned int line_length);

  int GetErrorLineCountF(int *id);

  void GetErrorLineF(int *id, int* n, char* line, unsigned int line_length);

  int GetWarningLineCountF(int *id);

  void GetWarningLineF(int *id, int* n, char* line, unsigned int line_length);

  int GetComponentCountF(int *id);

  void GetComponentF(int *id, int* n, char* line, unsigned int line_length);

  void OutputLastErrorF(int *id);

  void OutputLastWarningF(int *id);

  void OutputLinesF(int *id);


#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
