#ifndef __FWRAP__H
#define __FWRAP__H

#if defined(_WINDLL)
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#if defined(FC_FUNC)
#define CreateIPhreeqcF               FC_FUNC (createiphreeqcf,               CREATEIPHREEQCF)
#define DestroyIPhreeqcF              FC_FUNC (destroyiphreeqcf,              DESTROYIPHREEQCF)
#define LoadDatabaseF                 FC_FUNC (loaddatabasef,                 LOADDATABASEF)
#define LoadDatabaseStringF           FC_FUNC (loaddatabasestringf,           LOADDATABASESTRINGF)
#define UnLoadDatabaseF               FC_FUNC (unloaddatabasef,               UNLOADDATABASEF)
#define AccumulateLineF               FC_FUNC (accumulatelinef,               ACCUMULATELINEF)
#define RunAccumulatedF               FC_FUNC (runaccumulatedf,               RUNACCUMULATEDF)
#define RunFileF                      FC_FUNC (runfilef,                      RUNFILEF)
#define RunStringF                    FC_FUNC (runstringf,                    RUNSTRINGF)
#define GetSelectedOutputRowCountF    FC_FUNC (getselectedoutputrowcountf,    GETSELECTEDOUTPUTROWCOUNTF)
#define GetSelectedOutputColumnCountF FC_FUNC (getselectedoutputcolumncountf, GETSELECTEDOUTPUTCOLUMNCOUNTF)
#define GetSelectedOutputValueF       FC_FUNC (getselectedoutputvaluef,       GETSELECTEDOUTPUTVALUEF)
#define SetSelectedOutputOnF          FC_FUNC (setselectedoutputonf,          SETSELECTEDOUTPUTONF)
#define SetOutputOnF                  FC_FUNC (setoutputonf,                  SETOUTPUTONF)
#define SetErrorOnF                   FC_FUNC (seterroronf,                   SETERRORONF)
#define SetLogOnF                     FC_FUNC (setlogonf,                     SETLOGONF)
#define SetDumpOnF                    FC_FUNC (setdumponf,                    SETDUMPONF)
#define SetDumpStringOnF              FC_FUNC (setdumpstringonf,              SETDUMPSTRINGONF)
#define GetDumpLineCountF             FC_FUNC (getdumplinecountf,             GETDUMPLINECOUNTF)
#define GetDumpLineF                  FC_FUNC (getdumplinef,                  GETDUMPLINEF)
#define GetErrorLineCountF            FC_FUNC (geterrorlinecountf,            GETERRORLINECOUNTF)
#define GetErrorLineF                 FC_FUNC (geterrorlinef,                 GETERRORLINEF)
#define GetWarningLineCountF          FC_FUNC (getwarninglinecountf,          GETWARNINGLINECOUNTF)
#define GetWarningLineF               FC_FUNC (getwarninglinef,               GETWARNINGLINEF)
#define GetComponentCountF            FC_FUNC (getcomponentcountf,            GETCOMPONENTCOUNTF)
#define GetComponentF                 FC_FUNC (getcomponentf,                 GETCOMPONENTF)
#define OutputLastErrorF              FC_FUNC (outputlasterrorf,              OUTPUTLASTERRORF)
#define OutputLastWarningF            FC_FUNC (outputlastwarningf,            OUTPUTLASTWARNINGF)
#define OutputLinesF                  FC_FUNC (outputlinesf,                  OUTPUTLINESF)
#endif /* FC_FUNC */

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
