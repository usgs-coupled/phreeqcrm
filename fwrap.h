#ifndef __FWRAP__H
#define __FWRAP__H


#if defined(__cplusplus)
extern "C" {
#endif

  int CreateIPhreeqcF(void);

  int LoadDatabaseF(int *id, char* filename, unsigned int filename_length);

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

  int GetComponentCountF(int *id);

  void GetComponentF(int *id, int* n, char* line, unsigned int line_length);

  void OutputLastErrorF(int *id);

  void OutputLinesF(int *id);


#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
