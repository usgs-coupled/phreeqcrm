#ifndef __FWRAP__H
#define __FWRAP__H


#if defined(__cplusplus)
extern "C" {
#endif


  int LoadDatabaseF(char* filename, unsigned int filename_length);

  VRESULT AccumulateLineF(char *line, unsigned int line_length);

  int RunF(int* output_on, int* error_on, int* log_on, int* selected_output_on);

  int RunFileF(int* output_on, int* error_on, int* log_on, int* selected_output_on, char* filename, unsigned int filename_length);
  int GetSelectedOutputRowCountF(void);

  int GetSelectedOutputColumnCountF(void);

  VRESULT GetSelectedOutputValueF(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);

  void OutputLastErrorF(void);

  void OutputLinesF(void);


#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
