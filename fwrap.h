#ifndef __FWRAP__H
#define __FWRAP__H


#if defined(__cplusplus)
extern "C" {
#endif


  int LoadDatabaseF(char* filename, unsigned int filename_length);

  VRESULT AccumulateLineF(char *line, unsigned int line_length);

  int RunF(void);

  int RunFileF(char* filename, unsigned int filename_length);

  int RunStringF(char* input, unsigned int input_length);

  int GetSelectedOutputRowCountF(void);

  int GetSelectedOutputColumnCountF(void);

  VRESULT GetSelectedOutputValueF(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length);

  void SetSelectedOutputOnF(int* selected_output_on);

  void SetOutputOnF(int* output_on);

  void SetErrorOnF(int* error_on);

  void SetLogOnF(int* error_on);

  void OutputLastErrorF(void);

  void OutputLinesF(void);


#if defined(__cplusplus)
}
#endif

#endif  /* __FWRAP__H */
