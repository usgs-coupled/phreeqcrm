#ifndef _INC_IPHREEQC_HPP
#define _INC_IPHREEQC_HPP

#include "Phreeqc.h"                /* Phreeqc */
#include "IPhreeqcCallbacks.h"      /* PFN_PRERUN_CALLBACK, PFN_POSTRUN_CALLBACK, PFN_CATCH_CALLBACK */
#include "Var.h"                    /* VRESULT */
#include "SelectedOutput.hxx"


class IErrorReporter;

struct PhreeqcStop{};

class IPhreeqc : public Phreeqc
{
public:
	IPhreeqc(void);
	~IPhreeqc(void);

public:
	int LoadDatabase(const char* filename);
	int LoadDatabaseString(const char* input);

	void UnLoadDatabase(void);

	void OutputLastError(void);
	const char* GetLastErrorString(void);

	VRESULT AccumulateLine(const char *line);

	//{{
	void SetDumpOn(bool bValue);
	void SetErrorOn(bool bValue);
	void SetLogOn(bool bValue);
	void SetOutputOn(bool bValue);
	void SetSelectedOutputOn(bool bValue);
	//}}

	int Run(void);
	int RunFile(const char* filename);
	int RunString(const char* input);

	int GetSelectedOutputRowCount(void)const;
	int GetSelectedOutputColumnCount(void)const;
	VRESULT GetSelectedOutputValue(int row, int col, VAR* pVAR);

	void OutputLines(void);

	size_t AddError(const char* error_msg);

	const std::string& GetAccumulatedLines(void);
	void ClearAccumulatedLines(void);

	// Singleton for library
	static IPhreeqc* LibraryInstance();

	// Callbacks
	//

	// IPhreeqc.cpp
	static int handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);
	int output_handler(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);
	int open_handler(const int type, const char *file_name/*, void *cookie*/);

	// module_files.c
	static int module_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);
	int module_isopen_handler(const int type);
	int module_open_handler(const int type, const char *file_name);

	// module_output.c
	int output_isopen(const int type);

	virtual int EndRow(void);
	void AddSelectedOutput(const char* name, const char* format, va_list argptr);

	void check_database(const char* sz_routine);
	void do_run(const char* sz_routine, std::istream* pis, FILE* fp, int output_on, int error_on, int log_on, int selected_output_on, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie);

protected:
	void init(void);

protected:
	// Data
	IErrorReporter        *ErrorReporter;
	CSelectedOutput       *SelectedOutput;
	std::string            PunchFileName;
	bool                   DatabaseLoaded;
	std::string            StringInput;

	bool                   SelectedOutputOn;
	//{{
	bool                   OutputOn;
	bool                   LogOn;
	bool                   ErrorOn;
	bool                   DumpOn;
	bool                   DumpStringOn;
	//}}

private:
	static IPhreeqc* Instance;
};

#endif /* _INC_IPHREEQC_HPP */
