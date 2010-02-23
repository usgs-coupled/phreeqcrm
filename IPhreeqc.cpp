#if defined(_DEBUG)
#pragma warning(disable : 4786)   // disable truncation warning
#endif
#include <iostream>  // std::cout
#include <fstream>   // std::ifstream
#include <memory>    // std::auto_ptr
#include <sstream>   // std::istringstream
#include <cstdarg>

#if defined(_WIN32)
#include <windows.h> // ::OutputDebugString
#endif

#include "phreeqcns.hxx"
#include "ErrorReporter.hxx"
#include "SelectedOutput.hxx"
#include "Var.h"
#include "IPhreeqc.h"
#include "module_files.h"

#ifdef PHREEQC_CPP
extern int dump_entities(void);
extern int delete_entities(void);
extern int run_as_cells(void);
#endif


const char OUTPUT_FILENAME[] = "phreeqc.out";
const char ERROR_FILENAME[]  = "phreeqc.err";
const char LOG_FILENAME[]    = "phreeqc.log";
const char PUNCH_FILENAME[]  = "selected.out";


const std::string& GetAccumulatedLines(void);
void ClearAccumulatedLines(void);

void EndRow(void);
void AddSelectedOutput(const char* name, const char* format, va_list argptr);


struct PhreeqcStop{};

struct IPhreeqc
{
	// ctor
	IPhreeqc(void)
		: m_pIErrorReporter(0)
		, m_bDatabaseLoaded(false)
		, m_bSelectedOutputOn(false)
		, m_sPunchFileName(NULL)
	{
		// initialize
		UnLoadDatabase();
	};
	IPhreeqc(IErrorReporter *pIErrorReporter)
		: m_pIErrorReporter(0)
		, m_bDatabaseLoaded(false)
		, m_bSelectedOutputOn(false)
		, m_sPunchFileName(NULL)
	{
		// initialize
		UnLoadDatabase();
	};
	~IPhreeqc(void)
	{
		UnLoadDatabase();
		delete[] m_sPunchFileName;
	};
	IErrorReporter  *m_pIErrorReporter;
	bool             m_bDatabaseLoaded;
	bool             m_bSelectedOutputOn;
	char            *m_sPunchFileName;
};

static std::string s_string_input;
static std::ostringstream s_oss;
static CErrorReporter<std::ostringstream> s_errorReporter;
static IPhreeqc s_IPhreeqc(&s_errorReporter);

static void check_database(const char *sz_routine);

static void do_run(const char* sz_routine, std::istream* pis, FILE* fp, int output_on, int error_on, int log_on, int selected_output_on, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie);


int istream_getc(void *cookie);
int output_handler(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);
int open_handler(const int type, const char *file_name);
int handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);


int
istream_getc(void *cookie)
{
	if (cookie)
	{
		std::istream* is = (std::istream*)cookie;
		return is->get();
	}
	return EOF;
}

int
handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	int n = OK;
	switch (action)
	{
	case ACTION_OPEN:
		n = open_handler(type, err_str);
		break;

	case ACTION_OUTPUT:
		n = output_handler(type, err_str, stop, cookie, format, args);
		break;

	default:
		n = module_handler(action, type, err_str, stop, cookie, format, args);
		break;
	}

	if (stop == STOP)
	{
		throw PhreeqcStop();
	}
	return n;
}

int
output_handler(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	IPhreeqc* pIPhreeqc = (IPhreeqc*) cookie;

	switch (type)
	{
	case OUTPUT_ERROR:
		if (pIPhreeqc && pIPhreeqc->m_pIErrorReporter)
		{
			std::ostringstream oss;
			oss << "ERROR: " << err_str << "\n";
			if (stop == STOP) {
				oss << "Stopping.\n";
			}
			pIPhreeqc->m_pIErrorReporter->AddError(oss.str().c_str());
		}
		break;

	case OUTPUT_PUNCH:
		AddSelectedOutput(err_str, format, args);
		break;

	}
	return module_handler(ACTION_OUTPUT, type, err_str, stop, cookie, format, args);
}

int
open_handler(const int type, const char *file_name)
{
	int n = OK;
	switch (type)
	{
	case OUTPUT_PUNCH:
		if (file_name)
		{
			ASSERT(s_IPhreeqc.m_sPunchFileName != file_name);
			if (s_IPhreeqc.m_sPunchFileName != file_name)
			{
				delete[] s_IPhreeqc.m_sPunchFileName;
				s_IPhreeqc.m_sPunchFileName = NULL;
				s_IPhreeqc.m_sPunchFileName = new char[::strlen(file_name) + 1];
				::strcpy(s_IPhreeqc.m_sPunchFileName, file_name);
			}
		}
		if (s_IPhreeqc.m_bSelectedOutputOn)
		{
			n = module_handler(ACTION_OPEN, type, file_name, CONTINUE, NULL, NULL, NULL);
		}
		break;
	default:
		n = module_handler(ACTION_OPEN, type, file_name, CONTINUE, NULL, NULL, NULL);
		break;
	}
	return n;
}

int
LoadDatabase(const char* filename)
{
	try
	{
		// cleanup
		//
		UnLoadDatabase();

		CSelectedOutput::Instance()->Clear();

#if 0
		// open stream
		//
		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "LoadDatabase: Unable to open:" << "\"" << filename << "\".";
			error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
		}

		// read input
		//
		read_database(istream_getc, &ifs);
#else
		// open file
		//
		FILE* f = fopen(filename, "r");
		if (!f)
		{
			std::ostringstream oss;
			oss << "LoadDatabase: Unable to open:" << "\"" << filename << "\".";
			error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
		}

		// read input
		//
		read_database(getc_callback, f);
#endif
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "LoadDatabase: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	s_IPhreeqc.m_bDatabaseLoaded = (input_error == 0);
	return input_error;
}

int
LoadDatabaseString(const char* input)
{
	try
	{
		// cleanup
		//
		UnLoadDatabase();

		CSelectedOutput::Instance()->Clear();

		std::string s(input);
		std::istringstream iss(s);

		// read input
		//
		read_database(istream_getc, &iss);
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "LoadDatabaseString: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	s_IPhreeqc.m_bDatabaseLoaded = (input_error == 0);
	return input_error;
}

void
UnLoadDatabase(void)
{
	// init IPhreeqc
	//
	s_IPhreeqc.m_pIErrorReporter   = &s_errorReporter;
	s_IPhreeqc.m_bDatabaseLoaded   = false;
	s_IPhreeqc.m_bSelectedOutputOn = false;

	// clear error state
	//
	s_errorReporter.Clear();

	// free selectedoutput
	//
	CSelectedOutput::Release();

	// initialize phreeqc
	//
	clean_up();
	add_output_callback(handler, &s_IPhreeqc);
	do_initialize();
	input_error = 0;
}


void
OutputLastError(void)
{
	std::cout << s_errorReporter.GetOS()->str().c_str() << std::endl;
}

const char*
GetLastErrorString(void)
{
	static std::string str;
	str = s_errorReporter.GetOS()->str();
	return str.c_str();
}

// COMMENT: {11/27/2006 7:00:49 PM}#if defined(_WIN32)
// COMMENT: {11/27/2006 7:00:49 PM}void
// COMMENT: {11/27/2006 7:00:49 PM}DebugOutputLastError(void)
// COMMENT: {11/27/2006 7:00:49 PM}{
// COMMENT: {11/27/2006 7:00:49 PM}	std::istringstream iss(s_errorReporter.GetOS()->str());
// COMMENT: {11/27/2006 7:00:49 PM}	std::string line;
// COMMENT: {11/27/2006 7:00:49 PM}	while (std::getline(iss, line)) {
// COMMENT: {11/27/2006 7:00:49 PM}		::OutputDebugString(line.c_str());
// COMMENT: {11/27/2006 7:00:49 PM}		::OutputDebugString("\n");
// COMMENT: {11/27/2006 7:00:49 PM}	}
// COMMENT: {11/27/2006 7:00:49 PM}}
// COMMENT: {11/27/2006 7:00:49 PM}#endif


VRESULT
AccumulateLine(const char *line)
{
	try
	{
		s_errorReporter.Clear();
		s_string_input.append(line);
		s_string_input.append("\n");
		return VR_OK;
	}
	catch (...)
	{
		s_errorReporter.AddError("AccumulateLine: An unhandled exception occured.\n");
	}
	return VR_OUTOFMEMORY;
}

int
Run(int output_on, int error_on, int log_on, int selected_output_on)
{
	static const char *sz_routine = "Run";
	try
	{
		// this may throw
		check_database(sz_routine);

		input_error = 0;

		// create input stream
		std::istringstream iss(GetAccumulatedLines());

		// this may throw
		do_run(sz_routine, &iss, NULL, output_on, error_on, log_on, selected_output_on, NULL, NULL, NULL);
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "Run: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	::ClearAccumulatedLines();
	close_output_files();
	return input_error;
}

int
RunFile(const char* filename, int output_on, int error_on, int log_on, int selected_output_on)
{
	static const char *sz_routine = "RunFile";
	try
	{
		// this may throw
		check_database(sz_routine);

		input_error = 0;

#if 0
		// create input stream
		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			input_error = 1;
			error_msg(oss.str().c_str(), STOP); // throws
		}

		// this may throw
		do_run(sz_routine, &ifs, NULL, output_on, error_on, log_on, selected_output_on, NULL, NULL, NULL);
#else
		// open file
		//
		FILE* f = fopen(filename, "r");
		if (!f)
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
		}

		// this may throw
		do_run(sz_routine, NULL, f, output_on, error_on, log_on, selected_output_on, NULL, NULL, NULL);
#endif
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "RunFile: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	close_output_files();
	return input_error;
}

int
RunString(const char* input, int output_on, int error_on, int log_on, int selected_output_on)
{
	static const char *sz_routine = "RunString";
	try
	{
		// this may throw
		check_database(sz_routine);

		input_error = 0;

		// create input stream
		std::string s(input);
		std::istringstream iss(s);

		// this may throw
		do_run(sz_routine, &iss, NULL, output_on, error_on, log_on, selected_output_on, NULL, NULL, NULL);
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "RunString: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	close_output_files();
	return input_error;
}


int
GetSelectedOutputRowCount(void)
{
	return (int)CSelectedOutput::Instance()->GetRowCount();
}

int
GetSelectedOutputColumnCount(void)
{
	return (int)CSelectedOutput::Instance()->GetColCount();
}

VRESULT
GetSelectedOutputValue(int row, int col, VAR* pVAR)
{
	s_errorReporter.Clear();
	if (!pVAR) {
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDARG pVar is NULL.\n");
		return VR_INVALIDARG;
	}

	VRESULT v = CSelectedOutput::Instance()->Get(row, col, pVAR);
	switch (v) {
	case VR_OK:
		break;
	case VR_OUTOFMEMORY:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_OUTOFMEMORY Out of memory.\n");
		break;
	case VR_BADVARTYPE:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_BADVARTYPE pVar must be initialized(VarInit) and/or cleared(VarClear).\n");
		break;
	case VR_INVALIDARG:
		// not possible
		break;
	case VR_INVALIDROW:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDROW Row index out of range.\n");
		break;
	case VR_INVALIDCOL:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDCOL Column index out of range.\n");
		break;
	}
	return v;
}

size_t
AddError(const char* error_msg)
{
	return s_errorReporter.AddError(error_msg);
}

void
OutputLines(void)
{
	std::cout << s_string_input.c_str() << std::endl;
}

// COMMENT: {11/27/2006 7:01:16 PM}#if defined(WIN32)
// COMMENT: {11/27/2006 7:01:16 PM}void
// COMMENT: {11/27/2006 7:01:16 PM}DebugOutputLines(void)
// COMMENT: {11/27/2006 7:01:16 PM}{
// COMMENT: {11/27/2006 7:01:16 PM}	std::istringstream iss(s_string_input);
// COMMENT: {11/27/2006 7:01:16 PM}	std::string line;
// COMMENT: {11/27/2006 7:01:16 PM}	while (std::getline(iss, line)) {
// COMMENT: {11/27/2006 7:01:16 PM}		::OutputDebugString(line.c_str());
// COMMENT: {11/27/2006 7:01:16 PM}		::OutputDebugString("\n");
// COMMENT: {11/27/2006 7:01:16 PM}	}
// COMMENT: {11/27/2006 7:01:16 PM}}
// COMMENT: {11/27/2006 7:01:16 PM}#endif

const std::string&
GetAccumulatedLines(void)
{
	return s_string_input;
}

void
ClearAccumulatedLines(void)
{
	s_string_input.erase();
}

int
RunWithCallback(PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie, int output_on, int error_on, int log_on, int selected_output_on)
{
	static const char *sz_routine = "RunWithCallback";
	try
	{
		// this may throw
		check_database(sz_routine);

		input_error = 0;

		// this may throw
		do_run(sz_routine, NULL, NULL, output_on, error_on, log_on, selected_output_on, pfn_pre, pfn_post, cookie);
	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "RunWithCallback: An unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}

	ClearAccumulatedLines();
	close_output_files();
	output_close(OUTPUT_DUMP); // this should be closed in close_output_files
	return input_error;
}


int
CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie)
{
	int rvalue = OK;
	try
	{
		input_error = 0;

		if (pfn)
		{
			rvalue = pfn(cookie);
		}

	}
	catch (PhreeqcStop)
	{
		// do nothing
	}
	catch (...)
	{
		const char errmsg[] = "CatchErrors: Unhandled exception occured.\n";
		try
		{
			error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
		}
	}
	return rvalue;
}

static void
check_database(const char* sz_routine)
{
	s_errorReporter.Clear();
	CSelectedOutput::Instance()->Clear();

	if (!s_IPhreeqc.m_bDatabaseLoaded)
	{
		std::ostringstream oss;
		oss << sz_routine << ": No database is loaded";
		input_error = 1;
		error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
	}
}

static void
do_run(const char* sz_routine, std::istream* pis, FILE* fp, int output_on, int error_on, int log_on, int selected_output_on, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
	std::auto_ptr<std::istringstream> auto_iss(NULL);
	char token[MAX_LENGTH];

	if (output_on)
	{
		if (output_open(OUTPUT_MESSAGE, OUTPUT_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << OUTPUT_FILENAME << "\".\n";
			warning_msg(oss.str().c_str());
		}
	}
	if (error_on)
	{
		if (output_open(OUTPUT_ERROR, ERROR_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << ERROR_FILENAME << "\".\n";
			warning_msg(oss.str().c_str());
		}
	}
	if (log_on)
	{
		if (output_open(OUTPUT_LOG, LOG_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << LOG_FILENAME << "\".\n";
			warning_msg(oss.str().c_str());
		}
	}

	s_IPhreeqc.m_bSelectedOutputOn = (selected_output_on != 0);

/*
 *   call pre-run callback
 */
	if (pfn_pre)
	{
		pfn_pre(cookie);
	}

/*
 *   set read callback
 */
	if (!pis)
	{
		if (fp)
		{
			set_read_callback(getc_callback, fp, FALSE);
		}
		else
		{
			std::auto_ptr<std::istringstream> a_iss(new std::istringstream(GetAccumulatedLines()));
			auto_iss = a_iss;
			set_read_callback(istream_getc, auto_iss.get(), FALSE);
		}
	}
	else
	{
		set_read_callback(istream_getc, pis, FALSE);
	}


/*
 *   Read input data for simulation
 */
	for (simulation = 1; ; simulation++) {

#ifdef PHREEQ98
   		AddSeries = !connect_simulations;
#endif
		sprintf(token, "Reading input data for simulation %d.", simulation);

		output_msg(OUTPUT_GUI_ERROR, "\nSimulation %d\n", simulation);

#ifdef SWIG_SHARED_OBJ
		int save_punch_in = punch.in;
#endif // SWIG_SHARED_OBJ
		dup_print(token, TRUE);
		if (read_input() == EOF) break;

#ifdef SWIG_SHARED_OBJ
		if (simulation > 1 && save_punch_in == TRUE && punch.new_def == TRUE)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning SELECTED_OUTPUT has been redefined.\n";
			warning_msg(oss.str().c_str());

		}
		if (simulation > 1 && keyword[39].keycount > 0)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning USER_PUNCH has been redefined.\n";
			warning_msg(oss.str().c_str());
		}
#endif // SWIG_SHARED_OBJ

		if (title_x != NULL) {
			sprintf(token, "TITLE");
			dup_print(token, TRUE);
			if (pr.headings == TRUE) output_msg(OUTPUT_MESSAGE,"%s\n\n", title_x);
		}

#ifdef SWIG_SHARED_OBJ
		if (punch.in == TRUE)
		{
			//
			// (punch.in == TRUE) when any "RUN" has contained
			// a SELECTED_OUTPUT block since the last LoadDatabase call.
			//
			// Since LoadDatabase inititializes punch.in to FALSE
			// (via UnLoadDatabase...do_initialize)
			// and punch.in is set to TRUE in read_selected_output
			//
			// This causes the SELECTED_OUTPUT to contain the same headings
			// until another SELECTED_OUTPUT is defined which sets the variable
			// punch.new_def to TRUE
			//
			// WHAT IF A USER_PUNCH IS DEFINED?? IS punch.new_def SET TO
			// TRUE ???
			//
			//
			if (!selected_output_on) ASSERT(!::output_isopen(OUTPUT_PUNCH));

			if (pr.punch == FALSE)
			{
				// No selected_output for this simulation
				// this happens when
				//    PRINT;  -selected_output false
				// is given as input
				// Note: this also disables the CSelectedOutput object
			}
			else
			{
				if (punch.new_def == FALSE)
				{
					if (selected_output_on && !::output_isopen(OUTPUT_PUNCH))
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						std::string filename = s_IPhreeqc.m_sPunchFileName;
						output_open(OUTPUT_PUNCH, filename.c_str());
						if (!::output_isopen(OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							punch.new_def = TRUE;
							tidy_punch();
						}
					}
				}
				else
				{
					if (selected_output_on && !::output_isopen(OUTPUT_PUNCH))
					{
						// This is a special case which could not occur in
						// phreeqc
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run with SELECTED_OUTPUT
						//
						std::string filename = PUNCH_FILENAME;
						if (s_IPhreeqc.m_sPunchFileName && ::strlen(s_IPhreeqc.m_sPunchFileName))
						{
							filename = s_IPhreeqc.m_sPunchFileName;
						}
						output_open(OUTPUT_PUNCH, filename.c_str());
						if (!::output_isopen(OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							ASSERT(punch.new_def == TRUE);
							tidy_punch();
						}
					}
				}
			}
		}

		if (!selected_output_on) ASSERT(!::output_isopen(OUTPUT_PUNCH));
		/* the converse is not necessarily true */

		n_user_punch_index = -1;
#endif // SWIG_SHARED_OBJ
		tidy_model();
#ifdef PHREEQC_CPP
		//test_classes();
#endif
#ifdef PHREEQ98
                if (!phreeq98_debug) {
#endif


/*
 *   Calculate distribution of species for initial solutions
 */
		if (new_solution) initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (new_exchange) initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (new_surface) initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (new_gas_phase) initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
		reactions();
/*
 *   Calculate inverse models
 */
		inverse_models();
/*
 *   Calculate advection
 */
		if (use.advect_in == TRUE) {
			dup_print ("Beginning of advection calculations.", TRUE);
			advection();
		}
/*
 *   Calculate transport
 */
		if (use.trans_in == TRUE) {
			dup_print ("Beginning of transport calculations.", TRUE);
			transport();
		}

#ifdef PHREEQC_CPP
/*
 *   run
 */
			run_as_cells();
#endif

/*
 *   Copy
 */
		if (new_copy) copy_entities();
#ifdef PHREEQC_CPP

/*
 *   dump
 */
			dump_entities();
/*
 *   delete
 */
			delete_entities();
#endif

/*
 *   End of simulation
 */
		dup_print( "End of simulation.", TRUE);
#ifdef PHREEQ98
                } /* if (!phreeq98_debug) */
#endif
	}

/*
 *   Display successful status
 */
	do_status();

/*
 *   call post-run callback
 */
	if (pfn_post)
	{
		pfn_post(cookie);
	}

	if (input_error > 0)
	{
		std::ostringstream oss;
		oss << "<input>\n";
		oss << s_string_input.c_str();
		oss << "</input>\n";
		error_msg(oss.str().c_str(), CONTINUE);
	}
}
