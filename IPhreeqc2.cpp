#include "IPhreeqc2.h"
#include "Phreeqc.h"                // Phreeqc

#include <cassert>                  // assert
#include <memory>                   // auto_ptr
#include "ErrorReporter.hxx"        // CErrorReporter
#include "SelectedOutput.hxx"       // CSelectedOutput
#include "dumper.h"                 // dumper

int istream_getc(void *cookie);


typedef enum {
	ACTION_ISOPEN = -100,
} module_action_type;

const char OUTPUT_FILENAME[] = "phreeqc.out";
const char ERROR_FILENAME[]  = "phreeqc.err";
const char LOG_FILENAME[]    = "phreeqc.log";
const char PUNCH_FILENAME[]  = "selected.out";

int istream_getc(void *cookie)
{
	if (cookie)
	{
		std::istream* is = (std::istream*)cookie;
		int n = is->get();
		if (n == 13 && is->peek() == 10)
		{
			n = is->get();
		}
		return n;
	}
	return EOF;
}

IPhreeqc2::IPhreeqc2(void)
: DatabaseLoaded(false)
, SelectedOutputOn(false)
, OutputOn(false)
, LogOn(false)
, ErrorOn(false)
, DumpOn(false)
, DumpStringOn(false)
, ErrorReporter(0)
, WarningReporter(0)
, SelectedOutput(0)
, PhreeqcPtr(0)
{
	this->ErrorReporter   = new CErrorReporter<std::ostringstream>;
	this->WarningReporter = new CErrorReporter<std::ostringstream>;
	this->SelectedOutput  = new CSelectedOutput();
	this->PhreeqcPtr = new Phreeqc;
	ASSERT(this->PhreeqcPtr->phast == 0);
	this->UnLoadDatabase();
}

IPhreeqc2::~IPhreeqc2(void)
{
	delete this->PhreeqcPtr;
}

int IPhreeqc2::handler2(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	int n = OK;
	IPhreeqc2 *pThis = (IPhreeqc2*)cookie;
	switch (action)
	{
	case Phreeqc::ACTION_OPEN:
		n = pThis->open_handler2(type, err_str);
		break;

	case Phreeqc::ACTION_OUTPUT:
		n = pThis->output_handler2(type, err_str, stop, cookie, format, args);
		break;

	default:
		n = pThis->module_handler2(action, type, err_str, stop, cookie, format, args);
		break;
	}

	if (stop == STOP)
	{
		throw IPhreeqcStop();
	}
	return n;
}

int IPhreeqc2::output_handler2(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	assert(cookie == this);

	switch (type)
	{
	case Phreeqc::OUTPUT_ERROR:
		if (this)
		{
			this->AddError("ERROR: ");
			this->AddError(err_str);
			this->AddError("\n");
			if (stop == STOP)
			{
				this->AddError("Stopping.\n");
			}
		}
		break;

	case Phreeqc::OUTPUT_WARNING:
		if (this)
		{
			std::ostringstream oss;
			oss << "WARNING: " << err_str << "\n";
			this->AddWarning(oss.str().c_str());
		}
		break;

	case Phreeqc::OUTPUT_PUNCH:
		this->AddSelectedOutput(err_str, format, args);
		break;

	case Phreeqc::OUTPUT_PUNCH_END_ROW:
 		this->EndRow();
		break;

	}
	return module_handler2(Phreeqc::ACTION_OUTPUT, type, err_str, stop, cookie, format, args);
}

int IPhreeqc2::open_handler2(const int type, const char *file_name)
{
	int n = OK;
	switch (type)
	{
	case Phreeqc::OUTPUT_PUNCH:
		if (file_name)
		{
			if (this->PunchFileName.compare(file_name) != 0)
			{
				this->PunchFileName = file_name;
			}
		}
		if (this->SelectedOutputOn)
		{
			n = module_handler2(Phreeqc::ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		}
		break;
	default:
		n = module_handler2(Phreeqc::ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		break;
	}
	return n;
}

int IPhreeqc2::module_handler2(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	IPhreeqc2* pThis = (IPhreeqc2*) cookie;

	switch (action)
	{
	case Phreeqc::ACTION_OPEN:
		return pThis->module_open_handler2(type, err_str);
		break;
	case ACTION_ISOPEN:
		return pThis->module_isopen_handler2(type);
		break;
	default:
		return pThis->PhreeqcPtr->phreeqc_handler(action, type, err_str, stop, pThis->PhreeqcPtr, format, args);
		break;
	}
	return ERROR;
}

int IPhreeqc2::module_isopen_handler2(const int type)
{
	switch (type)
	{
	case Phreeqc::OUTPUT_PUNCH:
		if (this->PhreeqcPtr->punch_file) return 1;
		break;
	default:
		assert(0);
	}
	return 0;
}

int IPhreeqc2::module_open_handler2(const int type, const char *file_name)
{
	assert(file_name && ::strlen(file_name));
	switch (type)
	{
	case Phreeqc::OUTPUT_MESSAGE:
		if (this->PhreeqcPtr->output != NULL)
		{
			::fclose(this->PhreeqcPtr->output);
			this->PhreeqcPtr->output = NULL;
		}
		if ( (this->PhreeqcPtr->output = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case Phreeqc::OUTPUT_ERROR:
		assert(this->PhreeqcPtr->error_file != stderr);
		if (this->PhreeqcPtr->error_file != NULL)
		{
			::fclose(this->PhreeqcPtr->error_file);
			this->PhreeqcPtr->error_file = NULL;
		}
		if ( (this->PhreeqcPtr->error_file = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case Phreeqc::OUTPUT_LOG:
		if (this->PhreeqcPtr->log_file != NULL)
		{
			::fclose(this->PhreeqcPtr->log_file);
			this->PhreeqcPtr->log_file = NULL;
		}
		if ( (this->PhreeqcPtr->log_file = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	default:
		return this->PhreeqcPtr->open_handler(type, file_name);
		break;

	}
	return(OK);
}

int IPhreeqc2::output_isopen2(const int type)
{
	size_t i;
	int isopen;
	for (i = 0; i < this->PhreeqcPtr->count_output_callback; ++i)
	{
		isopen = (this->PhreeqcPtr->output_callbacks[i].callback)(ACTION_ISOPEN, type, NULL, CONTINUE, this->PhreeqcPtr->output_callbacks[i].cookie, NULL, NULL);
		if (isopen) return 1;
	}
	return 0;
}

int IPhreeqc2::EndRow(void)
{
	if (this->SelectedOutput->GetRowCount() <= 1)
	{
		// ensure all user_punch headings are included
		assert(this->PhreeqcPtr->n_user_punch_index >= 0);
		for (int i = this->PhreeqcPtr->n_user_punch_index; i < this->PhreeqcPtr->user_punch_count_headings; ++i)
		{
			this->SelectedOutput->PushBackEmpty(this->PhreeqcPtr->user_punch_headings[i]);
		}
	}
	return this->SelectedOutput->EndRow();
}

void IPhreeqc2::ClearAccumulatedLines(void)
{
	this->StringInput.erase();
}

size_t IPhreeqc2::AddError(const char* error_msg)
{
	return this->ErrorReporter->AddError(error_msg);
}

size_t IPhreeqc2::AddWarning(const char* error_msg)
{
	return this->WarningReporter->AddError(error_msg);
}

const std::string& IPhreeqc2::GetAccumulatedLines(void)
{
	return this->StringInput;
}

void IPhreeqc2::OutputLastError(void)
{
	std::cout << this->GetLastErrorString() << std::endl;
}

void IPhreeqc2::OutputLines(void)
{
	std::cout << this->StringInput.c_str() << std::endl;
}


int IPhreeqc2::LoadDatabase(const char* filename)
{
	try
	{
		// cleanup
		//
		this->UnLoadDatabase();
		this->SelectedOutput->Clear();

		// open file
		//
		std::ifstream ifs;
		ifs.open(filename);

		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "LoadDatabase: Unable to open:" << "\"" << filename << "\".";
			this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
		}

		// read input
		//
		this->PhreeqcPtr->read_database(istream_getc, &ifs);
	}
	catch (IPhreeqcStop)
	{
		this->PhreeqcPtr->close_input_files();
	}
	catch (...)
	{
		const char *errmsg = "LoadDatabase: An unhandled exception occured.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
	}

	this->DatabaseLoaded = (this->PhreeqcPtr->input_error == 0);
	return this->PhreeqcPtr->input_error;
}

int IPhreeqc2::LoadDatabaseString(const char* input)
{
	try
	{
		// cleanup
		//
		this->UnLoadDatabase();

		this->SelectedOutput->Clear();

		std::string s(input);
		std::istringstream iss(s);

		// read input
		//
		this->PhreeqcPtr->read_database(istream_getc, &iss);
	}
	catch (IPhreeqcStop)
	{
		this->PhreeqcPtr->close_input_files();
	}
	catch(...)
	{
		const char *errmsg = "LoadDatabaseString: An unhandled exception occured.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
	}

	this->DatabaseLoaded = (this->PhreeqcPtr->input_error == 0);
	return this->PhreeqcPtr->input_error;
}

void IPhreeqc2::UnLoadDatabase(void)
{
	// init IPhreeqc
	//
	this->DatabaseLoaded   = false;

	// clear error state
	//
	ASSERT(this->ErrorReporter);
	this->ErrorReporter->Clear();
	this->LastErrorString.clear();

	// clear warning state
	//
	ASSERT(this->WarningReporter);
	this->WarningReporter->Clear();
	this->LastWarningString.clear();

	// clear selectedoutput
	//
	ASSERT(this->SelectedOutput);
	this->SelectedOutput->Clear();

	// clear dump string
	//
	this->DumpString.clear();
	this->DumpLines.clear();

	// initialize phreeqc
	//
	this->PhreeqcPtr->clean_up();
	this->PhreeqcPtr->add_output_callback(IPhreeqc2::handler2, this);
	this->PhreeqcPtr->do_initialize();
	this->PhreeqcPtr->input_error = 0;
}

bool IPhreeqc2::GetOutputOn(void)const
{
	return this->OutputOn;
}

void IPhreeqc2::SetOutputOn(bool bValue)
{
	this->OutputOn = bValue;
}

bool IPhreeqc2::GetSelectedOutputOn(void)const
{
	return this->SelectedOutputOn;
}

void IPhreeqc2::SetSelectedOutputOn(bool bValue)
{
	this->SelectedOutputOn = bValue;
}

bool IPhreeqc2::GetLogOn(void)const
{
	return this->LogOn;
}

void IPhreeqc2::SetLogOn(bool bValue)
{
	this->LogOn = bValue;
}

bool IPhreeqc2::GetDumpOn(void)const
{
	return this->DumpOn;
}

void IPhreeqc2::SetDumpOn(bool bValue)
{
	this->DumpOn = bValue;
}

bool IPhreeqc2::GetDumpStringOn(void)const
{
	return this->DumpStringOn;
}

void IPhreeqc2::SetDumpStringOn(bool bValue)
{
	this->DumpStringOn = bValue;
}

bool IPhreeqc2::GetErrorOn(void)const
{
	return this->ErrorOn;
}

void IPhreeqc2::SetErrorOn(bool bValue)
{
	this->ErrorOn = bValue;
}

void IPhreeqc2::AddSelectedOutput(const char* name, const char* format, va_list argptr)
{
	int bInt;
	int bDouble;
	int bString;

	int state;
	int bLongDouble;
	char ch;


	/* state values
	0 Haven't found start(%)
	1 Just read start(%)
	2 Just read Flags(-0+ #) (zero or more)
	3 Just read Width
	4 Just read Precision start (.)
	5 Just read Size modifier
	6 Just read Type
	*/

	if (name == NULL) {
		return;
	}

	bDouble = 0;
	bInt = 0;
	bString = 0;

	bLongDouble = 0;

	state = 0;
	ch = *format++;
	while (ch != '\0')
	{
		switch (state)
		{
		case 0: /* looking for Start specification (%) */
			switch (ch)
			{
			case '%':
				state = 1;
				break;
			default:
				break;
			}
			ch = *format++;
			break;
		case 1: /* reading Flags (zero or more(-,+,0,# or space)) */
			switch (ch)
			{
				case '-': case '0': case '+': case ' ': case '#':
					ch = *format++;
					break;
				default:
					state = 2;
					break;
			}
			break;
		case 2: /* reading Minimum field width (decimal integer constant) */
			switch (ch)
			{
			case '.':
				state = 3;
				ch = *format++;
				break;
			case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
				ch = *format++;
				break;
			default:
				state = 4;
				break;
			}
			break;
		case 3: /* reading Precision specification (period already read) */
			switch (ch)
			{
			case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
				ch = *format++;
				break;
			default:
				state = 4;
				break;
			}
			break;
		case 4: /* reading Size modifier */
			switch (ch)
			{
			case 'l':
				ch = *format++;
				break;
			case 'L':
				bLongDouble = 1;
				ch = *format++;
				break;
			case 'h':
				ch = *format++;
				break;
			}
			state = 5;
			break;
		case 5: /* reading Conversion letter */
			switch (ch)
			{
			case 'c':
				break;
			case 'd':
			case 'i':
				bInt = 1;
				break;
			case 'n':
			case 'o':
			case 'p':
				break;
			case 's':
				bString = 1;
				break;
			case 'u':
			case 'x':
			case 'X':
			case '%':
				break;
			case 'f':
			case 'e':
			case 'E':
			case 'g':
			case 'G':
				bDouble = 1;
				break;
			default:
				ASSERT(false);
				break;
			}
			ch = '\0';  /* done */
			break;
		}
	}

	if (bDouble)
	{
		double valDouble;

		if (bLongDouble)
		{
			valDouble = (double)va_arg(argptr, long double);
		}
		else
		{
			valDouble = va_arg(argptr, double);
		}

		this->SelectedOutput->PushBackDouble(name, valDouble);
	}
	else if (bInt)
	{
		int valInt;
		valInt = va_arg(argptr, int);

		this->SelectedOutput->PushBackLong(name, (long)valInt);
	}
	else if (bString)
	{
		char* valString;
		valString = (char *)va_arg(argptr, char *);

		this->SelectedOutput->PushBackString(name, valString);
	}
	else
	{
		ASSERT(false);
		this->SelectedOutput->PushBackEmpty(name);
	}
}

const char* IPhreeqc2::GetLastErrorString(void)
{
	this->LastErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	return this->LastErrorString.c_str();
}

const char* IPhreeqc2::GetLastWarningString(void)
{
	this->LastWarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	return this->LastWarningString.c_str();
}

const char* IPhreeqc2::GetDumpString(void)
{
	return this->DumpString.c_str();
}


void IPhreeqc2::check_database(const char* sz_routine)
{
	this->ErrorReporter->Clear();
	this->SelectedOutput->Clear();

	if (!this->DatabaseLoaded)
	{
		std::ostringstream oss;
		oss << sz_routine << ": No database is loaded";
		this->PhreeqcPtr->input_error = 1;
		this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
	}
}

void IPhreeqc2::do_run(const char* sz_routine, std::istream* pis, FILE* fp, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
	std::auto_ptr<std::istringstream> auto_iss(NULL);
	char token[MAX_LENGTH];

	if (this->OutputOn)
	{
		if (this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_MESSAGE, OUTPUT_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << OUTPUT_FILENAME << "\".\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());
		}
	}
	if (this->ErrorOn)
	{
		if (this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_ERROR, ERROR_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << ERROR_FILENAME << "\".\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());
		}
	}
	if (this->LogOn)
	{
		if (this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_LOG, LOG_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << LOG_FILENAME << "\".\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());
		}
	}

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
			this->PhreeqcPtr->set_read_callback(Phreeqc::getc_callback, fp, FALSE);
		}
		else
		{
			std::auto_ptr<std::istringstream> a_iss(new std::istringstream(this->GetAccumulatedLines()));
			auto_iss = a_iss;
			this->PhreeqcPtr->set_read_callback(istream_getc, auto_iss.get(), FALSE);
		}
	}
	else
	{
		this->PhreeqcPtr->set_read_callback(istream_getc, pis, FALSE);
	}


/*
 *   Read input data for simulation
 */
	for (this->PhreeqcPtr->simulation = 1; ; this->PhreeqcPtr->simulation++) {

#ifdef PHREEQ98
   		AddSeries = !connect_simulations;
#endif
		::sprintf(token, "Reading input data for simulation %d.", this->PhreeqcPtr->simulation);

		this->PhreeqcPtr->output_msg(Phreeqc::OUTPUT_GUI_ERROR, "\nSimulation %d\n", this->PhreeqcPtr->simulation);

#ifdef SWIG_SHARED_OBJ
		int save_punch_in = this->PhreeqcPtr->punch.in;
#endif // SWIG_SHARED_OBJ
		this->PhreeqcPtr->dup_print(token, TRUE);
		if (this->PhreeqcPtr->read_input() == EOF) break;

#ifdef SWIG_SHARED_OBJ
		if (this->PhreeqcPtr->simulation > 1 && save_punch_in == TRUE && this->PhreeqcPtr->punch.new_def == TRUE)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning SELECTED_OUTPUT has been redefined.\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());

		}
		if (this->PhreeqcPtr->simulation > 1 && this->PhreeqcPtr->keyword[39].keycount > 0)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning USER_PUNCH has been redefined.\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());
		}
#endif // SWIG_SHARED_OBJ

		if (this->PhreeqcPtr->title_x != NULL) {
			::sprintf(token, "TITLE");
			this->PhreeqcPtr->dup_print(token, TRUE);
			if (this->PhreeqcPtr->pr.headings == TRUE) this->PhreeqcPtr->output_msg(Phreeqc::OUTPUT_MESSAGE, "%s\n\n", this->PhreeqcPtr->title_x);
		}

#ifdef SWIG_SHARED_OBJ
		if (this->PhreeqcPtr->punch.in == TRUE)
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
			if (!this->SelectedOutputOn) ASSERT(!this->output_isopen2(Phreeqc::OUTPUT_PUNCH));

			if (this->PhreeqcPtr->pr.punch == FALSE)
			{
				// No selected_output for this simulation
				// this happens when
				//    PRINT;  -selected_output false
				// is given as input
				// Note: this also disables the CSelectedOutput object
			}
			else
			{
				if (this->PhreeqcPtr->punch.new_def == FALSE)
				{
					if (this->SelectedOutputOn && !this->output_isopen2(Phreeqc::OUTPUT_PUNCH))
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						std::string filename = this->PunchFileName;
						this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_PUNCH, filename.c_str());
						if (!this->output_isopen2(Phreeqc::OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							this->PhreeqcPtr->warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							this->PhreeqcPtr->punch.new_def = TRUE;
							this->PhreeqcPtr->tidy_punch();
						}
					}
				}
				else
				{
					if (this->SelectedOutputOn && !this->output_isopen2(Phreeqc::OUTPUT_PUNCH))
					{
						// This is a special case which could not occur in
						// phreeqc
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run with SELECTED_OUTPUT
						//
						std::string filename = PUNCH_FILENAME;
						if (this->PunchFileName.size())
						{
							filename = this->PunchFileName;
						}
						this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_PUNCH, filename.c_str());
						if (!this->output_isopen2(Phreeqc::OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							this->PhreeqcPtr->warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							ASSERT(this->PhreeqcPtr->punch.new_def == TRUE);
							this->PhreeqcPtr->tidy_punch();
						}
					}
				}
			}
		}

		if (!this->SelectedOutputOn) ASSERT(!this->output_isopen2(Phreeqc::OUTPUT_PUNCH));
		/* the converse is not necessarily true */

		this->PhreeqcPtr->n_user_punch_index = -1;
#endif // SWIG_SHARED_OBJ
		this->PhreeqcPtr->tidy_model();
#ifdef PHREEQC_CPP
		//test_classes();
#endif
#ifdef PHREEQ98
                if (!phreeq98_debug) {
#endif


/*
 *   Calculate distribution of species for initial solutions
 */
		if (this->PhreeqcPtr->new_solution) this->PhreeqcPtr->initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (this->PhreeqcPtr->new_exchange) this->PhreeqcPtr->initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (this->PhreeqcPtr->new_surface) this->PhreeqcPtr->initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (this->PhreeqcPtr->new_gas_phase) this->PhreeqcPtr->initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
		this->PhreeqcPtr->reactions();
/*
 *   Calculate inverse models
 */
		this->PhreeqcPtr->inverse_models();
/*
 *   Calculate advection
 */
		if (this->PhreeqcPtr->use.advect_in == TRUE) {
			this->PhreeqcPtr->dup_print ("Beginning of advection calculations.", TRUE);
			this->PhreeqcPtr->advection();
		}
/*
 *   Calculate transport
 */
		if (this->PhreeqcPtr->use.trans_in == TRUE) {
			this->PhreeqcPtr->dup_print ("Beginning of transport calculations.", TRUE);
			this->PhreeqcPtr->transport();
		}

#ifdef PHREEQC_CPP
/*
 *   run
 */
			this->PhreeqcPtr->run_as_cells();
#endif

/*
 *   Copy
 */
		if (this->PhreeqcPtr->new_copy) this->PhreeqcPtr->copy_entities();
#ifdef PHREEQC_CPP

/*
 *   dump
 */
		dumper dump_info_save(this->PhreeqcPtr->dump_info);
		if (this->DumpOn)
		{
			this->PhreeqcPtr->dump_entities();
		}
		if (this->DumpStringOn)
		{
			this->PhreeqcPtr->dump_info = dump_info_save;
			if (this->PhreeqcPtr->dump_info.Get_bool_any())
			{
				std::ostringstream oss;
				this->PhreeqcPtr->dump_ostream(oss);
				if (this->PhreeqcPtr->dump_info.get_append())
				{
					this->DumpString += oss.str();
				}
				else
				{
					this->DumpString = oss.str();
				}

				/* Fill dump lines */
				this->DumpLines.clear();
				std::istringstream iss(this->DumpString);
				std::string line;
				while (std::getline(iss, line))
				{
					this->DumpLines.push_back(line);
				}
			}
		}
/*
 *   delete
 */
			this->PhreeqcPtr->delete_entities();
#endif

/*
 *   End of simulation
 */
		this->PhreeqcPtr->dup_print( "End of simulation.", TRUE);
#ifdef PHREEQ98
                } /* if (!phreeq98_debug) */
#endif
	}

/*
 *   Display successful status
 */
	this->PhreeqcPtr->do_status();

/*
 *   call post-run callback
 */
	if (pfn_post)
	{
		pfn_post(cookie);
	}

	if (this->PhreeqcPtr->input_error > 0)
	{
		std::ostringstream oss;
		oss << "<input>\n";
		oss << this->StringInput.c_str();
		oss << "</input>\n";
		this->PhreeqcPtr->error_msg(oss.str().c_str(), CONTINUE);
	}
	this->update_errors();
}

VRESULT IPhreeqc2::AccumulateLine(const char *line)
{
	try
	{
		this->ErrorReporter->Clear();
		this->WarningReporter->Clear();
		this->StringInput.append(line);
		this->StringInput.append("\n");
		return VR_OK;
	}
	catch (...)
	{
		this->AddError("AccumulateLine: An unhandled exception occured.\n");
	}
	return VR_OUTOFMEMORY;
}

int IPhreeqc2::RunAccumulated(void)
{
	static const char *sz_routine = "RunAccumulated";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;

		// create input stream
		std::istringstream iss(this->GetAccumulatedLines());

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL, NULL);
	}
	catch (IPhreeqcStop)
	{
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "RunAccumulated: An unhandled exception occured.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
	}

	this->ClearAccumulatedLines();
	this->PhreeqcPtr->close_output_files();
	this->update_errors();

	return this->PhreeqcPtr->input_error;
}

int IPhreeqc2::RunFile(const char* filename)
{
	static const char *sz_routine = "RunFile";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;

		// open file
		//
		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			this->PhreeqcPtr->error_msg(oss.str().c_str(), STOP); // throws
		}

		// this may throw
		this->do_run(sz_routine, &ifs, NULL, NULL, NULL, NULL);
	}
	catch (IPhreeqcStop)
	{
		this->PhreeqcPtr->close_input_files();
	}
	catch(...)
	{
		const char *errmsg = "RunFile: An unhandled exception occured.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
	}

	this->PhreeqcPtr->close_output_files();
	this->update_errors();

	return this->PhreeqcPtr->input_error;
}

int IPhreeqc2::RunString(const char* input)
{
	static const char *sz_routine = "RunString";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->PhreeqcPtr->input_error = 0;

		// create input stream
		std::string s(input);
		std::istringstream iss(s);

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL, NULL);
	}
	catch (IPhreeqcStop)
	{
		this->PhreeqcPtr->close_input_files();
	}
	catch(...)
	{
		const char *errmsg = "RunString: An unhandled exception occured.\n";
		try
		{
			this->PhreeqcPtr->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
	}

	this->PhreeqcPtr->close_output_files();
	this->update_errors();

	return this->PhreeqcPtr->input_error;
}

int IPhreeqc2::GetSelectedOutputRowCount(void)const
{
	return (int)this->SelectedOutput->GetRowCount();
}

int IPhreeqc2::GetSelectedOutputColumnCount(void)const
{
	return (int)this->SelectedOutput->GetColCount();
}

VRESULT IPhreeqc2::GetSelectedOutputValue(int row, int col, VAR* pVAR)
{
	this->ErrorReporter->Clear();
	if (!pVAR)
	{
		this->AddError("GetSelectedOutputValue: VR_INVALIDARG pVAR is NULL.\n");
		this->update_errors();
		return VR_INVALIDARG;
	}

	VRESULT v = this->SelectedOutput->Get(row, col, pVAR);
	switch (v)
	{
	case VR_OK:
		break;
	case VR_OUTOFMEMORY:
		this->AddError("GetSelectedOutputValue: VR_OUTOFMEMORY Out of memory.\n");
		break;
	case VR_BADVARTYPE:
		this->AddError("GetSelectedOutputValue: VR_BADVARTYPE pVar must be initialized(VarInit) and/or cleared(VarClear).\n");
		break;
	case VR_INVALIDARG:
		// not possible
		break;
	case VR_INVALIDROW:
		this->AddError("GetSelectedOutputValue: VR_INVALIDROW Row index out of range.\n");
		break;
	case VR_INVALIDCOL:
		this->AddError("GetSelectedOutputValue: VR_INVALIDCOL Column index out of range.\n");
		break;
	}
	this->update_errors();
	return v;
}

int IPhreeqc2::GetDumpLineCount(void)const
{
	return (int)this->DumpLines.size();
}

const char* IPhreeqc2::GetDumpLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetDumpLineCount())
	{
		return empty;
	}
	return this->DumpLines[n].c_str();
}

int IPhreeqc2::GetErrorLineCount(void)const
{
	return (int)this->ErrorLines.size();
}

const char* IPhreeqc2::GetErrorLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetErrorLineCount())
	{
		return empty;
	}
	return this->ErrorLines[n].c_str();
}

void IPhreeqc2::update_errors(void)
{
	this->LastErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();

	this->ErrorLines.clear();
	std::istringstream iss(this->LastErrorString);
	std::string line;
	while (std::getline(iss, line))
	{
		this->ErrorLines.push_back(line);
	}
}

std::list< std::string > IPhreeqc2::ListComponents(void)
{
	std::list< std::string > comps;
	this->PhreeqcPtr->list_components(comps);
	return comps;
}

size_t IPhreeqc2::GetComponentCount(void)
{
	std::list< std::string > comps;
	this->PhreeqcPtr->list_components(comps);
	return comps.size();
}

const char* IPhreeqc2::GetComponent(int n)
{
	static const char empty[] = "";
	this->Components = this->ListComponents();
	if (n < 0 || n >= (int)this->Components.size())
	{
		return empty;
	}
	std::list< std::string >::iterator it = this->Components.begin();
	for(int i = 0; i < n; ++i)
	{
		++it;
	}
	return (*it).c_str();
}
