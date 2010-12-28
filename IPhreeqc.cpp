#include <memory>                   // auto_ptr

#include "IPhreeqc.hpp"             // IPhreeqc
#include "Phreeqc.h"                // Phreeqc

#include "Debug.h"                  // ASSERT
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

IPhreeqc::IPhreeqc(void)
: DatabaseLoaded(false)
, ClearAccumulated(false)
, UpdateComponents(true)
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

IPhreeqc::~IPhreeqc(void)
{
	delete this->PhreeqcPtr;
	delete this->SelectedOutput;
	delete this->WarningReporter;
	delete this->ErrorReporter;
}

VRESULT IPhreeqc::AccumulateLine(const char *line)
{
	try
	{
		if (this->ClearAccumulated)
		{
			this->ClearAccumulatedLines();
			this->ClearAccumulated = false;
		}

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

size_t IPhreeqc::AddError(const char* error_msg)
{
	return this->ErrorReporter->AddError(error_msg);
}

size_t IPhreeqc::AddWarning(const char* warn_msg)
{
	return this->WarningReporter->AddError(warn_msg);
}

void IPhreeqc::ClearAccumulatedLines(void)
{
	this->StringInput.erase();
}

const std::string& IPhreeqc::GetAccumulatedLines(void)
{
	return this->StringInput;
}

const char* IPhreeqc::GetComponent(int n)
{
	static const char empty[] = "";
	this->ListComponents();
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

size_t IPhreeqc::GetComponentCount(void)
{
	this->ListComponents();
	return this->Components.size();
}

bool IPhreeqc::GetDumpFileOn(void)const
{
	return this->DumpOn;
}

const char* IPhreeqc::GetDumpString(void)
{
	static const char err_msg[] = "GetDumpString: DumpStringOn not set.\n";
	if (!this->DumpStringOn)
	{
		return err_msg;
	}
	return this->DumpString.c_str();
}

const char* IPhreeqc::GetDumpStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetDumpStringLineCount())
	{
		return empty;
	}
	return this->DumpLines[n].c_str();
}

int IPhreeqc::GetDumpStringLineCount(void)const
{
	return (int)this->DumpLines.size();
}

bool IPhreeqc::GetDumpStringOn(void)const
{
	return this->DumpStringOn;
}

bool IPhreeqc::GetErrorFileOn(void)const
{
	return this->ErrorOn;
}

const char* IPhreeqc::GetErrorString(void)
{
	this->ErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	return this->ErrorString.c_str();
}

const char* IPhreeqc::GetErrorStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetErrorStringLineCount())
	{
		return empty;
	}
	return this->ErrorLines[n].c_str();
}

int IPhreeqc::GetErrorStringLineCount(void)const
{
	return (int)this->ErrorLines.size();
}

bool IPhreeqc::GetLogFileOn(void)const
{
	return this->LogOn;
}

bool IPhreeqc::GetOutputFileOn(void)const
{
	return this->OutputOn;
}

int IPhreeqc::GetSelectedOutputColumnCount(void)const
{
	return (int)this->SelectedOutput->GetColCount();
}

bool IPhreeqc::GetSelectedOutputFileOn(void)const
{
	return this->SelectedOutputOn;
}

int IPhreeqc::GetSelectedOutputRowCount(void)const
{
	return (int)this->SelectedOutput->GetRowCount();
}

VRESULT IPhreeqc::GetSelectedOutputValue(int row, int col, VAR* pVAR)
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

const char* IPhreeqc::GetWarningString(void)
{
	this->WarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	return this->WarningString.c_str();
}

const char* IPhreeqc::GetWarningStringLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetWarningStringLineCount())
	{
		return empty;
	}
	return this->WarningLines[n].c_str();
}

int IPhreeqc::GetWarningStringLineCount(void)const
{
	return (int)this->WarningLines.size();
}

std::list< std::string > IPhreeqc::ListComponents(void)
{
	if (this->UpdateComponents)
	{
		this->Components.clear();
		this->PhreeqcPtr->list_components(this->Components);
		this->UpdateComponents = false;
	}
	return this->Components;
}

int IPhreeqc::LoadDatabase(const char* filename)
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

int IPhreeqc::LoadDatabaseString(const char* input)
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

void IPhreeqc::OutputAccumulatedLines(void)
{
	std::cout << this->StringInput.c_str() << std::endl;
}

void IPhreeqc::OutputErrorString(void)
{
	std::cout << this->GetErrorString() << std::endl;
}

void IPhreeqc::OutputWarningString(void)
{
	std::cout << this->GetWarningString() << std::endl;
}

int IPhreeqc::RunAccumulated(void)
{
	static const char *sz_routine = "RunAccumulated";
	try
	{
		// these may throw
		this->open_output_files(sz_routine);
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

	this->ClearAccumulated = true;
	this->PhreeqcPtr->close_output_files();
	this->update_errors();

	return this->PhreeqcPtr->input_error;
}

int IPhreeqc::RunFile(const char* filename)
{
	static const char *sz_routine = "RunFile";
	try
	{
		// these may throw
		this->open_output_files(sz_routine);
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

int IPhreeqc::RunString(const char* input)
{
	static const char *sz_routine = "RunString";
	try
	{
		// these may throw
		this->open_output_files(sz_routine);
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

void IPhreeqc::SetDumpFileOn(bool bValue)
{
	this->DumpOn = bValue;
}

void IPhreeqc::SetDumpStringOn(bool bValue)
{
	this->DumpStringOn = bValue;
}

void IPhreeqc::SetErrorFileOn(bool bValue)
{
	this->ErrorOn = bValue;
}

void IPhreeqc::SetLogFileOn(bool bValue)
{
	this->LogOn = bValue;
}

void IPhreeqc::SetOutputFileOn(bool bValue)
{
	this->OutputOn = bValue;
}

void IPhreeqc::SetSelectedOutputFileOn(bool bValue)
{
	this->SelectedOutputOn = bValue;
}

void IPhreeqc::UnLoadDatabase(void)
{
	// init IPhreeqc
	//
	this->DatabaseLoaded   = false;
	this->UpdateComponents = true;
	this->Components.clear();

	// clear error state
	//
	ASSERT(this->ErrorReporter);
	this->ErrorReporter->Clear();
	this->ErrorString.clear();

	// clear warning state
	//
	ASSERT(this->WarningReporter);
	this->WarningReporter->Clear();
	this->WarningString.clear();

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
	this->PhreeqcPtr->add_output_callback(IPhreeqc::handler, this);
	this->PhreeqcPtr->do_initialize();
	this->PhreeqcPtr->input_error = 0;
}

int IPhreeqc::handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	int n = OK;
	IPhreeqc *pThis = (IPhreeqc*)cookie;
	switch (action)
	{
	case Phreeqc::ACTION_OPEN:
		n = pThis->open_handler(type, err_str);
		break;

	case Phreeqc::ACTION_OUTPUT:
		n = pThis->output_handler(type, err_str, stop, cookie, format, args);
		break;

	default:
		n = pThis->module_handler(action, type, err_str, stop, cookie, format, args);
		break;
	}

	if (stop == STOP)
	{
		throw IPhreeqcStop();
	}
	return n;
}

int IPhreeqc::output_handler(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	ASSERT(cookie == this);

	switch (type)
	{
	case Phreeqc::OUTPUT_ERROR:
		if (this)
		{
			this->AddError("ERROR: ");
			this->AddError(err_str);
			this->AddError("\n");
#if 0
			if (stop == STOP)
			{
				this->AddError("Stopping.\n");
			}
#endif
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

// COMMENT: {5/19/2010 4:50:29 PM}	case Phreeqc::OUTPUT_LOG:
// COMMENT: {5/19/2010 4:50:29 PM}		if (this)
// COMMENT: {5/19/2010 4:50:29 PM}		{
// COMMENT: {5/19/2010 4:50:29 PM}			std::ostringstream oss;
// COMMENT: {5/19/2010 4:50:29 PM}			oss << "WARNING: " << err_str << "\n";
// COMMENT: {5/19/2010 4:50:29 PM}			this->AddWarning(oss.str().c_str());
// COMMENT: {5/19/2010 4:50:29 PM}		}
// COMMENT: {5/19/2010 4:50:29 PM}		break;
	}
	return module_handler(Phreeqc::ACTION_OUTPUT, type, err_str, stop, cookie, format, args);
}

int IPhreeqc::open_handler(const int type, const char *file_name)
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
			n = module_handler(Phreeqc::ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		}
		break;
	default:
		n = module_handler(Phreeqc::ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		break;
	}
	return n;
}

int IPhreeqc::module_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	IPhreeqc* pThis = (IPhreeqc*) cookie;

	switch (action)
	{
	case Phreeqc::ACTION_OPEN:
		return pThis->module_open_handler(type, err_str);
		break;
	case ACTION_ISOPEN:
		return pThis->module_isopen_handler(type);
		break;
	default:
		return pThis->PhreeqcPtr->phreeqc_handler(action, type, err_str, stop, pThis->PhreeqcPtr, format, args);
		break;
	}
	return ERROR;
}

int IPhreeqc::module_isopen_handler(const int type)
{
	switch (type)
	{
	case Phreeqc::OUTPUT_PUNCH:
		if (this->PhreeqcPtr->punch_file) return 1;
		break;
	default:
		ASSERT(0);
	}
	return 0;
}

int IPhreeqc::module_open_handler(const int type, const char *file_name)
{
	ASSERT(file_name && ::strlen(file_name));
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
		ASSERT(this->PhreeqcPtr->error_file != stderr);
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

int IPhreeqc::output_isopen(const int type)
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

int IPhreeqc::EndRow(void)
{
	if (this->SelectedOutput->GetRowCount() <= 1)
	{
		// ensure all user_punch headings are included
		ASSERT(this->PhreeqcPtr->n_user_punch_index >= 0);
		for (int i = this->PhreeqcPtr->n_user_punch_index; i < this->PhreeqcPtr->user_punch_count_headings; ++i)
		{
			this->SelectedOutput->PushBackEmpty(this->PhreeqcPtr->user_punch_headings[i]);
		}
	}
	return this->SelectedOutput->EndRow();
}

void IPhreeqc::AddSelectedOutput(const char* name, const char* format, va_list argptr)
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

void IPhreeqc::check_database(const char* sz_routine)
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

void IPhreeqc::do_run(const char* sz_routine, std::istream* pis, FILE* fp, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
	std::auto_ptr<std::istringstream> auto_iss(NULL);
	char token[MAX_LENGTH];

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
			if (!this->SelectedOutputOn) ASSERT(!this->output_isopen(Phreeqc::OUTPUT_PUNCH));

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
					if (this->SelectedOutputOn && !this->output_isopen(Phreeqc::OUTPUT_PUNCH))
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						std::string filename = this->PunchFileName;
						this->PhreeqcPtr->output_open(Phreeqc::OUTPUT_PUNCH, filename.c_str());
						if (!this->output_isopen(Phreeqc::OUTPUT_PUNCH))
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
					if (this->SelectedOutputOn && !this->output_isopen(Phreeqc::OUTPUT_PUNCH))
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
						if (!this->output_isopen(Phreeqc::OUTPUT_PUNCH))
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

		if (!this->SelectedOutputOn) ASSERT(!this->output_isopen(Phreeqc::OUTPUT_PUNCH));
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

	this->UpdateComponents = true;
	this->update_errors();
}

void IPhreeqc::update_errors(void)
{
	this->ErrorLines.clear();
	this->ErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	if (this->ErrorString.size())
	{
		std::istringstream iss(this->ErrorString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->ErrorLines.push_back(line);
		}
	}

	this->WarningLines.clear();
	this->WarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	if (this->WarningString.size())
	{
		std::istringstream iss(this->WarningString);
		std::string line;
		while (std::getline(iss, line))
		{
			this->WarningLines.push_back(line);
		}
	}
}

void IPhreeqc::open_output_files(const char* sz_routine)
{
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
}