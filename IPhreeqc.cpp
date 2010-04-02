#include <memory>
#include <cassert>
#include "IPhreeqc.h"
#include "IPhreeqc.hpp"
#include "ErrorReporter.hxx"

int istream_getc(void *cookie);

const char OUTPUT_FILENAME[] = "phreeqc.out";
const char ERROR_FILENAME[]  = "phreeqc.err";
const char LOG_FILENAME[]    = "phreeqc.log";
const char PUNCH_FILENAME[]  = "selected.out";

int
LoadDatabase(const char* filename)
{
	return IPhreeqc::LibraryInstance()->LoadDatabase(filename);
}

// COMMENT: {3/29/2010 10:18:52 PM}int
// COMMENT: {3/29/2010 10:18:52 PM}LoadDatabaseM(int n, const char* filename)
// COMMENT: {3/29/2010 10:18:52 PM}{
// COMMENT: {3/29/2010 10:18:52 PM}	IPhreeqc* IPhreeqcPtr = IPhreeqc::GetInstance(n);
// COMMENT: {3/29/2010 10:18:52 PM}	if (IPhreeqcPtr)
// COMMENT: {3/29/2010 10:18:52 PM}	{
// COMMENT: {3/29/2010 10:18:52 PM}		return IPhreeqcPtr->LoadDatabase(filename);
// COMMENT: {3/29/2010 10:18:52 PM}	}
// COMMENT: {3/29/2010 10:18:52 PM}	return PHR_BADINSTANCE;
// COMMENT: {3/29/2010 10:18:52 PM}}

int
LoadDatabaseString(const char* input)
{
	return IPhreeqc::LibraryInstance()->LoadDatabaseString(input);
}

void
UnLoadDatabase(void)
{
	IPhreeqc::LibraryInstance()->UnLoadDatabase();
}

void
OutputLastError(void)
{
	IPhreeqc::LibraryInstance()->OutputLastError();
}

const char*
GetLastErrorString(void)
{
	return IPhreeqc::LibraryInstance()->GetLastErrorString();
}

const char* 
GetLastWarningString(void)
{
	return IPhreeqc::LibraryInstance()->GetLastWarningString();
}

const char*
GetDumpString(void)
{
	return IPhreeqc::LibraryInstance()->GetDumpString();
}


VRESULT
AccumulateLine(const char *line)
{
	return IPhreeqc::LibraryInstance()->AccumulateLine(line);
}

void
SetSelectedOutputOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetSelectedOutputOn(value != 0);
}

void
SetOutputOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetOutputOn(value != 0);
}

void
SetErrorOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetErrorOn(value != 0);
}

void
SetLogOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetLogOn(value != 0);
}

void
SetDumpOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetDumpOn(value != 0);
}

void
SetDumpStringOn(int value)
{
	return IPhreeqc::LibraryInstance()->SetDumpStringOn(value != 0);
}

int
RunAccumulated(void)
{
	return IPhreeqc::LibraryInstance()->RunAccumulated();
}

int
RunFile(const char* filename)
{
	return IPhreeqc::LibraryInstance()->RunFile(filename);
}

int
RunString(const char* input)
{
	return IPhreeqc::LibraryInstance()->RunString(input);
}

int
GetSelectedOutputRowCount(void)
{
	return (int)IPhreeqc::LibraryInstance()->GetSelectedOutputRowCount();
}

int
GetSelectedOutputColumnCount(void)
{
	return (int)IPhreeqc::LibraryInstance()->GetSelectedOutputColumnCount();
}

VRESULT
GetSelectedOutputValue(int row, int col, VAR* pVAR)
{
	return IPhreeqc::LibraryInstance()->GetSelectedOutputValue(row, col, pVAR);
}

size_t
AddError(const char* error_msg)
{
	return IPhreeqc::LibraryInstance()->AddError(error_msg);
}


void
OutputLines(void)
{
	IPhreeqc::LibraryInstance()->OutputLines();
}

const std::string&
GetAccumulatedLines(void)
{
	return IPhreeqc::LibraryInstance()->GetAccumulatedLines();
}

void
ClearAccumulatedLines(void)
{
	IPhreeqc::LibraryInstance()->ClearAccumulatedLines();
}

int
GetDumpLineCount(void)
{
	return IPhreeqc::LibraryInstance()->GetDumpLineCount();
}

const char*
GetDumpLine(int n)
{
	return IPhreeqc::LibraryInstance()->GetDumpLine(n);
}

int
GetErrorLineCount(void)
{
	return IPhreeqc::LibraryInstance()->GetErrorLineCount();
}

const char*
GetErrorLine(int n)
{
	return IPhreeqc::LibraryInstance()->GetErrorLine(n);
}

int
GetComponentCount(void)
{
	return (int)IPhreeqc::LibraryInstance()->ListComponents().size();
}

const char*
GetComponent(int n)
{
	static const char empty[] = "";
	static std::string comp;
	std::list< std::string > comps = IPhreeqc::LibraryInstance()->ListComponents();
	if (n < 0 || n >= (int)comps.size())
	{
		return empty;
	}
	std::list< std::string >::iterator it = comps.begin();
	for(int i = 0; i < n; ++i)
	{
		++it;
	}
	comp = (*it);
	return comp.c_str();
}


IPhreeqc::IPhreeqc(void)
: Phreeqc()
, ErrorReporter(0)
, WarningReporter(0)
, SelectedOutput(0)
, DatabaseLoaded(false)
, SelectedOutputOn(false)
, OutputOn(false)
, LogOn(false)
, ErrorOn(false)
, DumpOn(false)
, DumpStringOn(false)
{
	ASSERT(this->phast == 0);
	this->ErrorReporter   = new CErrorReporter<std::ostringstream>;
	this->WarningReporter = new CErrorReporter<std::ostringstream>;
	this->SelectedOutput  = new CSelectedOutput();
// COMMENT: {3/25/2010 2:37:54 PM}	this->init();
	this->UnLoadDatabase();
}

IPhreeqc::~IPhreeqc(void)
{
	delete this->ErrorReporter;
	delete this->WarningReporter;
	delete this->SelectedOutput;
}

// COMMENT: {4/1/2010 5:52:55 PM}std::map<size_t, IPhreeqc*> IPhreeqc::Instances;
// COMMENT: {4/1/2010 5:52:55 PM}size_t IPhreeqc::InstancesIndex = 0;
// COMMENT: {4/1/2010 5:52:55 PM}
// COMMENT: {4/1/2010 5:52:55 PM}int
// COMMENT: {4/1/2010 5:52:55 PM}IPhreeqc::CreateIPhreeqc(void)
// COMMENT: {4/1/2010 5:52:55 PM}{
// COMMENT: {4/1/2010 5:52:55 PM}	int n = -1;
// COMMENT: {4/1/2010 5:52:55 PM}	try
// COMMENT: {4/1/2010 5:52:55 PM}	{
// COMMENT: {4/1/2010 5:52:55 PM}		IPhreeqc* IPhreeqcPtr = new IPhreeqc;
// COMMENT: {4/1/2010 5:52:55 PM}		if (IPhreeqcPtr)
// COMMENT: {4/1/2010 5:52:55 PM}		{
// COMMENT: {4/1/2010 5:52:55 PM}			std::map<size_t, IPhreeqc*>::value_type instance(IPhreeqc::InstancesIndex, IPhreeqcPtr);
// COMMENT: {4/1/2010 5:52:55 PM}			std::pair<std::map<size_t, IPhreeqc*>::iterator, bool> pr = IPhreeqc::Instances.insert(instance);
// COMMENT: {4/1/2010 5:52:55 PM}			if (pr.second)
// COMMENT: {4/1/2010 5:52:55 PM}			{
// COMMENT: {4/1/2010 5:52:55 PM}				n = (int) (*pr.first).first;
// COMMENT: {4/1/2010 5:52:55 PM}				++IPhreeqc::InstancesIndex;
// COMMENT: {4/1/2010 5:52:55 PM}			}
// COMMENT: {4/1/2010 5:52:55 PM}		}
// COMMENT: {4/1/2010 5:52:55 PM}	}
// COMMENT: {4/1/2010 5:52:55 PM}	catch(...)
// COMMENT: {4/1/2010 5:52:55 PM}	{
// COMMENT: {4/1/2010 5:52:55 PM}		return -1;
// COMMENT: {4/1/2010 5:52:55 PM}	}
// COMMENT: {4/1/2010 5:52:55 PM}	return n;
// COMMENT: {4/1/2010 5:52:55 PM}}
// COMMENT: {4/1/2010 5:52:55 PM}
// COMMENT: {4/1/2010 5:52:55 PM}int
// COMMENT: {4/1/2010 5:52:55 PM}IPhreeqc::DestroyIPhreeqc(int n)
// COMMENT: {4/1/2010 5:52:55 PM}{
// COMMENT: {4/1/2010 5:52:55 PM}	int retval = PHR_BADINSTANCE;
// COMMENT: {4/1/2010 5:52:55 PM}	if (n >= 0)
// COMMENT: {4/1/2010 5:52:55 PM}	{
// COMMENT: {4/1/2010 5:52:55 PM}		std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.find(size_t(n));
// COMMENT: {4/1/2010 5:52:55 PM}		if (it != IPhreeqc::Instances.end())
// COMMENT: {4/1/2010 5:52:55 PM}		{
// COMMENT: {4/1/2010 5:52:55 PM}			delete (*it).second;
// COMMENT: {4/1/2010 5:52:55 PM}			IPhreeqc::Instances.erase(it);
// COMMENT: {4/1/2010 5:52:55 PM}			retval = 0;
// COMMENT: {4/1/2010 5:52:55 PM}		}
// COMMENT: {4/1/2010 5:52:55 PM}	}
// COMMENT: {4/1/2010 5:52:55 PM}	return retval;
// COMMENT: {4/1/2010 5:52:55 PM}}
// COMMENT: {4/1/2010 5:52:55 PM}
// COMMENT: {4/1/2010 5:52:55 PM}IPhreeqc*
// COMMENT: {4/1/2010 5:52:55 PM}IPhreeqc::GetInstance(int n)
// COMMENT: {4/1/2010 5:52:55 PM}{
// COMMENT: {4/1/2010 5:52:55 PM}	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.find(size_t(n));
// COMMENT: {4/1/2010 5:52:55 PM}	if (it != IPhreeqc::Instances.end())
// COMMENT: {4/1/2010 5:52:55 PM}	{
// COMMENT: {4/1/2010 5:52:55 PM}		return (*it).second;
// COMMENT: {4/1/2010 5:52:55 PM}	}
// COMMENT: {4/1/2010 5:52:55 PM}	return 0;
// COMMENT: {4/1/2010 5:52:55 PM}}

// the library singleton
IPhreeqc* IPhreeqc::Instance = 0;

IPhreeqc* IPhreeqc::LibraryInstance()
{
	if (IPhreeqc::Instance == 0)
	{
		IPhreeqc::Instance = new IPhreeqc;
	}
	return IPhreeqc::Instance;
}

int Phreeqc::EndRow(void)
{
	return OK;
}

int IPhreeqc::EndRow(void)
{
	if (this->SelectedOutput->GetRowCount() <= 1)
	{
		// ensure all user_punch headings are included
		for (int i = this->n_user_punch_index; i < this->user_punch_count_headings; ++i)
		{
			this->SelectedOutput->PushBackEmpty(this->user_punch_headings[i]);
		}
	}
	return this->SelectedOutput->EndRow();
}

void IPhreeqc::ClearAccumulatedLines(void)
{
	this->StringInput.erase();
}

size_t IPhreeqc::AddError(const char* error_msg)
{
	return this->ErrorReporter->AddError(error_msg);
}

size_t IPhreeqc::AddWarning(const char* error_msg)
{
	return this->WarningReporter->AddError(error_msg);
}

const std::string& IPhreeqc::GetAccumulatedLines(void)
{
	return this->StringInput;
}

void IPhreeqc::OutputLines(void)
{
	std::cout << this->StringInput.c_str() << std::endl;
}

void IPhreeqc::UnLoadDatabase(void)
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
	this->clean_up();
	this->add_output_callback(IPhreeqc::handler, this);
	this->do_initialize();
	this->input_error = 0;
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
			this->error_msg(oss.str().c_str(), STOP); // throws
		}

		// read input
		//
		this->read_database(istream_getc, &ifs);
	}
	catch (PhreeqcStop)
	{
		assert(this->database_file == 0);
		assert(this->input_file == 0);
		assert(this->input_error > 0);
	}
	catch (...)
	{
		const char *errmsg = "LoadDatabase: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			assert(this->database_file == 0);
			assert(this->input_file == 0);
			assert(this->input_error > 0);
		}
	}

#if defined(_WIN32)
	int n = ::_fcloseall();
	assert(n == 0);
#endif

	this->DatabaseLoaded = (this->input_error == 0);
	return this->input_error;
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
		this->read_database(istream_getc, &iss);
	}
	catch (PhreeqcStop)
	{
		assert(this->database_file == 0);
		assert(this->input_file == 0);
		assert(this->input_error > 0);
	}
	catch(...)
	{
		const char *errmsg = "LoadDatabaseString: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			assert(this->database_file == 0);
			assert(this->input_file == 0);
			assert(this->input_error > 0);
		}
	}

	this->DatabaseLoaded = (this->input_error == 0);
	return this->input_error;
}

void IPhreeqc::OutputLastError(void)
{
	std::cout << this->GetLastErrorString() << std::endl;
}

const char* IPhreeqc::GetLastErrorString(void)
{
	this->LastErrorString = ((CErrorReporter<std::ostringstream>*)this->ErrorReporter)->GetOS()->str();
	return this->LastErrorString.c_str();
}

const char* IPhreeqc::GetLastWarningString(void)
{
	this->LastWarningString = ((CErrorReporter<std::ostringstream>*)this->WarningReporter)->GetOS()->str();
	return this->LastWarningString.c_str();
}


const char* IPhreeqc::GetDumpString(void)
{
	return this->DumpString.c_str();
}

void IPhreeqc::check_database(const char* sz_routine)
{
	this->ErrorReporter->Clear();
	this->SelectedOutput->Clear();

	if (!this->DatabaseLoaded)
	{
		std::ostringstream oss;
		oss << sz_routine << ": No database is loaded";
		this->input_error = 1;
		this->error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
	}
}

void IPhreeqc::do_run(const char* sz_routine, std::istream* pis, FILE* fp, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
	std::auto_ptr<std::istringstream> auto_iss(NULL);
	char token[MAX_LENGTH];

	if (this->OutputOn)
	{
		if (this->output_open(OUTPUT_MESSAGE, OUTPUT_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << OUTPUT_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
		}
	}
	if (this->ErrorOn)
	{
		if (this->output_open(OUTPUT_ERROR, ERROR_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << ERROR_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
		}
	}
	if (this->LogOn)
	{
		if (this->output_open(OUTPUT_LOG, LOG_FILENAME) != OK)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << LOG_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
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
			this->set_read_callback(getc_callback, fp, FALSE);
		}
		else
		{
			std::auto_ptr<std::istringstream> a_iss(new std::istringstream(this->GetAccumulatedLines()));
			auto_iss = a_iss;
			this->set_read_callback(istream_getc, auto_iss.get(), FALSE);
		}
	}
	else
	{
		this->set_read_callback(istream_getc, pis, FALSE);
	}


/*
 *   Read input data for simulation
 */
	for (this->simulation = 1; ; this->simulation++) {

#ifdef PHREEQ98
   		AddSeries = !connect_simulations;
#endif
		::sprintf(token, "Reading input data for simulation %d.", simulation);

		this->output_msg(OUTPUT_GUI_ERROR, "\nSimulation %d\n", simulation);

#ifdef SWIG_SHARED_OBJ
		int save_punch_in = this->punch.in;
#endif // SWIG_SHARED_OBJ
// COMMENT: {3/3/2010 10:46:12 PM}		dup_print(token, TRUE);
		if (this->read_input() == EOF) break;

#ifdef SWIG_SHARED_OBJ
		if (this->simulation > 1 && save_punch_in == TRUE && this->punch.new_def == TRUE)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning SELECTED_OUTPUT has been redefined.\n";
			this->warning_msg(oss.str().c_str());

		}
		if (this->simulation > 1 && this->keyword[39].keycount > 0)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning USER_PUNCH has been redefined.\n";
			this->warning_msg(oss.str().c_str());
		}
#endif // SWIG_SHARED_OBJ

		if (this->title_x != NULL) {
			::sprintf(token, "TITLE");
			this->dup_print(token, TRUE);
			if (this->pr.headings == TRUE) this->output_msg(OUTPUT_MESSAGE, "%s\n\n", this->title_x);
		}

#ifdef SWIG_SHARED_OBJ
		if (this->punch.in == TRUE)
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
			if (!this->SelectedOutputOn) ASSERT(!this->output_isopen(OUTPUT_PUNCH));

			if (this->pr.punch == FALSE)
			{
				// No selected_output for this simulation
				// this happens when
				//    PRINT;  -selected_output false
				// is given as input
				// Note: this also disables the CSelectedOutput object
			}
			else
			{
				if (this->punch.new_def == FALSE)
				{
					if (this->SelectedOutputOn && !this->output_isopen(OUTPUT_PUNCH))
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						std::string filename = this->PunchFileName;
						this->output_open(OUTPUT_PUNCH, filename.c_str());
						if (!this->output_isopen(OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							this->warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							this->punch.new_def = TRUE;
							this->tidy_punch();
						}
					}
				}
				else
				{
					if (this->SelectedOutputOn && !this->output_isopen(OUTPUT_PUNCH))
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
						this->output_open(OUTPUT_PUNCH, filename.c_str());
						if (!this->output_isopen(OUTPUT_PUNCH))
						{
							std::ostringstream oss;
							oss << sz_routine << ": Unable to open:" << "\"" << filename << "\".\n";
							this->warning_msg(oss.str().c_str());
						}
						else
						{
							// output selected_output headings
							ASSERT(punch.new_def == TRUE);
							this->tidy_punch();
						}
					}
				}
			}
		}

		if (!this->SelectedOutputOn) ASSERT(!this->output_isopen(OUTPUT_PUNCH));
		/* the converse is not necessarily true */

		this->n_user_punch_index = -1;
#endif // SWIG_SHARED_OBJ
		this->tidy_model();
#ifdef PHREEQC_CPP
		//test_classes();
#endif
#ifdef PHREEQ98
                if (!phreeq98_debug) {
#endif


/*
 *   Calculate distribution of species for initial solutions
 */
		if (this->new_solution) this->initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (this->new_exchange) this->initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (this->new_surface) this->initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (this->new_gas_phase) this->initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
		this->reactions();
/*
 *   Calculate inverse models
 */
		this->inverse_models();
/*
 *   Calculate advection
 */
		if (this->use.advect_in == TRUE) {
			this->dup_print ("Beginning of advection calculations.", TRUE);
			this->advection();
		}
/*
 *   Calculate transport
 */
		if (this->use.trans_in == TRUE) {
			this->dup_print ("Beginning of transport calculations.", TRUE);
			this->transport();
		}

#ifdef PHREEQC_CPP
/*
 *   run
 */
			this->run_as_cells();
#endif

/*
 *   Copy
 */
		if (this->new_copy) this->copy_entities();
#ifdef PHREEQC_CPP

/*
 *   dump
 */
		dumper dump_info_save(dump_info);
		if (this->DumpOn)
		{
			this->dump_entities();
		}
		if (this->DumpStringOn)
		{
			dump_info = dump_info_save;
			if (dump_info.Get_bool_any())
			{
				std::ostringstream oss;
				this->dump_ostream(oss);
				if (dump_info.get_append())
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
			this->delete_entities();
#endif

/*
 *   End of simulation
 */
		this->dup_print( "End of simulation.", TRUE);
#ifdef PHREEQ98
                } /* if (!phreeq98_debug) */
#endif
	}

/*
 *   Display successful status
 */
	this->do_status();

/*
 *   call post-run callback
 */
	if (pfn_post)
	{
		pfn_post(cookie);
	}

	if (this->input_error > 0)
	{
		std::ostringstream oss;
		oss << "<input>\n";
		oss << this->StringInput.c_str();
		oss << "</input>\n";
		this->error_msg(oss.str().c_str(), CONTINUE);
	}
	//{{
	this->update_errors();
	//}}
}



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

int IPhreeqc::handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	int n = OK;
	IPhreeqc *pThis = (IPhreeqc*)cookie;
	switch (action)
	{
	case ACTION_OPEN:
		n = pThis->open_handler(type, err_str);
		break;

	case ACTION_OUTPUT:
		n = pThis->output_handler(type, err_str, stop, cookie, format, args);
		break;

	default:
		n = pThis->module_handler(action, type, err_str, stop, cookie, format, args);
		break;
	}

	if (stop == STOP)
	{
		throw PhreeqcStop();
	}
	return n;
}

int IPhreeqc::output_handler(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	assert(cookie == this);

	switch (type)
	{
	case OUTPUT_ERROR:
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

	case OUTPUT_WARNING:
		if (this)
		{
			std::ostringstream oss;
			oss << "WARNING: " << err_str << "\n";
			this->AddWarning(oss.str().c_str());
		}
		break;

	case OUTPUT_PUNCH:
		this->AddSelectedOutput(err_str, format, args);
		break;

	}
	return module_handler(ACTION_OUTPUT, type, err_str, stop, cookie, format, args);
}

int IPhreeqc::open_handler(const int type, const char *file_name)
{
	int n = OK;
	switch (type)
	{
	case OUTPUT_PUNCH:
		if (file_name)
		{
			if (this->PunchFileName.compare(file_name) != 0)
			{
				this->PunchFileName = file_name;
			}
		}
		if (this->SelectedOutputOn)
		{
			n = module_handler(ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		}
		break;
	default:
		n = module_handler(ACTION_OPEN, type, file_name, CONTINUE, this, NULL, NULL);
		break;
	}
	return n;
}

// COMMENT: {3/25/2010 2:37:25 PM}void IPhreeqc::init(void)
// COMMENT: {3/25/2010 2:37:25 PM}{
// COMMENT: {3/25/2010 2:37:25 PM}	int i;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	moles_per_kilogram_string = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	pe_string = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	debug_model = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	debug_prep = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	debug_set = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	debug_diffuse_layer = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	debug_inverse = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	itmax = 100;
// COMMENT: {3/25/2010 2:37:25 PM}#ifdef USE_LONG_DOUBLE
// COMMENT: {3/25/2010 2:37:25 PM}	/* from float.h, sets tolerance for cl1 routine */
// COMMENT: {3/25/2010 2:37:25 PM}	ineq_tol = pow((long double) 10, (long double) -LDBL_DIG);
// COMMENT: {3/25/2010 2:37:25 PM}#else
// COMMENT: {3/25/2010 2:37:25 PM}	ineq_tol = pow((double) 10, (double) -DBL_DIG);
// COMMENT: {3/25/2010 2:37:25 PM}#endif
// COMMENT: {3/25/2010 2:37:25 PM}	convergence_tolerance = 1e-8;
// COMMENT: {3/25/2010 2:37:25 PM}#ifdef USE_LONG_DOUBLE
// COMMENT: {3/25/2010 2:37:25 PM}	/* from float.h, sets tolerance for cl1 routine */
// COMMENT: {3/25/2010 2:37:25 PM}	inv_tol_default = pow((long double) 10, (long double) -LDBL_DIG + 5);
// COMMENT: {3/25/2010 2:37:25 PM}#else
// COMMENT: {3/25/2010 2:37:25 PM}	inv_tol_default = pow((double) 10, (double) -DBL_DIG + 5);
// COMMENT: {3/25/2010 2:37:25 PM}#endif
// COMMENT: {3/25/2010 2:37:25 PM}	step_size               = 100.;
// COMMENT: {3/25/2010 2:37:25 PM}	pe_step_size            = 10.;
// COMMENT: {3/25/2010 2:37:25 PM}	pp_scale                = 1.0;
// COMMENT: {3/25/2010 2:37:25 PM}	pp_column_scale         = 1.0;
// COMMENT: {3/25/2010 2:37:25 PM}	diagonal_scale          = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	censor                  = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}	mass_water_switch       = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	delay_mass_water        = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	incremental_reactions   = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	aqueous_only            = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	negative_concentrations = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	LOG_10 = log(10.0);
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	max_solution       = MAX_SOLUTION;
// COMMENT: {3/25/2010 2:37:25 PM}	max_pp_assemblage  = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_exchange       = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_surface        = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_gas_phase      = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_kinetics       = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_s_s_assemblage = MAX_PP_ASSEMBLAGE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	max_elements       = MAX_ELEMENTS;
// COMMENT: {3/25/2010 2:37:25 PM}	max_elts           = MAX_ELTS;
// COMMENT: {3/25/2010 2:37:25 PM}	max_line           = MAX_LINE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_master         = MAX_MASTER;
// COMMENT: {3/25/2010 2:37:25 PM}	max_mb_unknowns    = MAX_TRXN;
// COMMENT: {3/25/2010 2:37:25 PM}	max_phases         = MAX_PHASES;
// COMMENT: {3/25/2010 2:37:25 PM}	max_s              = MAX_S;
// COMMENT: {3/25/2010 2:37:25 PM}	max_strings        = MAX_STRINGS;
// COMMENT: {3/25/2010 2:37:25 PM}	max_trxn           = MAX_TRXN;
// COMMENT: {3/25/2010 2:37:25 PM}	max_logk           = MAX_S;
// COMMENT: {3/25/2010 2:37:25 PM}	max_master_isotope = MAX_ELTS;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	count_solution       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_pp_assemblage  = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_exchange       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_surface        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_gas_phase      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_kinetics       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_s_s_assemblage = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	count_elements       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_irrev          = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_master         = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_mix            = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_phases         = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_s              = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_temperature    = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_logk           = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_master_isotope = 0;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   Initialize advection
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	count_ad_cells   = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	count_ad_shifts  = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	print_ad_modulus = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	punch_ad_modulus = 1;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	advection_punch            = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	advection_kin_time         = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}	advection_kin_time_defined = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	advection_print            = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	advection_warnings         = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   Initialize transport
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	count_cells      = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	count_shifts     = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	ishift           = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	bcon_first       = bcon_last = 3;
// COMMENT: {3/25/2010 2:37:25 PM}	diffc            = 0.3e-9;
// COMMENT: {3/25/2010 2:37:25 PM}	simul_tr         = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	tempr            = 2.0;
// COMMENT: {3/25/2010 2:37:25 PM}	heat_diffc       = -0.1;
// COMMENT: {3/25/2010 2:37:25 PM}	timest           = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}	multi_Dflag      = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	interlayer_Dflag = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	interlayer_tortf = 100.0;
// COMMENT: {3/25/2010 2:37:25 PM}	interlayer_Dpor  = 0.1;
// COMMENT: {3/25/2010 2:37:25 PM}/* !!!!        count_stag = 0; */
// COMMENT: {3/25/2010 2:37:25 PM}	mcd_substeps       = 1.0;
// COMMENT: {3/25/2010 2:37:25 PM}	print_modulus      = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	punch_modulus      = 1;
// COMMENT: {3/25/2010 2:37:25 PM}	dump_modulus       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	dump_in            = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	transport_warnings = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	pp_assemblage  = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	exchange       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	surface        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	gas_phase      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	kinetics       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	s_s_assemblage = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	cell_data      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	elements       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	elt_list       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	inverse       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_inverse = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	irrev = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	line = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	line_save = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	master = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	mb_unknowns = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	mix       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_mix = 0;
// COMMENT: {3/25/2010 2:37:25 PM}/* !!!! */
// COMMENT: {3/25/2010 2:37:25 PM}	stag_data = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	phases = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	trxn.token = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	s = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	logk = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	master_isotope = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	solution = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	temperature = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	title_x       = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	pe_x          = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	description_x = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	units_x       = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_x           = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	sum_mb1    = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	sum_mb2    = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	sum_jacob0 = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	sum_jacob1 = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	sum_jacob2 = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	sum_delta  = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	isotopes_x = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	x            = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	max_unknowns = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	array     = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	delta     = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	residual  = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_h2o     = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_hplus   = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_h3oplus = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_eminus  = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_co3     = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_h2      = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	s_o2      = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	logk_hash_table           = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	master_isotope_hash_table = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	strings_hash_table        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	elements_hash_table       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	species_hash_table        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	phases_hash_table         = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	keyword_hash_table        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *  Initialize use pointers
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	use.solution_in      = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	use.pp_assemblage_in = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	use.mix_in           = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	use.irrev_in         = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   Initialize punch
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	punch.in               = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_totals     = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.totals           = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_molalities = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.molalities       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_activities = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.activities        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_pure_phases = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.pure_phases = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_si    = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.si          = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_gases = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.gases     = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_s_s = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.s_s = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_kinetics = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.kinetics = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_isotopes = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.isotopes       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.count_calculate_values = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.calculate_values       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	count_save_values = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	save_values       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.inverse = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	punch.sim            = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.state          = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.soln           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.dist           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.time           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.step           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.rxn            = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.temp           = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.ph             = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.pe             = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.alk            = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.mu             = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.water          = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.high_precision = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.user_punch     = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.charge_balance = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	punch.percent_error  = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   last model
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.exchange       = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.gas_phase      = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.s_s_assemblage = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.kinetics       = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.pp_assemblage  = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.add_formula    = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.si             = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.surface_comp   = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	last_model.surface_charge = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   Update hash table
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	keyword_hash = 0;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   rates
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	rates = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	count_rates = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	initial_total_time = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_m = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_m0 = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_p = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_time = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_sim_time_start = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_sim_time_end = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_sim_time = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	rate_moles = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	initial_total_time = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   user_print, user_punch
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	user_print = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	user_punch = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	user_punch_headings = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	user_punch_count_headings = 0;
// COMMENT: {3/25/2010 2:37:25 PM}#ifdef PHREEQ98
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *   user_graph
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	user_graph                = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	user_graph_headings       = 0
// COMMENT: {3/25/2010 2:37:25 PM}	user_graph_count_headings = 0;
// COMMENT: {3/25/2010 2:37:25 PM}#endif
// COMMENT: {3/25/2010 2:37:25 PM}	/*
// COMMENT: {3/25/2010 2:37:25 PM}	   Initialize llnl aqueous model parameters
// COMMENT: {3/25/2010 2:37:25 PM}	 */
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_temp = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_count_temp = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_adh = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_count_adh = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_bdh = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_count_bdh = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_bdot = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_count_bdot = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_co2_coefs = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	llnl_count_co2_coefs = 0;
// COMMENT: {3/25/2010 2:37:25 PM}/*
// COMMENT: {3/25/2010 2:37:25 PM} *
// COMMENT: {3/25/2010 2:37:25 PM} */
// COMMENT: {3/25/2010 2:37:25 PM}	command_hash_table = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	change_surf       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	change_surf_count = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}#if defined(WINDOWS) || defined(_WINDOWS)
// COMMENT: {3/25/2010 2:37:25 PM}	/* SRC pr.status = FALSE; */
// COMMENT: {3/25/2010 2:37:25 PM}#endif
// COMMENT: {3/25/2010 2:37:25 PM}	/* Initialize print here, not in global.h */
// COMMENT: {3/25/2010 2:37:25 PM}	pr.all                = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.initial_solutions  = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.initial_exchangers = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.reactions          = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.gas_phase          = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.s_s_assemblage     = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.pp_assemblage      = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.surface            = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.exchange           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.kinetics           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.totals             = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.eh                 = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.species            = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.saturation_indices = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.irrev              = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.mix                = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.reaction           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.use                = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.logfile            = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.punch              = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	if (phast == TRUE)
// COMMENT: {3/25/2010 2:37:25 PM}	{
// COMMENT: {3/25/2010 2:37:25 PM}		pr.status = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	}
// COMMENT: {3/25/2010 2:37:25 PM}	else
// COMMENT: {3/25/2010 2:37:25 PM}	{
// COMMENT: {3/25/2010 2:37:25 PM}		pr.status = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	}
// COMMENT: {3/25/2010 2:37:25 PM}	pr.inverse            = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.dump               = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.user_print         = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.headings           = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.user_graph         = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.echo_input         = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	count_warnings = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.warnings           = 100;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.initial_isotopes   = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.isotope_ratios     = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.isotope_alphas     = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.hdf                = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	pr.alkalinity         = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	species_list = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	user_database             = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	first_read_input          = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	have_punch_name           = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	selected_output_file_name = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	dump_file_name            = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	/* calculate_value */
// COMMENT: {3/25/2010 2:37:25 PM}	max_calculate_value = MAX_ELTS;
// COMMENT: {3/25/2010 2:37:25 PM}	count_calculate_value = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	calculate_value = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	calculate_value_hash_table = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	/* isotope_ratio */
// COMMENT: {3/25/2010 2:37:25 PM}	max_isotope_ratio = MAX_ELTS;
// COMMENT: {3/25/2010 2:37:25 PM}	count_isotope_ratio = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	isotope_ratio = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	isotope_ratio_hash_table = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	/* isotope_value */
// COMMENT: {3/25/2010 2:37:25 PM}	max_isotope_alpha = MAX_ELTS;
// COMMENT: {3/25/2010 2:37:25 PM}	count_isotope_alpha = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	isotope_alpha = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	isotope_alpha_hash_table = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	phreeqc_mpi_myself = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	copy_solution.n_user       = copy_solution.start       = copy_solution.end       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_pp_assemblage.n_user  = copy_pp_assemblage.start  = copy_pp_assemblage.end  = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_exchange.n_user       = copy_exchange.start       = copy_exchange.end       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_surface.n_user        = copy_surface.start        = copy_surface.end        = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_s_s_assemblage.n_user = copy_s_s_assemblage.start = copy_s_s_assemblage.end = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_gas_phase.n_user      = copy_gas_phase.start      = copy_gas_phase.end      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_kinetics.n_user       = copy_kinetics.start       = copy_kinetics.end       = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_mix.n_user            = copy_mix.start            = copy_mix.end            = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_irrev.n_user          = copy_irrev.start          = copy_irrev.end          = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	copy_temperature.n_user    = copy_temperature.start    = copy_temperature.end    = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	set_forward_output_to_log(FALSE);
// COMMENT: {3/25/2010 2:37:25 PM}	simulation = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	/*
// COMMENT: {3/25/2010 2:37:25 PM}	 *  cvode
// COMMENT: {3/25/2010 2:37:25 PM}	 */
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	cvode_init();
// COMMENT: {3/25/2010 2:37:25 PM}	/*
// COMMENT: {3/25/2010 2:37:25 PM}	 *  Pitzer
// COMMENT: {3/25/2010 2:37:25 PM}	 */
// COMMENT: {3/25/2010 2:37:25 PM}	pitzer_model = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_pitz_param = 100;
// COMMENT: {3/25/2010 2:37:25 PM}	count_pitz_param = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	use_etheta = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	pitz_params = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	max_theta_param = 100;
// COMMENT: {3/25/2010 2:37:25 PM}	count_theta_param = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	theta_params = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	ICON = TRUE;
// COMMENT: {3/25/2010 2:37:25 PM}	OTEMP = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}	for (i = 0; i < 23; i++)
// COMMENT: {3/25/2010 2:37:25 PM}	{
// COMMENT: {3/25/2010 2:37:25 PM}		BK[i] = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}		DK[i] = 0.0;
// COMMENT: {3/25/2010 2:37:25 PM}	}
// COMMENT: {3/25/2010 2:37:25 PM}	pitzer_pe = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	/*
// COMMENT: {3/25/2010 2:37:25 PM}	 *  SIT
// COMMENT: {3/25/2010 2:37:25 PM}	 */
// COMMENT: {3/25/2010 2:37:25 PM}	sit_model       = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	max_sit_param   = 100;
// COMMENT: {3/25/2010 2:37:25 PM}	count_sit_param = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	sit_params      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	/*
// COMMENT: {3/25/2010 2:37:25 PM}	 * to facilitate debuging
// COMMENT: {3/25/2010 2:37:25 PM}	 */
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_use           = &use;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_solution      = solution;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_exchange      = exchange;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_surface       = surface;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_pp_assemblage = pp_assemblage;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_kinetics      = kinetics;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_irrev         = irrev;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_mix           = mix;
// COMMENT: {3/25/2010 2:37:25 PM}	dbg_master        = master;
// COMMENT: {3/25/2010 2:37:25 PM}	calculating_deriv = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}	numerical_deriv   = FALSE;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	zeros     = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	zeros_max = 1;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	cell_pore_volume = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	cell_volume      = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	cell_porosity    = 0;
// COMMENT: {3/25/2010 2:37:25 PM}	cell_saturation  = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	charge_group = NULL;
// COMMENT: {3/25/2010 2:37:25 PM}	print_density = 0;
// COMMENT: {3/25/2010 2:37:25 PM}
// COMMENT: {3/25/2010 2:37:25 PM}	return;
// COMMENT: {3/25/2010 2:37:25 PM}}

bool IPhreeqc::GetOutputOn(void)const
{
	return this->OutputOn;
}

void IPhreeqc::SetOutputOn(bool bValue)
{
	this->OutputOn = bValue;
}

bool IPhreeqc::GetSelectedOutputOn(void)const
{
	return this->SelectedOutputOn;
}

void IPhreeqc::SetSelectedOutputOn(bool bValue)
{
	this->SelectedOutputOn = bValue;
}

bool IPhreeqc::GetLogOn(void)const
{
	return this->LogOn;
}

void IPhreeqc::SetLogOn(bool bValue)
{
	this->LogOn = bValue;
}

bool IPhreeqc::GetDumpOn(void)const
{
	return this->DumpOn;
}

void IPhreeqc::SetDumpOn(bool bValue)
{
	this->DumpOn = bValue;
}

bool IPhreeqc::GetDumpStringOn(void)const
{
	return this->DumpStringOn;
}

void IPhreeqc::SetDumpStringOn(bool bValue)
{
	this->DumpStringOn = bValue;
}

bool IPhreeqc::GetErrorOn(void)const
{
	return this->ErrorOn;
}

void IPhreeqc::SetErrorOn(bool bValue)
{
	this->ErrorOn = bValue;
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

VRESULT IPhreeqc::AccumulateLine(const char *line)
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

int IPhreeqc::RunAccumulated(void)
{
	static const char *sz_routine = "RunAccumulated";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->input_error = 0;

		// create input stream
		std::istringstream iss(this->GetAccumulatedLines());

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL, NULL);
	}
	catch (PhreeqcStop)
	{
		assert(this->database_file == 0);
		assert(this->input_file == 0);
		assert(this->input_error > 0);
	}
	catch(...)
	{
		const char *errmsg = "RunAccumulated: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			assert(this->database_file == 0);
			assert(this->input_file == 0);
			assert(this->input_error > 0);
		}
	}

	this->ClearAccumulatedLines();
	this->close_output_files();
	this->update_errors();

	return this->input_error;
}

int IPhreeqc::RunFile(const char* filename)
{
	static const char *sz_routine = "RunFile";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->input_error = 0;

		// open file
		//
		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open())
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			this->error_msg(oss.str().c_str(), STOP); // throws
		}

		// this may throw
		this->do_run(sz_routine, &ifs, NULL, NULL, NULL, NULL);
	}
	catch (PhreeqcStop)
	{
		assert(this->database_file == 0);
		assert(this->input_file == 0);
		assert(this->input_error > 0);
	}
	catch(...)
	{
		const char *errmsg = "RunFile: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			assert(this->database_file == 0);
			assert(this->input_file == 0);
			assert(this->input_error > 0);
		}
	}

	this->close_output_files();
	this->update_errors();

	return this->input_error;
}

int IPhreeqc::RunString(const char* input)
{
	static const char *sz_routine = "RunString";
	try
	{
		// this may throw
		this->check_database(sz_routine);

		this->input_error = 0;

		// create input stream
		std::string s(input);
		std::istringstream iss(s);

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL, NULL);
	}
	catch (PhreeqcStop)
	{
		assert(this->database_file == 0);
		assert(this->input_file == 0);
		assert(this->input_error > 0);
	}
	catch(...)
	{
		const char *errmsg = "RunString: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			assert(this->database_file == 0);
			assert(this->input_file == 0);
			assert(this->input_error > 0);
		}
	}

	this->close_output_files();
	this->update_errors();

	return this->input_error;
}

int IPhreeqc::GetSelectedOutputRowCount(void)const
{
	return (int)this->SelectedOutput->GetRowCount();
}

int IPhreeqc::GetSelectedOutputColumnCount(void)const
{
	return (int)this->SelectedOutput->GetColCount();
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

int IPhreeqc::GetDumpLineCount(void)const
{
	return (int)this->DumpLines.size();
}

const char* IPhreeqc::GetDumpLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetDumpLineCount())
	{
		return empty;
	}
	return this->DumpLines[n].c_str();
}

int IPhreeqc::GetErrorLineCount(void)const
{
	return (int)this->ErrorLines.size();
}

const char* IPhreeqc::GetErrorLine(int n)
{
	static const char empty[] = "";
	if (n < 0 || n >= this->GetErrorLineCount())
	{
		return empty;
	}
	return this->ErrorLines[n].c_str();
}

void IPhreeqc::update_errors(void)
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
std::list< std::string > IPhreeqc::ListComponents(void)
{
	std::list< std::string > comps;
	this->list_components(comps);
	return comps;
}
