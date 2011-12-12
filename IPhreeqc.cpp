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
, input_file(0)
, database_file(0)
{
	this->ErrorReporter   = new CErrorReporter<std::ostringstream>;
	this->WarningReporter = new CErrorReporter<std::ostringstream>;
	this->SelectedOutput  = new CSelectedOutput();
	this->PhreeqcPtr      = new Phreeqc(this);

	ASSERT(this->PhreeqcPtr->phast == 0);
	this->UnLoadDatabase();
}

IPhreeqc::~IPhreeqc(void)
{
#ifdef _DEBUG
	this->OutputOn = false;
#endif
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

size_t IPhreeqc::AddError(const char* str)
{
	return this->ErrorReporter->AddError(str);
}

size_t IPhreeqc::AddWarning(const char* str)
{
	return this->WarningReporter->AddError(str);
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
	bool bSaveOutputOn = this->OutputOn;
	this->OutputOn = false;
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
		this->PhreeqcPtr->phrq_io->push_istream(&ifs, false);
		this->PhreeqcPtr->read_database();
	}
	catch (IPhreeqcStop)
	{
		this->close_input_files();
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
	this->PhreeqcPtr->phrq_io->clear_istream();
	this->DatabaseLoaded = (this->PhreeqcPtr->get_input_errors() == 0);
	this->OutputOn = bSaveOutputOn;
	return this->PhreeqcPtr->get_input_errors();
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
		ASSERT(this->PhreeqcPtr->phrq_io->get_istream() == NULL);
		this->PhreeqcPtr->phrq_io->push_istream(&iss, false);
		this->PhreeqcPtr->read_database();
	}
	catch (IPhreeqcStop)
	{
		this->close_input_files();
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

	this->PhreeqcPtr->phrq_io->clear_istream();
	this->DatabaseLoaded = (this->PhreeqcPtr->get_input_errors() == 0);
	return this->PhreeqcPtr->get_input_errors();
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
		this->io_error_count = 0;

		// create input stream
		std::istringstream iss(this->GetAccumulatedLines());

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL);
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
	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
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
		this->io_error_count = 0;

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
		this->do_run(sz_routine, &ifs, NULL, NULL, NULL);
	}
	catch (IPhreeqcStop)
	{
		this->close_input_files();
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

	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
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
		this->io_error_count = 0;

		// create input stream
		std::string s(input);
		std::istringstream iss(s);

		// this may throw
		this->do_run(sz_routine, &iss, NULL, NULL, NULL);
	}
	catch (IPhreeqcStop)
	{
		this->close_input_files();
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

	this->close_output_files();
	this->update_errors();
	this->PhreeqcPtr->phrq_io->clear_istream();

	return this->PhreeqcPtr->get_input_errors();
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
	this->PhreeqcPtr->do_initialize();
	this->PhreeqcPtr->input_error = 0;
	this->io_error_count = 0;
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

void IPhreeqc::do_run(const char* sz_routine, std::istream* pis, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
{
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
	std::auto_ptr<std::istringstream> auto_iss(0);
	if (!pis)
	{
		auto_iss.reset(new std::istringstream(this->GetAccumulatedLines()));
		this->PhreeqcPtr->phrq_io->push_istream(auto_iss.get(), false);
	}
	else
	{
		ASSERT(this->PhreeqcPtr->phrq_io->get_istream() == NULL);
		this->PhreeqcPtr->phrq_io->push_istream(pis, false);
	}


/*
 *   Read input data for simulation
 */
	for (this->PhreeqcPtr->simulation = 1; ; this->PhreeqcPtr->simulation++)
	{

#ifdef PHREEQ98
   		AddSeries = !connect_simulations;
#endif
		::sprintf(token, "Reading input data for simulation %d.", this->PhreeqcPtr->simulation);

		int save_punch_in = this->PhreeqcPtr->punch.in;

		this->PhreeqcPtr->dup_print(token, TRUE);
		if (this->PhreeqcPtr->read_input() == EOF)
			break;

		if (this->PhreeqcPtr->simulation > 1 && save_punch_in == TRUE && this->PhreeqcPtr->punch.new_def == TRUE)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning SELECTED_OUTPUT has been redefined.\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());

		}
		if (this->PhreeqcPtr->simulation > 1 && this->PhreeqcPtr->keycount[Keywords::KEY_USER_PUNCH] > 0)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Warning USER_PUNCH has been redefined.\n";
			this->PhreeqcPtr->warning_msg(oss.str().c_str());
		}

		if (this->PhreeqcPtr->title_x != NULL)
		{
			::sprintf(token, "TITLE");
			this->PhreeqcPtr->dup_print(token, TRUE);
			if (this->PhreeqcPtr->pr.headings == TRUE)
			{
				char *p = this->PhreeqcPtr->sformatf("%s\n\n", this->PhreeqcPtr->title_x);
				this->PhreeqcPtr->output_msg(p);
			}
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
			if (!this->SelectedOutputOn)
			{
				ASSERT(!this->punch_ostream);
			}

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
					if (this->SelectedOutputOn && !this->punch_ostream)
					{
						//
						// LoadDatabase
						// do_run -- containing SELECTED_OUTPUT ****TODO**** check -file option
						// another do_run without SELECTED_OUTPUT
						//
						std::string filename = this->PunchFileName;
						if (!this->punch_open(filename.c_str()))
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
					if (this->SelectedOutputOn && !this->punch_ostream)
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
						if (!this->punch_open(filename.c_str()))
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

		if (!this->SelectedOutputOn) ASSERT(!this->punch_ostream);
		/* the converse is not necessarily true */

		this->PhreeqcPtr->n_user_punch_index = -1;
#endif // SWIG_SHARED_OBJ
		this->PhreeqcPtr->tidy_model();
#ifdef PHREEQ98
                if (!phreeq98_debug)
				{
#endif

/*
 *   Calculate distribution of species for initial solutions
 */
		if (this->PhreeqcPtr->new_solution)
			this->PhreeqcPtr->initial_solutions(TRUE);
/*
 *   Calculate distribution for exchangers
 */
		if (this->PhreeqcPtr->new_exchange)
			this->PhreeqcPtr->initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
		if (this->PhreeqcPtr->new_surface)
			this->PhreeqcPtr->initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
		if (this->PhreeqcPtr->new_gas_phase)
			this->PhreeqcPtr->initial_gas_phases(TRUE);
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
		if (this->PhreeqcPtr->use.advect_in == TRUE)
		{
			this->PhreeqcPtr->dup_print("Beginning of advection calculations.", TRUE);
			this->PhreeqcPtr->advection();
		}
/*
 *   Calculate transport
 */
		if (this->PhreeqcPtr->use.trans_in == TRUE)
		{
			this->PhreeqcPtr->dup_print("Beginning of transport calculations.", TRUE);
			this->PhreeqcPtr->transport();
		}
/*
 *   run
 */
		this->PhreeqcPtr->run_as_cells();
/*
 *   Copy
 */
		if (this->PhreeqcPtr->new_copy) this->PhreeqcPtr->copy_entities();
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
				if (this->PhreeqcPtr->dump_info.Get_append())
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

void IPhreeqc::error_msg(const char *str, bool stop)
{
	ASSERT(!(this->ErrorOn ^ (this->error_ostream != 0)));
	this->PHRQ_io::error_msg(str);

	this->AddError(str);
	if (stop)
	{
		if (this->error_ostream)
		{
			(*this->error_ostream) << "Stopping.\n";
			this->error_ostream->flush();
		}
		throw IPhreeqcStop();
	}
}

void IPhreeqc::warning_msg(const char *str)
{
	ASSERT(!(this->ErrorOn ^ (this->error_ostream != 0)));
	this->PHRQ_io::warning_msg(str);

	std::ostringstream oss;
	oss << str << std::endl;
	this->AddWarning(oss.str().c_str());
}

void IPhreeqc::output_msg(const char * str)
{
	ASSERT(!(this->OutputOn ^ (this->output_ostream != 0)));
	this->PHRQ_io::output_msg(str);
}

void IPhreeqc::screen_msg(const char *err_str)
{
	// no-op
}

void IPhreeqc::punch_msg(const char *str)
{
	this->PHRQ_io::punch_msg(str);
}

void IPhreeqc::open_output_files(const char* sz_routine)
{
	if (this->OutputOn)
	{
		if (this->output_ostream != NULL)
		{
			delete this->output_ostream;
			this->output_ostream = NULL;
		}
		if ( (this->output_ostream = new std::ofstream(OUTPUT_FILENAME)) == NULL)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << OUTPUT_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
		}
	}
	if (this->ErrorOn)
	{
		if (this->error_ostream != NULL)
		{
			delete this->error_ostream;
			this->error_ostream = NULL;
		}
		if ( (this->error_ostream = new std::ofstream(ERROR_FILENAME)) == NULL)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << ERROR_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
		}
	}
	if (this->LogOn)
	{
		if (this->log_ostream != NULL)
		{
			delete this->log_ostream;
			this->log_ostream = NULL;
		}
		if ( (this->log_ostream = new std::ofstream(LOG_FILENAME)) == NULL)
		{
			std::ostringstream oss;
			oss << sz_routine << ": Unable to open:" << "\"" << LOG_FILENAME << "\".\n";
			this->warning_msg(oss.str().c_str());
		}
	}
}

int IPhreeqc::close_input_files(void)
{
	int i = 0;
	if (this->database_file)
	{
		i |= fclose(this->database_file);
	}
	if (this->input_file)
	{
		i |= fclose(this->input_file);
	}
	this->input_file = this->database_file = NULL;
	return (i);
}

int IPhreeqc::close_output_files(void)
{
	int ret = 0;

	if (this->output_ostream != NULL)
		delete this->output_ostream;
	if (this->log_ostream != NULL)
		delete this->log_ostream;
	if (this->punch_ostream != NULL)
		delete this->punch_ostream;
	if (this->dump_ostream != NULL)
		delete this->dump_ostream;
	if (this->error_ostream != NULL)
		delete this->error_ostream;
	this->error_ostream = NULL;
	this->output_ostream = this->log_ostream = this->punch_ostream = NULL;
	this->dump_ostream = NULL;
	return ret;
}

void IPhreeqc::fpunchf(const char *name, const char *format, double d)
{
	this->PHRQ_io::fpunchf(name, format, d);
	this->SelectedOutput->PushBackDouble(name, d);
}

void IPhreeqc::fpunchf(const char *name, const char *format, char *s)
{
	this->PHRQ_io::fpunchf(name, format, s);
	this->SelectedOutput->PushBackString(name, s);
}

void IPhreeqc::fpunchf(const char *name, const char *format, int i)
{
	this->PHRQ_io::fpunchf(name, format, i);
	this->SelectedOutput->PushBackLong(name, (long)i);
}

void IPhreeqc::fpunchf_end_row(const char *format)
{
	this->EndRow();
}

bool IPhreeqc::punch_open(const char *file_name, std::ios_base::openmode mode)
{
	if (file_name)
	{
		this->PunchFileName = file_name;
	}
	if (this->SelectedOutputOn)
	{
		return this->PHRQ_io::punch_open(file_name, mode);
	}
	return true;
}

bool IPhreeqc::output_open(const char *file_name, std::ios_base::openmode mode)
{
	if (file_name)
	{
		//this->PunchFileName = file_name;
	}
	if (this->OutputOn)
	{
		return this->PHRQ_io::output_open(file_name, mode);
	}
	return true;
}
