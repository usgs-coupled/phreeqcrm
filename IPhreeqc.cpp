#include <memory>
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
Run(void)
{
	return IPhreeqc::LibraryInstance()->Run();
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
	this->init();
	this->UnLoadDatabase();
}

IPhreeqc::~IPhreeqc(void)
{
	delete this->ErrorReporter;
	delete this->WarningReporter;
	delete this->SelectedOutput;
}

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
		FILE* f = fopen(filename, "r");
		if (!f)
		{
			std::ostringstream oss;
			oss << "LoadDatabase: Unable to open:" << "\"" << filename << "\".";
			this->error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
		}

		// read input
		//
		this->read_database(getc_callback, f);

	}
	catch (PhreeqcStop)
	{
		// do nothing
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
			// do nothing
		}
	}

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
		// do nothing
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
			// do nothing
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
	IPhreeqc* pIPhreeqc = (IPhreeqc*) cookie;

	switch (type)
	{
	case OUTPUT_ERROR:
		if (pIPhreeqc)
		{
			std::ostringstream oss;
			oss << "ERROR: " << err_str << "\n";
			if (stop == STOP)
			{
				oss << "Stopping.\n";
			}
			pIPhreeqc->AddError(oss.str().c_str());
		}
		break;

	case OUTPUT_WARNING:
		if (pIPhreeqc)
		{
			std::ostringstream oss;
			oss << "WARNING: " << err_str << "\n";
			pIPhreeqc->AddWarning(oss.str().c_str());
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

void IPhreeqc::init(void)
{
	int i;

	moles_per_kilogram_string = 0;
	pe_string = 0;

	debug_model = FALSE;
	debug_prep = FALSE;
	debug_set = FALSE;
	debug_diffuse_layer = FALSE;
	debug_inverse = FALSE;
	itmax = 100;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	ineq_tol = pow((long double) 10, (long double) -LDBL_DIG);
#else
	ineq_tol = pow((double) 10, (double) -DBL_DIG);
#endif
	convergence_tolerance = 1e-8;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	inv_tol_default = pow((long double) 10, (long double) -LDBL_DIG + 5);
#else
	inv_tol_default = pow((double) 10, (double) -DBL_DIG + 5);
#endif
	step_size               = 100.;
	pe_step_size            = 10.;
	pp_scale                = 1.0;
	pp_column_scale         = 1.0;
	diagonal_scale          = FALSE;
	censor                  = 0.0;
	mass_water_switch       = FALSE;
	delay_mass_water        = FALSE;
	incremental_reactions   = FALSE;
	aqueous_only            = 0;
	negative_concentrations = FALSE;

	LOG_10 = log(10.0);

	max_solution       = MAX_SOLUTION;
	max_pp_assemblage  = MAX_PP_ASSEMBLAGE;
	max_exchange       = MAX_PP_ASSEMBLAGE;
	max_surface        = MAX_PP_ASSEMBLAGE;
	max_gas_phase      = MAX_PP_ASSEMBLAGE;
	max_kinetics       = MAX_PP_ASSEMBLAGE;
	max_s_s_assemblage = MAX_PP_ASSEMBLAGE;

	max_elements       = MAX_ELEMENTS;
	max_elts           = MAX_ELTS;
	max_line           = MAX_LINE;
	max_master         = MAX_MASTER;
	max_mb_unknowns    = MAX_TRXN;
	max_phases         = MAX_PHASES;
	max_s              = MAX_S;
	max_strings        = MAX_STRINGS;
	max_trxn           = MAX_TRXN;
	max_logk           = MAX_S;
	max_master_isotope = MAX_ELTS;

	count_solution       = 0;
	count_pp_assemblage  = 0;
	count_exchange       = 0;
	count_surface        = 0;
	count_gas_phase      = 0;
	count_kinetics       = 0;
	count_s_s_assemblage = 0;

	count_elements       = 0;
	count_irrev          = 0;
	count_master         = 0;
	count_mix            = 0;
	count_phases         = 0;
	count_s              = 0;
	count_temperature    = 0;
	count_logk           = 0;
	count_master_isotope = 0;
/*
 *   Initialize advection
 */
	count_ad_cells   = 1;
	count_ad_shifts  = 1;
	print_ad_modulus = 1;
	punch_ad_modulus = 1;

	advection_punch            = 0;
	advection_kin_time         = 0.0;
	advection_kin_time_defined = FALSE;
	advection_print            = 0;
	advection_warnings         = TRUE;
/*
 *   Initialize transport
 */
	count_cells      = 1;
	count_shifts     = 1;
	ishift           = 1;
	bcon_first       = bcon_last = 3;
	diffc            = 0.3e-9;
	simul_tr         = 0;
	tempr            = 2.0;
	heat_diffc       = -0.1;
	timest           = 0.0;
	multi_Dflag      = FALSE;
	interlayer_Dflag = FALSE;
	interlayer_tortf = 100.0;
	interlayer_Dpor  = 0.1;
/* !!!!        count_stag = 0; */
	mcd_substeps       = 1.0;
	print_modulus      = 1;
	punch_modulus      = 1;
	dump_modulus       = 0;
	dump_in            = FALSE;
	transport_warnings = TRUE;

	pp_assemblage  = 0;
	exchange       = 0;
	surface        = 0;
	gas_phase      = 0;
	kinetics       = 0;
	s_s_assemblage = 0;
	cell_data      = 0;
	elements       = 0;
	elt_list       = 0;


	inverse       = 0;
	count_inverse = 0;

	irrev = 0;

	line = 0;
	line_save = 0;

	master = 0;

	mb_unknowns = 0;

	mix       = 0;
	count_mix = 0;
/* !!!! */
	stag_data = 0;

	phases = 0;

	trxn.token = 0;

	s = 0;

	logk = 0;

	master_isotope = 0;

	solution = 0;

	temperature = 0;

	title_x       = NULL;
	pe_x          = NULL;
	description_x = NULL;
	units_x       = NULL;
	s_x           = NULL;

	sum_mb1    = NULL;
	sum_mb2    = NULL;
	sum_jacob0 = NULL;
	sum_jacob1 = NULL;
	sum_jacob2 = NULL;
	sum_delta  = NULL;

	isotopes_x = 0;

	x            = NULL;
	max_unknowns = 0;

	array     = NULL;
	delta     = NULL;
	residual  = NULL;
	s_h2o     = NULL;
	s_hplus   = NULL;
	s_h3oplus = NULL;
	s_eminus  = NULL;
	s_co3     = NULL;
	s_h2      = NULL;
	s_o2      = NULL;

	logk_hash_table           = 0;
	master_isotope_hash_table = 0;
	strings_hash_table        = 0;
	elements_hash_table       = 0;
	species_hash_table        = 0;
	phases_hash_table         = 0;
	keyword_hash_table        = 0;
/*
 *  Initialize use pointers
 */
	use.solution_in      = FALSE;
	use.pp_assemblage_in = FALSE;
	use.mix_in           = FALSE;
	use.irrev_in         = FALSE;
/*
 *   Initialize punch
 */
	punch.in               = FALSE;
	punch.count_totals     = 0;
	punch.totals           = 0;
	punch.count_molalities = 0;

	punch.molalities       = 0;
	punch.count_activities = 0;

	punch.activities        = 0;
	punch.count_pure_phases = 0;

	punch.pure_phases = 0;
	punch.count_si    = 0;

	punch.si          = 0;
	punch.count_gases = 0;

	punch.gases     = 0;
	punch.count_s_s = 0;
	punch.s_s = 0;

	punch.count_kinetics = 0;
	punch.kinetics = 0;

	punch.count_isotopes = 0;
	punch.isotopes       = 0;

	punch.count_calculate_values = 0;
	punch.calculate_values       = 0;

	count_save_values = 0;
	save_values       = 0;


	punch.inverse = TRUE;

	punch.sim            = TRUE;
	punch.state          = TRUE;
	punch.soln           = TRUE;
	punch.dist           = TRUE;
	punch.time           = TRUE;
	punch.step           = TRUE;
	punch.rxn            = FALSE;
	punch.temp           = FALSE;
	punch.ph             = TRUE;
	punch.pe             = TRUE;
	punch.alk            = FALSE;
	punch.mu             = FALSE;
	punch.water          = FALSE;
	punch.high_precision = FALSE;
	punch.user_punch     = TRUE;
	punch.charge_balance = FALSE;
	punch.percent_error  = FALSE;
/*
 *   last model
 */
	last_model.exchange       = NULL;
	last_model.gas_phase      = NULL;
	last_model.s_s_assemblage = NULL;
	last_model.kinetics       = NULL;
	last_model.pp_assemblage  = NULL;
	last_model.add_formula    = NULL;
	last_model.si             = NULL;
	last_model.surface_comp   = NULL;
	last_model.surface_charge = NULL;
/*
 *   Update hash table
 */
	keyword_hash = 0;
/*
 *   rates
 */
	rates = 0;
	count_rates = 0;
	initial_total_time = 0;
	rate_m = 0;
	rate_m0 = 0;
	rate_p = NULL;
	rate_time = 0;
	rate_sim_time_start = 0;
	rate_sim_time_end = 0;
	rate_sim_time = 0;
	rate_moles = 0;
	initial_total_time = 0;

/*
 *   user_print, user_punch
 */
	user_print = 0;
	user_punch = 0;
	user_punch_headings = 0;
	user_punch_count_headings = 0;
#ifdef PHREEQ98
/*
 *   user_graph
 */
	user_graph                = 0;
	user_graph_headings       = 0
	user_graph_count_headings = 0;
#endif
	/*
	   Initialize llnl aqueous model parameters
	 */
	llnl_temp = 0;
	llnl_count_temp = 0;

	llnl_adh = 0;
	llnl_count_adh = 0;

	llnl_bdh = 0;
	llnl_count_bdh = 0;

	llnl_bdot = 0;
	llnl_count_bdot = 0;

	llnl_co2_coefs = 0;
	llnl_count_co2_coefs = 0;
/*
 *
 */
	command_hash_table = 0;

	change_surf       = 0;
	change_surf_count = 0;


#if defined(WINDOWS) || defined(_WINDOWS)
	/* SRC pr.status = FALSE; */
#endif
	/* Initialize print here, not in global.h */
	pr.all                = TRUE;
	pr.initial_solutions  = TRUE;
	pr.initial_exchangers = TRUE;
	pr.reactions          = TRUE;
	pr.gas_phase          = TRUE;
	pr.s_s_assemblage     = TRUE;
	pr.pp_assemblage      = TRUE;
	pr.surface            = TRUE;
	pr.exchange           = TRUE;
	pr.kinetics           = TRUE;
	pr.totals             = TRUE;
	pr.eh                 = TRUE;
	pr.species            = TRUE;
	pr.saturation_indices = TRUE;
	pr.irrev              = TRUE;
	pr.mix                = TRUE;
	pr.reaction           = TRUE;
	pr.use                = TRUE;
	pr.logfile            = FALSE;
	pr.punch              = TRUE;
	if (phast == TRUE)
	{
		pr.status = FALSE;
	}
	else
	{
		pr.status = TRUE;
	}
	pr.inverse            = TRUE;
	pr.dump               = TRUE;
	pr.user_print         = TRUE;
	pr.headings           = TRUE;
	pr.user_graph         = TRUE;
	pr.echo_input         = TRUE;
	count_warnings = 0;
	pr.warnings           = 100;
	pr.initial_isotopes   = TRUE;
	pr.isotope_ratios     = TRUE;
	pr.isotope_alphas     = TRUE;
	pr.hdf                = FALSE;
	pr.alkalinity         = FALSE;

	species_list = NULL;

	user_database             = NULL;
	first_read_input          = TRUE;
	have_punch_name           = FALSE;
	selected_output_file_name = NULL;
	dump_file_name            = NULL;

	/* calculate_value */
	max_calculate_value = MAX_ELTS;
	count_calculate_value = 0;

	calculate_value = 0;
	calculate_value_hash_table = 0;

	/* isotope_ratio */
	max_isotope_ratio = MAX_ELTS;
	count_isotope_ratio = 0;
	isotope_ratio = 0;
	isotope_ratio_hash_table = 0;

	/* isotope_value */
	max_isotope_alpha = MAX_ELTS;
	count_isotope_alpha = 0;
	isotope_alpha = 0;
	isotope_alpha_hash_table = 0;

	phreeqc_mpi_myself = 0;

	copy_solution.n_user       = copy_solution.start       = copy_solution.end       = 0;
	copy_pp_assemblage.n_user  = copy_pp_assemblage.start  = copy_pp_assemblage.end  = 0;
	copy_exchange.n_user       = copy_exchange.start       = copy_exchange.end       = 0;
	copy_surface.n_user        = copy_surface.start        = copy_surface.end        = 0;
	copy_s_s_assemblage.n_user = copy_s_s_assemblage.start = copy_s_s_assemblage.end = 0;
	copy_gas_phase.n_user      = copy_gas_phase.start      = copy_gas_phase.end      = 0;
	copy_kinetics.n_user       = copy_kinetics.start       = copy_kinetics.end       = 0;
	copy_mix.n_user            = copy_mix.start            = copy_mix.end            = 0;
	copy_irrev.n_user          = copy_irrev.start          = copy_irrev.end          = 0;
	copy_temperature.n_user    = copy_temperature.start    = copy_temperature.end    = 0;

	set_forward_output_to_log(FALSE);
	simulation = 0;
	/*
	 *  cvode
	 */

	cvode_init();
	/*
	 *  Pitzer
	 */
	pitzer_model = FALSE;
	max_pitz_param = 100;
	count_pitz_param = 0;
	use_etheta = TRUE;
	pitz_params = 0;

	max_theta_param = 100;
	count_theta_param = 0;
	theta_params = 0;

	ICON = TRUE;
	OTEMP = 0.0;
	for (i = 0; i < 23; i++)
	{
		BK[i] = 0.0;
		DK[i] = 0.0;
	}
	pitzer_pe = FALSE;


	/*
	 *  SIT
	 */
	sit_model       = FALSE;
	max_sit_param   = 100;
	count_sit_param = 0;
	sit_params      = 0;

	/*
	 * to facilitate debuging
	 */
	dbg_use           = &use;
	dbg_solution      = solution;
	dbg_exchange      = exchange;
	dbg_surface       = surface;
	dbg_pp_assemblage = pp_assemblage;
	dbg_kinetics      = kinetics;
	dbg_irrev         = irrev;
	dbg_mix           = mix;
	dbg_master        = master;
	calculating_deriv = FALSE;
	numerical_deriv   = FALSE;

	zeros     = 0;
	zeros_max = 1;

	cell_pore_volume = 0;
	cell_volume      = 0;
	cell_porosity    = 0;
	cell_saturation  = 0;

	charge_group = NULL;
	print_density = 0;

	return;
}

void IPhreeqc::SetOutputOn(bool bValue)
{
	this->OutputOn = bValue;
}

void IPhreeqc::SetSelectedOutputOn(bool bValue)
{
	this->SelectedOutputOn = bValue;
}

void IPhreeqc::SetLogOn(bool bValue)
{
	this->LogOn = bValue;
}

void IPhreeqc::SetDumpOn(bool bValue)
{
	this->DumpOn = bValue;
}

void IPhreeqc::SetDumpStringOn(bool bValue)
{
	this->DumpStringOn = bValue;
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

int IPhreeqc::Run(void)
{
	static const char *sz_routine = "Run";
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
		// do nothing
	}
	catch(...)
	{
		const char *errmsg = "Run: An unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
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

#if 0
		// create input stream
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
#else
		// open file
		//
		FILE* f = ::fopen(filename, "r");
		if (!f)
		{
			std::ostringstream oss;
			oss << "RunFile: Unable to open:" << "\"" << filename << "\".";
			this->error_msg(oss.str().c_str(), STOP); // throws PhreeqcStop
		}

		// this may throw
		this->do_run(sz_routine, NULL, f, NULL, NULL, NULL);
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
			this->error_msg(errmsg, STOP); // throws PhreeqcStop
		}
		catch (PhreeqcStop)
		{
			// do nothing
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
		// do nothing
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
			// do nothing
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
