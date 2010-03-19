#include "module_files.h"
#if defined(R_SO)
#include "phreeqc_files.inl"
#else
#include "phreeqcpp/phreeqc/phreeqc_files.c"
#endif

// COMMENT: {3/2/2010 4:06:35 PM}static int module_isopen_handler(const int type);
// COMMENT: {3/2/2010 4:06:35 PM}static int module_open_handler(const int type, const char *file_name);

#include "IPhreeqc.hpp"

int IPhreeqc::module_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	IPhreeqc* pThis = (IPhreeqc*) cookie;

	switch (action) {
	case ACTION_OPEN:
		return pThis->module_open_handler(type, err_str);
		break;
	case ACTION_ISOPEN:
		return pThis->module_isopen_handler(type);
		break;
	default:
		return pThis->phreeqc_handler(action, type, err_str, stop, cookie, format, args);
		break;
	}
	return ERROR;
}

int IPhreeqc::module_isopen_handler(const int type)
{
	switch (type)
	{
	case OUTPUT_PUNCH:
		if (this->punch_file) return 1;
		break;
	default:
		assert(0);
	}
	return 0;
}

int IPhreeqc::module_open_handler(const int type, const char *file_name)
{
	assert(file_name && ::strlen(file_name));
	switch (type)
	{
	case OUTPUT_MESSAGE:
		if (this->output != NULL)
		{
			::fclose(this->output);
			this->output = NULL;
		}
		if ( (this->output = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case OUTPUT_ERROR:
		assert(this->error_file != stderr);
		if (this->error_file != NULL)
		{
			::fclose(this->error_file);
			this->error_file = NULL;
		}
		if ( (this->error_file = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case OUTPUT_LOG:
		if (this->log_file != NULL)
		{
			::fclose(this->log_file);
			this->log_file = NULL;
		}
		if ( (this->log_file = ::fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	default:
		return this->Phreeqc::open_handler(type, file_name);
		break;

	}
	return(OK);
}
