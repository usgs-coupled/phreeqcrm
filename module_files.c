#include "module_files.h"
#if defined(R_SO)
#include "phreeqc_files.inl"
#else
#include "phreeqcpp/phreeqc/phreeqc_files.c"
#endif

static int module_isopen_handler(const int type);
static int module_open_handler(const int type, const char *file_name);

int
module_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
{
	switch (action) {
	case ACTION_OPEN:
		return module_open_handler(type, err_str);
		break;
	case ACTION_ISOPEN:
		return module_isopen_handler(type);
		break;
	default:
		return phreeqc_handler(action, type, err_str, stop, cookie, format, args);
		break;
	}
	return ERROR;
}

static int
module_isopen_handler(const int type)
{
	switch (type)
	{
	case OUTPUT_PUNCH:
		if (punch_file) return 1;
		break;
	default:
		assert(0);
	}
	return 0;
}

static int 
module_open_handler(const int type, const char *file_name)
{
	assert(file_name && strlen(file_name));
	switch (type)
	{
	case OUTPUT_MESSAGE:
		if (output != NULL)
		{
			fclose(output);
			output = NULL;
		}
		if ( (output = fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case OUTPUT_ERROR:
		assert(error_file != stderr);
		if (error_file != NULL)
		{
			fclose(error_file);
			error_file = NULL;
		}
		if ( (error_file = fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;

	case OUTPUT_LOG:
		if (log_file != NULL)
		{
			fclose(log_file);
			log_file = NULL;
		}
		if ( (log_file = fopen(file_name, "w")) == NULL)
		{
			return ERROR;
		}
		break;
	default:
		return open_handler(type, file_name);
		break;
	}
	return(OK);
}
