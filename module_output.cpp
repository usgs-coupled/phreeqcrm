#include "module_files.h"

#include "IPhreeqc.hpp"

#if defined(R_SO)
#include "output.inl"
#else
// COMMENT: {3/25/2010 12:36:24 PM}#include "phreeqcpp/phreeqc/output.c"
#endif


int IPhreeqc::output_isopen(const int type)
{
	size_t i;
	int isopen;
	for (i = 0; i < this->count_output_callback; ++i)
	{
		isopen = (this->output_callbacks[i].callback)(ACTION_ISOPEN, type, NULL, CONTINUE, this->output_callbacks[i].cookie, NULL, NULL);
		if (isopen) return 1;
	}
	return 0;
}
