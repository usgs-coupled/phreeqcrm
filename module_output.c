#include "module_files.h"

#if defined(R_SO)
#include "output.inl"
#else
#include "phreeqcpp/phreeqc/output.c"
#endif


/* ---------------------------------------------------------------------- */
int output_isopen(const int type)
/* ---------------------------------------------------------------------- */
{
	size_t i;
	int isopen;
	for (i = 0; i < count_output_callback; ++i) {
		isopen = (output_callbacks[i].callback)(ACTION_ISOPEN, type, NULL, CONTINUE, output_callbacks[i].cookie, NULL, NULL);
		if (isopen) return 1;
	}
	return 0;
}
