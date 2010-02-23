#ifndef __MODULE_FILES__H
#define __MODULE_FILES__H

#include <stdio.h>
#include <stdarg.h>


#if defined(__cplusplus)
extern "C" {
#endif

int module_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);
int output_isopen(const int type);


typedef enum {
	ACTION_ISOPEN = -100,
} module_action_type;


#if defined(__cplusplus)
}
#endif

#endif  /* __MODULE_FILES__H */
