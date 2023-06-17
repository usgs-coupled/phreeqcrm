#ifndef IRM_DLL_EXPORT
# if defined(_WIN32)
#  if defined(DLL_EXPORT)
#   define IRM_DLL_EXPORT __declspec(dllexport)
#  else
#   define IRM_DLL_EXPORT
#  endif
# elif defined(__GNUC__) && (__GNUC__ >= 4)
#  define IRM_DLL_EXPORT __attribute__ ((visibility ("default")))
# else
#  define IRM_DLL_EXPORT
# endif
#endif
