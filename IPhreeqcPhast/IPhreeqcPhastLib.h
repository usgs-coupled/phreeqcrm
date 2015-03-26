#if !defined(IPHREEQC_PHAST_LIB_H_INCLUDED)
#define IPHREEQC_PHAST_LIB_H_INCLUDED
#include "IPhreeqc.h"
#include "IPhreeqcPhast.h"
class IPhreeqcPhastLib
{
public:
	static void CleanupIPhreeqcPhast(void);
	static int CreateIPhreeqcPhast(void);
	static IPQ_RESULT DestroyIPhreeqcPhast(int n);
	static IPhreeqcPhast* GetInstance(int n);
};

#endif // !defined(IPHREEQC_PHAST_LIB_H_INCLUDED)
