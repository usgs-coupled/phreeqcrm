#include <cassert>
#include <iostream>
#include <map>

#include "IPhreeqc.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
#include "IPhreeqcPhastLib.h"
//#define CreateIPhreeqcPhast          createiphreeqcphast
//class IPhreeqcPhastLib
//{
//public:
//	static int CreateIPhreeqcPhast(void);
//	static IPQ_RESULT DestroyIPhreeqcPhast(int n);
//	static IPhreeqcPhast* GetInstance(int n);
//};


// helper functions
//

int
IPhreeqcPhastLib::CreateIPhreeqcPhast(void)
{
	int n = IPQ_OUTOFMEMORY;
	IPhreeqcPhast* IPhreeqcPhastPtr;
	try
	{
#ifdef USE_OPENMP
		#pragma omp critical(IPhreeqcPhastLib)
#endif
		{
			IPhreeqcPhastPtr = new IPhreeqcPhast;
		}
		n = (int) IPhreeqcPhastPtr->Index;
	}
	catch(...)
	{
		return IPQ_OUTOFMEMORY;
	}
	
	return n;
}

IPQ_RESULT
IPhreeqcPhastLib::DestroyIPhreeqcPhast(int id)
{
	IPQ_RESULT retval = IPQ_BADINSTANCE;
	if (id >= 0)
	{
		if (IPhreeqc *ptr = IPhreeqcPhastLib::GetInstance(id))
		{
#ifdef USE_OPENMP
			#pragma omp critical(IPhreeqcPhastLib)
#endif
			{
				delete ptr;
			}
				retval = IPQ_OK;
		}
	}
	return retval;
}

IPhreeqcPhast*
IPhreeqcPhastLib::GetInstance(int id)
{
	std::map<size_t, IPhreeqcPhast*>::iterator it;
	bool found=false;
#ifdef USE_OPENMP
	#pragma omp critical(IPhreeqcLib)
#endif
	{
		it = IPhreeqcPhast::PhastInstances.find(size_t(id));
		found = (it != IPhreeqcPhast::PhastInstances.end());
	}
	if (found)
	{
		return (*it).second;
	}
	return 0;
}
//// static method
void IPhreeqcPhastLib::CleanupIPhreeqcPhast(void)
{
	std::map<size_t, IPhreeqcPhast*>::iterator it = IPhreeqcPhast::PhastInstances.begin();
	std::vector<IPhreeqcPhast*> ipp_list;
	for ( ; it != IPhreeqcPhast::PhastInstances.end(); it++)
	{
		ipp_list.push_back(it->second);
	}
	for (size_t i = 0; i < ipp_list.size(); i++)
	{
		delete ipp_list[i];
	}
}
