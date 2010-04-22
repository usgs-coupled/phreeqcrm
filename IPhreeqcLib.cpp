#include <cassert>
#include <iostream>
#include <map>

#include "IPhreeqc.h"
#include "IPhreeqc.hpp"

class IPhreeqcLib
{
public:
	static int CreateIPhreeqc(void);
	static IPQ_RESULT DestroyIPhreeqc(int n);
	static IPhreeqc* GetInstance(int n);

private:
	static std::map<size_t, IPhreeqc*> Instances;
	static size_t InstancesIndex;
};

int
CreateIPhreeqc(void)
{
	return IPhreeqcLib::CreateIPhreeqc();
}

IPQ_RESULT
DestroyIPhreeqc(int id)
{
	return IPhreeqcLib::DestroyIPhreeqc(id);
}

int
LoadDatabase(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabase(filename);
	}
	return IPQ_BADINSTANCE;
}

int
LoadDatabaseString(int id, const char* input)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabaseString(input);
	}
	return IPQ_BADINSTANCE;
}

int
UnLoadDatabase(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->UnLoadDatabase();
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

void
OutputLastError(int id)
{
	static const char err_msg[] = "OutputLastError: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputLastError();
		return;
	}
	std::cout << err_msg << std::endl;
}

void
OutputLastWarning(int id)
{
	static const char err_msg[] = "OutputLastWarning: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputLastWarning();
		return;
	}
	std::cout << err_msg << std::endl;
}

const char*
GetLastErrorString(int id)
{
	static const char err_msg[] = "GetLastErrorString: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLastErrorString();
	}
	return err_msg;
}

const char*
GetLastWarningString(int id)
{
	static const char err_msg[] = "GetLastWarningString: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLastWarningString();
	}
	return err_msg;
}

const char*
GetDumpString(int id)
{
	static const char empty[] = "";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpString();
	}
	return empty;
}

int
GetDumpLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpLineCount();
	}
	return 0;
}

const char*
GetDumpLine(int id, int n)
{
	static const char err_msg[] = "GetDumpLine: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpLine(n);
	}
	return err_msg;
}

int
GetComponentCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->ListComponents().size();
	}
	return 0;
}

const char*
GetComponent(int id, int n)
{
	static const char err_msg[] = "GetComponent: Bad instance.\n";
	static const char empty[] = "";
	static std::string comp;

	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		std::list< std::string > comps = IPhreeqcPtr->ListComponents();
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
	return err_msg;
}

IPQ_RESULT
AccumulateLine(int id, const char *line)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->AccumulateLine(line))
		{
		case VR_OK:
			return IPQ_OK;
		case VR_OUTOFMEMORY:
			return IPQ_OUTOFMEMORY;
		default:
			assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetSelectedOutputOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetSelectedOutputOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
GetOutputOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetOutputOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetOutputOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
GetErrorOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetErrorOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetErrorOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
GetLogOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetLogOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetLogOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}


int
GetDumpOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetDumpOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetDumpOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
GetDumpStringOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetDumpStringOn())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetDumpStringOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpStringOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

int
RunAccumulated(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunAccumulated();
	}
	return IPQ_BADINSTANCE;
}

int
RunFile(int id, const char* filename)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunFile(filename);
	}
	return IPQ_BADINSTANCE;
}

int
RunString(int id, const char* input)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunString(input);
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputRowCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputRowCount();
	}
	return IPQ_BADINSTANCE;
}

int
GetSelectedOutputColumnCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputColumnCount();
	}
	return IPQ_BADINSTANCE;
}


IPQ_RESULT
GetSelectedOutputValue(int id, int row, int col, VAR* pVAR)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch(IPhreeqcPtr->GetSelectedOutputValue(row, col, pVAR))
		{
		case VR_OK:          return IPQ_OK;
		case VR_OUTOFMEMORY: return IPQ_OUTOFMEMORY;
		case VR_BADVARTYPE:  return IPQ_BADVARTYPE;
		case VR_INVALIDARG:  return IPQ_INVALIDARG;
		case VR_INVALIDROW:  return IPQ_INVALIDROW;
		case VR_INVALIDCOL:  return IPQ_INVALIDCOL;
		default:
			assert(false);
		}
	}
	return IPQ_BADINSTANCE;
}

int
AddError(int id, const char* error_msg)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->AddError(error_msg);
	}
	return IPQ_BADINSTANCE;
}

void
OutputLines(int id)
{
	static const char err_msg[] = "OutputLines: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputLines();
		return;
	}
	std::cout << err_msg << std::endl;
}

int
GetErrorLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetErrorLineCount();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetErrorLine(int id, int n)
{
	static const char err_msg[] = "GetErrorLine: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorLine(n);
	}
	return err_msg;
}

int
GetWarningLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetWarningLineCount();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetWarningLine(int id, int n)
{
	static const char err_msg[] = "GetWarningLine: Bad instance.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetWarningLine(n);
	}
	return err_msg;
}

std::map<size_t, IPhreeqc*> IPhreeqcLib::Instances;
size_t IPhreeqcLib::InstancesIndex = 0;

int
IPhreeqcLib::CreateIPhreeqc(void)
{
	int n = IPQ_OUTOFMEMORY;
	try
	{
		IPhreeqc* IPhreeqcPtr = new IPhreeqc;
		if (IPhreeqcPtr)
		{
			std::map<size_t, IPhreeqc*>::value_type instance(IPhreeqcLib::InstancesIndex, IPhreeqcPtr);
			std::pair<std::map<size_t, IPhreeqc*>::iterator, bool> pr = IPhreeqcLib::Instances.insert(instance);
			if (pr.second)
			{
				n = (int) (*pr.first).first;
				++IPhreeqcLib::InstancesIndex;
			}
		}
	}
	catch(...)
	{
		return IPQ_OUTOFMEMORY;
	}
	return n;
}

IPQ_RESULT
IPhreeqcLib::DestroyIPhreeqc(int id)
{
	IPQ_RESULT retval = IPQ_BADINSTANCE;
	if (id >= 0)
	{
		std::map<size_t, IPhreeqc*>::iterator it = IPhreeqcLib::Instances.find(size_t(id));
		if (it != IPhreeqcLib::Instances.end())
		{
			delete (*it).second;
			IPhreeqcLib::Instances.erase(it);
			retval = IPQ_OK;
		}
	}
	return retval;
}

IPhreeqc*
IPhreeqcLib::GetInstance(int id)
{
	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqcLib::Instances.find(size_t(id));
	if (it != IPhreeqcLib::Instances.end())
	{
		return (*it).second;
	}
	return 0;
}
