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
AddError(int id, const char* error_msg)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->AddError(error_msg);
	}
	return IPQ_BADINSTANCE;
}

// TODO AddWarning

int
CreateIPhreeqc(void)
{
	return IPhreeqcLib::CreateIPhreeqc();
}

IPQ_RESULT
ClearAccumulatedLines(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->ClearAccumulatedLines();
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
DestroyIPhreeqc(int id)
{
	return IPhreeqcLib::DestroyIPhreeqc(id);
}

// TODO Maybe GetAccumulatedLines

const char*
GetComponent(int id, int n)
{
	static const char err_msg[] = "GetComponent: Invalid instance id.\n";
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

int
GetComponentCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->ListComponents().size();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetDumpStringLine(int id, int n)
{
	static const char err_msg[] = "GetDumpStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpStringLine(n);
	}
	return err_msg;
}

int
GetDumpStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpStringLineCount();
	}
	return 0;
}

int
GetDumpFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetDumpFileOn())
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

const char*
GetErrorStringLine(int id, int n)
{
	static const char err_msg[] = "GetErrorStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorStringLine(n);
	}
	return err_msg;
}

int
GetErrorStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetErrorStringLineCount();
	}
	return IPQ_BADINSTANCE;
}

int
GetErrorFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetErrorFileOn())
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

const char*
GetErrorString(int id)
{
	static const char err_msg[] = "GetErrorString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetErrorString();
	}
	return err_msg;
}

int
GetLogFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetLogFileOn())
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

int
GetOutputFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetOutputFileOn())
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

int
GetSelectedOutputFileOn(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		if (IPhreeqcPtr->GetSelectedOutputFileOn())
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

const char*
GetWarningStringLine(int id, int n)
{
	static const char err_msg[] = "GetWarningStringLine: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetWarningStringLine(n);
	}
	return err_msg;
}

int
GetWarningStringLineCount(int id)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->GetWarningStringLineCount();
	}
	return IPQ_BADINSTANCE;
}

const char*
GetWarningString(int id)
{
	static const char err_msg[] = "GetWarningString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetWarningString();
	}
	return err_msg;
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

void
OutputErrorString(int id)
{
	static const char err_msg[] = "OutputErrorString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputErrorString();
		return;
	}
	std::cout << err_msg << std::endl;
}

void
OutputAccumulatedLines(int id)
{
	static const char err_msg[] = "OutputAccumulatedLines: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputAccumulatedLines();
		return;
	}
	std::cout << err_msg << std::endl;
}

void
OutputWarningString(int id)
{
	static const char err_msg[] = "OutputWarningString: Invalid instance id.\n";
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->OutputWarningString();
		return;
	}
	std::cout << err_msg << std::endl;
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

IPQ_RESULT
SetDumpFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpFileOn(value != 0);
		return IPQ_OK;
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

IPQ_RESULT
SetErrorFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetLogFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetOutputFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
SetSelectedOutputFileOn(int id, int value)
{
	IPhreeqc* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputFileOn(value != 0);
		return IPQ_OK;
	}
	return IPQ_BADINSTANCE;
}

IPQ_RESULT
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

// helper functions
//

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
