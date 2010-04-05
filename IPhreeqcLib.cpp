#include "IPhreeqcLib.h"
#include "IPhreeqc2.h"
#include <cassert>


int
CreateIPhreeqc(void)
{
	return IPhreeqcLib::CreateIPhreeqc();
}

IPL_RESULT
DestroyIPhreeqc(int id)
{
	return IPhreeqcLib::DestroyIPhreeqc(id);
}

int
LoadDatabaseM(int id, const char* filename)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabase(filename);
	}
	return IPL_BADINSTANCE;
}

int
LoadDatabaseStringM(int id, const char* input)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->LoadDatabaseString(input);
	}
	return IPL_BADINSTANCE;
}

int
UnLoadDatabaseM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->UnLoadDatabase();
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

const char*
GetLastErrorStringM(int id)
{
	static const char err_msg[] = "GetLastErrorString: Bad instance.\n";
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetLastErrorString();
	}
	return err_msg;
}

const char*
GetDumpStringM(int id)
{
	static const char empty[] = "";
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpString();
	}
	return empty;
}

int
GetDumpLineCountM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpLineCount();
	}
	return 0;
}

const char* 
GetDumpLineM(int id, int n)
{
	static const char err_msg[] = "GetDumpLine: Bad instance.\n";
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetDumpLine(n);
	}
	return err_msg;
}

int
GetComponentCountM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->ListComponents().size();
	}
	return 0;
}

const char*
GetComponentM(int id, int n)
{
	static const char err_msg[] = "GetComponent: Bad instance.\n";
	static const char empty[] = "";
	static std::string comp;

	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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

IPL_RESULT
AccumulateLineM(int id, const char *line)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch (IPhreeqcPtr->AccumulateLine(line))
		{
		case VR_OK:
			return IPL_OK;
		case VR_OUTOFMEMORY:
			return IPL_OUTOFMEMORY;
		default:
			assert(false);
		}
	}
	return IPL_BADINSTANCE;
}

int
GetSelectedOutputOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetSelectedOutputOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetSelectedOutputOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

int
GetOutputOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetOutputOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetOutputOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

int
GetErrorOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetErrorOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetErrorOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

int
GetLogOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetLogOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetLogOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}


int
GetDumpOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetDumpOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

int
GetDumpStringOnM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
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
	return IPL_BADINSTANCE;
}

IPL_RESULT
SetDumpStringOnM(int id, int value)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		IPhreeqcPtr->SetDumpStringOn(value != 0);
		return IPL_OK;
	}
	return IPL_BADINSTANCE;
}

int
RunAccumulatedM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunAccumulated();
	}
	return IPL_BADINSTANCE;
}

int
RunFileM(int id, const char* filename)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunFile(filename);
	}
	return IPL_BADINSTANCE;
}

int
RunStringM(int id, const char* input)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->RunString(input);
	}
	return IPL_BADINSTANCE;
}

int
GetSelectedOutputRowCountM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputRowCount();
	}
	return IPL_BADINSTANCE;
}

int
GetSelectedOutputColumnCountM(int id)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return IPhreeqcPtr->GetSelectedOutputColumnCount();
	}
	return IPL_BADINSTANCE;
}


IPL_RESULT
GetSelectedOutputValueM(int id, int row, int col, VAR* pVAR)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		switch(IPhreeqcPtr->GetSelectedOutputValue(row, col, pVAR))
		{
		case VR_OK:          return IPL_OK;
		case VR_OUTOFMEMORY: return IPL_OUTOFMEMORY;
		case VR_BADVARTYPE:  return IPL_BADVARTYPE;
		case VR_INVALIDARG:  return IPL_INVALIDARG;
		case VR_INVALIDROW:  return IPL_INVALIDROW;
		case VR_INVALIDCOL:  return IPL_INVALIDCOL;
		default:
			assert(false);
		}
	}
	return IPL_BADINSTANCE;
}

int
AddErrorM(int id, const char* error_msg)
{
	IPhreeqc2* IPhreeqcPtr = IPhreeqcLib::GetInstance(id);
	if (IPhreeqcPtr)
	{
		return (int)IPhreeqcPtr->AddError(error_msg);
	}
	return IPL_BADINSTANCE;
}

std::map<size_t, IPhreeqc2*> IPhreeqcLib::Instances;
size_t IPhreeqcLib::InstancesIndex = 0;

int
IPhreeqcLib::CreateIPhreeqc(void)
{
	int n = IPL_OUTOFMEMORY;
	try
	{
		IPhreeqc2* IPhreeqcPtr = new IPhreeqc2;
		if (IPhreeqcPtr)
		{
			std::map<size_t, IPhreeqc2*>::value_type instance(IPhreeqcLib::InstancesIndex, IPhreeqcPtr);
			std::pair<std::map<size_t, IPhreeqc2*>::iterator, bool> pr = IPhreeqcLib::Instances.insert(instance);
			if (pr.second)
			{
				n = (int) (*pr.first).first;
				++IPhreeqcLib::InstancesIndex;
			}
		}
	}
	catch(...)
	{
		return IPL_OUTOFMEMORY;
	}
	return n;
}

IPL_RESULT
IPhreeqcLib::DestroyIPhreeqc(int id)
{
	IPL_RESULT retval = IPL_BADINSTANCE;
	if (id >= 0)
	{
		std::map<size_t, IPhreeqc2*>::iterator it = IPhreeqcLib::Instances.find(size_t(id));
		if (it != IPhreeqcLib::Instances.end())
		{
			delete (*it).second;
			IPhreeqcLib::Instances.erase(it);
			retval = IPL_OK;
		}
	}
	return retval;
}

IPhreeqc2*
IPhreeqcLib::GetInstance(int id)
{
	std::map<size_t, IPhreeqc2*>::iterator it = IPhreeqcLib::Instances.find(size_t(id));
	if (it != IPhreeqcLib::Instances.end())
	{
		return (*it).second;
	}
	return 0;
}

