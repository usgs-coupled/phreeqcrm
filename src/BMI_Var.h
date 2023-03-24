#if !defined(BMI_VAR_H_INCLUDED)
#define BMI_VAR_H_INCLUDED
#include <iostream>
#include <map>
#include <string>

#if defined(_WINDLL)
#define IRM_DLL_EXPORT __declspec(dllexport)
#else
#define IRM_DLL_EXPORT
#endif

class IRM_DLL_EXPORT BMI_Var
{
private:
	std::string name;
	std::string type;
	std::string units;
	bool set;
	bool get;
	int Nbytes;
	int Itemsize;
	std::string ctype;
	std::string ftype;
	std::string ptype;
public:
	// methods
	BMI_Var()
	{
		this->set = false;
		this->get = false;
		this->Nbytes = 0;
		this->Itemsize = 0;
	}
	BMI_Var(std::string name_in, std::string units_in,
		bool set_in, bool get_in)
	{
		this->name = name_in;
		this->set = set_in;
		this->get = get_in;
		this->Nbytes = 0;
		this->Itemsize = 0;
	};
	BMI_Var(std::string name_in, std::string units_in,
		bool set_in, bool get_in, int Nbytes_in, int Itemsize_in)
	{
		this->name = name_in;
		this->units = units_in;
		this->set = set_in;
		this->get = get_in;
		this->Nbytes = Nbytes_in;
		this->Itemsize = Itemsize_in;
	};
	std::string GetName() { return this->name; }
	void SetName(std::string s) { this->name = s; }
	std::string GetUnits() { return this->units; }
	void SetUnits(std::string units_in) { this->units = units_in; }
	bool GetSet(void) { return this->set; }
	void SetSet(bool in_in) { this->set = in_in; }
	bool GetGet(void) { return this->get; }
	void SetGet(bool out_in) { this->get = out_in; }
	int GetNbytes(void) { return this->Nbytes; }
	void SetNbytes(int n) { this->Nbytes = n; }
	int GetItemsize(void) { return this->Itemsize; }
	void SetItemsize(int n) { this->Itemsize = n; }
	std::string GetCType() { return this->ctype; };
	void SetCType(std::string s) { this->ctype = s; };
	std::string GetFType() { return this->ftype; };
	void SetFType(std::string s) { this->ftype = s; };
	std::string GetPType() { return this->ptype; };
	void SetPType(std::string s) { this->ptype = s; };
	void SetTypes(std::string c, std::string f, std::string p)
	{
		this->ctype = c; this->ftype = f; this->ptype = p;
	}

};
#endif // BMI_VAR_H_INCLUDED