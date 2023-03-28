#if !defined(BMI_VAR_H_INCLUDED)
#define BMI_VAR_H_INCLUDED
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <assert.h>
//#include "VarManager.h"

#if defined(_WINDLL)
#define IRM_DLL_EXPORT __declspec(dllexport)
#else
#define IRM_DLL_EXPORT
#endif
class VarManager;
class PhreeqcRM;
class IRM_DLL_EXPORT Variant
{
	typedef void (VarManager::*VarFunction)(PhreeqcRM* rm_ptr); // function pointer type
private:
	std::string name;
	std::string type;
	std::string units;
	bool HasSetter;
	bool HasGetter;
	bool HasPtr;
	int Nbytes;
	int Itemsize;
	int dim;
	std::string ctype;
	std::string ftype;
	std::string ptype;
	bool							b_var;
	int								i_var;
	double							d_var;
	static std::string              string_var;
	static std::vector<int>         IntVector;
	static std::vector<double>      DoubleVector;
	static std::vector<std::string> StringVector;
	bool							NotImplemented;
	void*							VoidPtr;
	std::vector<const char*>        CharVector;
	VarFunction fn;
public:
	// methods
	Variant(VarFunction f);
	Variant()
	{
		HasSetter = false;
		HasGetter = false;
		HasPtr = false;
		Nbytes = 0;
		Itemsize = 0;
		dim = 0;
		b_var = false;
		i_var = 0;
		d_var = 0.0;
		NotImplemented = false;
		VoidPtr = NULL;
		fn = NULL;
	}
	void SetBasic(std::string name_in, std::string units_in,
		bool set_in, bool get_in, bool has_ptr_in)
	{
		this->name = name_in;
		this->units = units_in;
		this->HasSetter = set_in;
		this->HasGetter = get_in;
		this->HasPtr = has_ptr_in;
		//Nbytes = 0;
		//Itemsize = 0;
		//dim = 0;
		//b_var = false;
		//i_var = 0;
		//d_var = 0.0;
		//NotImplemented = false;
		//VoidPtr = NULL;
	};
	void SetSizes(int nbytes, int itemsize)
	{
		this->Nbytes = nbytes;
		this->Itemsize = itemsize;
		if (itemsize != 0) dim = nbytes / itemsize;
	}
	void SetTypes(std::string c, std::string f, std::string p)
	{
		this->ctype = c; this->ftype = f; this->ptype = p;
	}
	void SetBVar(bool v) { b_var = v; }
	void SetDVar(double v) { d_var = v; }
	void SetIVar(int v) { i_var = v; }
	void SetIVar(std::string v) { string_var = v; }
	void SetDoubleVector(std::vector<double> v) {
		assert(v.size() == DoubleVector.size());
		memcpy(DoubleVector.data(), v.data(), v.size() * sizeof(double));
	}
	void SetIntVector(std::vector<int> v) {
		assert(v.size() == IntVector.size());
		memcpy(DoubleVector.data(), v.data(), v.size() * sizeof(int));
	}
	void SetName(std::string s) { this->name = s; }
	void SetNotImplemented(bool v) { this->NotImplemented = v; }
	void SetVoidPtr(void* v) { this->VoidPtr = v; }
	void SetCharVector(std::vector<const char*> v) { this->CharVector = v; }
	void SetFn(VarFunction f) { this->fn = f; } // function pointer type

	std::string GetName() { return this->name; }
	std::string GetUnits() { return this->units; }
	bool GetHasGetter(void) { return this->HasGetter; }
	bool GetHasSetter(void) { return this->HasSetter; }
	bool GetHasPtr(void) { return this->HasPtr; }
	int GetNbytes(void) { return this->Nbytes; }
	int GetItemsize(void) { return this->Itemsize; }
	std::string GetCType() { return this->ctype; };;
	std::string GetFType() { return this->ftype; };
	std::string GetPType() { return this->ptype; };
	int GetDim() { return dim; }
	std::vector<const char*> GetCharVector() { return this->CharVector; }
	VarFunction GetFn() { 
		//((VarFunction*)fn)(rm_ptr);
		return this->fn; 
	}
};

class IRM_DLL_EXPORT BMI_Var
{
private:
	std::string name;
	std::string type;
	std::string units;
	bool set;
	bool get;
	bool has_ptr;
	int Nbytes;
	int Itemsize;
	int dim;
	std::string ctype;
	std::string ftype;
	std::string ptype;
public:
	// methods
	BMI_Var()
	{
		this->set = false;
		this->get = false;
		this->has_ptr = false;
		this->Nbytes = 0;
		this->Itemsize = 0;
		this->dim = 0;
	}
	BMI_Var(std::string name_in, std::string units_in,
		bool set_in, bool get_in, bool has_ptr_in)
	{
		this->name = name_in;
		this->set = set_in;
		this->get = get_in;
		this->has_ptr = has_ptr_in;
		this->Nbytes = 0;
		this->Itemsize = 0;
		this->dim = 0;
	};
	BMI_Var(std::string name_in, std::string units_in,
		bool set_in, bool get_in, bool has_ptr_in, int Nbytes_in, int Itemsize_in)
	{
		this->name = name_in;
		this->units = units_in;
		this->set = set_in;
		this->get = get_in;
		this->has_ptr = has_ptr_in;
		this->Nbytes = Nbytes_in;
		this->Itemsize = Itemsize_in;
		this->dim = 0;
		if (Itemsize != 0) dim = Nbytes / Itemsize;
	};
	std::string GetName() { return this->name; }
	void SetName(std::string s) { this->name = s; }
	std::string GetUnits() { return this->units; }
	void SetUnits(std::string units_in) { this->units = units_in; }
	bool GetSet(void) { return this->set; }
	void SetSet(bool in_in) { this->set = in_in; }
	bool GetGet(void) { return this->get; }
	void SetGet(bool out_in) { this->get = out_in; }
	bool GetHasPtr(void) { return this->has_ptr; }
	void SetHasPtr(bool tf) { this->has_ptr = tf; }
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
	int GetDim() { return dim; }

};
#endif // BMI_VAR_H_INCLUDED