#if !defined(BMI_VAR_H_INCLUDED)
#define BMI_VAR_H_INCLUDED
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>
//#include "VarManager.h"

#if defined(_WINDLL)
#define IRM_DLL_EXPORT __declspec(dllexport)
#else
#define IRM_DLL_EXPORT
#endif
class VarManager;
class PhreeqcRM;
class IRM_DLL_EXPORT BMIVariant
{
	//typedef void (VarManager::*VarFunction)(PhreeqcRM* rm_ptr); // function pointer type
	typedef void (VarManager::* VarFunction)();
private:
	bool Initialized;
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
	std::string              string_var;
	std::vector<int>         IntVector;
	std::vector<double>      DoubleVector;
	std::vector<std::string> StringVector;
	bool							NotImplemented;
	void*							VoidPtr;
	std::vector<const char*>        CharVector;
	VarFunction fn;
public:
	// methods
	BMIVariant(VarFunction f, std::string name);
	BMIVariant()
	{
		Initialized = false;
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
	void Clear();
	void SetBasic(std::string units_in,
		bool set_in, bool get_in, bool has_ptr_in, int nbytes, int itemsize)
	{;
		//d_var = 0.0
		//this->name = name_in;
		this->units = units_in;
		this->HasSetter = set_in;
		this->HasGetter = get_in;
		this->HasPtr = has_ptr_in;
		this->Nbytes = nbytes;
		this->Itemsize = itemsize;
		if (itemsize != 0) dim = nbytes / itemsize;
		//Nbytes = 0;
		//Itemsize = 0;
		//dim = 0;
		//b_var = false;
		//i_var = 0;
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
	void SetItemsize(int v) { this->Itemsize = v; }
	void SetInitialized(bool v) { Initialized = v; }
	void SetBVar(bool v) { b_var = v; }
	void SetDVar(double v) { d_var = v; }
	void SetIVar(int v) { i_var = v; }
	void SetStringVar(std::string v) { string_var = v; }
	void SetDoubleVector(const std::vector<double> &v) {
		assert(v.size() == DoubleVector.size());
		memcpy(DoubleVector.data(), v.data(), v.size() * sizeof(double));
	}
	void SetIntVector(const std::vector<int> &v) {
		assert(v.size() == IntVector.size());
		memcpy(DoubleVector.data(), v.data(), v.size() * sizeof(int));
	}
	void SetCType(std::string v) { this->ctype = v; }
	void SetFType(std::string v) { this->ftype = v; }
	void SetPType(std::string v) { this->ptype = v; }
	void SetName(std::string s) { this->name = s; }
	void SetNotImplemented(bool v) { this->NotImplemented = v; }
	void SetVoidPtr(void* v) { this->VoidPtr = v; }
	void SetCharVector(std::vector<const char*> v) { this->CharVector = v; }
	void SetFn(VarFunction f) { this->fn = f; } // function pointer type
	void SetStringVector(const std::vector<std::string>& v) { this->StringVector = v; }

	std::string GetName() { return this->name; }
	std::string GetUnits() { return this->units; }
	bool GetHasGetter(void) { return this->HasGetter; }
	bool GetHasSetter(void) { return this->HasSetter; }
	bool GetHasPtr(void) { return this->HasPtr; }
	bool GetInitialized(void) { return this->Initialized; }
	bool* GetBVarPtr() { return &this->b_var; }
	bool GetBVar() { return this->b_var; }
	double* GetDVarPtr() { return &this->d_var; }
	double GetDVar() { return this->d_var; }
	int* GetIVarPtr() { return &this->i_var; }
	int GetIVar() { return this->i_var; }
	std::string GetStringVar() { return this->string_var; }
	int GetNbytes(void) { return this->Nbytes; }
	int GetItemsize(void) { return this->Itemsize; }
	std::string& GetCType() { return this->ctype; };;
	std::string& GetFType() { return this->ftype; };
	std::string& GetPType() { return this->ptype; };
	double* GetDoubleVectorPtr() { return this->DoubleVector.data(); }
	int* GetIntVectorPtr() { return this->IntVector.data(); }
	std::vector<double>& GetDoubleVectorRef() { return this->DoubleVector; }
	std::vector<int>& GetIntVectorRef() { return this->IntVector; }
	std::vector<std::string>& GetStringVectorRef() { return this->StringVector; }
	std::string& GetStringRef() { return this->string_var; }
	int GetDim() { return dim; }
	void* GetVoidPtr() { return this->VoidPtr; }
	std::vector<const char*> GetCharVector() { return this->CharVector; }
	VarFunction GetFn() { 
		//((VarFunction*)fn)(rm_ptr);
		return this->fn; 
	}
	bool& GetNotImplementedRef() { return this->NotImplemented; }
	void CopyScalars(BMIVariant& bv);
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