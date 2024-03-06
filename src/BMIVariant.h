#if !defined(BMI_VAR_H_INCLUDED)
#define BMI_VAR_H_INCLUDED
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>

#if defined(WITH_PYBIND11)
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

#include "irm_dll_export.h"

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
	int column;
	std::string ctype;
	std::string clangtype;
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
#if defined(WITH_PYBIND11)
	bool hasPyArr;
	py::array pyArr;
#endif
public:
	// methods
	BMIVariant(VarFunction f, std::string name);
	BMIVariant(const std::string& name_in, const std::string& units, bool setter, bool getter, bool ptr, int bytes, int item_size);
	BMIVariant()
	{
		Initialized = false;
		HasSetter = false;
		HasGetter = false;
		HasPtr = false;
		Nbytes = 0;
		Itemsize = 0;
		dim = 0;
		column = -1;
		b_var = false;
		i_var = 0;
		d_var = 0.0;
		NotImplemented = false;
		VoidPtr = NULL;
		fn = NULL;
#if defined(WITH_PYBIND11)
		hasPyArr = false;
#endif
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
	void SetTypes(std::string c, std::string f, std::string p, std::string clang)
	{
		this->ctype = c; this->ftype = f; this->ptype = p; this->clangtype = clang;
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
	void SetClangType(std::string v) { this->clangtype = v; }
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
	std::string& GetCType() { return this->ctype; }
	std::string& GetClangType() { return this->clangtype; }
	std::string& GetFType() { return this->ftype; };
	std::string& GetPType() { return this->ptype; };
	double* GetDoubleVectorPtr() { return this->DoubleVector.data(); }
	int* GetIntVectorPtr() { return this->IntVector.data(); }
	std::vector<double>& GetDoubleVectorRef() { return this->DoubleVector; }
	std::vector<int>& GetIntVectorRef() { return this->IntVector; }
	std::vector<std::string>& GetStringVectorRef() { return this->StringVector; }
	std::string& GetStringRef() { return this->string_var; }
	int GetDim() { return dim; }
	int GetColumn() { return column; }
	void SetColumn(int v) { column = v; }
	void* GetVoidPtr() { return this->VoidPtr; }
	std::vector<const char*> GetCharVector() { return this->CharVector; }
	VarFunction GetFn() { 
		//((VarFunction*)fn)(rm_ptr);
		return this->fn; 
	}
	bool& GetNotImplementedRef() { return this->NotImplemented; }
	void CopyScalars(BMIVariant& bv);

#if defined(WITH_PYBIND11)
	bool HasPyArray()
	{
		return hasPyArr;
	}
	py::array GetPyArray()
	{
		assert(hasPyArr);
		return pyArr;
	}
	py::array SetPyArray(py::array arr)
	{
		assert(!hasPyArr);
		pyArr = arr;
		hasPyArr = true;
		return pyArr;
	}
#endif
};
#endif // BMI_VAR_H_INCLUDED