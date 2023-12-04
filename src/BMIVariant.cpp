#include "BMIVariant.h"
BMIVariant::BMIVariant(VarFunction f, std::string name_in)
{
	Initialized = false;
	this->name = name_in;
	HasSetter = false;
	HasGetter = false;
	HasPtr = false;
	Nbytes = 0;
	Itemsize = 0;
	dim = 0;
	this->column = -1;
	b_var = false;
	i_var = 0;
	d_var = 0.0;
	NotImplemented = false;
	VoidPtr = NULL;
	fn = f;
#if defined(WITH_PYBIND11)
	hasPyArr = false;
#endif
}
BMIVariant::BMIVariant(const std::string& name_in, const std::string& units_in, bool setter, bool getter, bool ptr, int bytes, int item_size)
{
	Initialized = false;
	this->name = name_in;
	this->units = units_in;
	HasSetter = setter;
	HasGetter = getter;
	HasPtr = ptr;
	Nbytes = bytes;
	Itemsize = item_size;
	dim = 0;
	if (item_size > 0) dim = Nbytes / Itemsize;
	this->column = -1;
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
void BMIVariant::CopyScalars(BMIVariant& bv)
{
	this->Initialized = bv.Initialized;
	this->HasSetter = bv.HasSetter;
	this->HasGetter = bv.HasGetter;
	this->HasPtr = bv.HasPtr;
	this->Nbytes = bv.Nbytes;
	this->Itemsize = bv.Itemsize;
	this->dim = bv.dim;
	this->column = bv.column;
	this->b_var = bv.b_var;
	this->i_var = bv.i_var;
	this->d_var = bv.d_var;
	this->NotImplemented = bv.NotImplemented;
	this->VoidPtr = bv.VoidPtr;
	this->fn = bv.fn;
	this->NotImplemented = bv.NotImplemented;
	this->ctype = bv.ctype;
	this->ftype = bv.ftype;
	this->ptype = bv.ptype;
	this->clangtype = bv.clangtype;
}
void BMIVariant::Clear()
{

	Initialized = false;
	this->name.clear();
	this->type.clear();
	this->units.clear();
	HasSetter = false;
	HasGetter = false;
	HasPtr = false;
	Nbytes = 0;
	Itemsize = 0;
	dim = 0;
	ctype.clear();
	ftype.clear();
	ptype.clear();
	b_var = false;
	i_var = 0;
	d_var = 0.0;
	string_var.clear();
	IntVector.clear();
	DoubleVector.clear();
	StringVector.clear();
	NotImplemented = false;
	VoidPtr = NULL;
	fn = NULL;
	CharVector.clear();
#if defined(WITH_PYBIND11)
	hasPyArr = false;
#endif
}