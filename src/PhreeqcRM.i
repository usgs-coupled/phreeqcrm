// %module(directors="1") phreeqcrm
%module phreeqcrm
%include "typemaps.i"
%include <std_vector.i>
%{
#define SWIG_FILE_WITH_INIT
#include "BMI_Var.h"
#include "IrmResult.h"
#include "PhreeqcRM.h"
#include "bmi.hxx"
#include "BMIPhreeqcRM.h"
%}

%ignore BMIVariant;
%ignore bmi::Bmi;

#if 0
%include "numpy.i"
%init %{
import_array();
%}
#endif

%include "std_string.i"
%include "std_vector.i"
%template(DoubleVector) std::vector<double>;
%template(IntVector)    std::vector<int>;
%template(StringVector) std::vector<std::string>;

#if defined(SWIGPYTHON)
%include "swig/python/PhreeqcRM_docstrings.i"
#if defined(USE_YAML)
%include "swig/python/YAMLPhreeqcRM.i"
#endif
#endif

// Adding (%include "IPhreeqcPhast.h") forces inclusion of the
// following classes cxxSolution, cxxExchange, cxxGasPhase,
// cxxKinetics, cxxPPassemblage, cxxSSassemblage, cxxSurface
// cxxMix, cxxReaction, cxxTemperature, cxxPressure
%include "BMI_Var.h"
%include "IrmResult.h"

// Ignore methods
%ignore PhreeqcRM::GetSpeciesStoichiometry(void);
// Switch argument to output variable
%apply std::vector < double > &OUTPUT { std::vector < double > &destination_c };
%include "PhreeqcRM.h"

%extend PhreeqcRM { %pythoncode 
%{ 
import phreeqcrm
def GetSpeciesStoichiometry(self):
	species = phreeqcrm.StringVector 
	nelt_in_species = phreeqcrm.IntVector 
	elts = phreeqcrm.StringVector 
	coefs = phreeqcrm.IntVector 
	self.GetSpeciesStoichiometrySerialized(species, nelt_in_species, elts, coefs)
	all_stoich = dict()
	j_tot = 0
	for i in range(len(species)):
		n = nelt_in_species[i]
		s_stoich = dict()
		for j in range(n):
			s_stoich[elts[j_tot]] = coefs[j_tot]
		all_stoich[species[i]] = s_stoich
	return all_stoich			
%} 
}


#include "bmi.hxx"

// Ignore methods
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< int,std::allocator< int > >);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< std::string,std::allocator< std::string > >);
// Rename overloaded methods
%rename(GetValue_string_vector) GetValue(const std::string name, std::vector<std::string>& OUTPUT);
%rename(GetValue_string) GetValue(const std::string name, std::string &OUTPUT);
%rename(GetValue_double_vector) GetValue(const std::string name, std::vector<double>& OUTPUT);
%rename(GetValue_double) GetValue(const std::string name, double &OUTPUT);
%rename(GetValue_int_vector) GetValue(const std::string name, std::vector<int>& OUTPUT);
%rename(GetValue_int) GetValue(const std::string name, int &OUTPUT);
%rename(GetValue_bool) GetValue(const std::string name, bool &OUTPUT);

%rename(SetValue_string_vector) SetValue(const std::string name, std::vector<std::string> src);
%rename(SetValue_string) SetValue(const std::string name, const std::string src);
%rename(SetValue_double_vector) SetValue(const std::string name, std::vector<double> src);
%rename(SetValue_double) SetValue(const std::string name, double src);
%rename(SetValue_int_vector) SetValue(const std::string name, std::vector<int> src);
%rename(SetValue_int) SetValue(const std::string name, int src);
%rename(SetValue_bool) SetValue(const std::string name, bool src);
%include "BMIPhreeqcRM.h"

// Write new python method GetValue with one argument
%extend BMIPhreeqcRM { %pythoncode 
%{ 
def GetValue(self, var): 
	Nbytes = self.GetVarNbytes(var)
	Itemsize = self.GetVarItemsize(var)
	dim = Nbytes / Itemsize
	type = self.GetVarType(var)
	#print(f"Type={type}")
	if type=="double":
		if dim==1:
			return self.GetValue_double(var)
		if dim>1:
			v = self.GetValue_double_vector(var)
			return list(v)
	if type=="int":
		if dim==1:
			return self.GetValue_int(var)
		if dim>1:
			v = self.GetValue_int_vector(var)
			return list(v)
	if type=="std::string":
		return self.GetValue_string(var)
	if type=="std::vector<std::string>":
		v = self.GetValue_string_vector(var)
		return list(v)
	if type=="bool":
		return self.GetValue_bool(var)
def SetValue(self, var, value):
	Nbytes = self.GetVarNbytes(var)
	Itemsize = self.GetVarItemsize(var)
	dim = Nbytes / Itemsize
	type = self.GetVarType(var)
	#print(f"Type={type}")
	if type=="double":
		if dim==1:
			self.SetValue_double(var, value)
		if dim>1:
			self.SetValue_double_vector(var, value)
	if type=="int":
		if dim==1:
			self.SetValue_int(var, value)
		if dim>1:
			self.SetValue_int_vector(var)
	if type=="std::string":
		self.GetValue_string(var, value)
	if type=="std::vector<std::string>":
		self.SetValue_string_vector(var, value)	
%} 
};


