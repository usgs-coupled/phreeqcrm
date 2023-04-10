// %module(directors="1") phreeqcrm
%module phreeqcrm

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
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< int,std::allocator< int > >);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< std::string,std::allocator< std::string > >);

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
%include "PhreeqcRM.h"
#include "bmi.hxx"
%include "BMIPhreeqcRM.h"
