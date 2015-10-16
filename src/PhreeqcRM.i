#ifdef SWIGJAVA
%module(directors="1") phreeqcrm_java
#elif SWIGRUBY
%module(directors="1") phreeqcrm_ruby
#else
%module(directors="1") phreeqcrm
#endif

%{
#include "IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp"
#include "IPhreeqcPhast/IPhreeqc/Var.h"
#include "IPhreeqcPhast/IPhreeqc/BasicCallback.h"
#include "IPhreeqcPhast/IPhreeqcPhast.h"
#include "IrmResult.h"
#include "PhreeqcRM.h"
%}

%feature("director") BasicCallback;

%include "stl.i"
%include "std_except.i"
%include "std_string.i"

/* ignore */
%ignore SetBasicFortranCallback;

/* instantiate the required template specializations */
namespace std {
    %template(IntVector)           vector<int>;
    %template(DoubleVector)        vector<double>;
    %template(BoolVector)          vector<bool>;
    %template(StringVector)        vector<string>;
    %template(IPhreeqcPhastVector) vector<IPhreeqcPhast*>;
    %template(IntVectorVector)     vector< vector<int> >;
}

%include "std_string.i"

%include "IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp"
%include "IPhreeqcPhast/IPhreeqc/Var.h"
%include "IPhreeqcPhast/IPhreeqc/BasicCallback.h"
%include "IPhreeqcPhast/IPhreeqcPhast.h"
%include "IrmResult.h"
%include "PhreeqcRM.h"
