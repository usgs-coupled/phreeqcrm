%module phreeqcrm
%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(vectori) vector<int>;
  %template(vectord) vector<double>;
  %template(vectors) vector<string>;
};

%include "../src/IrmResult.h"

%{
#include "../src/PhreeqcRM.h"
%}

// need to modify
%ignore PhreeqcRM::GetSelectedOutputHeading(int icol, std::string& heading);  // no passing by ref???

%ignore PhreeqcRM::ReturnHandler(IRM_RESULT result, const std::string & e_string);

%include "../src/PhreeqcRM.h"
