%module(directors="1") phreeqcrm

%{
#define SWIG_FILE_WITH_INIT
#include "IrmResult.h"
#include "PhreeqcRM.h"
%}

%ignore GetGridCellCountYAML(std::string YAML_file);  // TODO

%include "numpy.i"
%init %{
import_array();
%}

%include "std_string.i" // req'd for PhreeqcRM::LoadDatabase
%include "std_vector.i"
%template(DoubleVector) std::vector<double>;
%template(IntVector) std::vector<int>;
%template(StringVector) std::vector<std::string>;

// Adding (%include "IPhreeqcPhast.h") forces inclusion of the
// following classes cxxSolution, cxxExchange, cxxGasPhase,
// cxxKinetics, cxxPPassemblage, cxxSSassemblage, cxxSurface
// cxxMix, cxxReaction, cxxTemperature, cxxPressure
%include "IrmResult.h"
%include "PhreeqcRM.h"
