// %module(directors="1") phreeqcrm
%module phreeqcrm

%{
#define SWIG_FILE_WITH_INIT
#include "BMI_Var.h"
#include "IrmResult.h"
#include "PhreeqcRM.h"
%}

%ignore GetGridCellCountYAML(std::string YAML_file);  // TODO

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

%define SetPorosity_DOCSTRING
"Set the porosity for each reaction cell.

The volume of water in a reaction cell is the product of porosity, saturation
:meth:`SetSaturation`, and representative volume :meth:`SetRepresentativeVolume`.

Args:
    por (list, numpy.ndarray or tuple): A list, NumPy array, or tuple of floats.

Returns:
    IRM_RESULT:	0 is success, negative is failure (See :meth:`DecodeError`)."
%enddef

%feature("docstring") PhreeqcRM::SetPorosity SetPorosity_DOCSTRING

// Adding (%include "IPhreeqcPhast.h") forces inclusion of the
// following classes cxxSolution, cxxExchange, cxxGasPhase,
// cxxKinetics, cxxPPassemblage, cxxSSassemblage, cxxSurface
// cxxMix, cxxReaction, cxxTemperature, cxxPressure
%include "BMI_Var.h"
%include "IrmResult.h"
%include "PhreeqcRM.h"
