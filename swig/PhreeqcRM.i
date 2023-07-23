// %module(directors="1") phreeqcrm
%module phreeqcrm

%begin %{
#ifdef _MSC_VER
#define SWIG_PYTHON_INTERPRETER_NO_DEBUG
#endif
%}

%header %{
#if defined(SWIGPYTHON)
// This is reqd to use pyfragments.swg from numpy
// it fixes these type of errors:
// PhreeqcRMPYTHON_wrap.cxx: error C2065: 'Integer': undeclared identifier
// PhreeqcRMPYTHON_wrap.cxx: error C3861: 'PyArray_IsScalar': identifier not found
// see https://numpy.org/doc/stable/reference/swig.interface-file.html#numpy-array-scalars-and-swig
#include <numpy/arrayobject.h>
#endif
%}

%pythoncode 
%{ 
import numpy as np
import phreeqcrm
from enum import Enum, unique

@unique
class State(Enum):
    UNINITIALIZED = 1
    INITIALIZED = 2
%}
%include "typemaps.i"
%include <std_vector.i>
%{
#define SWIG_FILE_WITH_INIT
#include "BMIVariant.h"
#include "IrmResult.h"
#include "PhreeqcRM.h"
#include "bmi.hxx"
#include "BMIPhreeqcRM.h"
#if defined(USE_YAML)
#include "yaml-cpp/yaml.h"
#endif
%}
%ignore BMIVariant;
%ignore bmi::Bmi;

#if defined(SWIGPYTHON)
%include "numpy.i"
%init %{
import_array();
%}
%fragment("NumPy_Fragments");
#endif

%include "std_string.i"
%include "std_vector.i"
%template(BoolVector)   std::vector<bool>;
%template(DoubleVector) std::vector<double>;
%template(IntVector)    std::vector<int>;
%template(StringVector) std::vector<std::string>;

#if defined(SWIGPYTHON)
%include "python/BMIPhreeqcRM_docstrings.swg"
%include "python/PhreeqcRM_docstrings.swg"
#if defined(USE_YAML)
%include "python/YAMLPhreeqcRM.swg"
#endif
#endif

// Adding (%include "IPhreeqcPhast.h") forces inclusion of the
// following classes cxxSolution, cxxExchange, cxxGasPhase,
// cxxKinetics, cxxPPassemblage, cxxSSassemblage, cxxSurface
// cxxMix, cxxReaction, cxxTemperature, cxxPressure
%include "../src/BMIVariant.h"
%include "../src/IrmResult.h"

// Ignore methods
%ignore PhreeqcRM::GetIPhreeqcPointer(int i);
%ignore PhreeqcRM::GetSelectedOutputHeading(int icol, std::string &heading);
%ignore PhreeqcRM::GetWorkers();
%ignore PhreeqcRM::MpiAbort();
%ignore PhreeqcRM::SetMpiWorkerCallbackC(int (*fcn)(int *method, void * cookie));
%ignore PhreeqcRM::SetMpiWorkerCallbackCookie(void * cookie);
// Switch argument to output variable
%apply std::vector < double >     &OUTPUT { std::vector < double > &destination_c };
%apply std::vector<std::string>   &OUTPUT { std::vector<std::string> &species_output, std::vector<std::string> &elts_output };
%apply std::vector<int>           &OUTPUT { std::vector<int> &nelt_output }; 
%apply std::vector<double>        &OUTPUT { std::vector<double> &coef_output };
%apply std::vector<int>           &OUTPUT { std::vector<int>& nback_output, std::vector<int> &cellnumbers_output};

%apply std::vector< double >      &OUTPUT { std::vector< double > &c_output };
%apply std::vector< double >      &OUTPUT { std::vector< double > &d_output };
%apply std::vector< double >      &OUTPUT { std::vector< double > &gas_moles_output };
%apply std::vector< double >      &OUTPUT { std::vector< double > &gas_phi };
%apply std::vector< double >      &OUTPUT { std::vector< double > &gas_pressure };
%apply std::vector< double >      &OUTPUT { std::vector< double > &gas_volume_output };
%apply std::vector< double >      &OUTPUT { std::vector< double > &species_conc_output };
%apply std::vector< double >      &OUTPUT { std::vector< double > &species_log10gammas };
%apply std::vector< double >      &OUTPUT { std::vector< double > &species_log10molalities };
%apply std::vector< double >      &OUTPUT { std::vector< double > &s_output };
%apply std::vector< std::string > &OUTPUT { std::vector< std::string > &headings };
%apply std::vector< double >      &OUTPUT { std::vector< double > &sat_output };

// Rename method to avoid tuple

%rename(InitialPhreeqc2ConcentrationsSWIG)          InitialPhreeqc2Concentrations(
													std::vector < double > &destination_c, 
													const std::vector < int > &boundary_solution1);

%rename(InitialPhreeqc2ConcentrationsSWIG_mix)      InitialPhreeqc2Concentrations(
													std::vector < double > & destination_c,
													const std::vector < int >    & boundary_solution1,
													const std::vector < int >    & boundary_solution2,
													const std::vector < double > & fraction1);
%rename(InitialPhreeqc2ModuleSWIG)                  InitialPhreeqc2Module(
													const std::vector < int > & initial_conditions1);
%rename(InitialPhreeqc2ModuleSWIG_mix)              InitialPhreeqc2Module(
													const std::vector < int > & initial_conditions1,
													const std::vector < int > & initial_conditions2,	
													const std::vector < double > & fraction1);
%rename(InitialPhreeqc2SpeciesConcentrationsSWIG)   InitialPhreeqc2SpeciesConcentrations(
													std::vector < double > & destination_c,
													const std::vector < int >    & boundary_solution1);
%rename(InitialPhreeqc2SpeciesConcentrationsSWIG_mix) InitialPhreeqc2SpeciesConcentrations(
													std::vector < double > & destination_c,
													const std::vector < int >    & boundary_solution1,
													const std::vector < int >    & boundary_solution2,
													const std::vector < double > & fraction1);
%rename(InitialPhreeqcCell2ModuleSWIG)              InitialPhreeqcCell2Module(int n, 
		                                            const std::vector< int > &cell_numbers);
%rename(CreateMappingSWIG)                          CreateMapping(const std::vector< int > &grid2chem);
%rename(GetConcentrationsSWIG)                      GetConcentrations(std::vector< double > &c_output);
%rename(GetDensityCalculatedSWIG)                   GetDensityCalculated(std::vector< double > & d_output);
%rename(GetEndCellSWIG)                             GetEndCell();
%rename(GetForwardMappingSWIG)                      GetForwardMapping();
%rename(GetGasCompMolesSWIG)                        GetGasCompMoles(std::vector< double >& gas_moles);
%rename(GetGasCompPhiSWIG)                          GetGasCompPhi(std::vector< double >& gas_phi);
%rename(GetGasCompPressuresSWIG)                    GetGasCompPressures(std::vector< double >& gas_pressure);
%rename(GetGasPhaseVolumeSWIG)                      GetGasPhaseVolume(std::vector< double >& gas_volume);
%rename(GetGfwSWIG)                                 GetGfw();
%rename(GetIthConcentrationSWIG)                    GetIthConcentration(int i, std::vector< double >& c_output);  
%rename(GetIthSpeciesConcentrationSWIG)             GetIthSpeciesConcentration(int i, std::vector< double >& c_output);
%rename(GetPorositySWIG)                            GetPorosity();
%rename(GetPressureSWIG)                            GetPressure();
%rename(GetPrintChemistryMaskSWIG)                  GetPrintChemistryMask();
%rename(GetPrintChemistryOn_bool_vector)            GetPrintChemistryOn();
%rename(GetSaturationCalculatedSWIG)                GetSaturationCalculated(std::vector< double > & sat_output);
%rename(GetSelectedOutputSWIG)                      GetSelectedOutput(std::vector< double > &s_output);
%rename(GetSelectedOutputHeadingsSWIG)              GetSelectedOutputHeadings(std::vector< std::string >& headings);
%rename(GetSpeciesConcentrationsSWIG)               GetSpeciesConcentrations(std::vector< double > &species_conc_output);
%rename(GetSpeciesD25SWIG)                          GetSpeciesD25();
%rename(GetSpeciesLog10GammasSWIG)                  GetSpeciesLog10Gammas(std::vector< double > & species_log10gammas);
%rename(GetSpeciesLog10MolalitiesSWIG)              GetSpeciesLog10Molalities(std::vector< double >& species_log10molalities);
%rename(GetSpeciesZSWIG)                            GetSpeciesZ();
%rename(GetSpeciesConcentrationsSWIG)               GetSpeciesConcentrations(std::vector< double > &species_conc_output);
%rename(GetSolutionVolumeSWIG)                      GetSolutionVolume();
%rename(GetStartCellSWIG)                           GetStartCell();
%rename(GetTemperatureSWIG)                         GetTemperature();
%rename(GetViscositySWIG)                           GetViscosity();
%rename(InitialEquilibriumPhases2ModuleSWIG)        InitialEquilibriumPhases2Module(const std::vector < int >& equilibrium_phases);
%rename(InitialExchanges2ModuleSWIG)                InitialExchanges2Module(const std::vector < int >& exchangers);
%rename(InitialGasPhases2ModuleSWIG)                InitialGasPhases2Module(const std::vector < int >& gas_phases);
%rename(InitialKinetics2ModuleSWIG)                 InitialKinetics2Module(const std::vector < int >& kinetics);
%rename(InitialSolutions2ModuleSWIG)                InitialSolutions2Module(const std::vector < int >& solutions);
%rename(InitialSolidSolutions2ModuleSWIG)           InitialSolidSolutions2Module(const std::vector < int >& solid_solutions);
%rename(InitialSurfaces2ModuleSWIG)                 InitialSurfaces2Module(const std::vector < int >& surfaces);
%rename(SetPrintChemistryMaskSWIG)                  SetPrintChemistryMask(const std::vector<int> & cell_mask);
%rename(SetIthConcentrationSWIG)                    SetIthConcentration(int i, std::vector< double >& c); 
%rename(SetIthSpeciesConcentrationSWIG)             SetIthSpeciesConcentration(int i, std::vector< double >& c);

// Ignore methods
%ignore BMIPhreeqcRM::BMIPhreeqcRM(int,int);  // @todo
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
///%ignore BMIPhreeqcRM::GetValuePtr(std::string);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< int,std::allocator< int > >);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< std::string,std::allocator< std::string > >);

%ignore PhreeqcRM::Initializer;

%include "../src/PhreeqcRM.h"

%extend PhreeqcRM { %pythoncode 
%{ 
def GetSpeciesStoichiometry(self):
	v = self.GetSpeciesStoichiometrySWIG()
	species = v[0]
	nelt_in_species = v[1]
	elts = v[2]
	coefs = v[3]
	all_stoich = dict()
	j_tot = 0
	for i in range(len(species)):
		n = nelt_in_species[i]
		s_stoich = dict()
		for j in range(n):
			s_stoich[elts[j_tot]] = coefs[j_tot]
			j_tot += 1
		all_stoich[species[i]] = s_stoich
	return all_stoich	
def GetBackwardMapping(self):
	v = self.GetBackwardMappingSWIG()
	count = v[0]
	allcells = v[1]
	backward_mapping = dict()
	j_tot = 0
	for i in range(len(count)):
		n = count[i]
		back = []
		for j in range(n):
			back.append(allcells[j_tot]);
			j_tot += 1
		backward_mapping[i] = back
	return backward_mapping	
def InitialPhreeqc2Concentrations(self, bc1):
	if not isinstance(bc1, phreeqcrm.IntVector):
		bc1 = self.GetIntVector(bc1)
	return np.array(self.InitialPhreeqc2ConcentrationsSWIG(bc1)[1])
def InitialPhreeqc2Module(self, ic1):
	if not isinstance(ic1, phreeqcrm.IntVector):
		ic1 = self.GetIntVector(ic1)
	return self.InitialPhreeqc2ModuleSWIG(ic1)
def InitialPhreeqc2Module_mix(self, ic1, ic2, f1):
	if not isinstance(ic1, phreeqcrm.IntVector):
		ic1 = self.GetIntVector(ic1)
	if not isinstance(ic2, phreeqcrm.IntVector):
		ic2 = self.GetIntVector(ic2)
	return self.InitialPhreeqc2ModuleSWIG_mix(ic1,ic2,f1)
def InitialPhreeqc2Concentrations_mix(self, bc1, bc2, f1):
	if not isinstance(bc1, phreeqcrm.IntVector):
		bc1 = self.GetIntVector(bc1)
	if not isinstance(bc2, phreeqcrm.IntVector):
		bc2 = self.GetIntVector(bc2)
	return np.array(self.InitialPhreeqc2ConcentrationsSWIG_mix(bc1,bc2,f1)[1])
def InitialPhreeqc2SpeciesConcentrations(self, bc1):
	if not isinstance(bc1, phreeqcrm.IntVector):
		bc1 = self.GetIntVector(bc1)
	return np.array(self.InitialPhreeqc2SpeciesConcentrationsSWIG(bc1)[1])
def InitialPhreeqc2SpeciesConcentrations_mix(self, bc1, bc2, f1):
	if not isinstance(bc1, phreeqcrm.IntVector):
		bc1 = self.GetIntVector(bc1)
	if not isinstance(bc2, phreeqcrm.IntVector):
		bc2 = self.GetIntVector(bc2)
	return np.array(self.InitialPhreeqc2SpeciesConcentrationsSWIG_mix(bc1,bc2,f1)[1])
def InitialPhreeqcCell2Module(self, n, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialPhreeqcCell2ModuleSWIG(n, v)
def CreateMapping(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.CreateMappingSWIG(v)
def GetConcentrations(self):
	return np.array(self.GetConcentrationsSWIG()[1])
def GetDensityCalculated(self): 
	return np.array(self.GetDensityCalculatedSWIG()[1])     
def GetEndCell(self):
	return np.array(self.GetEndCellSWIG())    
def GetForwardMapping(self):
	return np.array(self.GetForwardMappingSWIG())  
def GetGasCompMoles(self):                 
	return np.array(self.GetGasCompMolesSWIG()[1])             
def GetGasCompPhi(self):                
	return np.array(self.GetGasCompPhiSWIG()[1])  
def GetGasCompPressures(self):      
	return np.array(self.GetGasCompPressuresSWIG()[1])           
def GetGasPhaseVolume(self):           
	return np.array(self.GetGasPhaseVolumeSWIG()[1])    
def GetGfw(self):
	return np.array(self.GetGfwSWIG())  
def GetIthConcentration(self, n):
	return np.array(self.GetIthConcentrationSWIG(n)[1])   
def GetIthSpeciesConcentration(self, n):
	return np.array(self.GetIthSpeciesConcentrationSWIG(n)[1])  
def GetPrintChemistryOn(self):
	return np.array(self.GetPrintChemistryOn_bool_vector())  
def GetPorosity(self):
	return np.array(self.GetPorositySWIG())
def GetPressure(self):
	return np.array(self.GetPressureSWIG())
def GetPrintChemistryMask(self):
	return np.array(self.GetPrintChemistryMaskSWIG())
def GetPrintChemistryOn(self):
	return np.array(self.GetPrintChemistryOn_bool_vector())
def GetSaturationCalculated(self):      
	return np.array(self.GetSaturationCalculatedSWIG()[1])   
def GetSelectedOutput(self):    
	return np.array(self.GetSelectedOutputSWIG()[1])     
def GetSelectedOutputHeadings(self):             
	return np.array(self.GetSelectedOutputHeadingsSWIG()[1])  
def GetSpeciesConcentrations(self):
	return np.array(self.GetSpeciesConcentrationsSWIG()[1])  
def GetSpeciesD25(self):
	return np.array(self.GetSpeciesD25SWIG())    
def GetSpeciesLog10Gammas(self):   
	return np.array(self.GetSpeciesLog10GammasSWIG()[1])     
def GetSpeciesLog10Molalities(self):        
	return np.array(self.GetSpeciesLog10MolalitiesSWIG()[1]) 
def GetSpeciesZ(self):
	return np.array(self.GetSpeciesZSWIG())
def GetSolutionVolume(self):
	return np.array(self.GetSolutionVolumeSWIG())
def GetStartCell(self):
	return np.array(self.GetStartCellSWIG())
def GetTemperature(self):
	return np.array(self.GetTemperatureSWIG())
def GetViscosity(self):
	return np.array(self.GetViscositySWIG())
def InitialEquilibriumPhases2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialEquilibriumPhases2ModuleSWIG(v)
def InitialExchanges2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialExchanges2ModuleSWIG(v)
def InitialGasPhases2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialGasPhases2ModuleSWIG(v)
def InitialKinetics2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialKinetics2ModuleSWIG(v)
def InitialSolutions2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialSolutions2ModuleSWIG(v)
def InitialSolidSolutions2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialSolidSolutions2ModuleSWIG(v)
def InitialSurfaces2Module(self, v):
	if not isinstance(v, phreeqcrm.IntVector):
		v = self.GetIntVector(v)
	return self.InitialSurfaces2ModuleSWIG(v)
def SetIthConcentration(self, i, c):
	if not isinstance(c, phreeqcrm.DoubleVector):
		c = self.GetDoubleVector(c)
	return self.SetIthConcentrationSWIG(i, c)
def SetIthSpeciesConcentration(self, i, c):
	if not isinstance(c, phreeqcrm.DoubleVector):
		c = self.GetDoubleVector(c)
	return self.SetIthSpeciesConcentrationSWIG(i, c)
def SetPrintChemistryMask(self, cell_mask):
	if not isinstance(cell_mask, phreeqcrm.IntVector):
		cell_mask = self.GetIntVector(cell_mask)
	return self.SetPrintChemistryMaskSWIG(cell_mask)
def GetIntVector(self, v):
	if isinstance(v, np.ndarray) and isinstance(v[0].item(), int):
		vv = phreeqcrm.IntVector()
		for i in range(len(v)):
			vv.push_back(v[i].item())
		return vv
	if (isinstance(v, tuple) or isinstance(v, list)) and isinstance(v[0], int):
		vv = phreeqcrm.IntVector()
		for i in range(len(v)):
			vv.push_back(v[i])
		return vv
	return v
def GetDoubleVector(self, v):
	if isinstance(v, np.ndarray) and isinstance(v[0].item(), float):
		vv = phreeqcrm.DoubleVector()
		for i in range(len(v)):
			vv.push_back(v[i].item())
		return vv
	if (isinstance(v, tuple) or isinstance(v, list)) and isinstance(v[0], float):
		vv = phreeqcrm.DoubleVector()
		for i in range(len(v)):
			vv.push_back(v[i])
		return vv
	return v

%} 
}


#include "bmi.hxx"

// Ignore methods
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
//%ignore BMIPhreeqcRM::GetValuePtr(std::string);
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

%rename(add_output_vars)             AddOutputVars(std::string option, std::string def);
%rename(initialize)                  Initialize(std::string config_file="");
%rename(update)                      Update();
%rename(update_until)                UpdateUntil(double end_time);
%rename(finalize)                    Finalize();
%rename(get_component_name)          GetComponentName();
%rename(get_input_item_count)        GetInputItemCount();
%rename(get_output_item_count)       GetOutputItemCount();
%rename(get_pointable_item_count)    GetPointableItemCount();
%rename(get_input_var_names)         GetInputVarNames();
%rename(get_output_var_names)        GetOutputVarNames();
%rename(get_pointable_var_names)     GetPointableVarNames();
%rename(get_readonly_var_names)      GetReadOnlyVarNames();
%rename(get_var_grid)                GetVarGrid(const std::string name);
%rename(get_var_type)                GetVarType(const std::string name);
%rename(get_var_units)               GetVarUnits(const std::string name);
%rename(get_var_itemsize)            GetVarItemsize(const std::string name);
%rename(get_var_nbytes)              GetVarNbytes(const std::string name);
%rename(get_var_location)            GetVarLocation(const std::string name);
%rename(get_current_time)            GetCurrentTime();
%rename(get_start_time)              GetStartTime();
%rename(get_end_time)                GetEndTime();
%rename(get_time_units)              GetTimeUnits();
%rename(get_time_step)               GetTimeStep();
//%rename(get_value_ptr)               GetValuePtr(std::string name);
//%rename(get_value_at_indices)        GetValueAtIndices(std::string name, void* dest, int* inds, int count);
//%rename(set_value_at_indices)        SetValueAtIndices(std::string name, int* inds, int count, void* src);
%rename(get_grid_rank)               GetGridRank(const int grid);
%rename(get_grid_size)               GetGridSize(const int grid);
%rename(get_grid_type)               GetGridType(const int grid);
%rename(get_grid_shape)              GetGridShape(const int grid, int* shape);
%rename(get_grid_spacing)            GetGridSpacing(const int grid, double* spacing);
%rename(get_grid_origin)             GetGridOrigin(const int grid, double* origin);
%rename(get_grid_x)                  GetGridX(const int grid, double* x);
%rename(get_grid_y)                  GetGridY(const int grid, double* y);
%rename(get_grid_z)                  GetGridZ(const int grid, double* z);
%rename(get_grid_node_count)         GetGridNodeCount(const int grid);
%rename(get_grid_edge_count)         GetGridEdgeCount(const int grid);
%rename(get_grid_face_count)         GetGridFaceCount(const int grid);
%rename(get_grid_edge_nodes)         GetGridEdgeNodes(const int grid, int* edge_nodes);
%rename(get_grid_face_edges)         GetGridFaceEdges(const int grid, int* face_edges);
%rename(get_grid_face_nodes)         GetGridFaceNodes(const int grid, int* face_nodes);
%rename(get_grid_nodes_per_face)     GetGridNodesPerFace(const int grid, int* nodes_per_face);

#if defined(USE_YAML)
#if defined(SWIGPYTHON)

%exception Initialize {
  try {
    $action
  } catch (YAML::BadFile &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

%exception Initialize {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

%exception GetGridEdgeCount {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridEdgeNodes {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridFaceCount {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridFaceEdges {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridFaceNodes {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridNodeCount {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridNodesPerFace {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridOrigin {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

// Should this be implemented? @todo
%exception GetGridShape {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridSpacing {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridX {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridY {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception GetGridZ {
  try {
    $action
  } catch (NotImplemented &e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    SWIG_fail;
  }
}

%exception get_value_ptr_int {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

%exception get_value_ptr_double {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

%exception get_value_ptr_vector_strings {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

%exception SetValue {
  try {
    $action
  } catch (std::runtime_error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}

#endif
#endif

// %numpy_typemaps(int,    NPY_INT32  , int)

%feature("pythonprepend") BMIPhreeqcRM::BMIPhreeqcRM() %{
    self._state = State.UNINITIALIZED
    self._values = {}
    self._pointables = {}
    self._readonlys = {}
%}

%feature("pythonappend") BMIPhreeqcRM::Initialize(std::string config_file="") %{
    self._pointables = {var.lower() for var in self.get_pointable_var_names()}
    self._readonlys = {var.lower() for var in self.get_readonly_var_names()}
    self._state = State.INITIALIZED
%}

%feature("pythonappend") BMIPhreeqcRM::Finalize() %{
    self._state = State.UNINITIALIZED
    self._values = {}
    self._pointables = {}
    self._readonlys = {}
%}

// initialize checks

%feature("pythonprepend") BMIPhreeqcRM::GetGridSize(const int grid) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetInputItemCount() %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetInputVarNames() %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetOutputItemCount() %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetOutputVarNames() %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetVarItemsize(const std::string name) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetVarNbytes(const std::string name) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetVarType(const std::string name) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::GetVarUnits(const std::string name) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::Update() %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%feature("pythonprepend") BMIPhreeqcRM::UpdateUntil(double end_time) %{
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")
%}

%include "../src/BMIPhreeqcRM.h"

%inline %{
void BMIPhreeqcRM::get_value_ptr_double(std::string var, double **ARGOUTVIEW_ARRAY1, int* DIM1) {
    *DIM1 = 0;
    *ARGOUTVIEW_ARRAY1 = nullptr;
    RMVARS v_enum = this->var_man->GetEnum(var);
    if (v_enum != RMVARS::NotFound)
    {
        BMIVariant& bv = this->var_man->VariantMap[v_enum];
        if (bv.GetVoidPtr() == nullptr)
        {
            this->var_man->task = VarManager::VAR_TASKS::GetPtr;
            ((*this->var_man).*bv.GetFn())();
        }
        *DIM1 = bv.GetDim();
        *ARGOUTVIEW_ARRAY1 = (double*)bv.GetVoidPtr();
    }
}

void BMIPhreeqcRM::get_value_ptr_int(std::string var, int** ARGOUTVIEW_ARRAY1, int *DIM1) {
    *DIM1 = 0;
    *ARGOUTVIEW_ARRAY1 = nullptr;
    RMVARS v_enum = this->var_man->GetEnum(var);
    if (v_enum != RMVARS::NotFound)
    {
        BMIVariant& bv = this->var_man->VariantMap[v_enum];
        if (bv.GetVoidPtr() == nullptr)
        {
            this->var_man->task = VarManager::VAR_TASKS::GetPtr;
            ((*this->var_man).*bv.GetFn())();
        }
        *DIM1 = bv.GetDim();
        *ARGOUTVIEW_ARRAY1 = (int*)bv.GetVoidPtr();
    }
}

std::vector<std::string>& BMIPhreeqcRM::get_value_ptr_vector_strings(std::string var)
{
    static std::vector<std::string> err = { "BAD Variable Name" };
    RMVARS v_enum = this->var_man->GetEnum(var);
    if (v_enum != RMVARS::NotFound)
    {
        BMIVariant& bv = this->var_man->VariantMap[v_enum];
        if (bv.GetVoidPtr() == nullptr)
        {
            this->var_man->task = VarManager::VAR_TASKS::GetPtr;
            ((*this->var_man).*bv.GetFn())();
        }
        return bv.GetStringVectorRef();
    }
    return err;
}
%}


// Write new python method GetValue with one argument
%extend BMIPhreeqcRM { %pythoncode 
%{ 
def get_value_ptr(self, var_name):

    if self._state != State.INITIALIZED:
      raise RuntimeError("must call initialize first")

    var_name_lower = var_name.lower()
    if var_name_lower in self._values:
        return self._values[var_name_lower]

    Nbytes = self.get_var_nbytes(var_name_lower)
    Itemsize = self.get_var_itemsize(var_name_lower)
    dim = 0
    if Itemsize != 0:
        dim = Nbytes / Itemsize
    vtype = self.get_var_type(var_name_lower)

    if vtype=="float64":
        if dim>=1:
            arry = self.get_value_ptr_double(var_name_lower)
            self._values[var_name_lower] = arry
            if var_name_lower in self._readonlys:
                arry.flags.writeable = False
            return arry
    if vtype=="int32":
        if dim>=1:
            arry = self.get_value_ptr_int(var_name_lower)
            self._values[var_name_lower] = arry
            if var_name_lower in self._readonlys:
                arry.flags.writeable = False
            return arry
    if vtype.startswith("<U"):
        arry = np.array(self.get_value_ptr_vector_strings(var_name_lower))
        self._values[var_name_lower] = arry
        # strings are always readonly
        arry.flags.writeable = False
        return arry
    return None

def get_value(self, var_name, dest):

    if self._state != State.INITIALIZED:
      raise RuntimeError("must call initialize first")

    var_name_lower = var_name.lower()
    if var_name_lower in self._pointables:
        dest[:] = self.get_value_ptr(var_name_lower).flatten()
        return dest

    Nbytes = self.get_var_nbytes(var_name_lower)
    Itemsize = self.get_var_itemsize(var_name_lower)
    dim = 0
    if Itemsize != 0:
        dim = Nbytes / Itemsize
    vtype = self.get_var_type(var_name_lower)

    if vtype=="float64":
        if dim==1:
            dest[:] = np.array(self.GetValue_double(var_name_lower)).flatten()
        if dim>1:
            dest[:] = np.array(self.GetValue_double_vector(var_name_lower)).flatten()
    if vtype=="int32":
        if dim==1:
            dest[:] = np.array(self.GetValue_int(var_name_lower)).flatten()
        if dim>1:
            dest[:] = np.array(self.GetValue_int_vector(var_name_lower)).flatten()
    if vtype.startswith("<U"):
        if dim==1:
            dest[:] = np.array([self.GetValue_string(var_name_lower)]).flatten()
        if dim>1:
            dest[:] = np.array(self.GetValue_string_vector(var_name_lower)).flatten()
    if vtype=="bool":
        dest[:] = np.array([self.GetValue_bool(var_name_lower)]).flatten()
    return dest

def get_value_at_indices(self, var_name, dest, indices):
    """Get values at particular indices.

    Parameters
    ----------
    var_name : str
        Name of variable as CSDMS Standard Name.
    dest : ndarray
        A numpy array into which to place the values.
    indices : array_like
        Array of indices.

    Returns
    -------
    array_like
        Values at indices.
    """
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")

    dest[:] = self.get_value_ptr(var_name).take(indices)
    return dest

def set_value_at_indices(self, var_name, inds, src):
    """Set model values at particular indices.

    Parameters
    ----------
    var_name : str
        Name of variable as CSDMS Standard Name.
    src : array_like
        Array of new values.
    indices : array_like
        Array of indices.
    """
    if self._state != State.INITIALIZED:
        raise RuntimeError("must call initialize first")

    val = self.get_value_ptr(var_name)
    val.flat[inds] = src

def set_value(self, var_name, value):
    '''@todo
    '''
    if self._state != State.INITIALIZED:
      raise RuntimeError("must call initialize first")

    Nbytes = self.get_var_nbytes(var_name)
    Itemsize = self.get_var_itemsize(var_name)
    dim = Nbytes / Itemsize
    vtype = self.get_var_type(var_name)
    if vtype=="float64":
        if dim==1:
            self.SetValue_double(var_name, value[0])
        if dim>1:
            self.SetValue_double_vector(var_name, value)
    if vtype=="int32":
        if dim==1:
            self.SetValue_int(var_name, value[0])
        if dim>1:
            # this needs testing
            self.SetValue_int_vector(var_name, value)
    if vtype.startswith("<U"):
        if dim==1:
            self.SetValue_string(var_name, value[0])
        if dim>1:
            self.SetValue_string_vector(var_name, value)
%}
};
