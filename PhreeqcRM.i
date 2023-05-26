// %module(directors="1") phreeqcrm
%module phreeqcrm
%pythoncode 
%{ 
import numpy as np
import phreeqcrm
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
%}
#%fragment("NumPy_Fragments"); 
%ignore BMIVariant;
%ignore bmi::Bmi;

#if 1
%include "numpy.i"
%init %{
import_array();
%}
#endif

%include "std_string.i"
%include "std_vector.i"
%template(BoolVector)   std::vector<bool>;
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
%include "BMIVariant.h"
%include "IrmResult.h"

// Ignore methods
%ignore PhreeqcRM::GetIPhreeqcPointer(int i);
%ignore PhreeqcRM::GetWorkers();
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

// Ignore methods
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
%ignore BMIPhreeqcRM::GetValuePtr(std::string);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< int,std::allocator< int > >);
%ignore BMIPhreeqcRM::SetValue(std::string,std::vector< std::string,std::allocator< std::string > >);


%include "PhreeqcRM.h"

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

%} 
}


#include "bmi.hxx"

// Ignore methods
%ignore BMIPhreeqcRM::GetValue(std::string const,bool *);
%ignore BMIPhreeqcRM::GetValue(std::string const,double *);
%ignore BMIPhreeqcRM::GetValue(std::string const,int *);
%ignore BMIPhreeqcRM::GetValuePtr(std::string);
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
	dim = 0
	if Itemsize != 0:
		dim = Nbytes / Itemsize
	type = self.GetVarType(var)
	#print(f"Type={type}")
	if type=="double":
		if dim==1:
			return self.GetValue_double(var)
		if dim>1:
			return np.array(self.GetValue_double_vector(var))
	if type=="int":
		if dim==1:
			return self.GetValue_int(var)
		if dim>1:
			return np.array(self.GetValue_int_vector(var))
	if type=="std::string":
		return self.GetValue_string(var)
	if type=="std::vector<std::string>":
		return np.array(self.GetValue_string_vector(var))
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
		self.SetValue_string(var, value)
	if type=="std::vector<std::string>":
		self.SetValue_string_vector(var, value)	
def GetValuePtr(self):
	return "Not Implemented."
%} 
};


