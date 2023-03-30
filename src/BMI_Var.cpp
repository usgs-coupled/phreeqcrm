#include "BMI_Var.h"
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
	b_var = false;
	i_var = 0;
	d_var = 0.0;
	NotImplemented = false;
	VoidPtr = NULL;
	fn = f;
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
	this->b_var = bv.b_var;
	this->i_var = bv.i_var;
	this->d_var = bv.d_var;
	this->NotImplemented = bv.NotImplemented;
	this->VoidPtr = bv.VoidPtr;
	this->fn = bv.fn;
	this->NotImplemented = bv.NotImplemented;
}
#ifdef SKIP
varfn_map["components"] = &Components_var;
varfn_map["componentcount"] = &ComponentCount_var;
varfn_map["concentrations"] = &Concentrations_var;
varfn_map["density"] = &Density_var;
varfn_map["errorstring"] = &ErrorString_var;
varfn_map["fileprefix"] = &FilePrefix_var;
varfn_map["gfw"] = &Gfw_var;
varfn_map["gridcellcount"] = &GridCellCount_var;
varfn_map["inputvarnames"] = &InputVarNames_var;
varfn_map["nthselectedoutput"] = &NthSelectedOutput_var;
varfn_map["outputvarnames"] = &OutputVarNames_var;
varfn_map["saturation"] = &Saturation_var;
varfn_map["selectedoutput"] = &SelectedOutput_var;
varfn_map["selectedoutputcolumncount"] = &SelectedOutputColumnCount_var;
varfn_map["selectedoutputcount"] = &SelectedOutputCount_var;
varfn_map["selectedoutputheadings"] = &SelectedOutputHeadings_var;
varfn_map["selectedoutputrowcount"] = &SelectedOutputRowCount_var;
varfn_map["solutionvolume"] = &SolutionVolume_var;
varfn_map["time"] = &Time_var;
varfn_map["timestep"] = &TimeStep_var;
varfn_map["currentselectedoutputusernumber"] = &CurrentSelectedOutputUserNumber_var;
varfn_map["porosity"] = &Porosity_var;
varfn_map["pressure"] = &Pressure_var;
varfn_map["selectedoutputon"] = &SelectedOutputOn_var;
varfn_map["temperature"] = &Temperature_var;
#endif