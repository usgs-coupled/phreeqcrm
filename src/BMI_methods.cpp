// PhreeqcRM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>
#include <string>
#include <iomanip>
#include <sstream>
#include "PhreeqcRM.h"
#include "IPhreeqc.h"
#include "IPhreeqcPhast.h"
#include "BMI_Var.h"
void PhreeqcRM::BMI_SetValue(std::string name, void* src)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        if (it->first == "Concentrations")
        {
            int ngrid = this->GetGridCellCount();
            int ncomps = this->GetComponentCount();
            std::vector<double> conc(ngrid * ncomps, INACTIVE_CELL_VALUE);
            memcpy(conc.data(), src, ngrid * ncomps * sizeof(double));
            this->SetConcentrations(conc);
            return;
        }
        if (it->first == "Density")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> density(ngrid, INACTIVE_CELL_VALUE);
            memcpy(density.data(), src, ngrid * sizeof(double));
            this->SetDensity(density);
            return;
        }
        if (it->first == "NthSelectedOutput")
        {
            int nth_so;
            memcpy(&nth_so, src, sizeof(int));
            int nuser = this->GetNthSelectedOutputUserNumber(nth_so);
            this->SetCurrentSelectedOutputUserNumber(nuser);
            return;
        }
        if (it->first == "Porosity")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> porosity(ngrid, INACTIVE_CELL_VALUE);
            memcpy(porosity.data(), src, ngrid * sizeof(double));
            this->SetPorosity(porosity);
            return;
        }
        if (it->first == "Pressure")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> pressure(ngrid, INACTIVE_CELL_VALUE);
            memcpy(pressure.data(), src, ngrid * sizeof(double));
            this->SetPressure(pressure);
            return;
        }
        if (it->first == "Saturation")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> sat(ngrid, INACTIVE_CELL_VALUE);
            memcpy(sat.data(), src, ngrid * sizeof(double));
            this->SetSaturation(sat);
            return;
        }
        if (it->first == "SelectedOutputOn")
        {
            int so_on;
            memcpy(&so_on, src, sizeof(int));
            bool so_on_bool = (bool)so_on;
            this->SetSelectedOutputOn(so_on_bool);
            return;
        }
        if (it->first == "Temperature")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> temp(ngrid, INACTIVE_CELL_VALUE);
            memcpy(temp.data(), src, ngrid * sizeof(double));
            this->SetTemperature(temp);
            return;
        }
    }
    throw LetItThrow("Item not found");
}
void PhreeqcRM::BMI_GetValue(std::string name, void* dest)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        if (it->first == "ComponentCount")
        {
            int count = this->GetComponentCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "Components")
        {
            int string_length = this->BMI_GetVarItemsize("Components");
            std::vector<std::string> comps = this->GetComponents();
            std::stringstream all_comps;
            for (size_t i = 0; i < comps.size(); i++)
            {
                all_comps << std::left << std::setfill(' ') << std::setw(string_length) << comps[i];
            }
            memcpy(dest, all_comps.str().c_str(), all_comps.str().size());
            return;
        }
        if (it->first == "Concentrations")
        {
            std::vector<double> rm_conc;
            this->GetConcentrations(rm_conc);
            memcpy(dest, rm_conc.data(), rm_conc.size()*sizeof(double));
            return;
        }
        if (it->first == "CurrentSelectedOutputUserNumber")
        {
            int count = this->GetIPhreeqcPointer(0)->GetCurrentSelectedOutputUserNumber();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "Density")
        {
            std::vector<double> density;
            this->GetDensity(density);
            memcpy(dest, density.data(), density.size()*sizeof(double));
            return;
        }
        if (it->first == "ErrorString")
        {
            std::string err = this->GetErrorString();
            memcpy(dest, err.c_str(), err.size());
            return;
        }
        if (it->first == "Gfw")
        {
            const std::vector<double> gfw = this->GetGfw();
            memcpy(dest, gfw.data(), gfw.size() * sizeof(double));
            return;
        }
        if (it->first == "GridCellCount")
        {
            int count = this->GetGridCellCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        //if (it->first == "NthSelectedOutputUserNumber")
        //{
        //    int num = this->GetNthSelectedOutputUserNumber();
        //    memcpy(dest, &num, sizeof(int));
        //    return;
        //}
        if (it->first == "SelectedOutput")
        {
            std::vector<double> so;
            IRM_RESULT status = this->GetSelectedOutput(so);
            memcpy(dest, so.data(), so.size() * sizeof(double));
            return;
        }
        if (it->first == "SelectedOutputColumnCount")
        {
            int count = this->GetSelectedOutputColumnCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "SelectedOutputCount")
        {
            int count = this->GetSelectedOutputCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "SelectedOutputHeadings")
        {
            int string_length = this->BMI_GetVarItemsize("SelectedOutputHeadings");
            int count = this->GetSelectedOutputColumnCount();
            std::stringstream all_headings;
            for (int i = 0; i < count; i++)
            {
                std::string heading;
                IRM_RESULT status = this->GetSelectedOutputHeading(i, heading);
                all_headings << std::left << std::setfill(' ') << std::setw(string_length) << heading;
            }
            memcpy(dest, all_headings.str().c_str(), all_headings.str().size());
            return;
        }
        if (it->first == "SelectedOutputOn")
        {
            int flag = (int)this->GetSelectedOutputOn();
            memcpy(dest, &flag, sizeof(int));
            return;
        }
        if (it->first == "SelectedOutputRowCount")
        {
            int count = this->GetSelectedOutputRowCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
    }
    throw LetItThrow("Item not found");
}
int PhreeqcRM::BMI_GetVarNbytes(std::string name)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        if (it->first == "ComponentCount")
        {
            return this->GetComponentCount();
        }
        if (it->first == "Components")
        {
            int string_size = this->BMI_GetVarItemsize("Components");
            int dim = this->GetComponentCount();
            return string_size * dim * (int)sizeof(char);
        }
        if (it->first == "Concentrations")
        {
            return (int)sizeof(double) * this->GetGridCellCount() * this->GetComponentCount();
        }
        if (it->first == "CurrentSelectedOutputUserNumber")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Density")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "ErrorString")
        {
            return (int)this->GetErrorString().size();
        }
        if (it->first == "Gfw")
        {
            return (int)sizeof(double) * this->GetComponentCount();
        }
        if (it->first == "GridCellCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "NthSelectedOutput")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Porosity")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "Pressure")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "Saturation")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "SelectedOutput")
        {
            return (int)sizeof(double) * this->GetSelectedOutputColumnCount() * this->GetSelectedOutputRowCount();
        }
        if (it->first == "SelectedOutputColumnCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputHeadings")
        {
            int string_size = this->BMI_GetVarItemsize("SelectedOutputHeadings");
            return string_size * this->GetSelectedOutputColumnCount() * (int)sizeof(char);
        }
        if (it->first == "SelectedOutputOn")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputRowCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Temperature")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
    }
	throw LetItThrow("Item not found");
}
std::string PhreeqcRM::BMI_GetVarUnits(std::string name)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        if (it->first == "Concentrations")
        {
            int units = this->GetUnitsSolution();
            switch (units)
            {
            case 1:
                return "mg L-1";
                break;
            case 2:
                return "mol L-1";
                    break;
            case 3:
                return "kg kgs-1";
				break;
			default:
				throw LetItThrow("Unknown units for GetUnitsSolution in GetVarUnits");
				break;
			}
        }
        return it->second.GetUnits();
    }
    throw LetItThrow("Item not found");
}
int PhreeqcRM::BMI_GetVarItemsize(std::string name)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        if (it->first == "ComponentCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Components")
        {
            const std::vector< std::string >& ref = this->GetComponents();
            size_t max = 0;
            for (size_t i = 0; i < this->GetComponents().size(); i++)
            {
                if (this->GetComponents()[i].size() > max) max = this->GetComponents()[i].size();
            }
            return (int)max;
        }
        if (it->first == "Concentrations")
        {
            return (int)sizeof(double);
        }
        if (it->first == "CurrentSelectedOutputUserNumber")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Density")
		{
            return (int)sizeof(double);
		}
        if (it->first == "ErrorString")
        {
            return (int)this->GetErrorString().size();
        }
        if (it->first == "Gfw")
		{
            return (int)sizeof(double);
		}
        if (it->first == "GridCellCount")
		{
            return (int)sizeof(int);
		}
		if (it->first == "NthSelectedOutput")
		{
            return (int)sizeof(int);
		}
		if (it->first == "Porosity")
		{
            return (int)sizeof(double);
		}
        if (it->first == "Pressure")
        {
            return (int)sizeof(double);
        }
        if (it->first == "Saturation")
        {
            return (int)sizeof(double);
        }
        if (it->first == "SelectedOutput")
        {
            return (int)sizeof(double);
        }
        if (it->first == "SelectedOutputColumnCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputHeadings")
        {
            int nhead = this->GetSelectedOutputColumnCount();
            size_t max = 0;
            for (int i = 0; i < nhead; i++)
            {
                std::string heading;
                this->GetSelectedOutputHeading(i, heading);
                if (heading.size() > max) max = heading.size();
            }
            return (int)max;
        }
        if (it->first == "SelectedOutputOn")
        {
            return (int)sizeof(int);
        }
        if (it->first == "SelectedOutputRowCount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "Temperature")
        {
            return (int)sizeof(double);
        }
    }
    else
    {
        throw LetItThrow("Item not found");
        return -1;
    }
    return -1;
}
std::string PhreeqcRM::BMI_GetVarType(std::string name)
{
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name);
    if (it != bmi_var_map.end())
    {
        return it->second.GetType();
    }
    else
    {
        throw LetItThrow("Item not found");
    }
}
void PhreeqcRM::BMI_MakeVarMap()
{
    //if (state >= 0)
    {
        //-1 before create
        //0  RM created
        //1  map created
        //2  PhreeqcInitial2module
        //var_map["ncells"] = Var_BMI("ncells", "int", "count", sizeof(int));
        //var_map["threads"] = Var_BMI("threads", "int", "count", sizeof(int));

        //var_map["BackwardMapping"] = Var_BMI("BackwardMapping", "int", "mapping", sizeof(int));
        //var_map["CreateMapping"] = Var_BMI("CreateMapping", "int", "mapping", sizeof(int));
        //var_map["ChemistryCellCount"] = Var_BMI("ChemistryCellCount", "int", "count", sizeof(int));
        bmi_var_map["Components"] = BMI_Var("Components", "string", "names", false, true);
        bmi_var_map["ComponentCount"] = BMI_Var("ComponentCount", "int", "names", false, true);
        bmi_var_map["Concentrations"] = BMI_Var("Concentrations", "double", "mol L-1", true, true);
        bmi_var_map["Density"] = BMI_Var("Density", "double", "kg L-1", true, true);
        //var_map["EndCell"] = Var_BMI("EndCell", "int", "cell numbers", sizeof(int));
        //var_map["EquilibriumPhasesNames"] = Var_BMI("EquilibriumPhasesNames", "string", "names", sizeof(char));
        bmi_var_map["ErrorString"] = BMI_Var("ErrorString", "string", "error", false, true);
        //var_map["ExchangeNames"] = Var_BMI("ExchangeNames", "string", "names", sizeof(char));
        //var_map["ExchangeSpeciesNames"] = Var_BMI("ExchangeSpeciesNames", "string", "names", sizeof(char));
        //var_map["FilePrefix"] = Var_BMI("FilePrefix", "string", "name", sizeof(char));
        //!var_map["GasComponentNames"] = Var_BMI("GasComponentsNames", "string", "names", sizeof(char));
        //!var_map["GasCompMoles"] = Var_BMI("GasCompMoles", "double", "mol", sizeof(double));
        //!var_map["GasCompPressures"] = Var_BMI("GasCompPressures", "double", "atm", sizeof(double));
        //!var_map["GasCompPhi"] = Var_BMI("GasCompPhi", "double", "atm-1", sizeof(double));
        //!var_map["GasPhaseVolume"] = Var_BMI("GasPhaseVolume", "double", "L", sizeof(double));
        bmi_var_map["Gfw"] = BMI_Var("Gfw", "double", "g mol-1", false, true);
        bmi_var_map["GridCellCount"] = BMI_Var("GridCellCount", "int", "count", false, true);
        //var_map["KineticReactions"] = Var_BMI("KineticReactions", "string", "names", sizeof(char));
        //var_map["MpiMyself"] = Var_BMI("MpiMyself", "int", "id", sizeof(int));
        //var_map["MpiTasks"] = Var_BMI("MpiTasks", "int", "count", sizeof(int));
        bmi_var_map["NthSelectedOutput"] = BMI_Var("NthSelectedOutput", "int", "id", true, false);
        bmi_var_map["Saturation"] = BMI_Var("Saturation", "double", "unitless", true, false);
        bmi_var_map["SelectedOutput"] = BMI_Var("SelectedOutput", "var", "user", false, true);
        bmi_var_map["SelectedOutputColumnCount"] = BMI_Var("SelectedOutputColumnCount", "int", "count", false, true);
        bmi_var_map["SelectedOutputCount"] = BMI_Var("SelectedOutputCount", "int", "count", false, true);
        bmi_var_map["SelectedOutputHeadings"] = BMI_Var("SelectedOutputHeadings", "string", "names", false, true);
        bmi_var_map["SelectedOutputRowCount"] = BMI_Var("SelectedOutputRowCount", "int", "count", false, true);
        //var_map["SINames"] = Var_BMI("SINames", "string", "names", 0);
        //var_map["SolidSolutionComponentsNames"] = Var_BMI("SolidSolutionComponentsNames", "string", "names", sizeof(char));
        //var_map["SolidSolutionNames"] = Var_BMI("SolidSolutionNames", "string", "names", sizeof(char));
        //var_map["SolutionVolume"] = Var_BMI("SolutionVolume", "double", "L", sizeof(double));
        //!var_map["SpeciesConcentrations"] = Var_BMI("SpeciesConcentrations", "double", "mg L-1", sizeof(double));
        //!var_map["SpeciesD25"] = Var_BMI("SpeciesD25", "double", "cm2 s-1", sizeof(double));
        //!var_map["SpeciesLog10Gammas"] = Var_BMI("SpeciesLog10Gammas", "double", "log L mol-1", sizeof(double));
        //!var_map["SpeciesLog10Molalities"] = Var_BMI("SpeciesLog10Molalities", "double", "log mol L-1", sizeof(double));
        //!var_map["SpeciesNames"] = Var_BMI("SpeciesNames", "string", "names", sizeof(char));
        //!var_map["SpeciesSaveOn"] = Var_BMI("SpeciesSaveOn", "int", "flag", sizeof(char));
        //!var_map["SpeciesZ"] = Var_BMI("SpeciesZ", "double", "charge number", sizeof(double));
        //var_map["StartCell"] = Var_BMI("StartCell", "int", "cell numbers", sizeof(int));
        //var_map["SurfaceNames"] = Var_BMI("SurfaceNames", "string", "names", sizeof(char));
        //var_map["SurfaceSpeciesNames"] = Var_BMI("SurfaceSpeciesNames", "string", "names", sizeof(char));
        //var_map["SurfaceTypes"] = Var_BMI("SurfaceTypes", "string", "names", sizeof(char));
        //var_map["ThreadCount"] = Var_BMI("ThreadCount", "int", "count", sizeof(int));
        //var_map["Time"] = Var_BMI("Time", "double", "s", sizeof(double));
        //var_map["TimeConversion"] = Var_BMI("TimeConversion", "double", "unitless", sizeof(double));
        //var_map["TimeStep"] = Var_BMI("TimeStep", "double", "s", sizeof(double));
        //var_map["MpiWorker"] = Var_BMI("MpiWorker", "int", "id", sizeof(int));
        //var_map["ComponentH2O"] = Var_BMI("ComponentH2O", "int", "flag", sizeof(int));
        bmi_var_map["CurrentSelectedOutputUserNumber"] = BMI_Var("CurrentSelectedOutputUserNumber", "int", "id", false, true);
        //var_map["DumpFileName"] = Var_BMI("DumpFileName", "string", "name", sizeof(char));
        //var_map["ErrorHandlerMode"] = Var_BMI("ErrorHandlerMode", "int", "flag", sizeof(int));
        //var_map["ErrorOn"] = Var_BMI("ErrorOn", "int", "flag", sizeof(int));
        //var_map["PartitionUZSolids"] = Var_BMI("PartitionUZSolids", "int", "flag", sizeof(int));
        bmi_var_map["Porosity"] = BMI_Var("Porosity", "double", "unitless", true, false);
        //var_map["PartitionUZSolids"] = Var_BMI("PartitionUZSolids", "int", "flag", sizeof(int));
        bmi_var_map["Pressure"] = BMI_Var("Pressure", "double", "atm", true, false);
        //var_map["PrintChemistryMask"] = Var_BMI("PrintChemistryMask", "int", "flags", sizeof(int));
        //var_map["PrintChemistryOn"] = Var_BMI("PrintChemistryOn", "int", "flag", sizeof(int));
        //var_map["RebalanceByCell"] = Var_BMI("RebalanceByCell", "int", "flag", sizeof(int));
        //var_map["RebalanceFraction"] = Var_BMI("RebalanceFraction", "double", "unitless", sizeof(double));
        //var_map["RepresentativeVolume"] = Var_BMI("RepresentativeVolume", "double", "L", sizeof(double));
        //var_map["RebalanceByCell"] = Var_BMI("RebalanceByCell", "int", "flag", sizeof(int));
        //var_map["ScreenOn"] = Var_BMI("ScreenOn", "int", "flag", sizeof(int));
        bmi_var_map["SelectedOutputOn"] = BMI_Var("SelectedOutputOn", "int", "flag", true, false);
        bmi_var_map["Temperature"] = BMI_Var("Temperature", "double", "C", true, false);
        //var_map["UnitsExchange"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsGasPhase"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsKinetics"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsPPassemblage"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSolution"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSSassemblage"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSurface"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UseSolutionDensityVolume"] = Var_BMI("UseSolutionDensityVolume", "int", "flag", sizeof(int));

    }
    this->bmi_input_vars.clear();
    this->bmi_output_vars.clear();
    std::map<std::string, class BMI_Var>::iterator it;
    for (it = bmi_var_map.begin(); it != bmi_var_map.end(); it++)
    {
        if (it->second.GetSet())
        {
            bmi_input_vars.push_back(it->first);
        }
        if (it->second.GetGet())
        {
            bmi_output_vars.push_back(it->first);
        }
    }
}

