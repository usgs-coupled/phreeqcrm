// PhreeqcRM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <iomanip>
#include <sstream>
#include "PhreeqcRM.h"
#include "IPhreeqc.h"
#include "IPhreeqcPhast.h"
#include "BMIVariant.h"
#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif

void PhreeqcRM::BMI_SetValue(std::string name, void* src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "concentrations")
        {
            int ngrid = this->GetGridCellCount();
            int ncomps = this->GetComponentCount();
            std::vector<double> conc(ngrid * ncomps, INACTIVE_CELL_VALUE);
            memcpy(conc.data(), src, ngrid * ncomps * sizeof(double));
            this->SetConcentrations(conc);
            return;
        }
        if (it->first == "density")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> density(ngrid, INACTIVE_CELL_VALUE);
            memcpy(density.data(), src, ngrid * sizeof(double));
            this->SetDensity(density);
            return;
        }
        if (it->first == "fileprefix")
        {
            std::string file = (char*)src;
            this->SetFilePrefix(file);
            return;
        }
        if (it->first == "nthselectedoutput")
        {
            int nth_so;
            memcpy(&nth_so, src, sizeof(int));
            int nuser = this->GetNthSelectedOutputUserNumber(nth_so - 1);
            this->SetCurrentSelectedOutputUserNumber(nuser);
            return;
        }
        if (it->first == "porosity")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> porosity(ngrid, INACTIVE_CELL_VALUE);
            memcpy(porosity.data(), src, ngrid * sizeof(double));
            this->SetPorosity(porosity);
            return;
        }
        if (it->first == "porosity")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> porosity(ngrid, INACTIVE_CELL_VALUE);
            memcpy(porosity.data(), src, ngrid * sizeof(double));
            this->SetPorosity(porosity);
            return;
        }
        if (it->first == "pressure")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> pressure(ngrid, INACTIVE_CELL_VALUE);
            memcpy(pressure.data(), src, ngrid * sizeof(double));
            this->SetPressure(pressure);
            return;
        }
        if (it->first == "saturation")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> sat(ngrid, INACTIVE_CELL_VALUE);
            memcpy(sat.data(), src, ngrid * sizeof(double));
            this->SetSaturation(sat);
            return;
        }
        if (it->first == "selectedoutputon")
        {
            int so_on;
            memcpy(&so_on, src, sizeof(int));
            bool so_on_bool = (bool)so_on;
            this->SetSelectedOutputOn(so_on_bool);
            return;
        }
        if (it->first == "temperature")
        {
            int ngrid = this->GetGridCellCount();
            std::vector<double> temp(ngrid, INACTIVE_CELL_VALUE);
            memcpy(temp.data(), src, ngrid * sizeof(double));
            this->SetTemperature(temp);
            return;
        }
        if (it->first == "time")
        {
            double time;
            memcpy(&time, src, sizeof(double));
            this->SetTime(time);
            return;
        }
        if (it->first == "timestep")
        {
            double timestep = 0;
            memcpy(&timestep, src, sizeof(double));
            this->SetTimeStep(timestep);
            return;
        }
    }
    ErrorMessage("Item not found");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_SetValue(std::string name, bool& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {

        if (it->first == "selectedoutputon")
        {
            this->SetSelectedOutputOn(src);
            return;
        }
    }
    ErrorMessage("Item not found for BMI_SetValue with bool argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_SetValue(std::string name, double& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {

        if (it->first == "time")
        {
            this->SetTime(src);
            return;
        }
        if (it->first == "timestep")
        {
            this->SetTimeStep(src);
            return;
        }
    }
    ErrorMessage("Item not found for BMI_SetValue with double argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_SetValue(const std::string name, int& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "nthselectedoutput")
        {
            this->SetNthSelectedOutput(src);
            return;
        }
    }
    ErrorMessage("Item not found for BMI_SetValue with int argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_SetValue(std::string name, const std::string& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "fileprefix")
        {
            this->SetFilePrefix(src);
            return;
        }
    }
    ErrorMessage("Item not found for BMI_SetValue with std::string argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_SetValue(std::string name, std::vector< double >& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "concentrations")
        {
            int ngrid = this->GetGridCellCount();
            int ncomps = this->GetComponentCount();
            if ((int)src.size() != ngrid * ncomps)
            {
                ErrorMessage("Dimension error for concentration vector.");
            }
            this->SetConcentrations(src);
            return;
        }
        if (it->first == "density")
        {
            int ngrid = this->GetGridCellCount();
            if ((int)src.size() != ngrid )
            {
                ErrorMessage("Dimension error for density vector.");
            }
            this->SetDensity(src);
            return;
        }

        if (it->first == "porosity")
        {
            int ngrid = this->GetGridCellCount();
            if ((int)src.size() != ngrid)
            {
                ErrorMessage("Dimension error for porosity vector.");
            }
            this->SetPorosity(src);
            return;
        }
        if (it->first == "pressure")
        {
            int ngrid = this->GetGridCellCount();
            if ((int)src.size() != ngrid)
            {
                ErrorMessage("Dimension error for pressure vector.");
            }
            this->SetPressure(src);
            return;
        }
        if (it->first == "saturation")
        {
            int ngrid = this->GetGridCellCount();
            if ((int)src.size() != ngrid)
            {
                ErrorMessage("Dimension error for saturation vector.");
            }
            this->SetSaturation(src);
            return;
        }
        if (it->first == "temperature")
        {
            int ngrid = this->GetGridCellCount();
            if ((int)src.size() != ngrid)
            {
                ErrorMessage("Dimension error for temperature vector.");
            }
            this->SetTemperature(src);
            return;
        }
    }
    ErrorMessage("Item not found for BMI_SetValue with std::vector < double > argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_SetValue(std::string name, std::vector< int>& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
    }
    ErrorMessage("Item not found for BMI_SetValue with std::vector < int > argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_SetValue(std::string name, std::vector<std::string>& src)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
    }
    ErrorMessage("Item not found for BMI_SetValue with std::vector < std::string > argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_GetValue(std::string name, bool& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {

        if (it->first == "selectedoutputon")
        {
            dest = this->GetSelectedOutputOn();
            return;
        }
    }
    ErrorMessage("Item not found for BMI_GetValue with bool argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_GetValue(std::string name, double& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "time")
        {
            dest = this->GetTime();
            return;
        }
        if (it->first == "timestep")
        {
            dest = this->GetTimeStep();
            return;
        }
    }
    //throw LetItThrow("Item not found");
    ErrorMessage("Item not found for BMI_GetValue with double argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_GetValue(std::string name, int& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "componentcount")
        {
            dest = this->GetComponentCount();
            return;
        }
        if (it->first == "currentselectedoutputusernumber")
        {
            dest = this->GetIPhreeqcPointer(0)->GetCurrentSelectedOutputUserNumber();
            return;
        }
        if (it->first == "gridcellcount")
        {
            dest = this->GetGridCellCount();
            return;
        }
        if (it->first == "selectedoutputcolumncount")
        {
            dest = this->GetSelectedOutputColumnCount();
            return;
        }
        if (it->first == "selectedoutputcount")
        {
            dest = this->GetSelectedOutputCount();
            return;
        }
        if (it->first == "selectedoutputrowcount")
        {
            dest = this->GetSelectedOutputRowCount();
            return;
        }
    }
    ErrorMessage("Item not found for BMI_GetValue with double argument.");
    throw PhreeqcRMStop();
}

void PhreeqcRM::BMI_GetValue(std::string name, std::string& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "fileprefix")
        {
            dest = this->GetFilePrefix();
            return;
        }
        if (it->first == "errorstring")
        {
            dest = this->GetErrorString();
            return;
        }
    }
    ErrorMessage("Item not found for BMI_GetValue with std::string argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_GetValue(std::string name, std::vector < double >& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "concentrations")
        {
            this->GetConcentrations(dest);
            return;
        }
        if (it->first == "density")
        {
            this->GetDensity(dest);
            return;
        }
        if (it->first == "gfw")
        {
            dest = this->GetGfw();
            return;
        }
        if (it->first == "porosity")
        {
            dest = this->GetPorosity();
            return;
        }
        if (it->first == "pressure")
        {
            dest = this->GetPressure();
            return;
        }
        if (it->first == "saturation")
        {
            this->GetSaturation(dest);
            return;
        }
        if (it->first == "selectedoutput")
        {
            IRM_RESULT status = this->GetSelectedOutput(dest);
            return;
        }
        if (it->first == "solutionvolume")
        {
            dest = this->GetSolutionVolume();
            return;
        }
        if (it->first == "temperature")
        {
           dest = this->GetTemperature();
           return;
        }

    }
    ErrorMessage("Item not found for BMI_GetValue with std::vector < double > argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_GetValue(std::string name, std::vector < std::string >& dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "components")
        {
            dest = this->GetComponents();
            return;
        }

        if (it->first == "inputvarnames")
        {
            dest = this->BMI_GetInputVarNames();
            return;
        }

        if (it->first == "outputvarnames")
        {
            dest = this->BMI_GetOutputVarNames();
            return;
        }
         if (it->first == "selectedoutputheadings")
        {
            int count = this->GetSelectedOutputColumnCount();
            dest.clear();
            for (int i = 0; i < count; i++)
            {
                std::string heading;
                IRM_RESULT status = this->GetSelectedOutputHeading(i, heading);
                dest.push_back(heading);
            }
            return;
        }
    }
    ErrorMessage("Item not found for BMI_GetValue with std::vector < std::string > argument.");
    throw PhreeqcRMStop();
}
void PhreeqcRM::BMI_GetValue(std::string name, void* dest)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "componentcount")
        {
            int count = this->GetComponentCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "components")
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
        if (it->first == "concentrations")
        {
            std::vector<double> rm_conc;
            this->GetConcentrations(rm_conc);
            memcpy(dest, rm_conc.data(), rm_conc.size() * sizeof(double));
            return;
        }
        if (it->first == "currentselectedoutputusernumber")
        {
            int count = this->GetIPhreeqcPointer(0)->GetCurrentSelectedOutputUserNumber();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "density")
        {
            std::vector<double> density;
            this->GetDensity(density);
            memcpy(dest, density.data(), density.size() * sizeof(double));
            return;
        }
        if (it->first == "fileprefix")
        {
            std::string filep = this->GetFilePrefix();
            memcpy(dest, filep.c_str(), filep.size());
            return;
        }
        if (it->first == "errorstring")
        {
            std::string err = this->GetErrorString();
            memcpy(dest, err.c_str(), err.size());
            return;
        }
        if (it->first == "gfw")
        {
            const std::vector<double> gfw = this->GetGfw();
            memcpy(dest, gfw.data(), gfw.size() * sizeof(double));
            return;
        }
        if (it->first == "gridcellcount")
        {
            int count = this->GetGridCellCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "inputvarnames")
        {
            int string_length = this->BMI_GetVarItemsize("inputvarnames");
            std::vector < std::string > varnames = this->BMI_GetInputVarNames();
            std::stringstream all_varnames;
            for (size_t i = 0; i < varnames.size(); i++)
            {
                all_varnames << std::left << std::setfill(' ') << std::setw(string_length) << varnames[i];
            }
            memcpy(dest, all_varnames.str().c_str(), all_varnames.str().size());
            return;
        }
        //if (it->first == "NthSelectedOutputUserNumber")
        //{
        //    int num = this->GetNthSelectedOutputUserNumber();
        //    memcpy(dest, &num, sizeof(int));
        //    return;
        //}
        if (it->first == "outputvarnames")
        {
            int string_length = this->BMI_GetVarItemsize("outputvarnames");
            std::vector < std::string > varnames = this->BMI_GetOutputVarNames();
            std::stringstream all_varnames;
            for (size_t i = 0; i < varnames.size(); i++)
            {
                all_varnames << std::left << std::setfill(' ') << std::setw(string_length) << varnames[i];
            }
            memcpy(dest, all_varnames.str().c_str(), all_varnames.str().size());
            return;
        }
        if (it->first == "porosity")
        {
            const std::vector<double>& porosity = this->GetPorosity();
            memcpy(dest, porosity.data(), porosity.size() * sizeof(double));
            return;
        }
        if (it->first == "pressure")
        {
            const std::vector<double>& pressure = this->GetPressure();
            memcpy(dest, pressure.data(), pressure.size() * sizeof(double));
            return;
        }
        if (it->first == "saturation")
        {
            std::vector<double> saturation;
            this->GetSaturation(saturation);
            memcpy(dest, saturation.data(), saturation.size() * sizeof(double));
            return;
        }
        if (it->first == "selectedoutput")
        {
            std::vector<double> so;
            IRM_RESULT status = this->GetSelectedOutput(so);
            memcpy(dest, so.data(), so.size() * sizeof(double));
            return;
        }
        if (it->first == "selectedoutputcolumncount")
        {
            int count = this->GetSelectedOutputColumnCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "selectedoutputcount")
        {
            int count = this->GetSelectedOutputCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "selectedoutputheadings")
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
        if (it->first == "selectedoutputon")
        {
            int tf = (int)this->GetSelectedOutputOn();
            memcpy(dest, &tf, sizeof(int));
            return;
        }
        if (it->first == "selectedoutputon")
        {
            int tf = (int)this->GetSelectedOutputOn();
            memcpy(dest, &tf, sizeof(bool));
            return;
        }
        if (it->first == "selectedoutputrowcount")
        {
            int count = this->GetSelectedOutputRowCount();
            memcpy(dest, &count, sizeof(int));
            return;
        }
        if (it->first == "solutionvolume")
        {
            const std::vector<double>& vol = this->GetSolutionVolume();
            memcpy(dest, vol.data(), vol.size() * sizeof(double));
            return;
        }
        if (it->first == "temperature")
        {
            const std::vector<double>& temperature = this->GetTemperature();
            memcpy(dest, temperature.data(), temperature.size() * sizeof(double));
            return;
        }
        if (it->first == "time")
        {
            double time = this->GetTime();
            memcpy(dest, &time, sizeof(double));
            return;
    }
        if (it->first == "timestep")
        {
            double timestep = this->GetTimeStep();
            memcpy(dest, &timestep, sizeof(double));
            return;
        }
    }
    ErrorMessage("Item not found");
    throw PhreeqcRMStop();
}
int PhreeqcRM::BMI_GetVarNbytes(std::string name)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "componentcount")
        {
            return this->GetComponentCount();
        }
        if (it->first == "components")
        {
            int string_size = this->BMI_GetVarItemsize("components");
            int dim = this->GetComponentCount();
            return string_size * dim * (int)sizeof(char);
        }
        if (it->first == "concentrations")
        {
            return (int)sizeof(double) * this->GetGridCellCount() * this->GetComponentCount();
        }
        if (it->first == "currentselectedoutputusernumber")
        {
            return (int)sizeof(int);
        }
        if (it->first == "density")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }

        if (it->first == "fileprefix")
        {
            return (int)this->GetFilePrefix().size();
        }
        if (it->first == "errorstring")
        {
            return (int)this->GetErrorString().size();
        }
        if (it->first == "gfw")
        {
            return (int)sizeof(double) * this->GetComponentCount();
        }
        if (it->first == "gridcellcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "inputvarnames")
        {
            int string_length = this->BMI_GetVarItemsize("inputvarnames");
            std::vector < std::string > varnames = this->BMI_GetInputVarNames();
            return string_length * (int)(varnames.size() * sizeof(char));
        }
        if (it->first == "nthselectedoutput")
        {
            return (int)sizeof(int);
        }
        if (it->first == "outputvarnames")
        {
            int string_length = this->BMI_GetVarItemsize("outputvarnames");
            std::vector < std::string > varnames = this->BMI_GetOutputVarNames();
            return string_length * (int)(varnames.size() * sizeof(char));
        }
        if (it->first == "porosity")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "pressure")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "saturation")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "selectedoutput")
        {
            return (int)sizeof(double) * this->GetSelectedOutputColumnCount() * this->GetSelectedOutputRowCount();
        }
        if (it->first == "selectedoutputcolumncount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputheadings")
        {
            int string_size = this->BMI_GetVarItemsize("selectedoutputheadings");
            return string_size * this->GetSelectedOutputColumnCount() * (int)sizeof(char);
        }
        if (it->first == "selectedoutputon")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputrowcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "solutionvolume")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "temperature")
        {
            return (int)sizeof(double) * this->GetGridCellCount();
        }
        if (it->first == "time")
        {
            return (int)sizeof(double);
        }
        if (it->first == "timestep")
        {
            return (int)sizeof(double);
        }
    }
	//throw LetItThrow("Item not found");
    ErrorMessage("Item not found");
    throw PhreeqcRMStop();
}
std::string PhreeqcRM::BMI_GetVarUnits(std::string name)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "concentrations")
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
				//throw LetItThrow("Unknown units for GetUnitsSolution in GetVarUnits");
                ErrorMessage("Unknown units for GetUnitsSolution in GetVarUnits");
                throw PhreeqcRMStop();
				break;
			}
        }
        return it->second.GetUnits();
    }
    //throw LetItThrow("Item not found");
    ErrorMessage("Item not found");
    throw PhreeqcRMStop();
}
int PhreeqcRM::BMI_GetVarItemsize(std::string name)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        if (it->first == "componentcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "components")
        {
            size_t max = 0;
            for (size_t i = 0; i < this->GetComponents().size(); i++)
            {
                if (this->GetComponents()[i].size() > max) max = this->GetComponents()[i].size();
            }
            return (int)max;
        }
        if (it->first == "concentrations")
        {
            return (int)sizeof(double);
        }
        if (it->first == "currentselectedoutputusernumber")
        {
            return (int)sizeof(int);
        }
        if (it->first == "density")
		{
            return (int)sizeof(double);
		}
        if (it->first == "fileprefix")
        {
            return (int)this->GetFilePrefix().size();
        }
        if (it->first == "errorstring")
        {
            return (int)this->GetErrorString().size();
        }
        if (it->first == "gfw")
		{
            return (int)sizeof(double);
		}
        if (it->first == "gridcellcount")
		{
            return (int)sizeof(int);
		}
        if (it->first == "inputvarnames")
        {
            std::vector < std::string > varnames = this->BMI_GetInputVarNames();
            size_t max = 0;
            for (size_t i = 0; i < varnames.size(); i++)
            {
                if (varnames[i].size() > max) max = varnames[i].size();
            }
            return (int)max;
        }
		if (it->first == "nthselectedoutput")
		{
            return (int)sizeof(int);
		}
        if (it->first == "outputvarnames")
        {
            std::vector < std::string > varnames = this->BMI_GetOutputVarNames();
            size_t max = 0;
            for (size_t i = 0; i < varnames.size(); i++)
            {
                if (varnames[i].size() > max) max = varnames[i].size();
            }
            return (int)max;
        }
		if (it->first == "porosity")
		{
            return (int)sizeof(double);
		}
        if (it->first == "pressure")
        {
            return (int)sizeof(double);
        }
        if (it->first == "saturation")
        {
            return (int)sizeof(double);
        }
        if (it->first == "selectedoutput")
        {
            return (int)sizeof(double);
        }
        if (it->first == "selectedoutputcolumncount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputheadings")
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
        if (it->first == "selectedoutputon")
        {
            return (int)sizeof(int);
        }
        if (it->first == "selectedoutputrowcount")
        {
            return (int)sizeof(int);
        }
        if (it->first == "solutionvolume")
        {
            return (int)sizeof(double);
        }
        if (it->first == "temperature")
        {
            return (int)sizeof(double);
        }
        if (it->first == "time")
        {
            return (int)sizeof(double);
        }
        if (it->first == "timestep")
        {
            return (int)sizeof(double);
        }
    }
    else
    {
        //throw LetItThrow("Item not found");
        ErrorMessage("Item not found");
        throw PhreeqcRMStop();
    }
    return -1;
}
std::string PhreeqcRM::BMI_GetVarType(std::string name)
{
    std::string name_lc = name;
    std::transform(name_lc.begin(), name_lc.end(), name_lc.begin(), tolower);
    std::map < std::string, BMI_Var >::iterator it = this->bmi_var_map.find(name_lc);
    if (it != bmi_var_map.end())
    {
        return it->second.GetType();
    }
    else
    {
        //throw LetItThrow("Item not found");
        ErrorMessage("Item not found");
        throw PhreeqcRMStop();
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
        bmi_var_map["components"] = BMI_Var("Components", "character,1d", "names", false, true);
        bmi_var_map["componentcount"] = BMI_Var("ComponentCount", "integer", "names", false, true);
        bmi_var_map["concentrations"] = BMI_Var("Concentrations", "double,2d", "mol L-1", true, true);
        bmi_var_map["density"] = BMI_Var("Density", "double,1d", "kg L-1", true, true);
        //var_map["EndCell"] = Var_BMI("EndCell", "int", "cell numbers", sizeof(int));
        //var_map["EquilibriumPhasesNames"] = Var_BMI("EquilibriumPhasesNames", "string", "names", sizeof(char));
        bmi_var_map["errorstring"] = BMI_Var("ErrorString", "character", "error", false, true);
        //var_map["ExchangeNames"] = Var_BMI("ExchangeNames", "string", "names", sizeof(char));
        //var_map["ExchangeSpeciesNames"] = Var_BMI("ExchangeSpeciesNames", "string", "names", sizeof(char));
        bmi_var_map["fileprefix"] = BMI_Var("FilePrefix", "character", "name", true, true);
        //!var_map["GasComponentNames"] = Var_BMI("GasComponentsNames", "string", "names", sizeof(char));
        //!var_map["GasCompMoles"] = Var_BMI("GasCompMoles", "double", "mol", sizeof(double));
        //!var_map["GasCompPressures"] = Var_BMI("GasCompPressures", "double", "atm", sizeof(double));
        //!var_map["GasCompPhi"] = Var_BMI("GasCompPhi", "double", "atm-1", sizeof(double));
        //!var_map["GasPhaseVolume"] = Var_BMI("GasPhaseVolume", "double", "L", sizeof(double));
        bmi_var_map["gfw"] = BMI_Var("Gfw", "double,1d", "g mol-1", false, true);
        bmi_var_map["gridcellcount"] = BMI_Var("GridCellCount", "integer", "count", false, true);
        //var_map["KineticReactions"] = Var_BMI("KineticReactions", "string", "names", sizeof(char));
        //var_map["MpiMyself"] = Var_BMI("MpiMyself", "int", "id", sizeof(int));
        //var_map["MpiTasks"] = Var_BMI("MpiTasks", "int", "count", sizeof(int));
        bmi_var_map["nthselectedoutput"] = BMI_Var("NthSelectedOutput", "integer", "id", true, false);
        bmi_var_map["saturation"] = BMI_Var("Saturation", "double,1d", "unitless", true, true);
        bmi_var_map["selectedoutput"] = BMI_Var("SelectedOutput", "double,2d", "user", false, true);
        bmi_var_map["selectedoutputcolumncount"] = BMI_Var("SelectedOutputColumnCount", "integer", "count", false, true);
        bmi_var_map["selectedoutputcount"] = BMI_Var("SelectedOutputCount", "integer", "count", false, true);
        bmi_var_map["selectedoutputheadings"] = BMI_Var("SelectedOutputHeadings", "character,1d", "names", false, true);
        bmi_var_map["selectedoutputrowcount"] = BMI_Var("SelectedOutputRowCount", "integer", "count", false, true);
        //var_map["SINames"] = Var_BMI("SINames", "string", "names", 0);
        //var_map["SolidSolutionComponentsNames"] = Var_BMI("SolidSolutionComponentsNames", "string", "names", sizeof(char));
        //var_map["SolidSolutionNames"] = Var_BMI("SolidSolutionNames", "string", "names", sizeof(char));
        bmi_var_map["solutionvolume"] = BMI_Var("SolutionVolume", "double,1d", "L", false, true);
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
        bmi_var_map["time"] = BMI_Var("Time", "double", "s", true, true);
        //var_map["TimeConversion"] = Var_BMI("TimeConversion", "double", "unitless", sizeof(double));
        bmi_var_map["timestep"] = BMI_Var("TimeStep", "double", "s", true, true);
        //var_map["MpiWorker"] = Var_BMI("MpiWorker", "int", "id", sizeof(int));
        //var_map["ComponentH2O"] = Var_BMI("ComponentH2O", "int", "flag", sizeof(int));
        bmi_var_map["currentselectedoutputusernumber"] = BMI_Var("CurrentSelectedOutputUserNumber", "integer", "id", false, true);
        //var_map["DumpFileName"] = Var_BMI("DumpFileName", "string", "name", sizeof(char));
        //var_map["ErrorHandlerMode"] = Var_BMI("ErrorHandlerMode", "int", "flag", sizeof(int));
        //var_map["ErrorOn"] = Var_BMI("ErrorOn", "int", "flag", sizeof(int));
        //var_map["PartitionUZSolids"] = Var_BMI("PartitionUZSolids", "int", "flag", sizeof(int));
        bmi_var_map["porosity"] = BMI_Var("Porosity", "double,1d", "unitless", true, true);
        //var_map["PartitionUZSolids"] = Var_BMI("PartitionUZSolids", "int", "flag", sizeof(int));
        bmi_var_map["pressure"] = BMI_Var("Pressure", "double,1d", "atm", true, true);
        //var_map["PrintChemistryMask"] = Var_BMI("PrintChemistryMask", "int", "flags", sizeof(int));
        //var_map["PrintChemistryOn"] = Var_BMI("PrintChemistryOn", "int", "flag", sizeof(int));
        //var_map["RebalanceByCell"] = Var_BMI("RebalanceByCell", "int", "flag", sizeof(int));
        //var_map["RebalanceFraction"] = Var_BMI("RebalanceFraction", "double", "unitless", sizeof(double));
        //var_map["RepresentativeVolume"] = Var_BMI("RepresentativeVolume", "double", "L", sizeof(double));
        //var_map["RebalanceByCell"] = Var_BMI("RebalanceByCell", "int", "flag", sizeof(int));
        //var_map["ScreenOn"] = Var_BMI("ScreenOn", "int", "flag", sizeof(int));
        bmi_var_map["selectedoutputon"] = BMI_Var("SelectedOutputOn", "logical", "flag", true, true);
        bmi_var_map["temperature"] = BMI_Var("Temperature", "double,1d", "C", true, true);
        //var_map["UnitsExchange"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsGasPhase"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsKinetics"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsPPassemblage"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSolution"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSSassemblage"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UnitsSurface"] = Var_BMI("UnitsExchange", "int", "flag", sizeof(int));
        //var_map["UseSolutionDensityVolume"] = Var_BMI("UseSolutionDensityVolume", "int", "flag", sizeof(int));
        bmi_var_map["inputvarnames"] = BMI_Var("InputVarNames", "character,1d", "string", false, true);
        bmi_var_map["outputvarnames"] = BMI_Var("OutputVarNames", "character,1d", "string", false, true);
    }
    this->bmi_input_vars.clear();
    this->bmi_output_vars.clear();
    std::map<std::string, class BMI_Var>::iterator it;
    for (it = bmi_var_map.begin(); it != bmi_var_map.end(); it++)
    {
        if (it->second.GetSet())
        {
            bmi_input_vars.push_back(it->second.GetName());
        }
        if (it->second.GetGet())
        {
            bmi_output_vars.push_back(it->second.GetName());
        }
    }
}




