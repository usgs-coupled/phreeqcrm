#include <map>
#include "PhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
class BMI_Variant
{
public:
    BMI_Var                  bmi_var;
    size_t                   Nbytes;
    size_t                   Itemsize;
    bool                     b_var;
    int                      i_var;
    double                   d_var;
    std::string              string_var;
    std::vector<int>         IntVector;
    std::vector<double>      DoubleVector;
    std::vector<std::string> StringVector;
    bool                     NotImplemented;
    std::string GetType() { return this->bmi_var.GetType(); }
    std::string GetUnits() { return this->bmi_var.GetUnits(); }
    void Clear();
};
class BMIPhreeqcRM : public bmi::Bmi, public PhreeqcRM
{
public:
    BMIPhreeqcRM(int nxyz, int nthreads);
    enum class BMI_TASKS {
        count_output,
        count_input,
        Nbytes,
        Itemsize,
        Type,
        Units,
        GetVar,
        SetVar,
        no_op
    };
    // Model information functions.
    std::string GetComponentName() {return "BMI PhreeqcRM";};
    int GetInputItemCount();
    std::string BMIPhreeqcRM::GetVarType(const std::string name);
    std::string BMIPhreeqcRM::GetVarUnits(const std::string name);
    // data
    BMI_TASKS task;
    BMI_Variant bmi_variant;
    typedef void (*VarFunction)(BMIPhreeqcRM& bmi_rm_ref); // function pointer type
    typedef std::map<std::string, VarFunction> VarFunction_map;
    VarFunction_map varfn_map;
};

