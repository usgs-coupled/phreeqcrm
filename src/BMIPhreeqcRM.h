#if !defined(BMIPHREEQCRM_H_INCLUDED)
#define BMIPHREEQCRM_H_INCLUDED
#include <map>
#include "PhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
class BMI_Variant
{
public:
    BMI_Var                  bmi_var;
    bool                     b_var;
    int                      i_var;
    double                   d_var;
    std::string              string_var;
    std::vector<int>         IntVector;
    std::vector<double>      DoubleVector;
    std::vector<std::string> StringVector;
    bool                     NotImplemented;
    bool GetGet() { return this->bmi_var.GetGet(); }
    void SetGet(bool tf) { this->bmi_var.SetGet(tf); }
    std::string GetName() { return this->bmi_var.GetName(); }
    void SetName(std::string s) { this->bmi_var.SetName(s); }
    bool GetSet() { return this->bmi_var.GetSet(); }
    void SetSet(bool tf) { this->bmi_var.SetSet(tf); }
    std::string GetUnits() { return this->bmi_var.GetUnits(); }
    void SetUnits(std::string s) { this->bmi_var.SetUnits(s); }
    int GetNbytes() { return (int)this->bmi_var.GetNbytes(); }
    void SetNbytes(int n) { this->bmi_var.SetNbytes(n); }
    int GetItemsize() { return this->bmi_var.GetItemsize(); }
    void SetItemsize(int n) { this->bmi_var.SetItemsize(n); }
    std::string GetCType() { return this->bmi_var.GetCType(); }
    void SetCType(std::string s) { this->bmi_var.SetCType(s); }
    std::string GetFType() { return this->bmi_var.GetFType(); }
    void SetFType(std::string s) { this->bmi_var.SetFType(s); }
    std::string GetPType() { return this->bmi_var.GetPType(); }
    void SetPType(std::string s) { this->bmi_var.SetPType(s); }
    void SetTypes(std::string c, std::string f, std::string p)
    {
        this->bmi_var.SetTypes(c, f, p);
    }

    void Clear();
};
class IRM_DLL_EXPORT BMIPhreeqcRM : /*public bmi::Bmi,*/ public PhreeqcRM
{
public:
    static void             CleanupBMIModuleInstances(void);
    static int              CreateBMIModule(int nxyz, MP_TYPE nthreads);
    static IRM_RESULT       DestroyBMIModule(int n);
    static BMIPhreeqcRM*    GetInstance(int n);

    BMIPhreeqcRM(int nxyz, int nthreads);
    enum class BMI_TASKS {
        Info,
        GetVar,
        SetVar,
        no_op
    };
    // Model control functions.
    void Initialize(std::string config_file);
    void Update();
    void UpdateUntil(double time);
    void Finalize();

    // Model information functions.
    std::string GetComponentName() {return "BMI PhreeqcRM";};
    int GetInputItemCount();
    int GetOutputItemCount();
    std::vector<std::string> GetInputVarNames();
    std::vector<std::string> GetOutputVarNames();
    // Variable information functions
    // Not applicable
    //virtual int GetVarGrid(const std::string name) = 0;             
    std::string GetVarType(const std::string name);
    std::string GetVarUnits(const std::string name);
    int GetVarItemsize(const std::string name);
    int GetVarNbytes(const std::string name);
    // Not implemented
    //virtual std::string GetVarLocation(const std::string name) = 0; 

    // Time functions
    double GetCurrentTime();
    double GetStartTime();
    double GetEndTime();
    std::string GetTimeUnits() { return "seconds"; };
    //double GetTimeStep();

    // Variable getters
    void GetValue(const std::string name, void* dest);
    void GetValue(const std::string name, bool& dest);
    void GetValue(const std::string name, double& dest);
    void GetValue(const std::string name, int& dest);
    void GetValue(const std::string name, std::string& dest);
    void GetValue(const std::string name, std::vector<double>& dest);
    void GetValue(const std::string name, std::vector<int>& dest);
    void GetValue(const std::string name, std::vector<std::string>& dest);
    // Not implemented
    //virtual void* GetValuePtr(std::string name) = 0;                  
    // Not implemented
    //virtual void GetValueAtIndices(std::string name, void* dest, int* inds, int count) = 0;
    
    // Variable setters
    void SetValue(std::string name, void* src);
    void SetValue(std::string name, bool src);
    void SetValue(std::string name, char* src);
    void SetValue(std::string name, double src);
    void SetValue(std::string name, int src);
    void SetValue(std::string name, std::string src);
    void SetValue(std::string name, std::vector<double> src);
    void SetValue(std::string name, std::vector<int>  src);
    void SetValue(std::string name, std::vector<std::string>  src);
    // Not implemented
    //virtual void SetValueAtIndices(std::string name, int* inds, int count, void* src) = 0;

    // Grid information functions
    // Not implemented
    // PhreeqcRM has no grid 
    //virtual int GetGridRank(const int grid) = 0;
    //virtual int GetGridSize(const int grid) = 0;
    //virtual std::string GetGridType(const int grid) = 0;

    //virtual void GetGridShape(const int grid, int* shape) = 0;
    //virtual void GetGridSpacing(const int grid, double* spacing) = 0;
    //virtual void GetGridOrigin(const int grid, double* origin) = 0;

    //virtual void GetGridX(const int grid, double* x) = 0;
    //virtual void GetGridY(const int grid, double* y) = 0;
    //virtual void GetGridZ(const int grid, double* z) = 0;

    //virtual int GetGridNodeCount(const int grid) = 0;
    //virtual int GetGridEdgeCount(const int grid) = 0;
    //virtual int GetGridFaceCount(const int grid) = 0;

    //virtual void GetGridEdgeNodes(const int grid, int* edge_nodes) = 0;
    //virtual void GetGridFaceEdges(const int grid, int* face_edges) = 0;
    //virtual void GetGridFaceNodes(const int grid, int* face_nodes) = 0;
    //virtual void GetGridNodesPerFace(const int grid, int* nodes_per_face) = 0;
    // data
    BMI_TASKS task;
    BMI_Variant bmi_variant;
    std::string language;
    typedef void (*VarFunction)(BMIPhreeqcRM* brm_ptr); // function pointer type
    typedef std::map<std::string, VarFunction> VarFunction_map;
    VarFunction_map varfn_map;
    VarFunction GetFn(const std::string name);

private:
    //friend class RM_interface;
    static std::map<size_t, BMIPhreeqcRM*> Instances;
    static size_t InstancesIndex;
};
#endif //BMIPHREEQCRM_H_INCLUDED
