/*! @file BMIPhreeqcRM.h
*	@brief C++ header file for BMIPhreeqcRM
*/
#if !defined(BMIPHREEQCRM_H_INCLUDED)
#define BMIPHREEQCRM_H_INCLUDED
#include <map>

#if defined(WITH_PYBIND11)
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

class NotIntialized : public std::runtime_error {
public:
    NotIntialized() : std::runtime_error("must call initialize first") { };
};
#endif

#include "PhreeqcRM.h"
#include "BMIVariant.h"
#include "bmi.hxx"
#include "VarManager.h"
/**
 * @class NotImplemented
 *
 * @brief Throws an exception for Basic Model Interface methods that are
 * not implemented in BMIPhreeqcRM
 */
class NotImplemented : public std::logic_error {
public:
    NotImplemented() : std::logic_error("Not Implemented") { };
};
/**
 * @class BMIPhreeqcRM
 *
 * @brief Basic Model Interface implementation of the 
 * geochemical reaction module PhreeqcRM
 */

class IRM_DLL_EXPORT BMIPhreeqcRM : public bmi::Bmi, public PhreeqcRM
{
public:
    static void             CleanupBMIModuleInstances(void);
    static int              CreateBMIModule();
    static int              CreateBMIModule(int nxyz, MP_TYPE nthreads);
    static IRM_RESULT       DestroyBMIModule(int n);
    static BMIPhreeqcRM*    GetInstance(int n);
    /**
    Default constructor for the BMIPhreeqcRM subclass of PhreeqcRM.
    Definition of the number of cells and threads (or MPI communicator) is deferred.
    @ref Initialize must be called to initialize the BMIPhreeqcRM instance.
    */
    BMIPhreeqcRM();
    /**
    Constructor for the BMIPhreeqcRM subclass of PhreeqcRM. A BMIPhreeqcRM
    instance has the BMI methods plus all of the PhreeqcRM methods. The
    constructor requires two arguments: the number of cells in the user's
    model, and either (a) the number of threads for OpenMP parallelization, or
    (b) an MPI communicator.
    @param ngrid Number of cells in the user's model. 
    @param nthreads Number of threads for parallelization with OpenMP or
    an MPI communicator if PhreeqcRM is compiled with MPI. With OpenMP,
    a value of zero causes the program to set nthreads to the number 
    of logical processors of the computer.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    int nthreads = 0;
    std::string yaml_file = "myfile.yaml";
    int nxyz = GetGridCellCountYAML(yaml_file);
    BMIPhreeqcRM brm(nxyz, nthreads);
    brm.Initialize(yaml_file);
    </PRE>
    </CODE>
    @endhtmlonly
    */
    BMIPhreeqcRM(int ngrid, MP_TYPE nthreads);

    ~BMIPhreeqcRM() override;

    // Model control functions.
    /**
    @a Initialize must be called to initialize a BMIPhreeqcRM instance. 
    A YAML file is normally used for initialization; however, an empty string can be used for the 
    file name when initializing without use of a YAML file.
    <p>
    The YAML file contains a YAML map of PhreeqcRM
    methods and data corresponding to each PhreeqcRM method. 
    For example,
    </p>
    @htmlonly
    <CODE>
    <PRE>
        - key: LoadDatabase
          database: phreeqc.dat
        - key: RunFile
          workers: true
          initial_phreeqc: true
          utility: true
          chemistry_name: advect.pqi
    </PRE>
    </CODE>
    @endhtmlonly

    @a Initialize will read the YAML file and execute the specified methods with
    the specified arguments. Using YAML terminology, the argument(s) for a method
    may be a scalar, a sequence, or a map, depending if the argument is
    a single item, a single vector, or there are multiple arguments. In the case
    of a map, the name associated with each argument (for example "chemistry_name"
    above) is arbitrary. The names of the map keys for map arguments are not used
    in parsing the YAML file; only the order of the arguments is important.

    The PhreeqcRM methods that can be specified in a YAML file include:
    @htmlonly
    <CODE>
    <PRE>
    CloseFiles();
    CreateMapping(std::vector< int >& grid2chem);
    DumpModule();
    FindComponents();
    InitialEquilibriumPhases2Module(std::vector< int > equilibrium_phases);
    InitialExchanges2Module(std::vector< int > exchanges);
    InitialGasPhases2Module(std::vector< int > gas_phases);
    InitialKineticss2Module(std::vector< int > kinetics);
    InitialSolidSolutions2Module(std::vector< int > solid_solutions);
    InitialSolutions2Module(std::vector< int > solutions);
    InitialSurfaces2Module(std::vector< int > surfaces);
    InitialPhreeqc2Module(std::vector< int > initial_conditions1);
    InitialPhreeqc2Module(std::vector< int > initial_conditions1,
        std::vector< int > initial_conditions2, std::vector< double > fraction1);
    InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
    LoadDatabase(std::string database);
    OpenFiles();
    OutputMessage(std::string str);
    RunCells();
    RunFile(bool workers, bool initial_phreeqc,
         bool utility, std::string chemistry_name);
    RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    ScreenMessage(std::string str);
    SetComponentH2O(bool tf);
    SetConcentrations(std::vector< double > c);
    SetCurrentSelectedOutputUserNumber(int n_user);
    SetDensityUser(std::vector< double > density);
    SetDumpFileName(std::string dump_name);
    SetErrorHandlerMode(int mode);
    SetErrorOn(bool tf);
    SetFilePrefix(std::string prefix);
    SetGasCompMoles(std::vector< double > gas_moles);
    SetGasPhaseVolume(std::vector< double > gas_volume);
    SetPartitionUZSolids(bool tf);
    SetPorosity(std::vector< double > por);
    SetPressure(std::vector< double > p);
    SetPrintChemistryMask(std::vector< int > cell_mask);
    SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
    SetRebalanceByCell(bool tf);
    SetRebalanceFraction(double f);
    SetRepresentativeVolume(std::vector< double > rv);
    SetSaturationUser(std::vector< double > sat);
    SetScreenOn(bool tf);
    SetSelectedOutputOn(bool tf);
    SetSpeciesSaveOn(bool save_on);
    SetTemperature(std::vector< double > t);
    SetTime(double time);
    SetTimeConversion(double conv_factor);
    SetTimeStep(double time_step);
    SetUnitsExchange(int option);
    SetUnitsGasPhase(int option);
    SetUnitsKinetics(int option);
    SetUnitsPPassemblage(int option);
    SetUnitsSolution(int option);
    SetUnitsSSassemblage(int option);
    SetUnitsSurface(int option);
    SpeciesConcentrations2Module(std::vector< double > species_conc);
    StateSave(int istate);
    StateApply(int istate);
    StateDelete(int istate);
    UseSolutionDensityVolume(bool tf);
    WarningMessage(std::string warnstr);
    </PRE>
    </CODE>
    @endhtmlonly
    @see
    @ref Update.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    int nthreads = 0;
    std::string yaml_file = "myfile.yaml";
    int nxyz = GetGridCellCountYAML(yaml_file);
    BMIPhreeqcRM brm(nxyz, nthreads);
    brm.Initialize(yaml_file);
    int ncomps;
    brm.GetValue("ComponentCount", &ncomps);
    int ngrid;
    brm.GetValue("GridCellCount", &ngrid);
    std::vector< double > c;
    brm.GetValue("Concentrations", c);
    brm.SetValue("TimeStep", 86400);
    for(double time = 0; time < 864000; time+=86400)
    {
        // Take a transport time step here and update the vector c.
        brm.SetValue("Time", time);
        brm.SetValue("Concentrations", c);
        brm.Update();
        brm.GetValue("Concentrations", c);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
     */
    void Initialize(std::string config_file = "") override;
    /**
    @a Update runs PhreeqcRM for one time step. PhreeqcRM will equilibrate 
    the solutions with all equilibrium
    reactants (EQUILIBRIUM_PHASES, EXCHANGE, GAS_PHASE, SOLID_SOLUTIONS, and SURFACE)
    and integrate KINETICS reactions for the specified time step
    (@ref SetValue "TimeStep" or @ref SetTimeStep).
    @see
    @ref Initialize,
    @ref UpdateUntil.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    BMIPhreeqcRM brm(nxyz, nthreads);
    brm.Initialize("myfile.yaml");
    int ncomps;
    brm.GetValue("ComponentCount", &ncomps);
    int ngrid;
    brm.GetValue("GridCellCount", &ngrid);
    std::vector< double > c;
    brm.GetValue("Concentrations", c);
    brm.SetValue("TimeStep", 86400);
    for(double time = 0; time < 864000; time+=86400)
    {
        // Take a transport time step here and update the vector c.
        brm.SetValue("Time", time);
        brm.SetValue("Concentrations", c);
        brm.Update();
        brm.GetValue("Concentrations", c);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
     */
    void Update() override;

    /**
    @a UpdateUntil is the same as @ref Update, except the time step is calculated
    from the argument @a end_time. The time step is calculated to be @a end_time minus
    the current time (@ref GetCurrentTime).
    @param end_time Time at the end of the time step.
    @see
    @ref Initialize,
    @ref Update.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    BMIPhreeqcRM brm(nxyz, nthreads);
    brm.Initialize("myfile.yaml");
    int ncomps;
    brm.GetValue("ComponentCount", &ncomps);
    int ngrid;
    brm.GetValue("GridCellCount", &ngrid);
    std::vector< double > c;
    brm.GetValue("Concentrations", c);
    for(double time = 0; time < 864000; time+=86400)
    {
        // Take a transport time step here and update the vector c.
        brm.SetValue("Time", time);
        brm.SetValue("Concentrations", c);
        brm.UpdateUntil(time + 86400.0);
        brm.GetValue("Concentrations", c);
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
     */
    void UpdateUntil(double end_time) override;
    /**
    @a Finalize closes any files open in the BMIPhreeqcRM instance.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    brm.Finalize();
     </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
    */
    void Finalize() override;

    // Model information functions.

    /**
    @a GetComponentName returns the component name--"BMI PhreeqcRM".
    BMI PhreeqcRM is a partial interface to PhreeqcRM, and provides
    the methods to implement chemical reactions in a
    multicomponent transport model. All of the native PhreeqcRM methods
    (non BMI methods) are also available, which provides a complete
    interface to PhreeqcRM.
    @retval The string "BMI PhreeqcRM".
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::string comp_name = brm.GetComponentName();
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
    */
#if defined(_MSC_VER)
#if defined(_WIN64)
    std::string GetComponentName() override {
        char buffer[400];
        sprintf(buffer, "BMI PhreeqcRM [MSC v.%d 64 bit (AMD64)]", _MSC_VER);
        return buffer;
    };
#elif defined(_WIN32)
    std::string GetComponentName() override {
        char buffer[400];
        sprintf(buffer, "BMI PhreeqcRM [MSC v.%d 32 bit (Intel)]", _MSC_VER);
        return buffer;
    };
#else
    std::string GetComponentName() override {
        char buffer[400];
        sprintf(buffer, "BMI PhreeqcRM [MSC v.%d Unknown]", _MSC_VER);
        return buffer;
    };
#endif
#else
    std::string GetComponentName() override { return "BMI PhreeqcRM"; };
#endif

    /**
    @a GetInputVarNames returns the count of input variables that can
    be set with @ref SetValue.
    @retval  Count of input variables that can be set with @ref SetValue.
    @see
    @ref GetInputVarNames,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits,
    @ref SetValue.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    int count = brm.GetInputItemCount();
    std::vector< std::string > InputVarNames = brm.GetInputVarNames();
    oss << "SetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << InputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(InputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(InputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(InputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(InputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(InputVarNames[i]) /
               brm.GetVarItemsize(InputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    int GetInputItemCount() override;

    /**
    @a GetOutputItemCount returns the count of output variables that can
    be retrieved with @ref GetValue.
    @retval  Count of output variables that can be retrieved with @ref GetValue.

    @see
    @ref GetOutputVarNames,
    @ref GetValue,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    int count = brm.GetOutputItemCount();
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
               brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    int GetOutputItemCount() override;
    /**
    @a GetPointableItemCount returns the count of variables for which
    pointers can be obtained with @ref GetValuePtr. The pointers point to
    current copies of the variables. Setting a value with one of the pointers
    will have no effect on the simulation, but will corrupt the copy of the variable.
    @retval  Count of pointers to variables that can be accessed with
    @ref GetValuePtr.

    @see
    @ref GetPointableVarNames,
    @ref GetValue,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    int count = brm.GetPointableItemCount();
    std::vector< std::string > PointableVarNames = brm.GetPointableVarNames();
    oss << "GetValuePtr variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << PointableVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(PointableVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(PointableVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(PointableVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(PointableVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(PointableVarNames[i]) /
               brm.GetVarItemsize(PointableVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    int GetPointableItemCount();

    /**
    @a GetInputVarNames returns a list of the variable names that can be
    set with @ref SetValue.
    @retval  A std::vector of the names of variables that can be set with
    @ref SetValue.

    @see
    @ref GetInputItemCount,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits,
    @ref SetValue.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > InputVarNames = brm.GetInputVarNames();
    int count = brm.GetInputItemCount();
    oss << "SetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << InputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(InputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(InputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(InputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(InputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(InputVarNames[i]) /
               brm.GetVarItemsize(InputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::vector<std::string> GetInputVarNames() override;

    /**
    @a GetOutputVarNames returns a list of the variable names that can
    be retrieved with @ref GetValue.
    @retval  A list of the names of variable that can be retrieved with
    @ref GetValue.

    @see
    @ref GetOutputItemCount,
    @ref GetValue,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    int count = brm.GetOutputItemCount();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
               brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::vector<std::string> GetOutputVarNames() override;

    /**
    @a GetPointableVarNames returns a list of the names of variables
    for which pointers can be retrieved with @ref GetValuePtr.
    @retval  A list of the names of variables for which pointers can
    be retieved with @ref GetValuePtr.

    @see
    @ref GetPointableItemCount,
    @ref GetValue,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > PointableVarNames = brm.GetPointableVarNames();
    int count = brm.GetPointableItemCount();
    oss << "GetValuePtr variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << PointableVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(PointableVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(PointableVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(PointableVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(PointableVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(PointableVarNames[i]) /
               brm.GetVarItemsize(PointableVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::vector<std::string> GetPointableVarNames();

    /**
    @todo
    */
    std::vector<std::string> GetReadOnlyVarNames();

    // Variable information functions
    /**
    @a GetVarGrid returns a value of 1, indicating points.
    BMIPhreeqcRM does not have a grid of its own. The cells
    of BMIPhreeqcRM are associated with the user's model grid,
    and all spatial characterists are assigned by the user's
    model.
    @param name Varaiable name. (Return value is the same regardless of @a name.)
    @retval 1 BMIPhreeqcRM cells derive meaning from the user's
    model.
    */
    int GetVarGrid(const std::string name) override { return 0; }

    /**
    @a GetVarType retrieves the type of a variable
    that can be set with @ref SetValue, retrieved with @ref GetValue,
    or pointed to by @ref GetValuePtr.
    Types are "int", "double", "std::string", or "std::vector<std::string>".

    @param name Name of the variable to retrieve type.
    @retval Character string of variable type.

    @see
    @ref GetVarNbytes,
    @ref GetVarItemsize,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    int count = brm.GetOutputItemCount();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
               brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::string GetVarType(const std::string name) override;
    /**
    @a GetVarUnits retrieves the units of a variable
    that can be set with @ref SetValue, retrieved with @ref GetValue,
    or pointed to by @ref GetValuePtr.
    @param name Name of the variable to retrieve units.
    @retval Character string of units for variable.

    @see
    @ref GetVarNbytes,
    @ref GetVarItemsize,
    @ref GetVarType.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    int count = brm.GetOutputItemCount();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
               brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::string GetVarUnits(const std::string name) override;

    /**
    @a GetVarItemsize retrieves size of an individual item that
    can be set or retrived. Sizes may be sizeof(int), sizeof(double),
    or a character length for string variables.
    @param name Name of the variable to retrieve size.
    @retval Size of one element of the variable.

    @see
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    int count = brm.GetOutputItemCount();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
               brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    int GetVarItemsize(const std::string name) override;

    /**
    @a GetVarNbytes retrieves the total number of bytes that are
    set for a variable with @ref SetValue, retrieved for a variable with
    @ref GetValue, or pointed to by @ref GetValuePtr.
    @param name Name of the variable to retrieve total bytes.
    @retval Total number of bytes set, retrieved, or pointed to for variable.
    @see
    @ref GetVarItemsize,
    @ref GetVarType,
    @ref GetVarUnits.

    @par C++ Example:
    @htmlonly
        <CODE>
        <PRE>
    std::vector< std::string > OutputVarNames = brm.GetOutputVarNames();
    int count = brm.GetOutputItemCount();
    oss << "GetValue variables:\n";
    for (size_t i = 0; i < count; i++)
    {
        oss << "  " << i << "  " << OutputVarNames[i] << "\n";
        oss << "     Type:        " << brm.GetVarType(OutputVarNames[i]) << "\n";
        oss << "     Units:       " << brm.GetVarUnits(OutputVarNames[i]) << "\n";
        oss << "     Total bytes: " << brm.GetVarNbytes(OutputVarNames[i]) << "\n";
        oss << "     Item bytes:  " << brm.GetVarItemsize(OutputVarNames[i]) << "\n";
        oss << "     Dim:         " << brm.GetVarNbytes(OutputVarNames[i]) /
            brm.GetVarItemsize(OutputVarNames[i]) << "\n";
    }
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI :
    Called by root.
     */
    int GetVarNbytes(const std::string name) override;

    /**
    @a GetVarLocation has no explicit meaning in BMIPhreeqcRM. All
    grid-related information derives from the user's model.
    @param name Name of the variable, but not used.
    @retval The string "Unknown" is returned.
    */
    std::string GetVarLocation(const std::string name) override { return "Unknown"; }

    // Time functions

    /**
    @a GetCurrentTime returns the current simulation time, in seconds.
    (Same as @ref GetTime.)
    @retval                 The current simulation time, in seconds.
    @see
    @ref GetEndTime,
    @ref GetTimeStep,
    @ref SetValue,
    @ref GetTime,
    @ref GetTimeStep,
    @ref SetTime,
    @ref SetTimeStep.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::cout << "Current time: " << GetCurrentTime() << " seconds\n";
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    double GetCurrentTime() override;

    /**
    @a GetStartTime returns the current simulation time, in seconds.
    (Same as @ref GetCurrentTime or @ref GetTime.)
    @retval    The current simulation time, in seconds.
    */
    double GetStartTime() override;

    /**
    @a GetEndTime returns @ref GetCurrentTime plus @ref GetTimeStep, in seconds.
    @retval    The end of the time step, in seconds.
    @see
    @ref GetCurrentTime,
    @ref GetTimeStep,
    @ref SetValue,
    @ref GetTime,
    @ref GetTimeStep,
    @ref SetTime,
    @ref SetTimeStep.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::cout << "End of time step " << GetEndTime() << " seconds\n";
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    double GetEndTime() override;

    /**
    @a GetTimeUnits returns the time units of PhreeqcRM.
    All time units are seconds.
    @retval                 Returns the string "seconds".
    @see
    @ref GetCurrentTime,
    @ref GetEndTime,
    @ref GetTimeStep,
    @ref SetValue,
    @ref GetTime,
    @ref GetTimeStep,
    @ref SetTime,
    @ref SetTimeStep.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::cout << "BMIPhreeqcRM time units are " << GetTimeUnits() << ".\n";
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    std::string GetTimeUnits() override { return "seconds"; };
    double GetTimeStep() override
    {
        return PhreeqcRM::GetTimeStep();
    };

    // Variable getters
    /**
    @a GetValue retrieves model variables. Only variables in the list
    provided by @ref GetOutputVarNames can be retrieved.

    @param name Name of the variable to retrieve.
    @param dest Variable in which to place results. The @a dest variable
    can be of different types depending on the variable retrieved. 
    @n@n
    The
    following list gives the name in the first argument and the
    corresponding data type of the @a dest argument:

    @n "ComponentCount", @a dest: int;
    @n "Components", @a dest: std::vector< std::string >;
    @n "Concentrations", @a dest: std::vector< double >;
    @n "CurrentSelectedOutputUserNumber", @a dest: int;
    @n "DensityCalculated", @a dest: std::vector< double >;
    @n "ErrorString", @a dest: std::string;
    @n "FilePrefix", @a dest: std::string;
    @n "Gfw", @a dest: std::vector< double >;
    @n "GridCellCount", @a dest: int;
    @n "Porosity", @a dest: std::vector< double >;
    @n "Pressure", @a dest: std::vector< double >;
    @n "SaturationCalculated", @a dest: std::vector< double >;
    @n "SelectedOutput", @a dest: std::vector< double >;
    @n "SelectedOutputColumnCount", @a dest: int;
    @n "SelectedOutputCount", @a dest: int;
    @n "SelectedOutputHeadings", @a dest: std::vector< std::string >;
    @n "SelectedOutputOn", @a dest: bool;
    @n "SelectedOutputRowCount", @a dest: int;
    @n "SolutionVolume", @a dest: std::vector< double >;
    @n "Temperature", @a dest: std::vector< double >;
    @n "Time",	@a dest: double;
    @n "TimeStep",	@a dest: double,
    @n "Viscosity", @a dest: std::vector< double >.

    @see
    @ref GetOutputVarNames,
    @ref GetOutputItemCount,
    @ref GetValuePtr;
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
        std::vector< double > density;
        phreeqc_rm.GetValue("DensityCalculated", density);
        std::vector< std::string > comps;
        phreeqc_rm.GetValue("Components", comps);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
    */
    //Add NEW_VARIABLE to GetValue Documentation
    void GetValue(const std::string name, void* dest) override;
    /*!
    * \overload void GetValue(std::string name, bool& OUTPUT);
    */
    void GetValue(const std::string name, bool& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, bool* OUTPUT);
    */
    void GetValue(const std::string name, bool* OUTPUT);
    /*!
    * \overload void GetValue(std::string name, double& OUTPUT);
    */
    void GetValue(const std::string name, double& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, double* OUTPUT);
    */
    void GetValue(const std::string name, double* OUTPUT);
    /*!
    * \overload void GetValue(std::string name, int& OUTPUT);
    */
    void GetValue(const std::string name, int& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, int* OUTPUT);
    */
    void GetValue(const std::string name, int* OUTPUT);
    /*!
    * \overload void GetValue(std::string name, std::string& OUTPUT);
    */
    void GetValue(const std::string name, std::string& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, std::vector<double>& OUTPUT);
    */
    void GetValue(const std::string name, std::vector<double>& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, std::vector<int>& OUTPUT);
    */
    void GetValue(const std::string name, std::vector<int>& OUTPUT);
    /*!
    * \overload void GetValue(std::string name, std::vector<int>& OUTPUT);
    */
    void GetValue(const std::string name, std::vector<std::string>& OUTPUT);
    /**
    @a GetValuePtr takes a variable name and returns a 
    pointer to a current copy of the variable values. Unlike the buffer 
    returned from @ref GetValue, 
    the reference always points to the current values of the variable, 
    even if the model's state has changed.
    @param name Name of the variable to retrieve.
    @retval Void pointer to data. The
    following list gives the name in the argument and the
    data type the void pointer should be cast to:
    @n "ComponentCount": int*;
    @n "Concentrations": double*;
    @n "DensityCalculated": double*;
    @n "Gfw": double*;
    @n "GridCellCount": int*;
    @n "Porosity": double*;
    @n "Pressure": double*;
    @n "SaturationCalculated": double*;
    @n "SelectedOutputOn": bool*;
    @n "SolutionVolume": double*;
    @n "Temperature": double*;
    @n "Time": double*;
    @n "TimeStep": double*;
    @n "Viscosity": double*;
    */
    //Add NEW_VARIABLE to GetValuePtr Documentation
    void* GetValuePtr(std::string name) override;  
    /**
    @a GetValueAtIndices is not implemented
    */
    void GetValueAtIndices(std::string name, void* dest, int* inds, int count) override
    {
        throw NotImplemented();
    };
    
    // Variable setters
    /**
    @a SetValue sets model variables. Only variables in the list
    provided by @ref GetInputVarNames can be set. 
    @param name Name of the variable to retrieve.
    @param src Data used to set the variable. The @a src data type
    can vary depending on the variable retrieved. 
    @n@n
    The following list gives the name in the first argument and the 
    corresponding data type and size of the @a src argument:

    @htmlonly
    <CODE>
    <PRE>
    "Concentrations", std::vector<double>, [GridCellCount*ComponentCount];
    "DensityUser", std::vector<double>, [GridCellCount];
    "FilePrefix", std::string;
    "NthSelectedOutput", int;
    "Porosity", std::vector<double>, [GridCellCount];
    "Pressure", std::vector<double>, [GridCellCount];
    "SaturationUser", std::vector<double>, [GridCellCount];
    "SelectedOutputOn", bool;
    "Temperature", std::vector<double>, [GridCellCount];
    "Time", double;
    "TimeStep", double;
    </PRE>
    </CODE>
    @endhtmlonly
    @see
    @ref GetInputVarNames,
    @ref GetInputItemCount,
    @ref GetValue,
    @ref GetVarItemsize,
    @ref GetVarNbytes,
    @ref GetVarType,
    @ref GetVarUnits.
    @par C++ Example:
    @htmlonly
    <CODE>
    <PRE>
    std::vector< double > temperature(ngrid, 28.0);
    phreeqc_rm.SetValue("Temperature", temperature);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
     */
    void SetValue(const std::string name, void* src) override;
    /*!
    * \overload void SetValue(const std::string name, bool src);
    */
    void SetValue(const std::string name, bool src);
    /*!
    * \overload void SetValue(const std::string name, const char* src);
    */
    void SetValue(const std::string name, const char* src);
    /*!
    * \overload void SetValue(const std::string name, double src);
    */
    void SetValue(const std::string name, double src);
    /*!
    * \overload void SetValue(const std::string name, int src);
    */
    void SetValue(const std::string name, int src);
    /*!
    * \overload void SetValue(const std::string name, std::string src);
    */
    void SetValue(const std::string name, const std::string src);
    /*!
    * \overload void SetValue(const std::string name, std::vector<double> src);
    */
    void SetValue(const std::string name, std::vector<double> src);
    /*!
    * \overload void SetValue(const std::string name, std::vector<int>  src);
    */
    void SetValue(const std::string name, std::vector<int>  src);
    /*!
    * \overload void SetValue(const std::string name, std::vector<std::string>  src);
    */
    void SetValue(const std::string name, std::vector<std::string>  src);
    /**
    @a SetValueAtIndices is not implemented.
    */
    void SetValueAtIndices(std::string name, int* inds, int count, void* src) override
    {
        throw NotImplemented();
    };
    // Grid information functions
    /**
    @a GetGridRank returns a rank of 1 for grid 0. 
    BMIPhreeqcRM only has a 1D series of
    cells; any grid or spatial information must
    be found in the user's model.
    @param grid Grid number, only grid 0 is considered.
    @retval Rank of 1 is returned for grid 0; 0 for
    all other values of @a grid.
    */
    int GetGridRank(const int grid) override;
    /**
    @a GetGridSize returns the number of cells specified
    at creation of the BMIPhreeqcRM instance. 
    @param grid Grid number, only grid 0 is considered.
    @retval Number of cells in the user's modle (same value as GetGridCellCount) is returned for grid 0; 
    0 for all other values of @a grid.
    */
    int GetGridSize(const int grid) override;
    /**
    @a GetGridType is considered to be points. No grid
    information is available in BMIPhreeqcRM; all grid 
    information must be found in the user's model.
    @param grid Grid number, only grid 0 is considered.
    @retval  "Points" is returned for grid 0;
    "Undefined grid identifier" is returned for all other 
    values of @a grid.
    */
    std::string GetGridType(const int grid) override;
    /**
    @a GetGridShape is not implemented.
    */
    void GetGridShape(const int grid, int* shape) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridSpacing is not implemented.
    */
    void GetGridSpacing(const int grid, double* spacing) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridOrigin is not implemented.
    */
    void GetGridOrigin(const int grid, double* origin) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridX is not implemented.
    */
    void GetGridX(const int grid, double* x) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridY is not implemented.
    */
    void GetGridY(const int grid, double* y) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridZ is not implemented.
    */
    void GetGridZ(const int grid, double* z) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridNodeCount is not implemented.
    */
    int GetGridNodeCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridEdgeCount is not implemented.
    */
    int GetGridEdgeCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridFaceCount is not implemented.
    */
    int GetGridFaceCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridEdgeNodes is not implemented.
    */
    void GetGridEdgeNodes(const int grid, int* edge_nodes) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridFaceEdges is not implemented.
    */
    void GetGridFaceEdges(const int grid, int* face_edges) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridFaceNodes is not implemented.
    */
    void GetGridFaceNodes(const int grid, int* face_nodes) override
    {
        throw NotImplemented();
    }
    /**
    @a GetGridNodesPerFace is not implemented.
    */
    void GetGridNodesPerFace(const int grid, int* nodes_per_face) override
    {
        throw NotImplemented();
    }
    /**
    @a AddOutputVars allows selection of sets of variables that can be retieved
    by the @ref GetValue method. Sets of variables can be included or excluded with
    multiple calls to this method. All calls must precede the final call to
    @ref FindComponents. @ref FindComponents generates SELECTED_OUTPUT 333 and
    USER_PUNCH 333 data blocks that make the variables accessible. Variables will
    only be accessible if the system includes the given reactant; for example, no
    gas variables will be created if there are no GAS_PHASEs in the model.

    @param option A string value, among those listed below, that includes or
    excludes variables from @ref GetOutputVarNames, @ref GetValue, and other
    BMI methods.
    @param def A string value that can be "false", "true", or a list of items to be included as
    accessible variables. A value of "false", excludes all variables of the given type; a
    value of "true" includes all variables of the given type for the current system; a list
    specifies a subset of items of the given type.

    Values for the the parameter @a option:
    @n@n
	@a AddOutputVars: 
    False excludes all variables; True causes the settings for each variable group
    to determine the variables that will be defined. Default True;
    @n@n
	@a SolutionProperties: 
    False excludes all solution property variables; True includes variables pH, pe,
    alkalinity, ionic strength, water mass, charge balance, percent error, and specific conductance.
    Default True.
    @n@n
	@a SolutionTotalMolalities: 
    False excludes all total element and element redox state variables;
    True includes all elements and element redox state variables for the system defined for the
    calculation; list restricts variables to the specified elements and redox states.
    Default True.
    @n@n
	@a ExchangeMolalities: 
    False excludes all variables related to exchange; True includes all
    variables related to exchange; list includes variables for the specified exchange species.
    Default True.
    @n@n
	@a SurfaceMolalities: 
    False excludes all variables related to surfaces; True includes all
    variables related to surfaces; list includes variables for the specified surface species.
    Default True.
    @n@n
	@a EquilibriumPhases: 
    False excludes all variables related to equilibrium phases; True includes all
    variables related to equilibrium phases; list includes variables for the specified
    equilibiurm phases. Default True.
    @n@n
	@a Gases: 
    False excludes all variables related to gases; True includes all
    variables related to gases; list includes variables for the specified gas components. Default True.
    @n@n
	@a KineticReactants: 
    False excludes all variables related to kinetic reactants; True includes all
    variables related to kinetic reactants; list includes variables for the specified kinetic
    reactants. Default True.
    @n@n
	@a SolidSolutions: 
    False excludes all variables related to solid solutions; True includes all
    variables related to solid solutions; list includes variables for the specified solid solutions
    components. Default True.
    @n@n
	@a CalculateValues: 
    False excludes all calculate values; True includes all
    calculate values; list includes the specified calculate values. CALCLUATE_VALUES can be
    used to calculate geochemical quantities not available in the other sets of variables.
    Default True.
    @n@n
	@a SolutionActivities: 
    False excludes all aqueous species; True includes all
    aqueous species; list includes only the specified aqueous species. Default False.
    @n@n
	@a SolutionMolalities: 
    False excludes all aqueous species; True includes all
    aqueous species; list includes only the specified aqueous species. Default False.
    @n@n
	@a SaturationIndices: 
    False excludes all saturation indices; True includes all
    saturation indices; list includes only the specified saturation indices. Default False.
    */
    void AddOutputVars(std::string option, std::string def) override;

    IRM_RESULT SetLanguage(const char* string) { this->language = string; return IRM_OK; };
    // data
    std::string language;
    // typedef void (*VarFunction)(BMIPhreeqcRM* brm_ptr); // function pointer type
    // typedef std::map<std::string, VarFunction> VarFunction_map;
    // VarFunction_map varfn_map;
    // VarFunction GetFn(const std::string name);
    //  std::set<std::string> UpdateMap;
    // std::set<std::string>& GetUpdateMap() { return UpdateMap; }

#if defined(WITH_PYBIND11)

    py::array BMIPhreeqcRM::get_value(std::string name, py::array arr);

    //py::array get_value_test(std::string arg, py::array dest/* = py::none()*/);
    //py::array BMIPhreeqcRM::get_value_test(std::string name, py::array_t<double> dest = py::none());
    py::array get_value_ptr(std::string name);

    void set_value(std::string name, py::array src);

    py::array get_value_at_indices(std::string name, py::array dest, py::array indices);

    void set_value_at_indices(std::string name, py::array indices, py::array src);

    py::sequence process_sequence(py::sequence seq);

    bool _initialized;   // { var_man != nullptr }
#endif

#if defined(SWIG) || defined(swig_python_EXPORTS)
    void get_value_ptr_double(std::string var, double** ARGOUTVIEW_ARRAY1, int* DIM1);
    void get_value_ptr_int(std::string var, int** ARGOUTVIEW_ARRAY1, int* DIM1);
    std::vector<std::string>& get_value_ptr_vector_strings(std::string var);
#endif

protected:
    void Construct(Initializer initializer) override;

private:
    //friend class RM_interface;
    VarManager* var_man;
	bool constructed;

    void ClearBMISelectedOutput() override;
    void GenerateAutoOutputVars() override;
    void UpdateBMI(RMVARS v_enum) override;
    void UpdateVariables();
    RMVARS GetEnum(const std::string name);
};
#endif //BMIPHREEQCRM_H_INCLUDED
