#if !defined(BMIPHREEQCRM_H_INCLUDED)
#define BMIPHREEQCRM_H_INCLUDED
#include <map>
#include "PhreeqcRM.h"
#include "BMI_Var.h"
#include "bmi.hxx"
class NotImplemented : public std::logic_error {
public:
    NotImplemented() : std::logic_error("Not Implemented") { };
};

class IRM_DLL_EXPORT BMIPhreeqcRM : public bmi::Bmi, public PhreeqcRM
{
public:
    static void             CleanupBMIModuleInstances(void);
    static int              CreateBMIModule(int nxyz, MP_TYPE nthreads);
    static IRM_RESULT       DestroyBMIModule(int n);
    static BMIPhreeqcRM*    GetInstance(int n);
    /**
    Constructor for the BMIPhreeqcRM subclass of PhreeqcRM. A BMIPhreeqcRM 
    instance has the BMI methods plus all of the PhreeqcRM methods. The 
    constructor requires two arguments: the number of cells in the user's
    model, and either (a) the number of threads for OpenMP parallelization, or
    (b) an MPI communicator.
    @param ngrid Number of cells in the user's model. A value of zero causes
    the program to set nthreads to the number of logical processors of the
    computer.
    @param nthreads Number of threads for parallelization with OpenMP or
    an MPI communicator if PhreeqcRM is compiled with MPI.
    @retval A BMIPhreeqcRM instance.
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
    BMIPhreeqcRM(int ngrid, int nthreads);
    // Model control functions.
    /**
    @ref Initialize is used to initialize a PhreeqcRM instance. This method is equivalent to
    @ref InitializeYAML. A YAML file used for initialization contains a YAML map of PhreeqcRM 
    methods and the arguments corresponding to the method. For example,
    @htmlonly
    <CODE>
    <PRE>
    LoadDatabase: phreeqc.dat
    RunFile:
      workers: true
      initial_phreeqc: true
      utility: true
      chemistry_name: advect.pqi
    </PRE>
    </CODE>
    @endhtmlonly

    @ref Initialize will read the YAML file and execute the specified methods with 
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
    @n CloseFiles();
    @n CreateMapping(std::vector< int >& grid2chem);
    @n DumpModule();
    @n FindComponents();
    @n InitialPhreeqc2Module(std::vector< int > initial_conditions1);
    @n InitialPhreeqc2Module(std::vector< int > initial_conditions1, 
    @n     std::vector< int > initial_conditions2, std::vector< double > fraction1);
    @n InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
    @n LoadDatabase(std::string database);
    @n OpenFiles(void);
    @n OutputMessage(std::string str);
    @n RunCells(void);
    @n RunFile(bool workers, bool initial_phreeqc, 
    @n      bool utility, std::string chemistry_name);
    @n RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    @n ScreenMessage(std::string str);
    @n SetComponentH2O(bool tf);
    @n SetConcentrations(std::vector< double > c);
    @n SetCurrentSelectedOutputUserNumber(int n_user);
    @n SetDensity(std::vector< double > density);
    @n SetDumpFileName(std::string dump_name);
    @n SetErrorHandlerMode(int mode);
    @n SetErrorOn(bool tf);
    SetFilePrefix(std::string prefix);
    @n SetGasCompMoles(std::vector< double > gas_moles);
    @n SetGasPhaseVolume(std::vector< double > gas_volume);
    @n SetPartitionUZSolids(bool tf);
    @n SetPorosity(std::vector< double > por);
    @n SetPressure(std::vector< double > p);
    @n SetPrintChemistryMask(std::vector< int > cell_mask);
    @n SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
    @n SetRebalanceByCell(bool tf);
    @n SetRebalanceFraction(double f);
    @n SetRepresentativeVolume(std::vector< double > rv);
    @n SetSaturation(std::vector< double > sat);
    @n SetScreenOn(bool tf);
    @n SetSelectedOutputOn(bool tf);
    @n SetSpeciesSaveOn(bool save_on);
    @n SetTemperature(std::vector< double > t);
    @n SetTime(double time);
    @n SetTimeConversion(double conv_factor);
    @n SetTimeStep(double time_step);
    @n SetUnitsExchange(int option);
    @n SetUnitsGasPhase(int option);
    @n SetUnitsKinetics(int option);
    @n SetUnitsPPassemblage(int option);
    @n SetUnitsSolution(int option);
    @n SetUnitsSSassemblage(int option);
    @n SetUnitsSurface(int option);
    @n SpeciesConcentrations2Module(std::vector< double > species_conc);
    @n StateSave(int istate);
    @n StateApply(int istate);
    @n StateDelete(int istate);
    @n UseSolutionDensityVolume(bool tf);
    @n WarningMessage(std::string warnstr);
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
    void Initialize(std::string config_file) override;
    /**
    @ref Update runs PhreeqcRM for one time step. This method is equivalent to
    @ref RunCells. PhreeqcRM will equilibrate the solutions with all equilibrium 
    reactants (EQUILIBRIUM_PHASES, EXCHANGE, GAS_PHASE, SOLID_SOLUTIONS, and SURFACE) 
    and integrate KINETICS reactions for the specified time step (@ref SetTimeStep).
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
    @ref UpdateUntil is the same as @ref Update, except the time step is calculated
    from the argument @end_time. The time step is calculated to be @a end_time minus 
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
    @ref Finalize closes any files open in the BMIPhreeqcRM instance.
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
    @ref GetComponentName returns the component name--"BMI PhreeqcRM". 
    BMI PhreeqcRM is a partial interface to PhreeqcRM, and provides 
    the most commonly used methods to implement chemical reactions in a
    multicomponent transport model. All of the native PhreeqcRM methods 
    (non BMI methods) provide are available, which provides a complete 
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
    std::string GetComponentName() override {return "BMI PhreeqcRM";};

    /**
    @ref GetInputVarNames returns the count of input variables that can 
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
    @ref GetOutputItemCount returns the count of output variables that can 
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
    @ref GetPointableItemCount returns the count of output variables for which
    pointers can be obtained with @ref GetValuePtr. The pointers point to 
    current copies of the variables. Setting a value with one of the pointers
    will have no effect on the simulation but will corrupt the copy of the variable.
    @retval  Count of pointers to variables that can be retrieved with 
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
    std::vector< std::string > OutputVarNames = brm.GetPointableVarNames();
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
    @ref GetInputVarNames returns a list of the variable names that can be 
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
    @ref GetOutputVarNames returns a list of the variable names that can 
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
    @ref GetPointableVarNames returns a list of the names of variables
    for which pointers can be retrieved with @ref GetValuePt.
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
    int count = brm.GetOutputItemCount();
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

    // Variable information functions
    /**
    @ref GetVarGrid returns a value of 1, indicating points.
    BMIPhreeqcRM does not have a grid of its own. The cells
    of BMIPhreeqcRM are associated with the user's model grid,
    and all spatial characterists are assigned by the user's
    model.
    @retval 1 BMIPhreeqcRM cells derive meaning from the user's
    model. 
    */
    int GetVarGrid(const std::string name) override {return 1;}

    /**
    Basic Model Interface method that retrieves the type of a variable that 
    can be set with @ref SetValue or retrieved with @ref GetValue. Types are 
    "int", "double", "std::string", or "std::vector<std::string>".

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
    @ref GetVarUnits retrieves the units of a variable 
    that can be set with @ref SetValue, retrieved with @ref GetValue,
    or pointed to by @ref GetValuePtr.
    @param name Name of the variable to retrieve type.
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
    @ref GetVarItemsize retrieves size of an individual item that 
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
    @ref GetVarNbytes retrieves the total number of bytes that are 
    set for a variable with @ref SetValue or retrieved for a variable with 
    @ref GetValue.
    @param name Name of the variable to retrieve total bytes.
    @retval Total number of bytes set or retrieved for variable.
    @see
    @ref GetVarItemsize,
    @ref GetVarType,
    @ref GetVarUnits.
    .
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
    @ref GetVarLocation has no explicit meaning in BMIPhreeqcRM. All
    grid-related information derives from the user's model.
    @param name Name of the variable to retrieve total bytes.
    @retval The string "Unknown" is returned. 
    */
    std::string GetVarLocation(const std::string name) override { return "Unknown"; }

    // Time functions

    /**
    @ref GetCurrentTime returns the current simulation time, in seconds. 
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
    std::cout << "Current time: "
         << GetCurrentTime()
         << " seconds\n";
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    double GetCurrentTime() override;

    /**
    @ref GetStartTime returns the current simulation time, in seconds.
    (Same as @ref GetCurrentTime or @ref GetTime.)
    @retval                 The current simulation time, in seconds.
    */
    double GetStartTime() override;

    /**
    @ref GetEndTime returns @ref GetCurrentTime plus @ref GetTimeStep, in seconds.
    @retval                 The end of the time step, in seconds.
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
    std::cout << "End of time step "
         << GetEndTime()
         << " seconds\n";
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root.
     */
    double GetEndTime() override;

    /**
    @ref GetTimeUnits returns the time units of PhreeqcRM.
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
    std::cout << "BMIPhreeqcRM time units are "
         << GetTimeUnits() << ".\n";
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
    @ref GetValue retrieves model variables. Only variables in the list
    provided by @ref GetOutputVarNames can be retrieved. 

    @param name Name of the variable to retrieve.
    @param dest Variable in which to place results. The @a dest variable
    can be of different types depending on the variable retrieved. The
    following list gives the name in the first argument and the 
    corresponding data type of the @a dest argument:

    @n "ComponentCount", @a dest: int;
    @n "Components", @a dest: std::vector< std::string >;
    @n "Concentrations", @a dest: std::vector< double >;
    @n "CurrentSelectedOutputUserNumber", @a dest: int;
    @n "Density", @a dest: std::vector< double >;
    @n "ErrorString", @a dest: std::string;
    @n "FilePrefix", @a dest: std::string;
    @n "Gfw", @a dest: std::vector< double >;
    @n "GridCellCount", @a dest: int;
    @n "Porosity", @a dest: std::vector< double >;
    @n "Pressure", @a dest: std::vector< double >;
    @n "Saturation", @a dest: std::vector< double >;
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
        phreeqc_rm.GetValue("Density", density);
        std::vector< std::string > comps;
        phreeqc_rm.GetValue("Components", comps);
    </PRE>
    </CODE>
    @endhtmlonly
    @par MPI:
    Called by root, workers must be in the loop of @ref MpiWorker.
    */
    void GetValue(const std::string name, void* dest) override;
    /*!
    * \overload void GetValue(std::string name, bool& dest);
    */
    void GetValue(const std::string name, bool& dest);
    /*!
    * \overload void GetValue(std::string name, bool* dest);
    */
    void GetValue(const std::string name, bool* dest);
    /*!
    * \overload void GetValue(std::string name, double& dest);
    */
    void GetValue(const std::string name, double& dest);
    /*!
    * \overload void GetValue(std::string name, double* dest);
    */
    void GetValue(const std::string name, double* dest);
    /*!
    * \overload void GetValue(std::string name, int& dest);
    */
    void GetValue(const std::string name, int& dest);
    /*!
    * \overload void GetValue(std::string name, int* dest);
    */
    void GetValue(const std::string name, int* dest);
    /*!
    * \overload void GetValue(std::string name, std::string& dest);
    */
    void GetValue(const std::string name, std::string& dest);
    /*!
    * \overload void GetValue(std::string name, std::vector<double>& dest);
    */
    void GetValue(const std::string name, std::vector<double>& dest);
    /*!
    * \overload void GetValue(std::string name, std::vector<int>& dest);
    */
    void GetValue(const std::string name, std::vector<int>& dest);
    /*!
    * \overload void GetValue(std::string name, std::vector<int>& dest);
    */
    void GetValue(const std::string name, std::vector<std::string>& dest);
    /**
    @ref GetValuePtr takes a variable name and returns a 
    reference to a variable. Unlike the buffer returned from @ref GetValue, 
    the reference always points to the current values of the variable, 
    even if the model's state has changed.
    @param name Name of the variable to retrieve.
    @retval Void pointer to data.
    */
    void* GetValuePtr(std::string name) override;  
    /**
    @ref GetValueAtIndices is not implemented
    */
    void GetValueAtIndices(std::string name, void* dest, int* inds, int count) override
    {
        throw NotImplemented();
    };
    
    // Variable setters
    /**
    @ref SetValue sets model variables. Only variables in the list
    provided by @ref GetInputVarNames can be set. 
    @param name Name of the variable to retrieve.
    @param src Data used to set the variable. The @a src data type
    can vary depending on the variable retrieved. The
    following list gives the name in the first argument and the 
    corresponding data type and size of the @a src argument:

    "Concentrations", std::vector<double>, [GridCellCount*ComponentCount];
    "Density", std::vector<double>, [GridCellCount];
    "FilePrefix", std::string;
    "NthSelectedOutput", int;
    "Porosity", std::vector<double>, [GridCellCount];
    "Pressure", std::vector<double>, [GridCellCount];
    "Saturation", std::vector<double>, [GridCellCount];
    "SelectedOutputOn", bool;
    "Temperature", std::vector<double>, [GridCellCount];
    "Time", double;
    "TimeStep", double;

    @see
    @ref GetInputVarNames,
    @ref GetInputItemCount,,
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
    void SetValue(std::string name, void* src) override;
    /*!
    * \overload void SetValue(std::string name, bool src);
    */
    void SetValue(std::string name, bool src);
    /*!
    * \overload void SetValue(std::string name, char* src);
    */
    void SetValue(std::string name, char* src);
    /*!
    * \overload void SetValue(std::string name, double src);
    */
    void SetValue(std::string name, double src);
    /*!
    * \overload void SetValue(std::string name, int src);
    */
    void SetValue(std::string name, int src);
    /*!
    * \overload void SetValue(std::string name, std::string src);
    */
    void SetValue(std::string name, std::string src);
    /*!
    * \overload void SetValue(std::string name, std::vector<double> src);
    */
    void SetValue(std::string name, std::vector<double> src);
    /*!
    * \overload void SetValue(std::string name, std::vector<int>  src);
    */
    void SetValue(std::string name, std::vector<int>  src);
    /*!
    * \overload void SetValue(std::string name, std::vector<std::string>  src);
    */
    void SetValue(std::string name, std::vector<std::string>  src);
    /**
    @ref SetValueAtIndices is not implemented.
    */
    void SetValueAtIndices(std::string name, int* inds, int count, void* src) override
    {
        throw NotImplemented();
    };
    // Grid information functions
    /**
    @ref GetGridRank returns a rank of 1 for grid 0. 
    BMIPhreeqcRM only has a 1D series of
    cells; any grid or spatial information must
    be found in the user's model.
    @param grid Grid number, only grid 0 is considered.
    @retval Rank of 1 is returned for grid 0; 0 for
    all other values of @a grid.
    */
    int GetGridRank(const int grid) override;
    /**
    @ref GetGridSize returns the number of cells specified
    at creation of the BMIPhreeqcRM instance. 
    @param grid Grid number, only grid 0 is considered.
    @retval Same value as GetGridCellCount is returned for grid 0; 
    0 for all other values of @a grid.
    */
    int GetGridSize(const int grid) override;
    /**
    @ref GetGridType is considered to be points. No grid
    information is available in BMIPhreeqcRM; all grid 
    information must be found in the user's model.
    @param grid Grid number, only grid 0 is considered.
    @retval  "Points" is returned for grid 0;
    "Undefined grid identifier" is returned for all other 
    values of @a grid.
    */
    std::string GetGridType(const int grid) override;
    /**
    @ref GetGridShape is not implemented.
    */
    void GetGridShape(const int grid, int* shape) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridSpacing is not implemented.
    */
    void GetGridSpacing(const int grid, double* spacing) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridOrigin is not implemented.
    */
    void GetGridOrigin(const int grid, double* origin) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridX is not implemented.
    */
    void GetGridX(const int grid, double* x) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridY is not implemented.
    */
    void GetGridY(const int grid, double* y) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridZ is not implemented.
    */
    void GetGridZ(const int grid, double* z) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridNodeCount is not implemented.
    */
    int GetGridNodeCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridEdgeCount is not implemented.
    */
    int GetGridEdgeCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridFaceCount is not implemented.
    */
    int GetGridFaceCount(const int grid) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridEdgeNodes is not implemented.
    */
    void GetGridEdgeNodes(const int grid, int* edge_nodes) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridFaceEdges is not implemented.
    */
    void GetGridFaceEdges(const int grid, int* face_edges) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridFaceNodes is not implemented.
    */
    void GetGridFaceNodes(const int grid, int* face_nodes) override
    {
        throw NotImplemented();
    }
    /**
    @ref GetGridNodesPerFace is not implemented.
    */
    void GetGridNodesPerFace(const int grid, int* nodes_per_face) override
    {
        throw NotImplemented();
    }
    // data
    std::string language;
   // typedef void (*VarFunction)(BMIPhreeqcRM* brm_ptr); // function pointer type
   // typedef std::map<std::string, VarFunction> VarFunction_map;
   // VarFunction_map varfn_map;
   // VarFunction GetFn(const std::string name);
   //  std::set<std::string> UpdateMap;
   // std::set<std::string>& GetUpdateMap() { return UpdateMap; }

private:
    //friend class RM_interface;
    static std::map<size_t, BMIPhreeqcRM*> Instances;
    static size_t InstancesIndex;
    void UpdateVariables();
};
#endif //BMIPHREEQCRM_H_INCLUDED
