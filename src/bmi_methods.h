// BMI data and methods
private:
	std::map<std::string, class BMI_Var> bmi_var_map;
	std::vector< std::string > bmi_input_vars;
	std::vector< std::string > bmi_output_vars;
public:

	//void BMI_Finalize();
/**
Basic Model Interface method that returns the component name--PhreeqcRM. The BMI interface to PhreeqcRM is
only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
prefix) provide a complete interface, and it is expected that the native methods will be used in preference to the BMI_
methods.

@retval The string "PhreeqcRM".
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::string comp_name = phreeqc_rm.BMI_GetComponentName();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
	std::string BMI_GetComponentName() { return "PhreeqcRM"; }

	//--------------------------	

	/**
	Basic Model Interface method that returns the current simulation time, in seconds. (Same as @ref GetTime.)
	The reaction module does not change the time value, so the
	returned value is equal to the default (0.0) or the last time set by
	@ref BMI_SetValue("Time", time) or @ref SetTime.
	@retval                 The current simulation time, in seconds.
	@see
	@ref BMI_GetEndTime,
	@ref BMI_GetTimeStep,
	@ref BMI_SetValue,
	@ref GetTime,
	@ref GetTimeStep,
	@ref SetTime,
	@ref SetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::cout << "Current time: "
		 << BMI_GetCurrentTime()
		 << " seconds\n";
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	double BMI_GetCurrentTime() { return this->GetTime(); }

	//--------------------------	

	/**
	Basic Model Interface method that returns @ref BMI_GetCurrentTime plus @ref BMI_GetTimeStep, in seconds.
	@retval                 The end of the time step, in seconds.
	@see
	@ref BMI_GetCurrentTime,
	@ref BMI_GetTimeStep,
	@ref BMI_SetValue,
	@ref GetTime,
	@ref GetTimeStep,
	@ref SetTime,
	@ref SetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::cout << "End of time step "
		 << BMI_GetEndTime()
		 << " seconds\n";
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	double BMI_GetEndTime() { return this->GetTime() + this->GetTimeStep(); }

	//--------------------------	

	/**
	Basic Model Interface method that returns count of input variables that can be set with @ref BMI_SetValue.
	@retval  Count of input variables that can be set with @ref BMI_SetValue.

	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits,
	@ref BMI_SetValue.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			std::vector< std::string > InputVarNames = phreeqc_rm.BMI_GetInputVarNames();
			int count = phreeqc_rm.BMI_GetInputItemCount();
			oss << "BMI_SetValue variables:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << InputVarNames[i] << "\n";
				oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(InputVarNames[i]) << "\n";
				oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(InputVarNames[i]) << "\n";
				oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) << "\n";
				oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
				oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) /
											   phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	int BMI_GetInputItemCount() { return (int)this->bmi_input_vars.size(); }

	//--------------------------	

	/**
	Basic Model Interface method that returns a list of the variable names that can be set with @ref BMI_SetValue.
	@retval  A list of the variable names that can be set with @ref BMI_SetValue.

	@see
	@ref BMI_GetInputItemCount,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits,
	@ref BMI_SetValue.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			std::vector< std::string > InputVarNames = phreeqc_rm.BMI_GetInputVarNames();
			int count = phreeqc_rm.BMI_GetInputItemCount();
			oss << "BMI_SetValue variables:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << InputVarNames[i] << "\n";
				oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(InputVarNames[i]) << "\n";
				oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(InputVarNames[i]) << "\n";
				oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) << "\n";
				oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
				oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) /
											   phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	std::vector< std::string > BMI_GetInputVarNames() { return this->bmi_input_vars; }

	//--------------------------	

	/**
	Basic Model Interface method that returns count of output variables that can be retrieved with @ref BMI_GetValue.
	@retval  Count of output variables that can be retrieved with @ref BMI_GetValue.

	@see
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetValue,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			std::vector< std::string > OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
			int count = phreeqc_rm.BMI_GetOutputItemCount();
			oss << "BMI_GetValue variables:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << OutputVarNames[i] << "\n";
				oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
				oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
				oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) << "\n";
				oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
				oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) /
											   phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	int BMI_GetOutputItemCount() { return (int)this->bmi_output_vars.size(); }

	//--------------------------	

	/**
	Basic Model Interface method that returns a list of the variable names that can be retrieved with @ref BMI_GetValue.
	@retval  A list of the variable names that can be retrieved with @ref BMI_GetValue.

	@see
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetValue,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			std::vector< std::string > OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
			int count = phreeqc_rm.BMI_GetOutputItemCount();
			oss << "BMI_GetValue variables:\n";
			for (size_t i = 0; i < count; i++)
			{
				oss << "  " << i << "  " << OutputVarNames[i] << "\n";
				oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
				oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
				oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) << "\n";
				oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
				oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) /
											   phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	std::vector< std::string > BMI_GetOutputVarNames() { return this->bmi_output_vars; };

	//--------------------------	

	/**
	Basic Model Interface method that returns the current simulation time step, in seconds. (Same as @ref GetTimeStep.)
	The reaction module does not change the time-step value, so the
	returned value is equal to the last time step set by
	@ref BMI_SetValue("TimeStep", time_step) or @ref SetTimeStep.
	@retval                 The current simulation time step, in seconds.
	@see
	@ref BMI_GetCurrentTime,
	@ref BMI_GetEndTime,
	@ref BMI_SetValue,
	@ref GetTime,
	@ref GetTimeStep,
	@ref SetTime,
	@ref SetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::cout << "Current time step: "
		 << BMI_GetTimeStep()
		 << " seconds\n";
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	double BMI_GetTimeStep() { return this->GetTimeStep(); }

	//--------------------------	

	/**
	Basic Model Interface method that returns the time units of PhreeqcRM.
	All time units are seconds for PhreeqcRM.
	@retval                 Returns the string "seconds".
	@see
	@ref BMI_GetCurrentTime,
	@ref BMI_GetEndTime,
	@ref BMI_GetTimeStep,
	@ref BMI_SetValue,
	@ref GetTime,
	@ref GetTimeStep,
	@ref SetTime,
	@ref SetTimeStep.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::cout << "PhreeqcRM time units are "
		 << BMI_GetTimeUnits() << ".\n";
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	std::string BMI_GetTimeUnits() { return "seconds"; };

	//--------------------------	

	/**
	Basic Model Interface method that retrieves model variables. Only variables in the list
	provided by @ref BMI_GetOutputVarNames can be retrieved. The BMI interface to PhreeqcRM is
	only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
	prefix) provide a complete interface.
	@param name Name of the variable to retrieve.
	@param dest Variable in which to place results.

	Variable names for the first argument (@a name) and variable type of the
	second argument (@a dest).
	@n "ComponentCount", @a dest: int;
	@n "Components", @a dest: std::vector< std::string >;
	@n "Concentrations", @a dest: std::vector< double >;
	@n "CurrentSelectedOutputUserNumber", @a dest: int;
	@n "Density", @a dest: std::vector< double >;
	@n "ErrorString", @a dest: std::string;
	@n "FilePrefix", @a dest: std::string;
	@n "Gfw", @a dest: std::vector< double >;
	@n "GridCellCount", @a dest: int;
	@n "InputVarNames", @a dest: std::vector< std::string >;
	@n "OutputVarNames", @a dest: std::vector< std::string >;
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
	@n "TimeStep",	@a dest: double.

	@see
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits,
	@ref BMI_SetValue.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
		std::vector< double > bmi_density;
		phreeqc_rm.BMI_GetValue("Density", bmi_density);
		std::vector< std::string > bmi_comps;
		phreeqc_rm.BMI_GetValue("Components", bmi_comps);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root, workers must be in the loop of @ref MpiWorker.
	*/
	void BMI_GetValue(std::string name, bool& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, double& dest)
 	*/
	void BMI_GetValue(std::string name, double& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, int& dest)
 	*/
	void BMI_GetValue(std::string name, int& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, std::string& dest)
 	*/
	void BMI_GetValue(std::string name, std::string& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, std::vector < double >& dest)
 	*/
	void BMI_GetValue(std::string name, std::vector < double >& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, std::vector < std::string >& dest)
 	*/
	void BMI_GetValue(std::string name, std::vector < std::string >& dest);
	/*!
 	* \overload void BMI_GetValue(std::string name, void* dest)
 	*/
	void BMI_GetValue(std::string name, void* dest);

	//--------------------------	

	/**
	Basic Model Interface method that retrieves size of an 
	individual item that can be set or retrived.
	Sizes may be sizeof(int), sizeof(double), 
	or a character length for string variables. Only variables in the list
	provided by @ref BMI_GetInputVarNames can be set. 
	Only variables in the list
	provided by @ref BMI_GetOutputVarNames can be retrieved. 
	@param name Name of the variable to retrieve size.
	@retval Size of one element of variable.

	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetInputItemCount,
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetValue,
	@ref BMI_GetVarNbytes,
	@ref BMI_SetValue.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
		int nbytes = phreeqc_rm.BMI_GetVarNbytes("Temperature");
		int item_size = phreeqc_rm.BMI_GetVarItemSize("Temperature");
		int dim = nbytes/item_size;
		std::vector< double > bmi_temperature(dim, 25.0);
		phreeqc_rm.BMI_SetValue("Temperature", bmi_temperature);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	int BMI_GetVarItemsize(std::string name);

	//--------------------------	

	/**
	Basic Model Interface method that retrieves the total number of bytes that are set for a variable with
	@ref BMI_SetValue or retrieved for a variable with @ref BMI_GetValue.
	Only variables in the list
	provided by @ref BMI_GetInputVarNames can be set. 
	Only variables in the list
	provided by @ref BMI_GetOutputVarNames can be retrieved. 
	@param name Name of the variable to retrieve total bytes.
	@retval Total number of bytes set or retrieved for variable.
	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetInputItemCount,
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetValue,
	@ref BMI_GetVarItemsize,
	@ref BMI_SetValue.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
		int nbytes = phreeqc_rm.BMI_GetVarNbytes("Temperature");
		int item_size = phreeqc_rm.BMI_GetVarItemSize("Temperature");
		int dim = nbytes/item_size;
		std::vector< double > bmi_temperature(dim, 25.0);
		phreeqc_rm.BMI_SetValue("Temperature", bmi_temperature);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	int BMI_GetVarNbytes(std::string name);

	//--------------------------	

	/**
	Basic Model Interface method that retrieves the type of a variable that can be set with
	@ref BMI_SetValue or retrieved with @ref BMI_GetValue. Types are "int", "double", or "string".
	Only variables in the list
	provided by @ref BMI_GetInputVarNames can be set. 
	Only variables in the list
	provided by @ref BMI_GetOutputVarNames can be retrieved. 

	@param name Name of the variable to retrieve type.
	@retval Character string of variable type.

	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetInputItemCount,
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetVarUnits.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector< std::string > OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
	int count = phreeqc_rm.BMI_GetOutputItemCount();
	oss << "BMI_GetValue variables:\n";
	for (size_t i = 0; i < count; i++)
	{
		oss << "  " << i << "  " << OutputVarNames[i] << "\n";
		oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
		oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
		oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) << "\n";
		oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
		oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) /
										phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
	}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	std::string BMI_GetVarType(std::string name);

	//--------------------------	

	/**
	Basic Model Interface method that retrieves the units of a variable that can be set with
	@ref BMI_SetValue or retrieved with @ref BMI_GetValue.
	Only variables in the list
	provided by @ref BMI_GetInputVarNames can be set. 
	Only variables in the list
	provided by @ref BMI_GetOutputVarNames can be retrieved. 

	@param name Name of the variable to retrieve type.
	@retval Character string of units for variable.
	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetInputItemCount,
	@ref BMI_GetOutputVarNames,
	@ref BMI_GetOutputItemCount,
	@ref BMI_GetVarType.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	std::vector< std::string > OutputVarNames = phreeqc_rm.BMI_GetOutputVarNames();
	int count = phreeqc_rm.BMI_GetOutputItemCount();
	oss << "BMI_GetValue variables:\n";
	for (size_t i = 0; i < count; i++)
	{
		oss << "  " << i << "  " << OutputVarNames[i] << "\n";
		oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(OutputVarNames[i]) << "\n";
		oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(OutputVarNames[i]) << "\n";
		oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) << "\n";
		oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
		oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(OutputVarNames[i]) /
										phreeqc_rm.BMI_GetVarItemsize(OutputVarNames[i]) << "\n";
	}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root.
	 */
	std::string BMI_GetVarUnits(std::string name);

	//--------------------------
#ifdef USE_YAML
	/**
	Basic Model Interface method that can be used to initialize a PhreeqcRM instance. This method is equivalent to
	@ref InitializeYAML. A YAML file can be used in initialization. The file contains a YAML map of PhreeqcRM methods
	and the arguments corresponding to the method. For example,
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

	BMI_Initialize will read the YAML file and execute the specified methods with the specified arguments. Using YAML
	terminology, the argument(s) for a method may be a scalar, a sequence, or a map, depending if the argument is
	a single item, a single vector, or there are multiple arguments. In the case of a map, the name associated
	with each argument (for example "chemistry_name" above) is arbitrary. The names of the map keys for map
	arguments are not used in parsing the YAML file; only the order of the arguments is important.

	The PhreeqcRM methods that can be specified in a YAML file include:
	@htmlonly
	<CODE>
	<PRE>
	CloseFiles(void);
	CreateMapping(std::vector< int >& grid2chem);
	DumpModule();
	FindComponents();
	InitialPhreeqc2Module(std::vector< int > initial_conditions1);
	InitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
	InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
	LoadDatabase(std::string database);
	OpenFiles(void);
	OutputMessage(std::string str);
	RunCells(void);
	RunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
	RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
	ScreenMessage(std::string str);
	SetComponentH2O(bool tf);
	SetConcentrations(std::vector< double > c);
	SetCurrentSelectedOutputUserNumber(int n_user);
	SetDensity(std::vector< double > density);
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
	SetSaturation(std::vector< double > sat);
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
	@ref BMI_Update.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			int nthreads = 0;
			std::string yaml_file = "myfile.yaml";
			int nxyz = GetGridCellCountYAML(yaml_file);
			PhreeqcRM phreeqc_rm(nxyz, nthreads);
			phreeqc_rm.BMI_Initialize(yaml_file);
			int ncomps;
			phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
			int ngrid;
			phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
			std::vector< double > c;
			phreeqc_rm.BMI_GetValue("Concentrations", c);
			phreeqc_rm.BMI_SetValue("TimeStep", 86400);
			for(double time = 0; time < 864000; time+=86400)
			{
				// Take a transport time step here and update the vector c.
				phreeqc_rm.BMI_SetValue("Time", time);
				phreeqc_rm.BMI_SetValue("Concentrations", c);
				phreeqc_rm.BMI_Update();
				phreeqc_rm.BMI_GetValue("Concentrations", c);
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root, workers must be in the loop of @ref MpiWorker.
	 */
	void BMI_Initialize(std::string config_file) { InitializeYAML(config_file); };
#endif
	//--------------------------	

	/**
	Basic Model Interface method that sets model variables. Only variables in the list
	provided by @ref BMI_GetInputVarNames can be set. The BMI interface to PhreeqcRM is
	only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
	prefix) provide a complete interface, and it is expected that the native methods will be used in preference to the BMI_
	methods.

	Variable names for the first argument
	of BMI_SetValue and the equivalent PhreeqcRM method are as follows:
	"Concentrations", @ref SetConcentrations;
	"Density", @ref SetDensity;
	"FilePrefix", @ref SetFilePrefix;
	"NthSelectedOutput", @ref SetNthSelectedOutput;
	"Porosity", @ref SetPorosity;
	"Pressure", @ref SetPressure;
	"Saturation", @ref SetSaturation;
	"SelectedOutputOn", @ref SetSelectedOutputOn;
	"Temperature", @ref SetTemperature;
	"Time", @ref SetTime;
	"TimeStep", @ref SetTimeStep.

	@see
	@ref BMI_GetInputVarNames,
	@ref BMI_GetInputItemCount,,
	@ref BMI_GetValue,
	@ref BMI_GetVarItemsize,
	@ref BMI_GetVarNbytes,
	@ref BMI_GetVarType,
	@ref BMI_GetVarUnits.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			std::vector< double > bmi_temperature(ngrid, 28.0);
			phreeqc_rm.BMI_SetValue("Temperature", bmi_temperature);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root, workers must be in the loop of @ref MpiWorker.
	 */
	void BMI_SetValue(std::string name, void* src);
	void BMI_SetValue(std::string name, bool& src);
	void BMI_SetValue(std::string name, double& src);
	void BMI_SetValue(std::string name, int& src);
	void BMI_SetValue(std::string name, const std::string& src);
	void BMI_SetValue(std::string name, std::vector < double >& src);
	void BMI_SetValue(std::string name, std::vector < int >& src);
	void BMI_SetValue(std::string name, std::vector < std::string >& src);

	//--------------------------	
	/**
	Basic Model Interface method that runs PhreeqcRM for one time step. This method is equivalent to
	@ref RunCells. PhreeqcRM will equilibrate the solutions with all equilibrium reactants (EQUILIBRIUM_PHASES,
	EXCHANGE, GAS_PHASE, SOLID_SOLUTIONS, and SURFACE) and
	integrate KINETICS reactions for the specified time step (@ref SetTimeStep).
	@see
	@ref BMI_Initialize.
	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
			PhreeqcRM phreeqc_rm(nxyz, nthreads);
			phreeqc_rm.BMI_Initialize("myfile.yaml");
			int ncomps;
			phreeqc_rm.BMI_GetValue("ComponentCount", &ncomps);
			int ngrid;
			phreeqc_rm.BMI_GetValue("GridCellCount", ngrid);
			std::vector< double > c;
			phreeqc_rm.BMI_GetValue("Concentrations", c);
			phreeqc_rm.BMI_SetValue("TimeStep", 86400);
			for(double time = 0; time < 864000; time+=86400)
			{
				// Take a transport time step here and update the vector c.
				phreeqc_rm.BMI_SetValue("Time", time);
				phreeqc_rm.BMI_SetValue("Concentrations", c);
				phreeqc_rm.BMI_Update();
				phreeqc_rm.BMI_GetValue("Concentrations", c);
			}
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root, workers must be in the loop of @ref MpiWorker.
	 */
	void BMI_Update(void) { this->RunCells(); };
#ifdef NOT_IMPLEMENTED
	void BMI_UpdateUntil(double time)
	{
		//throw LetItThrow("Not implemented");
		ErrorMessage("Not implemented");
		throw PhreeqcRMStop();
	}
	int BMI_GetVarGrid(std::string name)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	std::string BMI_GetVarLocation(std::string name)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	double BMI_GetStartTime()
	{
		//throw LetItThrow("Not implemented");
		ErrorMessage("Not implemented");
		throw PhreeqcRMStop();
	};
	void* BMI_GetValuePtr(std::string name)
	{
		//throw LetItThrow("Not implemented");
		ErrorMessage("Not implemented");
		throw PhreeqcRMStop();
	}
	void BMI_GetValueAtIndices(std::string name, void* dest, int* inds, int count)
	{
		//throw LetItThrow("Not implemented");
		ErrorMessage("Not implemented");
		throw PhreeqcRMStop();
	};
	void BMI_SetValueAtIndices(std::string name, int* inds, int len, void* src)
	{
		//throw LetItThrow("Not implemented");
		ErrorMessage("Not implemented");
		throw PhreeqcRMStop();
	};
	int BMI_GetGridRank(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	int BMI_GetGridSize(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	std::string BMI_GetGridType(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridShape(const int grid, int* shape)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridSpacing(const int grid, double* spacing)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void GetGridOrigin(const int grid, double* origin)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridX(const int grid, double* x)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridY(const int grid, double* y)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridZ(const int grid, double* z)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};

	int BMI_GetGridNodeCount(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	int BMI_GetGridEdgeCount(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	int BMI_GetGridFaceCount(const int grid)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridEdgeNodes(const int grid, int* edge_nodes)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridFaceEdges(const int grid, int* face_edges)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridFaceNodes(const int grid, int* face_nodes)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
	void BMI_GetGridNodesPerFace(const int grid, int* nodes_per_face)
	{
		//throw LetItThrow("Not applicable");
		ErrorMessage("Not applicable");
		throw PhreeqcRMStop();
	};
#endif
private:
	void BMI_MakeVarMap();
// End BMI data and methods
