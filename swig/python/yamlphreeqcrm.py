import numpy as np
import yaml as yaml
import sys
"""
    Helper module for generating a YAML file for 
    initializing a PhreeqcRM instance. The YAML
    file is likely to be written by a GUI or preprocessor,
    and the PhreeqcRM method InitializeYAML can be 
    used to execute the methods recorded in the YAML
    file to set properties and initial conditions in
    the PhreeqcRM instance.

    Functions that accept an array should be callable
    with any of three Python data types, 
    np.ndarray, list, or tuple.

    For example, the porosity can be initialized:

    # numpy.ndarray
    porosity = np.full((nxyz,), 0.2)

    # list (using list repetition)
    porosity = [0.2] * nxyz

    # list (using list comprehension)
    porosity = [0.2 for i in range(nxyz)]

    # tuple (using a generator expression)
    porosity = tuple(0.2 for i in range(nxyz))

    # And used:
    yrm = yamlphreeqcrm.YAMLPhreeqcRM()
    yrm.YAMLSetPorosity(por)

    # likewise for integer arrays

"""

class YAMLPhreeqcRM(object):
    def __init__(self):
        self.yaml_doc = []
    def Clear(self):
        """
        Clears all definitions from the YAML document.
        """
        self.yaml_doc.clear()
    def WriteYAMLDoc(self, file_name):
        """
        Writes the YAML document to file.
	
        @param file_name        Name of file where YAML document will be written.
        """
        doc = yaml.dump_all(self.yaml_doc)
        file = open(file_name, "w")
        yaml.dump(self.yaml_doc, file, sort_keys=False)
        file.close()
    def YAMLAddOutputVars(self, option, definition):
        """
        Inserts data into the YAML document for the method AddOutputVars.

        When the YAML document is written to file it can be processed by the BMI method Initialize,
        or the PhreeqcRM method InitializeYAML to initialize a PhreeqcRM instance. 
        AddOutputVars allows selection of sets of output variables that are available through 
        the BMI method GetValue. 
        @param option A string value, among those listed below, that includes or
        excludes variables from @ref GetOutputVarNames, @ref GetValue, and other
        BMI methods.
        @param def A string value that can be "false", "true", or a list of items to be included as
        accessible variables.A value of "false", excludes all variables of the given type; a
        value of "true" includes all variables of the given type for the current system; a list
        specifies a subset of items of the given type.

        Values for the the parameter @a option :
        @n AddOutputVars : False excludes all variables; True causes the settings for each variable group
        to determine the variables that will be defined.Default True;
        @n SolutionProperties : False excludes all solution property variables; True includes variables pH, pe,
        alkalinity, ionic strength, water mass, charge balance, percent error, and specific conductance.
        Default True.
        @n SolutionTotalMolalities : False excludes all total element and element redox state variables;
        True includes all elements and element redox state variables for the system defined for the
        calculation; list restricts variables to the specified elements and redox states.
        Default True.
        @n ExchangeMolalities : False excludes all variables related to exchange; True includes all
        variables related to exchange; list includes variables for the specified exchange species.
        Default True.
        @n SurfaceMolalities : False excludes all variables related to surfaces; True includes all
        variables related to surfaces; list includes variables for the specified surface species.
        Default True.
        @n EquilibriumPhases : False excludes all variables related to equilibrium phases; True includes all
        variables related to equilibrium phases; list includes variables for the specified
        equilibiurm phases.Default True.
        @n Gases : False excludes all variables related to gases; True includes all
        variables related to gases; list includes variables for the specified gas components.Default True.
        @n KineticReactants : False excludes all variables related to kinetic reactants; True includes all
        variables related to kinetic reactants; list includes variables for the specified kinetic
        reactants.Default True.
        @n SolidSolutions : False excludes all variables related to solid solutions; True includes all
        variables related to solid solutions; list includes variables for the specified solid solutions
        components.Default True.
        @n CalculateValues : False excludes all calculate values; True includes all
        calculate values; list includes the specified calculate values.CALCLUATE_VALUES can be
        used to calculate geochemical quantities not available in the other sets of variables.
        Default True.
        @n SolutionActivities : False excludes all aqueous species; True includes all
        aqueous species; list includes only the specified aqueous species.Default False.
        @n SolutionMolalities : False excludes all aqueous species; True includes all
        aqueous species; list includes only the specified aqueous species.Default False.
        @n SaturationIndices : False excludes all saturation indices; True includes all
        saturation indices; list includes only the specified saturation indices. Default False.
        """
        node = dict()
        node["key"] = "AddOutputVars"
        node["option"] = option
        node["definition"] = definition
        self.yaml_doc.append(node)
    def YAMLCloseFiles(self, void):
        """
        Inserts data into the YAML document for the PhreeqcRM method CloseFiles.

        When the YAML document is written to file it can be processed by the method InitializeYAML
        to initialize a PhreeqcRM instance. CloseFiles closes the output and log files.
        """
        node = dict()
        node["key"] = "CloseFiles"
        self.yaml_doc.append(node)
    def YAMLCreateMapping(self, grid2chem):
        """
        Inserts data into the YAML document for the PhreeqcRM method CreateMapping.
	
        When the YAML document is written to file it can be processed by the method 
        InitializeYAML to 	initialize a PhreeqcRM instance. CreateMapping 	
        provides a mapping from grid cells in the user's model to reaction cells 
        for which chemistry needs to be run. The mapping is used to eliminate 
        inactive cells and to use symmetry to decrease the number of cells 	for 
        which chemistry must be run. The array @a grid2chem of size @a nxyz (the 
        number of grid cells) 	must contain the set of all integers 0 <= @a i < @a 
        count_chemistry, where @a count_chemistry is a number less than or equal 
        to @a nxyz. Inactive cells are assigned a negative integer. The 
        mapping may be many-to-one to account for symmetry. Default is a 
        one-to-one mapping--all user grid cells are reaction cells 	(equivalent 
        to @a grid2chem values of 0,1,2,3,...,nxyz-1).
	
        @param grid2chem A list, tuple, or numpy array of integers: 
        Nonnegative is a reaction-cell number (0 based), negative is an inactive 
        cell. A list, tuple, or numpy array is of size @a nxyz (number of grid cells).
        """
        node = dict()
        node["key"] = "CreateMapping"
        node["grid2chem"] = self.GetList(grid2chem)
        self.yaml_doc.append(node)
    def YAMLDumpModule(self, dump_on, append):
        """
        Inserts data into the YAML document for the PhreeqcRM method DumpModule.
	
        When the YAML document is written to file it can be processed by the method 
        InitializeYAML to 	initialize a PhreeqcRM instance. DumpModule writes 
        the contents of all workers to file in _RAW formats (see appendix of 
        PHREEQC version 3 manual), 	including SOLUTIONs and all reactants.
        @param dump_on Signal for writing the dump file, true or false.
        @param append Signal to append to the contents of the dump file, true 
        or false.
	    """
        node = dict()
        node["key"] = "DumpModule"
        node["dump_on"] = dump_on
        node["append"] = append
        self.yaml_doc.append(node)
    def YAMLFindComponents(self):
        """
        Inserts data into the YAML document for the PhreeqcRM method FindComponents.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. FindComponents accumulates a list of elements. Elements 
        are those that have been defined in a solution or any other reactant
        (EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
        This method can be called multiple times and the list that is created is cummulative.
        The list is the set of components that needs to be transported. By default the list
        includes water, excess H and excess O (the H and O not contained in water);
        alternatively, the list may be set to contain total H and total O (YAMLSetComponentH2O),
        which requires transport results to be accurate to eight or nine significant digits.
        If multicomponent diffusion (MCD) is to be modeled, there is a capability to retrieve 
        aqueous species concentrations and to set new solution concentrations after
        MCD by using individual species concentrations (YAMLSpeciesConcentrations2Module).
        To use these methods, the save-species property needs to be turned on (YAMLSetSpeciesSaveOn).
        If the save-species property is on, FindComponents will generate
        a list of aqueous species, their diffusion coefficients at 25 C, and their charge.

        The FindComponents method also generates lists of reactants--equilibrium phases,
        exchangers, gas components, kinetic reactants, solid solution components, and surfaces.
        The lists are cumulative, including all reactants that were
        defined in the initial phreeqc instance at any time FindComponents was called.
        In addition, a list of phases is generated for which saturation indices may be calculated from the
        cumulative list of components.
        """
        node = dict()
        node["key"] = "FindComponents"
        self.yaml_doc.append(node)
    def YAMLInitialSolutions2Module(self, solutions):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialSolutions2Module.

        When the YAML document is written to file it can be processed by the method 
        InitializeYAML to initialize a PhreeqcRM instance. InitialSolutions2Module 
        transfers SOLUTION definitions from the InitialPhreeqc instance to the reaction-module 
        workers.
        @param solutions A list, tuple, or numpy array of index numbers that is 
        dimensioned nxyz, where nxyz is the number of grid cells in the user's model. 
        """
        node = dict()
        node["key"] = "InitialSolutions2Module"
        node["solutions"] = self.GetList(solutions)
        self.yaml_doc.append(node)
    def YAMLInitialEquilibriumPhases2Module(self, equilibrium_phases):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialEquilibriumPhases2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialEquilibriumPhases2Module 
        transfers EQUILIBRIUM_PHASES definitions from the InitialPhreeqc instance to the 
        reaction-module workers.
        @param equilibrium_phases A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialEquilibriumPhases2Module"
        node["equilibrium_phases"] = self.GetList(equilibrium_phases)
        self.yaml_doc.append(node)
    def YAMLInitialExchanges2Module(self, exchanges):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialExchanges2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialExchanges2Module 
        transfers EXCHANGE definitions from the InitialPhreeqc instance to the reaction-module 
        workers. 
        @param exchanges A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialExchanges2Module"
        node["exchanges"] = self.GetList(exchanges)
        self.yaml_doc.append(node)
    def YAMLInitialSurfaces2Module(self, surfaces):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialSurfaces2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialSurfaces2Module transfers 
        SURFACE definitions from the InitialPhreeqc instance to the reaction-module workers.
        @param surfaces A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialSurfaces2Module"
        node["surfaces"] = self.GetList(surfaces)
        self.yaml_doc.append(node)
    def YAMLInitialGasPhases2Module(self, gas_phases):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialGasPhases2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialGasPhases2Module 
        transfers GAS_PHASE definitions from the InitialPhreeqc instance to the 
        reaction-module workers.
        @param gas_phases A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialGasPhases2Module"
        node["gas_phases"] = self.GetList(gas_phases)
        self.yaml_doc.append(node)
    def YAMLInitialSolidSolutions2Module(self, solid_solutions):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialSolidSolutions2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialSolidSolutions2Module 
        transfers SOLID_SOLUTIONS definitions from the InitialPhreeqc instance to the 
        reaction-module workers.
        @param solid_solutions A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialSolidSolutions2Module"
        node["solid_solutions"] = self.GetList(solid_solutions)
        self.yaml_doc.append(node)
    def YAMLInitialKinetics2Module(self, kinetics):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialKinetics2Module.

        When the YAML document is written to file it can be processed by the method
        InitializeYAML to initialize a PhreeqcRM instance. InitialKinetics2Module 
        transfers KINETICS definitions from the InitialPhreeqc instance to the 
        reaction-module workers.
        @param kinetics A list, tuple, or numpy array of index numbers that is dimensioned nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "InitialKinetics2Module"
        node["kinetics"] = self.GetList(kinetics)
        self.yaml_doc.append(node)
    def YAMLInitialPhreeqc2Module(self, initial_conditions1):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. InitialPhreeqc2Module transfers solutions and reactants from 
        the InitialPhreeqc instance to the reaction-module workers.
        Initial_conditions1 is used to select initial conditions, including solutions and reactants,
        for each cell of the model, without mixing. Initial_conditions1 is dimensioned 7 times nxyz, 
        where nxyz is the number of grid cells in the user's model. The dimension of 7 refers to 
        solutions and reactants in the following order: (0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, 
        (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE, (5) SOLID_SOLUTIONS, and (6) KINETICS.
        The definition initial_solution1[3*nxyz + 99] = 2, indicates that cell 99 (0 based) contains 
        the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance.
        @param initial_conditions1 A list, tuple, or numpy array of solution and reactant index numbers that refer to
        definitions in the InitialPhreeqc instance. Negative values are ignored, resulting in no 
        definition of that entity for that cell.
        """
        node = dict()
        node["key"] = "InitialPhreeqc2Module"
        node["ic"] = self.GetList(initial_conditions1)
        self.yaml_doc.append(node)
    def YAMLInitialPhreeqc2Module_mix(self, ic1,  ic2,  f1):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqc2Module with three arguments.
        
        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. InitialPhreeqc2Module transfers solutions and reactants from the 
        InitialPhreeqc instance to the reaction-module workers, possibly with mixing.
        In its simplest form, initial_conditions1 is used to select initial conditions, 
        including solutions and reactants, for each cell of the model, without mixing.
        Initial_conditions1 is dimensioned 7 times nxyz, where nxyz is the number of grid cells 
        in the user's model. The dimension of 7 refers to solutions and reactants in the following 
        order: (0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
        (5) SOLID_SOLUTIONS, and (6) KINETICS.The definition initial_solution1[3*nxyz + 99] = 2, 
        indicates that cell 99 (0 based) contains the SURFACE definition (index 3) defined by 
        SURFACE 2 in the InitialPhreeqc instance (either by RunFile or RunString).

        It is also possible to mix solutions and reactants to obtain the initial conditions for cells. 
        For mixing, initials_conditions2 contains numbers for a second entity that mixes with the 
        entity defined in initial_conditions1. Fraction1 contains the mixing fraction for 
        initial_conditions1, whereas (1 - fraction1) is the mixing fraction for initial_conditions2.
        The definitions initial_conditions1[3*nxyz + 99] = 2, initial_conditions2[3*nxyz + 99] = 3,
        fraction1[3*nxyz + 99] = 0.25 indicates that cell 99 (0 based) contains a mixture of 
        0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface compositions have been defined in the 
        InitialPhreeqc instance. If the user number in initial_conditions2 is negative, no mixing occurs.
        @param initial_conditions1 A list, tuple, or numpy array of solution and reactant index numbers that refer to
        definitions in the InitialPhreeqc instance.
        @param initial_conditions2  A list, tuple, or numpy array of solution and reactant index numbers that refer to
        definitions in the InitialPhreeqc instance. Nonnegative values of @a initial_conditions2 
        result in mixing with the entities defined in @a initial_conditions1. Negative values result 
        in no mixing.
        @param fraction1 A list, tuple, or numpy array of fraction of initial_conditions1 that mixes with (1 - fraction1)
        of initial_conditions2.
        """
        node = dict()
        node["key"] = "InitialPhreeqc2Module_mix"
        node["ic1"] = self.GetList(ic1)
        node["ic2"] = self.GetList(ic2)
        node["f1"] = self.GetList(f1)
        self.yaml_doc.append(node)
    def YAMLInitialPhreeqcCell2Module(self, n, cell_numbers):
        """
        Inserts data into the YAML document for the PhreeqcRM method InitialPhreeqcCell2Module.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. InitialPhreeqcCell2Module uses a cell numbered n in the 
        InitialPhreeqc instance to populate a series of transport cells. All reactants with the 
        number n are transferred along with the solution. If MIX n exists, it is used for the 
        definition of the solution. If n is negative, n is redefined to be the largest solution 
        or MIX number in the InitialPhreeqc instance. All reactants for each cell in the list 
        cell_numbers are removed before the cell definition is copied from the InitialPhreeqc 
        instance to the workers.
        @param n Number that refers to a solution or MIX and associated reactants in the InitialPhreeqc 
        instance.
        @param cell_numbers A list, tuple, or numpy array of grid-cell numbers (user's grid-cell numbering system) that
        will be populated with cell n from the InitialPhreeqc instance.
        """
        node = dict()
        node["key"] = "InitialPhreeqcCell2Module"
        node["n"] = n
        node["cell_numbers"] = self.GetList(cell_numbers)
        self.yaml_doc.append(node)
    def YAMLLoadDatabase(self, database):
        """
        Inserts data into the YAML document for the PhreeqcRM method LoadDatabase.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. LoadDatabase loads a database for all IPhreeqc 
        instances--workers, InitialPhreeqc, and Utility. All definitions of the reaction module 
        are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
        @param database String containing the database name.
        """
        node = dict()
        node["key"] = "LoadDatabase"
        node["database"] = database
        self.yaml_doc.append(node)
    def YAMLLogMessage(self, str):
        """
        Inserts data into the YAML document for the PhreeqcRM method LogMessage.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. LogMessage prints a message to the log file.
        @param str              String to be printed.
        """
        node = dict()
        node["key"] = "LogMessage"
        node["str"] = str
        self.yaml_doc.append(node)
    def YAMLOpenFiles(self):
        """
        Inserts data into the YAML document for the PhreeqcRM method OpenFiles.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. OpenFiles opens the output and log files. Files are named 
        prefix.chem.txt and prefix.log.txt based on the prefix defined by YAMLSetFilePrefix.
        """
        node = dict()
        node["key"] = "OpenFiles"
        self.yaml_doc.append(node)
    def YAMLOutputMessage(self, str):
        """
        Inserts data into the YAML document for the PhreeqcRM method OutputMessage.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. OutputMessage prints a message to the output file.
        @param str String to be printed.
        """
        node = dict()
        node["key"] = "OutputMessage"
        node["str"] = str
        self.yaml_doc.append(node)
    def YAMLRunCells(self):
        """
        Inserts data into the YAML document for the PhreeqcRM method RunCells.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. RunCells runs reactions for all cells in the reaction module.
        During initialization, RunCells can be used to equilibrate each solution with all reactants 
        in a cell while using a time step of zero (@ref YAMLSetTimeStep) to avoid kinetic reactions.
        Other properties that may need to be initialized before RunCells is invoked
        include porosity (YAMLSetPorosity), saturation (YAMLSetSaturationUser),
        temperature (YAMLSetTemperature), and pressure (YAMLSetPressure).
        """
        node = dict()
        node["key"] = "RunCells"
        self.yaml_doc.append(node)
    def YAMLRunFile(self, workers, initial_phreeqc, utility, chemistry_name):
        """
        Inserts data into the YAML document for the PhreeqcRM method RunFile.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. RunFile runs a PHREEQC input file. The first three 
        arguments determine which IPhreeqc instances will run the file--the workers, the 
        InitialPhreeqc instance, and (or) the Utility instance. Input files that modify the thermodynamic 
        database should be run by all three sets of instances. Files with SELECTED_OUTPUT definitions 
        that will be used during the time-stepping loop need to be run by the workers. Files that contain 
        initial conditions or boundary conditions should be run by the InitialPhreeqc instance.
        @param workers True, the workers will run the file; False, the workers will not run the file.
        @param initial_phreeqc True, the InitialPhreeqc instance will run the file; False, the 
        InitialPhreeqc will not run the file.
        @param utility True, the Utility instance will run the file; False, the Utility instance 
        will not run the file.
        @param chemistry_name Name of the file to run.
        """
        node = dict()
        node["key"] = "RunFile"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["chemistry_name"] = chemistry_name
        self.yaml_doc.append(node)
    def YAMLRunString(self, workers, initial_phreeqc, utility, input_string):
        """
        Inserts data into the YAML document for the PhreeqcRM method RunString.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. RunString runs a PHREEQC input string. The first three 
        arguments determine which IPhreeqc instances will run the string--the workers, the 
        InitialPhreeqc instance, and (or) the Utility instance. Input strings that modify the 
        thermodynamic database should be run by all three sets of instances. Strings with 
        SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
        be run by the workers. Strings that contain initial conditions or boundary conditions 
        should be run by the InitialPhreeqc instance.
        @param workers True, the workers will run the string; False, the workers will not run the string.
        @param initial_phreeqc True, the InitialPhreeqc instance will run the string; False, the 
        InitialPhreeqc will not run the string.
        @param utility True, the Utility instance will run the string; False, the Utility instance 
        will not run the string.
        @param input_string String containing PHREEQC input.
        """
        node = dict()
        node["key"] = "RunString"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["input_string"] = input_string
        self.yaml_doc.append(node)
    def YAMLScreenMessage(self, str):
        """
        Inserts data into the YAML document for the PhreeqcRM method ScreenMessage.

        When the YAML document is written to file it can be processed by the method 
        InitializeYAML to initialize a PhreeqcRM instance. ScreenMessage prints a message 
        to the screen.
        @param str String to be printed.
        """
        node = dict()
        node["key"] = "ScreenMessage"
        node["str"] = str
        self.yaml_doc.append(node)
    def YAMLSetComponentH2O(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetComponentH2O.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetComponentH2O selects whether to include H2O in the 
        component list. The concentrations of H and O must be known accurately (8 to 10 significant 
        digits) for the numerical method of PHREEQC to produce accurate pH and pe values.
        Because most of the H and O are in the water species, it may be more robust (require less 
        accuracy in transport) to transport the excess H and O (the H and O not in water) and water.
        The default setting (true) is to include water, excess H, and excess O as components.
        A setting of false will include total H and total O as components. YAMLSetComponentH2O must 
        be called before @ref YAMLFindComponents.
        @param tf True (default), excess H, excess O, and water are included in the component list;
        False, total H and O are included in the component list.
        """
        node = dict()
        node["key"] = "SetComponentH2O"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetConcentrations(self, c):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetConcentrations.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. The only way to use this method is to have pre-calculated 
        PHREEQC solution concentrations, which is not common. Concentrations are normally initialized
        with YAMLInitialPhreeqc2Module, YAMLInitialPhreeqcCell2Module, or the methods 
        YAMLInitialSolutions2Module, YAMLInitialEquilibriumPhases2Module, and other similar methods
        for other reactants.
        @param c A list, tuple, or numpy array of component concentrations. Size of vector is ncomps times nxyz,
        where ncomps is the number of components as determined by FindComponents or GetComponentCount and
        nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetConcentrations"
        node["c"] = self.GetList(c)
        self.yaml_doc.append(node)
    def YAMLSetCurrentSelectedOutputUserNumber(self, n_user):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetCurrentSelectedOutputUserNumber.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetCurrentSelectedOutputUserNumber selects the current selected 
        output by user number. The user may define multiple SELECTED_OUTPUT data blocks for the workers. 
        A user number is specified for each data block. The value of the argument @a n_user selects 
        which of the SELECTED_OUTPUT definitions will be used for selected-output operations.
        @param n_user User number of the SELECTED_OUTPUT data block that is to be used.
        """
        node = dict()
        node["key"] = "SetCurrentSelectedOutputUserNumber"
        node["n_user"] = n_user
        self.yaml_doc.append(node)
    def YAMLSetDensityUser(self, density):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetDensityUser.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetDensityUser sets the density for each reaction cell. These 
        density values are used when converting from transported mass-fraction concentrations 
        (YAMLSetUnitsSolution) to produce per liter concentrations during a call to SetConcentrations
        when converting from reaction-cell concentrations to transport concentrations,
        if UseSolutionDensityVolume is set to false.
        @param density A list, tuple, or numpy array of densities. Size of vector is nxyz, where nxyz is the number
        of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetDensityUser"
        node["density"] = self.GetList(density)
        self.yaml_doc.append(node)
    def YAMLSetDumpFileName(self, dump_name):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetDumpFileName.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetDumpFileName sets the name of the dump file. It is the
        name used by the method DumpModule.
        @param dump_name Name of dump file.
        """
        node = dict()
        node["key"] = "SetDumpFileName"
        node["dump_name"] = dump_name
        self.yaml_doc.append(node)
    def YAMLSetErrorHandlerMode(self, mode):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetErrorHandlerMode.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetErrorHandlerMode sets the action to be taken when the 
        reaction module encounters an error. ptions are 0, return to calling program with an error 
        return code (default); 1, throw an exception, in C++, the exception can be caught, for C and 
        Fortran, the program will exit; or 2, attempt to exit gracefully.
        @param mode Error handling mode: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetErrorHandlerMode"
        node["mode"] = mode
        self.yaml_doc.append(node)
    def YAMLSetErrorOn(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetErrorOn.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetErrorOn sets the property that controls whether error 
        messages are generated and displayed. Messages include PHREEQC "ERROR" messages, and
        any messages written with the method ErrorMessage.
        @param tf True, enable error messages; False, disable error messages. Default is true.
        """
        node = dict()
        node["key"] = "SetErrorOn"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetFilePrefix(self, prefix):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetFilePrefix.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetFilePrefix sets the prefix for the output (prefix.chem.txt) 
        and log (prefix.log.txt) files. These files are opened by the method OpenFiles.
        @param prefix Prefix used when opening the output and log files.
        """
        node = dict()
        node["key"] = "SetFilePrefix"
        node["prefix"] = prefix
        self.yaml_doc.append(node)
    def YAMLSetGasCompMoles(self,  gas_moles):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetGasCompMoles.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetGasCompMoles transfers moles of gas components from
        the vector given in the argument list (gas_moles) to each reaction cell.
        @param  gas_moles A list, tuple, or numpy array of moles of gas components. Dimension of the vector is  
        ngas_comps times nxyz, where, ngas_comps is the result of GetGasComponentsCount,
        and nxyz is the number of user grid cells. If the number of moles is set to a negative number, 
        the gas component will not be defined for the GAS_PHASE of the reaction cell.
        """
        node = dict()
        node["key"] = "SetGasCompMoles"
        node["gas_moles"] = self.GetList(gas_moles)
        self.yaml_doc.append(node)
    def YAMLSetGasPhaseVolume(self,  gas_volume):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetGasPhaseVolume.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetGasPhaseVolume transfers volumes of gas phases from
        the vector given in the argument list (gas_volume) to each reaction cell. The gas-phase volume 
        affects the gas-component pressures calculated for fixed-volume gas phases. If a gas-phase volume 
        is defined with this methood for a GAS_PHASE in a cell, the gas phase is forced to be a 
        fixed-volume gas phase.
        @param  gas_volume A list, tuple, or numpy array of volumes for each gas phase. Dimension of the vector is nxyz, where 
        nxyz is the number of user grid cells. If the volume is set to a negative number for a cell, 
        the gas-phase volume for that cell is not changed.
        """
        node = dict()
        node["key"] = "SetGasPhaseVolume"
        node["gas_volume"] = self.GetList(gas_volume)
        self.yaml_doc.append(node)
    def YAMLSetGridCellCount(self, count):
        """
        Inserts data into the YAML document to define the number of cells in the user's model.

        Once the YAML document is written, the number of model cells is used to initialize the
        BMIPhreeqcRM instance if the default constructor (the constructor with no arguments)
        was used. Alternatively, the value can be extracted with the method GetGridCellCountYAML
        and used in the constructor with two arguments. SetGridCellCount will be ignored once the 
        BMIPhreeqcRM instance is initialized.
        @param n Number of cells for the PhreeqcRM instance. 
        """
        node = dict()
        node["key"] = "SetGridCellCount"
        node["count"] = count
        self.yaml_doc.append(node)
    def YAMLSetNthSelectedOutput(self, n):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetNthSelectedOutput.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetNthSelectedOutput specifies the current selected output 
        by sequence number. The user may define multiple SELECTED_OUTPUT data blocks for the workers. 
        A user number is specified for each data block, and the blocks are stored in user-number order. 
        The value of the argument n selects the sequence number of the SELECTED_OUTPUT definition that 
        will be used for selected-output operations.
        @param n Sequence number of the SELECTED_OUTPUT data block that is to be used.
        """
        node = dict()
        node["key"] = "SetNthSelectedOutput"
        node["n"] = n
        self.yaml_doc.append(node)
    def YAMLSetPartitionUZSolids(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetPartitionUZSolids.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetPartitionUZSolids sets the property for partitioning solids 
        between the saturated and unsaturated parts of a partially saturated cell. The option is intended 
        to be used by saturated-only flow codes that allow a variable water table. The value has meaning 
        only when saturations less than 1.0 are encountered. The partially saturated cells may have a 
        small water-to-rock ratio that causes reactions to proceed differently relative to fully saturated 
        cells. By setting SetPartitionUZSolids to true, the amounts of solids and gases are partioned 
        according to the saturation. If a cell has a saturation of 0.5, then the water interacts with 
        only half of the solids and gases; the other half is unreactive until the water table rises. 
        As the saturation in a cell varies, solids and gases are transferred between the saturated and 
        unsaturated (unreactive) reservoirs of the cell. Unsaturated-zone flow and transport codes will 
        probably use the default (false), which assumes all gases and solids are reactive regardless of 
        saturation.
        @param tf True, the fraction of solids and gases available for reaction is equal to the 
        saturation; False (default), all solids and gases are reactive regardless of saturation.
        """
        node = dict()
        node["key"] = "SetPartitionUZSolids"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetPorosity(self,  por):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetPorosity.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetPorosity sets the porosity for each reaction cell.
        The volume of water in a reaction cell is the product of porosity, saturation
        (SetSaturationUser), and representative volume (SetRepresentativeVolume).
        @param por A list, tuple, or numpy array of porosities, unitless. Size of vector is nxyz, where 
        nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetPorosity"
        node["por"] = self.GetList(por)
        self.yaml_doc.append(node)
    def YAMLSetPressure(self,  p):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetPressure.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetPressure sets the pressure for each reaction cell. 
        Pressure effects are considered only in three of the databases distributed with PhreeqcRM: 
        phreeqc.dat, Amm.dat, and pitzer.dat.
        @param p A list, tuple, or numpy array of pressures, in atm. Size of vector is nxyz, where nxyz is the 
        number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetPressure"
        node["p"] = self.GetList(p)
        self.yaml_doc.append(node)
    def YAMLSetPrintChemistryMask(self, cell_mask):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetPrintChemistryMask.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetPrintChemistryMask enables or disables detailed output 
        for each reaction cell. Printing for a reaction cell will occur only when the
        printing is enabled with SetPrintChemistryOn and the cell_mask value is 1.
        @param cell_mask A list, tuple, or numpy array of integers. Size of vector is nxyz, where nxyz is the number
        of grid cells in the user's model. A value of 0 will disable printing detailed output for 
        the cell; a value of 1 will enable printing detailed output for a cell.
        """
        node = dict()
        node["key"] = "SetPrintChemistryMask"
        node["cell_mask"] = self.GetList(cell_mask)
        self.yaml_doc.append(node)
    def YAMLSetPrintChemistryOn(self, workers, initial_phreeqc, utility):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetPrintChemistryOn.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetPrintChemistryOn sets the property that enables or disables 
        printing detailed output from reaction calculations to the output file for a set of cells defined 
        by SetPrintChemistryMask. The detailed output prints all of the output typical of a PHREEQC 
        reaction calculation, which includes solution descriptions and the compositions of all other 
        reactants. The output can be several hundred lines per cell, which can lead to a very
        large output file (prefix.chem.txt opened by the method OpenFiles).For the worker instances, 
        the output can be limited to a set of cells (method SetPrintChemistryMask) and, in general, the
        amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
        (applied by using methods RunFile or RunString). Printing the detailed output for the workers is 
        generally used only for debugging, and PhreeqcRM will run significantly faster when printing 
        detailed output for the workers is disabled.

        @param workers True, enable detailed printing in the worker instances; False, disable detailed 
        printing in the worker instances.
        @param initial_phreeqc True, enable detailed printing in the InitialPhreeqc instance; False, 
        disable detailed printing in the InitialPhreeqc instance.
        @param utility True, enable detailed printing in the Utility instance; False, disable detailed 
        printing in the Utility instance.
        """
        node = dict()
        node["key"] = "SetPrintChemistryOn"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        self.yaml_doc.append(node)
    def YAMLSetRebalanceByCell(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetRebalanceByCell.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetRebalanceByCell sets the load-balancing algorithm.
        PhreeqcRM attempts to rebalance the load of each thread or process such that each thread 
        or process takes the same amount of time to run its part of a RunCells calculation. 
        Two algorithms are available; one uses individual times for each cell and accounts for cells 
        that were not run because saturation was zero (default), and the other assigns an average time 
        to all cells. The methods are similar, but limited testing indicates the default method performs 
        better.
        @param tf True, indicates individual cell times are used in rebalancing (default); False, 
        indicates average times are used in rebalancing.
        """
        node = dict()
        node["key"] = "SetRebalanceByCell"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetRebalanceFraction(self, f):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetRebalanceFraction.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetRebalanceFraction sets the fraction of cells that are 
        transferred among threads or processes when rebalancing. PhreeqcRM attempts to rebalance 
        the load of each thread or process such that each thread or process takes the same amount of 
        time to run its part of a RunCells calculation. The rebalancing transfers cell calculations 
        among threads or processes to try to achieve an optimum balance. SetRebalanceFraction
        adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
        determine the actual number of cell transfers. A value of zero eliminates load rebalancing. 
        A value less than 1.0 is suggested to slow the approach to the optimum cell distribution and 
        avoid possible oscillations when too many cells are transferred at one iteration, requiring 
        reverse transfers at the next iteration. Default is 0.5.
        @param f Fraction from 0.0 to 1.0.
        """
        node = dict()
        node["key"] = "SetRebalanceFraction"
        node["f"] = f
        self.yaml_doc.append(node)
    def YAMLSetRepresentativeVolume(self,  rv):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetRepresentativeVolume.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetRepresentativeVolume sets the representative volume of 
        each reaction cell. By default the representative volume of each reaction cell is 1 liter.
        The volume of water in a reaction cell is determined by the product of the representative 
        volume, the porosity (SetPorosity), and the saturation (SetSaturationUser). The numerical method 
        of PHREEQC is more robust if the water volume for a reaction cell is within a couple orders 
        of magnitude of 1.0. Small water volumes caused by small porosities and (or) small saturations 
        (and (or) small representative volumes) may cause non-convergence of the numerical method.
        In these cases, a larger representative volume may help. Note that increasing the representative 
        volume also increases the number of moles of the reactants in the reaction cell (minerals, 
        surfaces, exchangers, and others), which are defined as moles per representative volume.
        SetRepresentativeVolume should be called before initial conditions are defined for the 
        reaction cells.
        @param rv A list, tuple, or numpy array of representative volumes, in liters. 
        Size of array is nxyz, where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetRepresentativeVolume"
        node["rv"] = self.GetList(rv)
        self.yaml_doc.append(node)
    def YAMLSetSaturationUser(self,  sat):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetSaturationUser.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetSaturationUser sets the saturation of each reaction cell. 
        Saturation is a fraction ranging from 0 to 1. The volume of water in a cell is the product 
        of porosity (SetPorosity), saturation (SetSaturationUser), and representative volume 
        (SetRepresentativeVolume). As a result of a reaction calculation, solution properties (density 
        and volume) will change; the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar 
        volume data to calculate these changes. The methods GetDensity, GetSolutionVolume, and 
        GetSaturation can be used to account for these changes in the succeeding transport calculation.
        @param sat A list, tuple, or numpy array of saturations, unitless. Size of vector is nxyz,
        where nxyz is the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetSaturationUser"
        node["sat"] = self.GetList(sat)
        self.yaml_doc.append(node)
    def YAMLSetScreenOn(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetScreenOn.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetScreenOn sets the property that controls whether messages 
        are written to the screen. Messages include information about rebalancing during RunCells, and
        any messages written with ScreenMessage.
        @param tf True, enable screen messages; False, disable screen messages. Default is true.
        """
        node = dict()
        node["key"] = "SetScreenOn"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetSelectedOutputOn(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetSelectedOutputOn.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetSelectedOutputOn sets the property that controls whether 
        selected-output results are available to be retrieved with GetSelectedOutput. 
        @param tf True indicates that selected-output results will be accumulated during RunCells 
        and can be retrieved with GetSelectedOutput; False indicates that selected-output results 
        will not be accumulated during RunCells.
        """
        node = dict()
        node["key"] = "SetSelectedOutputOn"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLSetSpeciesSaveOn(self, save_on):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetSpeciesSaveOn.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetSpeciesSaveOn sets the value of the species-save property.
        This method enables or disables use of PhreeqcRM with multicomponent-diffusion transport 
        calculations. By default, concentrations of aqueous species are not saved. Setting the 
        species-save property to @a true allows aqueous species concentrations to be retrieved
        with GetSpeciesConcentrations, and solution compositions to be set with 
        SpeciesConcentrations2Module. SetSpeciesSaveOn must be called before calls to FindComponents.
        @param save_on True indicates species concentrations are saved; False indicates species 
        concentrations are not saved.
        """
        node = dict()
        node["key"] = "SetSpeciesSaveOn"
        node["save_on"] = save_on
        self.yaml_doc.append(node)
    def YAMLSetTemperature(self,  t):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetTemperature.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetTemperature sets the temperature for each reaction cell. 
        If SetTemperature is not called, worker solutions will have temperatures as defined by initial 
        conditions (InitialPhreeqc2Module and InitialPhreeqcCell2Module).
        @param t A list, tuple, or numpy array of temperatures, in degrees C. Size of vector is @a nxyz, where @a nxyz is 
        the number of grid cells in the user's model.
        """
        node = dict()
        node["key"] = "SetTemperature"
        node["t"] = self.GetList(t)
        self.yaml_doc.append(node)
    def YAMLSetTime(self, time):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetTime.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetTime sets current simulation time for the reaction module.
        @param time Current simulation time, in seconds.
        """
        node = dict()
        node["key"] = "SetTime"
        node["time"] = time
        self.yaml_doc.append(node)
    def YAMLSetTimeConversion(self, conv_factor):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetTimeConversion.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetTimeConversion sets a factor to convert from seconds to 
        user time units. Factor times seconds produces user time units, which are used in some PhreeqcRM 
        printing.
        @param conv_factor Factor to convert seconds to user time units.
        """
        node = dict()
        node["key"] = "SetTimeConversion"
        node["conv_factor"] = conv_factor
        self.yaml_doc.append(node)
    def YAMLSetTimeStep(self, time_step):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetTimeStep.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetTimeStep sets current time step for the reaction module. 
        This is the length of time over which kinetic reactions are integrated.
        @param time_step Time step, in seconds.
        """
        node = dict()
        node["key"] = "SetTimeStep"
        node["time_step"] = time_step
        self.yaml_doc.append(node)
    def YAMLSetUnitsExchange(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsExchange.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsExchange sets input units for exchangers.
        In PHREEQC input, exchangers are defined by moles of exchange sites (Mp). SetUnitsExchange 
        specifies how the number of moles of exchange sites in a reaction cell (Mc) is calculated from 
        the input value (Mp). Options are 0, Mp is mol/L of RV (default), Mc = a Mp*RV, where RV is 
        the representative volume (SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, 
        Mc = Mp*P*RV, where P is porosity (SetPorosity); or 2, Mp is mol/L of rock in the RV,  
        Mc = Mp*(1-P)*RV.

        If a single EXCHANGE definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of exchangers will 
        be the same regardless of porosity. For option 1, the number of moles of exchangers will be vary 
        directly with porosity and inversely with rock volume. For option 2, the number of moles of 
        exchangers will vary directly with rock volume and inversely with porosity.

        @param option Units option for exchangers: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsExchange"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsGasPhase(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsGasPhase.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsGasPhase sets input units for gas phases.
        In PHREEQC input, gas phases are defined by moles of component gases (Mp). SetUnitsGasPhase 
        specifies how the number of moles of component gases in a reaction cell (Mc)
        is calculated from the input value (Mp). Options are 0, Mp is mol/L of RV (default), 
        Mc = Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
        1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (SetPorosity); or
        2, Mp is mol/L of rock in the RV,  Mc = Mp*(1-P)*RV.

        If a single GAS_PHASE definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of a gas component 
        will be the same regardless of porosity. For option 1, the number of moles of a gas component 
        will vary directly with porosity and inversely with rock volume. For option 2, the number of 
        moles of a gas component will vary directly with rock volume and inversely with porosity.

        @param option Units option for gas phases: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsGasPhase"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsKinetics(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsKinetics.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsKinetics sets input units for kinetic reactants.
        In PHREEQC input, kinetics are defined by moles of kinetic reactants (Mp). SetUnitsKinetics 
        specifies how the number of moles of kinetic reactants in a reaction cell (Mc) is calculated 
        from the input value (Mp). Options are 0, Mp is mol/L of RV (default), Mc = @a Mp*RV, where 
        RV is the representative volume (SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, 
        Mc = Mp*P*RV, where P is porosity (SetPorosity); or 2,  is mol/L of rock in the RV,  
        Mc = Mp*(1-P)*RV.

        If a single KINETICS definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of kinetic 
        reactants will be the same regardless of porosity. For option 1, the number of moles of kinetic 
        reactants will vary directly with porosity and inversely with rock volume. For option 2, the 
        number of moles of kinetic reactants will vary directly with rock volume and inversely with 
        porosity.

        Note that the volume of water in a cell in the reaction module is equal to the product of
        porosity (SetPorosity), the saturation (SetSaturationUser), and representative volume 
        (SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the RATES
        definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
        water, often by calculating the rate of reaction per liter of water and multiplying by the volume
        of water (Basic function SOLN_VOL).

        Rates that depend on surface area of solids, are not dependent on the volume of water. However, 
        it is important to get the correct surface area for the kinetic reaction. To scale the surface 
        area with the number of moles, the specific area (m^2 per mole of reactant) can be defined as a 
        parameter (KINETICS; -parm), which is multiplied by the number of moles of reactant (Basic 
        function M) in RATES to obtain the surface area.

        @param option Units option for kinetic reactants: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsKinetics"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsPPassemblage(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsPPassemblage.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsPPassemblage sets input units for pure phase 
        assemblages (equilibrium phases). In PHREEQC input, equilibrium phases are defined by moles 
        of each phase (Mp). SetUnitsPPassemblage specifies how the number of moles of phases in a 
        reaction cell (Mc) is calculated from the input value (Mp). Options are 0, @a Mp is mol/L of 
        RV (default), Mc = Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
        1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (SetPorosity); or
        2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P)*RV.

        If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of a mineral will 
        be the same regardless of porosity. For option 1, the number of moles of a mineral will vary 
        directly with porosity and inversely with rock volume. For option 2, the number of moles of a 
        mineral will vary directly with rock volume and inversely with porosity.

        @param option Units option for equilibrium phases: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsPPassemblage"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsSolution(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsSolution.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsSolution sets solution concentration units used by 
        the transport model. Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs. PHREEQC defines 
        solutions by the number of moles of each element in the solution. To convert from mg/L to moles
        of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
        multiplied by the solution volume, which is the product of porosity (SetPorosity), saturation 
        (SetSaturationUser), and representative volume (SetRepresentativeVolume). To convert from mol/L 
        to moles of element in the representative volume of a reaction cell, mol/L is multiplied by the 
        solution volume. To convert from mass fraction to moles of element in the representative volume 
        of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density (SetDensityUser) and
        multiplied by the solution volume.

        To convert from moles of element in the representative volume of a reaction cell to mg/L, the 
        number of moles of an element is divided by the solution volume resulting in mol/L, and then 
        converted to mg/L. To convert from moles of element in a cell to mol/L,  the number of moles of 
        an element is divided by the solution volume resulting in mol/L. To convert from moles
        of element in a cell to mass fraction, the number of moles of an element is converted to 
        kg and divided by the total mass of the solution.

        Two options are available for the volume and mass of solution that are used in converting 
        to transport concentrations: (1) the volume and mass of solution are calculated by PHREEQC, 
        or (2) the volume of solution is the product of porosity (SetPorosity), saturation 
        (SetSaturationUser), and representative volume (SetRepresentativeVolume), and the mass of solution 
        is volume times density as defined by SetDensityUser. Which option is used is determined by 
        UseSolutionDensityVolume.

        @param option Units option for solutions: 1, 2, or 3, default is 1, mg/L.
        """
        node = dict()
        node["key"] = "SetUnitsSolution"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsSSassemblage(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsSSassemblage.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsSSassemblage sets input units for solid-solution 
        assemblages. In PHREEQC, solid solutions are defined by moles of each component (Mp).
        SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a 
        reaction cell (Mc) is calculated from the input value (Mp). Options are 0, Mp is mol/L of 
        RV (default), Mc = Mp*RV, where RV is the representative volume (SetRepresentativeVolume);
        1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity (SetPorosity); or
        2, Mp is mol/L of rock in the RV,  Mc = Mp*(1-P)*RV.

        If a single SOLID_SOLUTION definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of a solid-solution 
        component will be the same regardless of porosity. For option 1, the number of moles of a 
        solid-solution component will vary directly with porosity and inversely with rock volume.
        For option 2, the number of moles of a solid-solution component will vary directly with rock 
        volume and inversely with porosity.

        @param option Units option for solid solutions: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsSSassemblage"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSetUnitsSurface(self, option):
        """
        Inserts data into the YAML document for the PhreeqcRM method SetUnitsSurface.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SetUnitsSurface sets input units for surfaces. In PHREEQC input, 
        surfaces are defined by moles of surface sites (Mp). SetUnitsSurface specifies how the number 
        of moles of surface sites in a reaction cell (@a Mc) is calculated from the input value (@a Mp).
        Options are 0, Mp is mol/L of RV (default), Mc = Mp*RV, where RV is the representative volume 
        (SetRepresentativeVolume); 1, Mp is mol/L of water in the RV, Mc = Mp*P*RV, where P is porosity 
        (SetPorosity); or 2, Mp is mol/L of rock in the RV, Mc = Mp*(1-P)*RV.

        If a single SURFACE definition is used for cells with different initial porosity,
	    the three options scale quite differently. For option 0, the number of moles of surface sites 
        will be the same regardless of porosity. For option 1, the number of moles of surface sites 
        will vary directly with porosity and inversely with rock volume. For option 2, the number of 
        moles of surface sites will vary directly with rock volume and inversely with porosity.
        @param option Units option for surfaces: 0, 1, or 2.
        """
        node = dict()
        node["key"] = "SetUnitsSurface"
        node["option"] = option
        self.yaml_doc.append(node)
    def YAMLSpeciesConcentrations2Module(self,  species_conc):
        """
        Inserts data into the YAML document for the PhreeqcRM method SpeciesConcentrations2Module.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. SpeciesConcentrations2Module sets solution concentrations in 
        the reaction cells based on the vector of aqueous species concentrations (@a species_conc).
        This method is intended for use with multicomponent-diffusion transport calculations, and 
        SetSpeciesSaveOn must be set to true.  The list of aqueous species is determined by 
        FindComponents and includes all aqueous species that can be made from the set of components.
        The method determines the total concentration of a component by summing the molarities of the 
        individual species times the stoichiometric coefficient of the element in each species.
        Solution compositions in the reaction cells are updated with these component concentrations.
        Usually, accurate concentrations will not be known to use YAMLSpeciesConcentrations2Module during
        initialization.

        @param species_conc A list, tuple, or numpy array of aqueous species concentrations. Dimension of the array is nspecies 
        times nxyz, where  nspecies is the number of aqueous species, and nxyz is the number of user 
        grid cells. Concentrations are moles per liter.
        """
        node = dict()
        node["key"] = "SpeciesConcentrations2Module"
        node["species_conc"] = self.GetList(species_conc)
        self.yaml_doc.append(node)
    def YAMLStateSave(self, istate):
        """
        Inserts data into the YAML document for the PhreeqcRM method StateSave.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. StateSave saves the state of the chemistry in all model cells, 
        including SOLUTIONs, EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and 
        SURFACEs. Although not generally used, MIXes, REACTIONs, REACTION_PRESSUREs, and 
        REACTION_TEMPERATUREs will be saved for each cell, if they have been defined in the worker 
        IPhreeqc instances. The distribution of cells among the workers and the chemistry of fully or 
        partially unsaturated cells are also saved. The state is saved in memory; use DumpModule to 
        save the state to file. PhreeqcRM can be reset to this state by using StateApply. A state is 
        identified by an integer, and multiple states can be saved.

        @param istate Integer identifying the state that is saved.
        """
        node = dict()
        node["key"] = "StateSave"
        node["istate"] = istate
        self.yaml_doc.append(node)
    def YAMLStateApply(self, istate):
        """
        Inserts data into the YAML document for the PhreeqcRM method StateApply.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. StateApply resets the state of the module to a state previously 
        saved with StateSave. The chemistry of all model cells are reset, including SOLUTIONs,
        EQUILIBRIUM_PHASES, EXCHANGEs, GAS_PHASEs, KINETICS, SOLID_SOLUTIONs, and SURFACEs.
        MIXes, REACTIONs, REACTION_PRESSUREs, and REACTION_TEMPERATUREs will be reset for each cell, 
        if they were defined in the worker IPhreeqc instances at the time the state was saved.
        The distribution of cells among the workers and the chemistry of fully or partially unsaturated 
        cells are also reset to the saved state. The state to be applied is identified by an integer.

        @param istate Integer identifying the state that is to be applied.
        """
        node = dict()
        node["key"] = "StateApply"
        node["istate"] = istate
        self.yaml_doc.append(node)
    def YAMLStateDelete(self, istate):
        """
        Inserts data into the YAML document for the PhreeqcRM method StateDelete.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. StateDelete deletes a state previously saved with StateSave.

        @param istate     Integer identifying the state that is to be deleted.
        """
        node = dict()
        node["key"] = "StateDelete"
        node["istate"] = istate
        self.yaml_doc.append(node)
    def YAMLThreadCount(self, nthreads):
        """
        Inserts data into the YAML document for the method ThreadCount.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance.ThreadCount provides the number of threads to use in 
        OpenMP multiprocessing when used to initialize a BMIPhreeqcRM instance, provided the 
        BMIPhreeqcRM instance was created with the default constructor--the constructor with 
        no arguments.
        """
        node = dict()
        node["key"] = "ThreadCount"
        node["nthreads"] = nthreads
        self.yaml_doc.append(node)
    def YAMLUseSolutionDensityVolume(self, tf):
        """
        Inserts data into the YAML document for the PhreeqcRM method UseSolutionDensityVolume.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. UseSolutionDensityVolume determines the volume and density 
        to use when converting from the reaction-cell concentrations to transport concentrations 
        (GetConcentrations). Two options are available to convert concentration units: (1) the density 
        and solution volume calculated by PHREEQC are used, or (2) the specified density (SetDensityUser)
        and solution volume are determined by the product of saturation (SetSaturationUser), porosity 
        (SetPorosity), and representative volume (SetRepresentativeVolume). Transport models that 
        consider density-dependent flow will probably use the PHREEQC-calculated density and solution 
        volume (default), whereas transport models that assume constant-density flow will probably use
        specified values of density and solution volume. Only the following databases distributed with 
        PhreeqcRM have molar-volume information needed to accurately calculate density and solution 
        volume: phreeqc.dat, Amm.dat, and pitzer.dat. Density is only used when converting to or from 
        transport units of mass fraction.

        @param tf True indicates that the solution density and volume as calculated by PHREEQC will 
        be used to calculate concentrations. False indicates that the solution density set by SetDensityUser 
        and the volume determined by the product of  SetSaturationUser, SetPorosity, and 
        SetRepresentativeVolume, will be used to calculate concentrations retrieved by GetConcentrations.
        """
        node = dict()
        node["key"] = "UseSolutionDensityVolume"
        node["tf"] = tf
        self.yaml_doc.append(node)
    def YAMLWarningMessage(self, str):
        """
        Inserts data into the YAML document for the PhreeqcRM method WarningMessage.

        When the YAML document is written to file it can be processed by the method InitializeYAML to
        initialize a PhreeqcRM instance. WarningMessage prints a warning message to the screen and 
        the log file.

        @param str String to be printed.
        """
        node = dict()
        node["key"] = "WarningMessage"
        node["str"] = str
        self.yaml_doc.append(node)
    def GetList(self, v):
        if isinstance(v, np.ndarray) and (isinstance(v[0].item(), int) or isinstance(v[0], float)):
            return v.tolist()
        elif isinstance(v, tuple):
            return list(v)
        return v
