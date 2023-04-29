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

    Functions that accept a constant vector reference
    (const std::vector<double>&) should be callable
    with any of the three Python data types, namely
    np.ndarray, list, and tuple.

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
    nxyz = 20
    nthreads = 3
    phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)
    status = phreeqc_rm.SetPorosity(porosity)

    # likewise for integers (const std::vector<int>&)

    Currently functions that take non-const vectors (ie OUT or INOUT)
    must use the SWIG wrapped vector class like this:

    c_dbl_vect = phreeqcrm.DoubleVector(nxyz * len(components))
    status = phreeqc_rm.GetConcentrations(c_dbl_vect)
"""

class YAMLPhreeqcRM(object):
    yaml_doc = dict()

    def Clear(self):
        yaml_doc = " "
    def WriteYAMLDoc(self, file_name):
        file = open(file_name, "w")
        yaml.dump(self.yaml_doc, file, sort_keys=False)
        file.close()
    def YAMLAddOutputVars(option, definition):
        node = dict()
        node["key"] = "AddOutputVars"
        node["definition"] = "definition"
        self.yaml_doc.push_back(node)
    def YAMLCloseFiles(self, void):
        node = dict()
        node["key"] = "CloseFiles"
        self.yaml_doc.push_back(node)
    def YAMLCreateMapping(self, g2c):
        node = dict()
        node["key"] = "CreateMapping"
        node["CreateMapping"] = g2c
        self.yaml_doc.push_back(node)
    def YAMLDumpModule(self, dump_on, append):
        node = dict()
        node["key"] = "DumpModule"
        node["dump_on"] = dump_on
        node["append"] = append
        self.yaml_doc.push_back(node)
    def YAMLFindComponents(self):
        node = dict()
        node["key"] = "FindComponents"
        self.yaml_doc.push_back(node)
    def YAMLInitialSolutions2Module(self, solutions):
        node = dict()
        node["key"] = "InitialSolutions2Module"
        node["solutions"] = solutions
        self.yaml_doc.push_back(node)
    def YAMLInitialEquilibriumPhases2Module(self, equilibrium_phases):
        node = dict()
        node["key"] = "InitialEquilibriumPhases2Module"
        node["equilibrium_phases"] = equilibrium_phases
        self.yaml_doc.push_back(node)
    def YAMLInitialExchanges2Module(self, exchanges):
        node = dict()
        node["key"] = "InitialExchanges2Module"
        node["exchanges"] = exchanges
        self.yaml_doc.push_back(node)
    def YAMLInitialSurfaces2Module(self, surfaces):
        node = dict()
        node["key"] = "InitialSurfaces2Module"
        node["surfaces"] = surfaces
        self.yaml_doc.push_back(node)
    def YAMLInitialGasPhases2Module(self, gas_phases):
        node = dict()
        node["key"] = "InitialGasPhases2Module"
        node["gas_phases"] = gas_phases
        self.yaml_doc.push_back(node)
    def YAMLInitialSolidSolutions2Module(self, solid_solutions):
        node = dict()
        node["key"] = "InitialSolidSolutions2Module"
        node["solid_solutions"] = solid_solutions
        self.yaml_doc.push_back(node)
    def YAMLInitialKinetics2Module(self, kinetics):
        node = dict()
        node["key"] = "InitialKinetics2Module"
        node["kinetics"] = kinetics
        self.yaml_doc.push_back(node)
    def YAMLInitialPhreeqc2Module(self, initial_conditions1):
        node = dict()
        node["key"] = "InitialPhreeqc2Module"
        node["ic"] = initial_conditions1
        self.yaml_doc.push_back(node)
    def YAMLInitialPhreeqc2Module_mix(self, ic1,  ic2,  f1):
        node = dict()
        node["key"] = "YAMLInitialPhreeqc2Module_mix"
        node["ic1"] = ic1
        node["ic2"] = ic2
        node["f1"] = f1
        self.yaml_doc.push_back(node)
    def YAMLInitialPhreeqcCell2Module(self, n,  cell_numbers):
        node = dict()
        node["key"] = "InitialPhreeqcCell2Module"
        node["n"] = n
        node["cell_numbers"] = cell_numbers
        self.yaml_doc.push_back(node)
    def YAMLLoadDatabase(self, database):
        node = dict()
        node["key"] = "LoadDatabase"
        node["database"] = database
        self.yaml_doc.push_back(node)
    def YAMLLogMessage(self, str):
        node = dict()
        node["key"] = "LogMessage"
        node["str"] = str
        self.yaml_doc.push_back(node)
    def YAMLOpenFiles(self):
        node = dict()
        node["key"] = "OpenFiles"
        self.yaml_doc.push_back(node)
    def YAMLOutputMessage(self, str):
        node = dict()
        node["key"] = "OutputMessage"
        node["str"] = str
        self.yaml_doc.push_back(node)
    def YAMLRunCells(self):
        node = dict()
        node["key"] = "RunCells"
        self.yaml_doc.push_back(node)
    def YAMLRunFile(self, workers, initial_phreeqc, utility, chemistry_name):
        node = dict()
        node["key"] = "RunFile"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["chemistry_name"] = chemistry_name
        self.yaml_doc.push_back(node)
    def YAMLRunString(self, workers, initial_phreeqc, utility, input_string):
        node = dict()
        node["key"] = "RunString"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["input_string"] = input_string
        self.yaml_doc.push_back(node)
    def YAMLScreenMessage(self, str):
        node = dict()
        node["key"] = "ScreenMessage"
        node["str"] = str
        self.yaml_doc.push_back(node)
    def YAMLSetComponentH2O(self, tf):
        node = dict()
        node["key"] = "SetComponentH2O"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetConcentrations(self, c):
        node = dict()
        node["key"] = "SetConcentrations"
        node["c"] = c
        self.yaml_doc.push_back(node)
    def YAMLSetCurrentSelectedOutputUserNumber(self, n_user):
        node = dict()
        node["key"] = "SetCurrentSelectedOutputUserNumber"
        node["n_user"] = n_user
        self.yaml_doc.push_back(node)
    def YAMLSetDensity(self, density):
        node = dict()
        node["key"] = "SetDensity"
        node["density"] = density
        self.yaml_doc.push_back(node)
    def YAMLSetDumpFileName(self, dump_name):
        node = dict()
        node["key"] = "SetDumpFileName"
        node["dump_name"] = dump_name
        self.yaml_doc.push_back(node)
    def YAMLSetErrorHandlerMode(self, mode):
        node = dict()
        node["key"] = "SetErrorHandlerMode"
        node["mode"] = mode
        self.yaml_doc.push_back(node)
    def YAMLSetErrorOn(self, tf):
        node = dict()
        node["key"] = "SetErrorOn"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetFilePrefix(self, prefix):
        node = dict()
        node["key"] = "SetFilePrefix"
        node["prefix"] = prefix
        self.yaml_doc.push_back(node)
    def YAMLSetGasCompMoles(self,  gas_moles):
        node = dict()
        node["key"] = "SetGasCompMoles"
        node["gas_moles"] = gas_moles
        self.yaml_doc.push_back(node)
    def YAMLSetGasPhaseVolume(self,  gas_volume):
        node = dict()
        node["key"] = "SetGasPhaseVolume"
        node["gas_volume"] = gas_volume
        self.yaml_doc.push_back(node)
    def YAMLSetGridCellCount(self, n):
        node = dict()
        node["key"] = "SetGridCellCount"
        node["n"] = n
        self.yaml_doc.push_back(node)
    def YAMLSetNthSelectedOutput(self, n):
        node = dict()
        node["key"] = "SetNthSelectedOutput"
        node["n"] = n
        self.yaml_doc.push_back(node)
    def YAMLSetPartitionUZSolids(self, tf):
        node = dict()
        node["key"] = "SetPartitionUZSolids"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetPorosity(self,  por):
        node = dict()
        node["key"] = "SetPorosity"
        node["por"] = por
        self.yaml_doc.push_back(node)
    def YAMLSetPressure(self,  p):
        node = dict()
        node["key"] = "SetPressure"
        node["p"] = p
        self.yaml_doc.push_back(node)
    def YAMLSetPrintChemistryMask(self, cell_mask):
        node = dict()
        node["key"] = "SetPrintChemistryMask"
        node["cell_mask"] = cell_mask
        self.yaml_doc.push_back(node)
    def YAMLSetPrintChemistryOn(self, workers, initial_phreeqc, utility):
        node = dict()
        node["key"] = "SetPrintChemistryOn"
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        self.yaml_doc.push_back(node)
    def YAMLSetRebalanceByCell(self, tf):
        node = dict()
        node["key"] = "SetRebalanceByCell"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetRebalanceFraction(self, f):
        node = dict()
        node["key"] = "SetRebalanceFraction"
        node["f"] = f
        self.yaml_doc.push_back(node)
    def YAMLSetRepresentativeVolume(self,  rv):
        node = dict()
        node["key"] = "SetRepresentativeVolume"
        node["rv"] = rv
        self.yaml_doc.push_back(node)
    def YAMLSetSaturation(self,  sat):
        node = dict()
        node["key"] = "SetSaturation"
        node["sat"] = sat
        self.yaml_doc.push_back(node)
    def YAMLSetScreenOn(self, tf):
        node = dict()
        node["key"] = "SetScreenOn"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetSelectedOutputOn(self, tf):
        node = dict()
        node["key"] = "SetSelectedOutputOn"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLSetSpeciesSaveOn(self, save_on):
        node = dict()
        node["key"] = "SetSpeciesSaveOn"
        node["save_on"] = save_on
        self.yaml_doc.push_back(node)
    def YAMLSetTemperature(self,  t):
        node = dict()
        node["key"] = "SetTemperature"
        node["t"] = t
        self.yaml_doc.push_back(node)
    def YAMLSetTime(self, time):
        node = dict()
        node["key"] = "SetTime"
        node["time"] = time
        self.yaml_doc.push_back(node)
    def YAMLSetTimeConversion(self, conv_factor):
        node = dict()
        node["key"] = "SetTimeConversion"
        node["conv_factor"] = conv_factor
        self.yaml_doc.push_back(node)
    def YAMLSetTimeStep(self, time_step):
        node = dict()
        node["key"] = "SetTimeStep"
        node["time_step"] = time_step
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsExchange(self, option):
        node = dict()
        node["key"] = "SetUnitsExchange"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsGasPhase(self, option):
        node = dict()
        node["key"] = "SetUnitsGasPhase"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsKinetics(self, option):
        node = dict()
        node["key"] = "SetUnitsKinetics"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsPPassemblage(self, option):
        node = dict()
        node["key"] = "SetUnitsPPassemblage"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsSolution(self, option):
        node = dict()
        node["key"] = "SetUnitsSolution"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsSSassemblage(self, option):
        node = dict()
        node["key"] = "SetUnitsSSassemblage"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSetUnitsSurface(self, option):
        node = dict()
        node["key"] = "SetUnitsSurface"
        node["option"] = option
        self.yaml_doc.push_back(node)
    def YAMLSpeciesConcentrations2Module(self,  species_conc):
        node = dict()
        node["key"] = "SpeciesConcentrations2Module"
        node["species_conc"] = species_conc
        self.yaml_doc.push_back(node)
    def YAMLStateSave(self, istate):
        node = dict()
        node["key"] = "StateSave"
        node["istate"] = istate
        self.yaml_doc.push_back(node)
    def YAMLStateApply(self, istate):
        node = dict()
        node["key"] = "StateApply"
        node["istate"] = istate
        self.yaml_doc.push_back(node)
    def YAMLStateDelete(self, istate):
        node = dict()
        node["key"] = "StateDelete"
        node["istate"] = istate
        self.yaml_doc.push_back(node)
    def YAMLThreadCount(self, nthreads):
        node = dict()
        node["key"] = "ThreadCount"
        node["nthreads"] = nthreads
        self.yaml_doc.push_back(node)
    def YAMLUseSolutionDensityVolume(self, tf):
        node = dict()
        node["key"] = "UseSolutionDensityVolume"
        node["tf"] = tf
        self.yaml_doc.push_back(node)
    def YAMLWarningMessage(self, str):
        node = dict()
        node["key"] = "WarningMessage"
        node["str"] = str
        self.yaml_doc.push_back(node)


#nxyz = 40
#yrm = YAMLPhreeqcRM()
#grid2cell = [i for i in range(nxyz)]
#yrm.YAMLCreateMapping(grid2cell)
#yrm.YAMLRunFile(True, True, True, "myrunfile")
#yrm.YAMLSetFilePrefix("myfile")
#print(yrm.yaml_doc)
#yaml.dump(yrm.yaml_doc, sys.stdout, sort_keys=False)
#yrm.WriteYAMLDoc("dumpfile")