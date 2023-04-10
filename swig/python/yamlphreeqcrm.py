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
    def YAMLCloseFiles(self, void):
        self.yaml_doc["CloseFiles"] = ""
    def YAMLCreateMapping(self, g2c):
        self.yaml_doc["CreateMapping"] = g2c
    def YAMLDumpModule(self, dump_on, append):
        node = dict()
        node["dump_on"] = dump_on
        node["append"] = append
        self.yaml_doc["DumpModule"] = node
    def YAMLFindComponents(self):
        self.yaml_doc["FindComponents"] = ""
    def YAMLInitialPhreeqc2Module(self, initial_conditions1):
        self.yaml_doc["InitialPhreeqc2Module"] = 1
    def YAMLInitialPhreeqc2Module_mix(self, ic1,  ic2,  f1):
        node = dict()
        node["ic1"] = ic1
        node["ic2"] = ic2
        node["f1"] = f1
        self.yaml_doc["InitialPhreeqc2Module_mix"] = node
    def YAMLInitialPhreeqcCell2Module(self, n,  cell_numbers):
        node = dict()
        node["n"] = n
        node["cell_numbers"] = cell_numbers
        self.yaml_doc["InitialPhreeqcCell2Module"] = node
    def YAMLLoadDatabase(self, database):
        self.yaml_doc["LoadDatabase"] = database
    def YAMLLogMessage(self, str):
        self.yaml_doc["LogMessage"] = str
    def YAMLOpenFiles(self):
        self.yaml_doc["OpenFiles"] = ""
    def YAMLOutputMessage(self, str):
        self.yaml_doc["OutputMessage"] = str
    def YAMLRunCells(self):
        self.yaml_doc["RunCells"] = ""
    def YAMLRunFile(self, workers, initial_phreeqc, utility, chemistry_name):
        node = dict()
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["chemistry_name"] = chemistry_name
        self.yaml_doc["RunFile"] = node
    def YAMLRunString(self, workers, initial_phreeqc, utility, input_string):
        node = dict()
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        node["input_string"] = input_string
        self.yaml_doc["RunString"] = node
    def YAMLScreenMessage(self, str):
        self.yaml_doc["ScreenMessage"] = str
    def YAMLSetComponentH2O(self, tf):
        self.yaml_doc["SetComponentH2O"] = tf
    def YAMLSetConcentrations(self, c):
        self.yaml_doc["SetConcentrations"] = c
    def YAMLSetCurrentSelectedOutputUserNumber(self, n_user):
        self.yaml_doc["SetCurrentSelectedOutputUserNumber"] = n_user
    def YAMLSetDensity(self, density):
        self.yaml_doc["SetDensity"] = density
    def YAMLSetDumpFileName(self, dump_name):
        self.yaml_doc["SetDumpFileName"] = dump_name
    def YAMLSetErrorHandlerMode(self, mode):
        self.yaml_doc["SetErrorHandlerMode"] = mode
    def YAMLSetErrorOn(self, tf):
        self.yaml_doc["SetErrorOn"] = tf
    def YAMLSetFilePrefix(self, prefix):
        self.yaml_doc["SetFilePrefix"] = prefix
    def YAMLSetGasCompMoles(self,  gas_moles):
        self.yaml_doc["SetGasCompMoles"] = gas_moles
    def YAMLSetGasPhaseVolume(self,  gas_volume):
        self.yaml_doc["SetGasPhaseVolume"] = gas_volume
    def YAMLSetGridCellCount(self, n):
        self.yaml_doc["SetGridCellCount"] = n
    def YAMLSetNthSelectedOutput(self, n):
        self.yaml_doc["SetNthSelectedOutput"] = n
    def YAMLSetPartitionUZSolids(self, tf):
        self.yaml_doc["SetPartitionUZSolids"] = tf
    def YAMLSetPorosity(self,  por):
        self.yaml_doc["SetPorosity"] = por
    def YAMLSetPressure(self,  p):
        self.yaml_doc["SetPressure"] = p
    def YAMLSetPrintChemistryMask(self,  cell_mask):
        self.yaml_doc["SetPrintChemistryMask"] = cell_mask
    def YAMLSetPrintChemistryOn(self, workers, initial_phreeqc, utility):
        node = dict()
        node["workers"] = workers
        node["initial_phreeqc"] = initial_phreeqc
        node["utility"] = utility
        self.yaml_doc["SetPrintChemistryOn"] = node
    def YAMLSetRebalanceByCell(self, tf):
        self.yaml_doc["SetRebalanceByCell"] = tf
    def YAMLSetRebalanceFraction(self, f):
        self.yaml_doc["SetRebalanceFraction"] = f
    def YAMLSetRepresentativeVolume(self,  rv):
        self.yaml_doc["SetRepresentativeVolume"] = rv
    def YAMLSetSaturation(self,  sat):
        self.yaml_doc["SetSaturation"] = sat
    def YAMLSetScreenOn(self, tf):
        self.yaml_doc["SetScreenOn"] = tf
    def YAMLSetSelectedOutputOn(self, tf):
        self.yaml_doc["SetSelectedOutputOn"] = tf
    def YAMLSetSpeciesSaveOn(self, save_on):
        self.yaml_doc["SetSpeciesSaveOn"] = save_on
    def YAMLSetTemperature(self,  t):
        self.yaml_doc["SetTemperature"] = t
    def YAMLSetTime(self, time):
        self.yaml_doc["SetTime"] = time
    def YAMLSetTimeConversion(self, conv_factor):
        self.yaml_doc["SetTimeConversion"] = conv_factor
    def YAMLSetTimeStep(self, time_step):
        self.yaml_doc["SetTimeStep"] = time_step
    def YAMLSetUnitsExchange(self, option):
        self.yaml_doc["SetUnitsExchange"] = option
    def YAMLSetUnitsGasPhase(self, option):
        self.yaml_doc["SetUnitsGasPhase"] = option
    def YAMLSetUnitsKinetics(self, option):
        self.yaml_doc["SetUnitsKinetics"] = option
    def YAMLSetUnitsPPassemblage(self, option):
        self.yaml_doc["SetUnitsPPassemblage"] = option
    def YAMLSetUnitsSolution(self, option):
        self.yaml_doc["SetUnitsSolution"] = option
    def YAMLSetUnitsSSassemblage(self, option):
        self.yaml_doc["SetUnitsSSassemblage"] = option
    def YAMLSetUnitsSurface(self, option):
        self.yaml_doc["SetUnitsSurface"] = option
    def YAMLSpeciesConcentrations2Module(self,  species_conc):
        self.yaml_doc["SpeciesConcentrations2Module"] = species_conc
    def YAMLStateSave(self, istate):
        self.yaml_doc["StateSave"] = istate
    def YAMLStateApply(self, istate):
        self.yaml_doc["StateApply"] = istate
    def YAMLStateDelete(self, istate):
        self.yaml_doc["StateDelete"] = istate
    def YAMLUseSolutionDensityVolume(self, tf):
        self.yaml_doc["UseSolutionDensityVolume"] = tf
    def YAMLWarningMessage(self, warnstr):
        self.yaml_doc["WarningMessage"] = warnstr


#nxyz = 40
#yrm = YAMLPhreeqcRM()
#grid2cell = [i for i in range(nxyz)]
#yrm.YAMLCreateMapping(grid2cell)
#yrm.YAMLRunFile(True, True, True, "myrunfile")
#yrm.YAMLSetFilePrefix("myfile")
#print(yrm.yaml_doc)
#yaml.dump(yrm.yaml_doc, sys.stdout, sort_keys=False)
#yrm.WriteYAMLDoc("dumpfile")