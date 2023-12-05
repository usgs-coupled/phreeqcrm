#ifdef USE_YAML
#ifndef INC_YAML_interface_C_H
#define INC_YAML_interface_C_H

#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif

IRM_DLL_EXPORT int CreateYAMLPhreeqcRM(void);
IRM_DLL_EXPORT int DestroyYAMLPhreeqcRM(int id);
IRM_DLL_EXPORT IRM_RESULT YAMLClear(int id);
IRM_DLL_EXPORT IRM_RESULT WriteYAMLDoc(int id, const char* file_name);
IRM_DLL_EXPORT IRM_RESULT YAMLAddOutputVars(int id, char* option_in, char* def_in);
IRM_DLL_EXPORT IRM_RESULT YAMLCloseFiles(int id);
IRM_DLL_EXPORT IRM_RESULT YAMLCreateMapping(int id, int* grid2chem, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLDumpModule(int id, int dump_on, int append);
IRM_DLL_EXPORT IRM_RESULT YAMLFindComponents(int id);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolutions2Module(int id, int* solutions, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialEquilibriumPhases2Module(int id, int* equilibrium_phases, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialExchanges2Module(int id, int* exchanges, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSurfaces2Module(int id, int* surfaces, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialGasPhases2Module(int id, int* gas_phases, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolidSolutions2Module(int id, int* solid_solutions, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialKinetics2Module(int id, int* kinetics, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module(int id, int* ic1, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module_mix(int id, int* ic1, int* ic2, double* f1, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqcCell2Module(int id,
							int n, int* cell_numbers, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLLoadDatabase(int id, const char* database);
IRM_DLL_EXPORT IRM_RESULT YAMLLogMessage(int id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLOpenFiles(int id);
IRM_DLL_EXPORT IRM_RESULT YAMLOutputMessage(int id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLRunCells(int id);
IRM_DLL_EXPORT IRM_RESULT YAMLRunFile(int id, int workers, int initial_phreeqc,
							int utility, const char* file_name);
IRM_DLL_EXPORT IRM_RESULT YAMLRunString(int id, int workers, int initial_phreeqc,
							int utility, const char* input_string);
IRM_DLL_EXPORT IRM_RESULT YAMLScreenMessage(int id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLSetComponentH2O(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetConcentrations(int id, double* c, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetCurrentSelectedOutputUserNumber(int id, int n_user);
IRM_DLL_EXPORT IRM_RESULT YAMLSetDensityUser(int id, double* density, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetDumpFileName(int id, const char* dump_name);
IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorHandlerMode(int id, int mode);
IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorOn(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetFilePrefix(int id, const char* prefix);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGasCompMoles(int id, double* gas_moles, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGasPhaseVolume(int id, double* gas_volume, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGridCellCount(int id, int count);
IRM_DLL_EXPORT IRM_RESULT YAMLSetNthSelectedOutput(int id, int n);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPartitionUZSolids(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPorosity(int id, double* por, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPressure(int id, double* p, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryMask(int id, int* cell_mask, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryOn(int id, int workers, int initial_phreeqc,
							int utility);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceByCell(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceFraction(int id, double f);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRepresentativeVolume(int id, double* rv, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSaturationUser(int id, double* sat, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetScreenOn(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSelectedOutputOn(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSpeciesSaveOn(int id, int save_on);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTemperature(int id, double* t, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTime(int id, double time);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeConversion(int id, double conv_factor);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeStep(int id, double time_step);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsExchange(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsGasPhase(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsKinetics(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsPPassemblage(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSolution(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSSassemblage(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSurface(int id, int option);
IRM_DLL_EXPORT IRM_RESULT YAMLSpeciesConcentrations2Module(int id, double* species_conc, int dim);
IRM_DLL_EXPORT IRM_RESULT YAMLStateSave(int id, int istate);
IRM_DLL_EXPORT IRM_RESULT YAMLStateApply(int id, int istate);
IRM_DLL_EXPORT IRM_RESULT YAMLStateDelete(int id, int istate);
IRM_DLL_EXPORT IRM_RESULT YAMLThreadCount(int id, int nthreads);
IRM_DLL_EXPORT IRM_RESULT YAMLUseSolutionDensityVolume(int id, int tf);
IRM_DLL_EXPORT IRM_RESULT YAMLWarningMessage(int id, const char* warnstr);

#if defined(__cplusplus)
}
#endif

#endif // INC_YAML_interface_C_H
#endif // USE_YAML
