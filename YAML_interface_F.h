#ifdef USE_YAML
#ifndef INC_YAML_interface_F_H
#define INC_YAML_interface_F_H

#include "irm_dll_export.h"

#if defined(__cplusplus)
extern "C" {
#endif

IRM_DLL_EXPORT int CreateYAMLPhreeqcRM_F(void);
IRM_DLL_EXPORT int DestroyYAMLPhreeqcRM_F(int* id);
IRM_DLL_EXPORT IRM_RESULT YAMLClear_F(int* id);
IRM_DLL_EXPORT IRM_RESULT WriteYAMLDoc_F(int* id, const char* file_name);
IRM_DLL_EXPORT IRM_RESULT YAMLAddOutputVars_F(int* id, char* option_in, char* def_in);
IRM_DLL_EXPORT IRM_RESULT YAMLCloseFiles_F(int* id);
IRM_DLL_EXPORT IRM_RESULT YAMLCreateMapping_F(int* id, int* grid2chem, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLDumpModule_F(int* id, int* dump_on, int* append);
IRM_DLL_EXPORT IRM_RESULT YAMLFindComponents_F(int* id);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolutions2Module_F(int* id, int* solutions, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialEquilibriumPhases2Module_F(int* id, int* equilibrium_phases, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialExchanges2Module_F(int* id, int* exchanges, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSurfaces2Module_F(int* id, int* surfaces, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialGasPhases2Module_F(int* id, int* gas_phases, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialSolidSolutions2Module_F(int* id, int* solid_solutions, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialKinetics2Module_F(int* id, int* kinetics, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module_F(int* id, int* ic1, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqc2Module_mix_F(int* id, int* ic1, int* ic2, double* f1, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLInitialPhreeqcCell2Module_F(int* id,
							int* n, int* cell_numbers, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLLoadDatabase_F(int* id, const char* database);
IRM_DLL_EXPORT IRM_RESULT YAMLLogMessage_F(int* id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLOpenFiles_F(int* id);
IRM_DLL_EXPORT IRM_RESULT YAMLOutputMessage_F(int* id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLRunCells_F(int* id);
IRM_DLL_EXPORT IRM_RESULT YAMLRunFile_F(int* id, int* workers, int* initial_phreeqc,
							int* utility, const char* file_name);
IRM_DLL_EXPORT IRM_RESULT YAMLRunString_F(int* id, int* workers, int* initial_phreeqc,
							int* utility, const char* input_string);
IRM_DLL_EXPORT IRM_RESULT YAMLScreenMessage_F(int* id, const char* str);
IRM_DLL_EXPORT IRM_RESULT YAMLSetComponentH2O_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetConcentrations_F(int* id, double* c, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetCurrentSelectedOutputUserNumber_F(int* id, int* n_user);
IRM_DLL_EXPORT IRM_RESULT YAMLSetDensityUser_F(int* id, double* density, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetDumpFileName_F(int* id, const char* dump_name);
IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorHandlerMode_F(int* id, int* mode);
IRM_DLL_EXPORT IRM_RESULT YAMLSetErrorOn_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetFilePrefix_F(int* id, const char* prefix);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGasCompMoles_F(int* id, double* gas_moles, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGasPhaseVolume_F(int* id, double* gas_volume, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetGridCellCount_F(int* id, int* count);
IRM_DLL_EXPORT IRM_RESULT YAMLSetNthSelectedOutput_F(int* id, int* n);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPartitionUZSolids_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPorosity_F(int* id, double* por, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPressure_F(int* id, double* p, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryMask_F(int* id, int* cell_mask, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetPrintChemistryOn_F(int* id, int* workers, int* initial_phreeqc,
							int* utility);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceByCell_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRebalanceFraction_F(int* id, double* f);
IRM_DLL_EXPORT IRM_RESULT YAMLSetRepresentativeVolume_F(int* id, double* rv, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSaturationUser_F(int* id, double* sat, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetScreenOn_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSelectedOutputOn_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLSetSpeciesSaveOn_F(int* id, int* save_on);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTemperature_F(int* id, double* t, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTime_F(int* id, double* time);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeConversion_F(int* id, double* conv_factor);
IRM_DLL_EXPORT IRM_RESULT YAMLSetTimeStep_F(int* id, double* time_step);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsExchange_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsGasPhase_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsKinetics_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsPPassemblage_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSolution_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSSassemblage_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSetUnitsSurface_F(int* id, int* option);
IRM_DLL_EXPORT IRM_RESULT YAMLSpeciesConcentrations2Module_F(int* id, double* species_conc, int* dim);
IRM_DLL_EXPORT IRM_RESULT YAMLStateSave_F(int* id, int* istate);
IRM_DLL_EXPORT IRM_RESULT YAMLStateApply_F(int* id, int* istate);
IRM_DLL_EXPORT IRM_RESULT YAMLStateDelete_F(int* id, int* istate);
IRM_DLL_EXPORT IRM_RESULT YAMLThreadCount_F(int* id, int* nthreads);
IRM_DLL_EXPORT IRM_RESULT YAMLUseSolutionDensityVolume_F(int* id, int* tf);
IRM_DLL_EXPORT IRM_RESULT YAMLWarningMessage_F(int* id, const char* warnstr);

#if defined(__cplusplus)
}
#endif

#endif // INC_YAML_interface_F_H
#endif // USE_YAML
