#!/bin/sh
for f in html/*; do
  sed -i "s/getgridcellcountyaml/GetGridCellCountYAML/g" $f 
done
for f in html/*; do
  sed -i "s/setgridcellcountyaml/SetGridCellCountYAML/g" $f 
done
for f in html/*; do
  sed -i "s/rm_abort/RM_Abort/g" $f 
done
for f in html/*; do
  sed -i "s/rm_closefiles/RM_CloseFiles/g" $f 
done
for f in html/*; do
  sed -i "s/rm_concentrations2utility/RM_Concentrations2Utility/g" $f 
done
for f in html/*; do
  sed -i "s/rm_create/RM_Create/g" $f 
done
for f in html/*; do
  sed -i "s/rm_createmapping/RM_CreateMapping/g" $f 
done
for f in html/*; do
  sed -i "s/rm_decodeerror/RM_DecodeError/g" $f 
done
for f in html/*; do
  sed -i "s/rm_destroy/RM_Destroy/g" $f 
done
for f in html/*; do
  sed -i "s/rm_dumpmodule/RM_DumpModule/g" $f 
done
for f in html/*; do
  sed -i "s/rm_errormessage/RM_ErrorMessage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_findcomponents/RM_FindComponents/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getbackwardmapping/RM_GetBackwardMapping/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getchemistrycellcount/RM_GetChemistryCellCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getcomponentcount/RM_GetComponentCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getcomponent/RM_GetComponent/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getconcentrations/RM_GetConcentrations/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getcurrentselectedoutputusernumber/RM_GetCurrentSelectedOutputUserNumber/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getdensitycalculated/RM_GetDensityCalculated/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getendcell/RM_GetEndCell/g" $f 
done
for f in html/*; do
  sed -i "s/rm_geterrorstringlength/RM_GetErrorStringLength/g" $f 
done
for f in html/*; do
  sed -i "s/rm_geterrorstring/RM_GetErrorString/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getfileprefix/RM_GetFilePrefix/g" $f 
done

for f in html/*; do
  sed -i "s/rm_getgascompmoles/RM_GetGasCompMoles/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getgascomppressures/RM_GetGasCompPressures/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getgascompphi/RM_GetGasCompPhi/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getgasphasevolume/RM_GetGasPhaseVolume/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setgascompmoles/RM_SetGasCompMoles/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setgasphasevolume/RM_SetGasPhaseVolume/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspecieslog10molalities/RM_GetSpeciesLog10Molalities/g" $f 
done

for f in html/*; do
  sed -i "s/rm_getgfw/RM_GetGfw/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getgridcellcount/RM_GetGridCellCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getiphreeqcid/RM_GetIPhreeqcId/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getmpimyself/RM_GetMpiMyself/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getmpitasks/RM_GetMpiTasks/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getnthselectedoutputusernumber/RM_GetNthSelectedOutputUserNumber/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getporosity/RM_GetPorosity/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getpressure/RM_GetPressure/g" $f  
done
for f in html/*; do
  sed -i "s/rm_getsaturationcalculated/RM_GetSaturationCalculated/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getselectedoutputcolumncount/RM_GetSelectedOutputColumnCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getselectedoutputcount/RM_GetSelectedOutputCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getselectedoutputheading/RM_GetSelectedOutputHeading/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getselectedoutputrowcount/RM_GetSelectedOutputRowCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getselectedoutput/RM_GetSelectedOutput/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getsolutionvolume/RM_GetSolutionVolume/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciesconcentrations/RM_GetSpeciesConcentrations/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciescount/RM_GetSpeciesCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciesd25/RM_GetSpeciesD25/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciesname/RM_GetSpeciesName/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciessaveon/RM_GetSpeciesSaveOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspeciesz/RM_GetSpeciesZ/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getstartcell/RM_GetStartCell/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getthreadcount/RM_GetThreadCount/g" $f 
done
for f in html/*; do
  sed -i "s/rm_gettime/RM_GetTime/g" $f 
done
for f in html/*; do
  sed -i "s/rm_gettimeconversion/RM_GetTimeConversion/g" $f 
done
for f in html/*; do
  sed -i "s/rm_gettimestep/RM_GetTimeStep/g" $f 
done
for f in html/*; do
  sed -i "s/rm_initializeyaml/RM_InitializeYAML/g" $f 
done
for f in html/*; do
  sed -i "s/rm_initialphreeqc2concentrations/RM_InitialPhreeqc2Concentrations/g" $f 
done
for f in html/*; do
  sed -i "s/rm_initialphreeqc2module/RM_InitialPhreeqc2Module/g" $f 
done
for f in html/*; do
  sed -i "s/rm_initialphreeqccell2module/RM_InitialPhreeqcCell2Module/g" $f 
done
for f in html/*; do
  sed -i "s/rm_initialphreeqc2speciesconcentrations/RM_InitialPhreeqc2SpeciesConcentrations/g" $f 
done
for f in html/*; do
  sed -i "s/rm_loaddatabase/RM_LoadDatabase/g" $f 
done
for f in html/*; do
  sed -i "s/rm_logmessage/RM_LogMessage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_mpiworkerbreak/RM_MpiWorkerBreak/g" $f 
done
for f in html/*; do
  sed -i "s/rm_mpiworker/RM_MpiWorker/g" $f 
done
for f in html/*; do
  sed -i "s/rm_openfiles/RM_OpenFiles/g" $f 
done
for f in html/*; do
  sed -i "s/rm_outputmessage/RM_OutputMessage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_runcells/RM_RunCells/g" $f 
done
for f in html/*; do
  sed -i "s/rm_runfile/RM_RunFile/g" $f 
done
for f in html/*; do
  sed -i "s/rm_runstring/RM_RunString/g" $f 
done
for f in html/*; do
  sed -i "s/rm_screenmessage/RM_ScreenMessage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setcomponenth2o/RM_SetComponentH2O/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setconcentrations/RM_SetConcentrations/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setcurrentselectedoutputusernumber/RM_SetCurrentSelectedOutputUserNumber/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setdensityuser/RM_SetDensityUser/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setdumpfilename/RM_SetDumpFileName/g" $f 
done
for f in html/*; do
  sed -i "s/rm_seterrorhandlermode/RM_SetErrorHandlerMode/g" $f 
done
for f in html/*; do
  sed -i "s/rm_seterroron/RM_SetErrorOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setfileprefix/RM_SetFilePrefix/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setmpiworkercallback/RM_SetMpiWorkerCallback/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setnthselectedoutput/RM_SetNthSelectedOutput/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setpartitionuzsolids/RM_SetPartitionUZSolids/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setporosity/RM_SetPorosity/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setprintchemistrymask/RM_SetPrintChemistryMask/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setprintchemistryon/RM_SetPrintChemistryOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setpressure/RM_SetPressure/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setrebalancefraction/RM_SetRebalanceFraction/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setrebalancebycell/RM_SetRebalanceByCell/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setrepresentativevolume/RM_SetRepresentativeVolume/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setsaturationuser/RM_SetSaturationUser/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setscreenon/RM_SetScreenOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setselectedoutputon/RM_SetSelectedOutputOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setspeciessaveon/RM_SetSpeciesSaveOn/g" $f 
done
for f in html/*; do
  sed -i "s/rm_settemperature/RM_SetTemperature/g" $f 
done
for f in html/*; do
  sed -i "s/rm_settimeconversion/RM_SetTimeConversion/g" $f 
done
for f in html/*; do
  sed -i "s/rm_settimestep/RM_SetTimeStep/g" $f 
done
for f in html/*; do
  sed -i "s/rm_settime/RM_SetTime/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitsexchange/RM_SetUnitsExchange/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitsgasphase/RM_SetUnitsGasPhase/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitskinetics/RM_SetUnitsKinetics/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitsppassemblage/RM_SetUnitsPPassemblage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitssolution/RM_SetUnitsSolution/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_setunitssurface/RM_SetUnitsSurface/g" $f 
done
for f in html/*; do
  sed -i "s/rm_speciesconcentrations2module/RM_SpeciesConcentrations2Module/g" $f 
done
for f in html/*; do
  sed -i "s/rm_usesolutiondensityvolume/RM_UseSolutionDensityVolume/g" $f 
done
for f in html/*; do
  sed -i "s/rm_warningmessage/RM_WarningMessage/g" $f 
done
for f in html/*; do
  sed -i "s/rm_getspecieslog10gammas/RM_GetSpeciesLog10Gammas/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getexchangespeciescount/RM_GetExchangeSpeciesCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getexchangespeciesname/RM_GetExchangeSpeciesName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getexchangename/RM_GetExchangeName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsurfacespeciescount/RM_GetSurfaceSpeciesCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsurfacespeciesname/RM_GetSurfaceSpeciesName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsurfacetype/RM_GetSurfaceType/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsurfacename/RM_GetSurfaceName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getequilibriumphasescount/RM_GetEquilibriumPhasesCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getequilibriumphasesname/RM_GetEquilibriumPhasesName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getgascomponentscount/RM_GetGasComponentsCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getgascomponentsname/RM_GetGasComponentsName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getkineticreactionscount/RM_GetKineticReactionsCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getkineticreactionsname/RM_GetKineticReactionsName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsolidsolutioncomponentscount/RM_GetSolidSolutionComponentsCount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsolidsolutioncomponentsname/RM_GetSolidSolutionComponentsName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsolidsolutionname/RM_GetSolidSolutionName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsicount/RM_GetSICount/g"  $f
done
for f in html/*; do
  sed -i "s/rm_getsiname/RM_GetSIName/g"  $f
done
for f in html/*; do
  sed -i "s/rm_statesave/RM_StateSave/g"  $f
done
for f in html/*; do
  sed -i "s/rm_stateapply/RM_StateApply/g"  $f
done
for f in html/*; do
  sed -i "s/rm_statedelete/RM_StateDelete/g"  $f
done
#for f in html/*; do
#  sed -i "s/rm_bmi_finalize/RM_BMI_Finalize/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getcomponentname/RM_BMI_GetComponentName/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getcurrenttime/RM_BMI_GetCurrentTime/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getendtime/RM_BMI_GetEndTime/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getinputitemcount/RM_BMI_GetInputItemCount/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getinputvarnames/RM_BMI_GetInputVarNames/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getoutputitemcount/RM_BMI_GetOutputItemCount/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getoutputvarnames/RM_BMI_GetOutputVarNames/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_gettimestep/RM_BMI_GetTimeStep/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_gettimeunits/RM_BMI_GetTimeUnits/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getvalue/RM_BMI_GetValue/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getvaritemsize/RM_BMI_GetVarItemsize/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getvarnbytes/RM_BMI_GetVarNbytes/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getvartype/RM_BMI_GetVarType/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_getvarunits/RM_BMI_GetVarUnits/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_initialize/RM_BMI_Initialize/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_setvalue/RM_BMI_SetValue/g" $f
#done
#for f in html/*; do
#  sed -i "s/rm_bmi_update/RM_BMI_Update/g" $f
#done
for f in html/*; do
  sed -i "s/rm_getcomponents/RM_GetComponents/g" $f
done
for f in html/*; do
  sed -i "s/rm_getgridcellcountyaml/RM_GetGridCellCountYAML/g" $f
done
for f in html/*; do
  sed -i "s/rm_getselectedoutputheadings/RM_GetSelectedOutputHeadings/g" $f
done
for f in html/*; do
  sed -i "s/rm_gettemperature/RM_GetTemperature/g" $f
done
for f in html/*; do
  sed -i "s/rm_getviscosity/RM_GetViscosity/g" $f
done
for f in html/*; do
	sed -i "s/rm_getithconcentration/RM_GetIthConcentration/g" $f
done
for f in html/*; do
	sed -i "s/rm_getithspeciesconcentration/RM_GetIthSpeciesConcentration/g" $f
done
for f in html/*; do
	sed -i "s/rm_setithconcentration/RM_SetIthConcentration/g" $f
done
for f in html/*; do
	sed -i "s/rm_setithspeciesconcentration/RM_SetIthSpeciesConcentration/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialsolutions2module/RM_InitialSolutions2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialequilibriumphases2module/RM_InitialEquilibriumPhases2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialexchanges2module/RM_InitialExchanges2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialsurfaces2module/RM_InitialSurfaces2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialgasphases2module/RM_InitialGasPhases2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialsolidsolutions2module/RM_InitialSolidSolutions2Module/g" $f
done
for f in html/*; do
	sed -i "s/rm_initialkinetics2module/RM_InitialKinetics2Module/g" $f
done


for f in html/*; do
	sed -i "s/destroyyamlphreeqcrm/DestroyYAMLPhreeqcRM/g"  $f
done
for f in html/*; do
	sed -i "s/createyamlphreeqcrm/CreateYAMLPhreeqcRM/g"  $f
done
for f in html/*; do
	sed -i "s/writeyamldoc/WriteYAMLDoc/g"  $f
done
for f in html/*; do
	sed -i "s/yamladdoutputvars/YAMLAddOutputVars/g" $f
done
for f in html/*; do
	sed -i "s/yamlclear/YAMLClear/g" $f
done
for f in html/*; do
	sed -i "s/yamlclosefiles/YAMLCloseFiles/g" $f
done
for f in html/*; do
	sed -i "s/yamlcreatemapping/YAMLCreateMapping/g" $f
done
for f in html/*; do
	sed -i "s/yamldumpmodule/YAMLDumpModule/g" $f
done
for f in html/*; do
	sed -i "s/yamlfindcomponents/YAMLFindComponents/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialsolutions2module/YAMLInitialSolutions2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialequilibriumphases2module/YAMLInitialEquilibriumPhases2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialexchanges2module/YAMLInitialExchanges2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialsurfaces2module/YAMLInitialSurfaces2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialgasphases2module/YAMLInitialGasPhases2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialsolidsolutions2module/YAMLInitialSolidSolutions2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialkinetics2module/YAMLInitialKinetics2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialphreeqc2module/YAMLInitialPhreeqc2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialphreeqc2module_mix/YAMLInitialPhreeqc2Module_mix/g" $f
done
for f in html/*; do
	sed -i "s/yamlinitialphreeqccell2module/YAMLInitialPhreeqcCell2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlloaddatabase/YAMLLoadDatabase/g" $f
done
for f in html/*; do
	sed -i "s/yamllogmessage/YAMLLogMessage/g" $f
done
for f in html/*; do
	sed -i "s/yamlopenfiles/YAMLOpenFiles/g" $f
done
for f in html/*; do
	sed -i "s/yamloutputmessage/YAMLOutputMessage/g" $f
done
for f in html/*; do
	sed -i "s/yamlruncells/YAMLRunCells/g" $f
done
for f in html/*; do
	sed -i "s/yamlrunfile/YAMLRunFile/g" $f
done
for f in html/*; do
	sed -i "s/yamlrunstring/YAMLRunString/g" $f
done
for f in html/*; do
	sed -i "s/yamlscreenmessage/YAMLScreenMessage/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetcomponenth2o/YAMLSetComponentH2O/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetconcentrations/YAMLSetConcentrations/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetcurrentselectedoutputusernumber/YAMLSetCurrentSelectedOutputUserNumber/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetdensityuser/YAMLSetDensityUser/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetdumpfilename/YAMLSetDumpFileName/g" $f
done
for f in html/*; do
	sed -i "s/yamlseterrorhandlermode/YAMLSetErrorHandlerMode/g" $f
done
for f in html/*; do
	sed -i "s/yamlseterroron/YAMLSetErrorOn/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetfileprefix/YAMLSetFilePrefix/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetgascompmoles/YAMLSetGasCompMoles/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetgasphasevolume/YAMLSetGasPhaseVolume/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetgridcellcount/YAMLSetGridCellCount/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetnthselectedoutput/YAMLSetNthSelectedOutput/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetpartitionuzsolids/YAMLSetPartitionUZSolids/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetporosity/YAMLSetPorosity/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetpressure/YAMLSetPressure/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetprintchemistrymask/YAMLSetPrintChemistryMask/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetprintchemistryon/YAMLSetPrintChemistryOn/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetrebalancebycell/YAMLSetRebalanceByCell/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetrebalancefraction/YAMLSetRebalanceFraction/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetrepresentativevolume/YAMLSetRepresentativeVolume/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetsaturationuser/YAMLSetSaturationUser/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetscreenon/YAMLSetScreenOn/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetselectedoutputon/YAMLSetSelectedOutputOn/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetspeciessaveon/YAMLSetSpeciesSaveOn/g" $f
done
for f in html/*; do
	sed -i "s/yamlsettemperature/YAMLSetTemperature/g" $f
done
for f in html/*; do
	sed -i "s/yamlsettimeconversion/YAMLSetTimeConversion/g" $f
done
for f in html/*; do
	sed -i "s/yamlsettimestep/YAMLSetTimeStep/g" $f
done
for f in html/*; do
	sed -i "s/yamlsettime/YAMLSetTime/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitsexchange/YAMLSetUnitsExchange/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitsgasphase/YAMLSetUnitsGasPhase/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitskinetics/YAMLSetUnitsKinetics/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitsppassemblage/YAMLSetUnitsPPassemblage/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitssolution/YAMLSetUnitsSolution/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitsssassemblage/YAMLSetUnitsSSassemblage/g" $f
done
for f in html/*; do
	sed -i "s/yamlsetunitssurface/YAMLSetUnitsSurface/g" $f
done
for f in html/*; do
	sed -i "s/yamlspeciesconcentrations2module/YAMLSpeciesConcentrations2Module/g" $f
done
for f in html/*; do
	sed -i "s/yamlstatesave/YAMLStateSave/g" $f
done
for f in html/*; do
	sed -i "s/yamlstateapply/YAMLStateApply/g" $f
done
for f in html/*; do
	sed -i "s/yamlstatedelete/YAMLStateDelete/g" $f
done
for f in html/*; do
	sed -i "s/yamlthreadcount/YAMLThreadCount/g" $f
done
for f in html/*; do
	sed -i "s/yamlusesolutiondensityvolume/YAMLUseSolutionDensityVolume/g" $f
done
for f in html/*; do
	sed -i "s/yamlwarningmessage/YAMLWarningMessage/g" $f
done

# Rename to mixed case in all files
for f in html/*; do
sed -i "s/yamlphreeqcrm/YAMLPhreeqcRM/g" $f
done
for f in html/*; do
sed -i "s/bmiphreeqcrm/BMIPhreeqcRM/g" $f
done
for f in html/*; do
sed -i "s/phreeqcrm/PhreeqcRM/g" $f
done
for f in html/*; do
sed -i "s/yaml_interface/YAML_interface/g" $f
done
#for f in html/*; do
#  sed -i "s/\bphreeqcrm\b/PhreeqcRM/g" $f
#done

# Change Titles
sed -i "s/PhreeqcRM Module Reference/PhreeqcRM Fortran Module Reference/g" html/namespacephreeqcrm.html
sed -i "s/RM_interface_C.h File Reference/PhreeqcRM C Reference/g" html/RM__interface__C_8h.html
sed -i "s/PhreeqcRM Class Reference/PhreeqcRM C++ Class Reference/g" html/classPhreeqcRM.html

echo finished
