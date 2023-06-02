@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/getgridcellcountyaml/GetGridCellCountYAML/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/setgridcellcountyaml/SetGridCellCountYAML/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_abort/RM_Abort/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_closefiles/RM_CloseFiles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_concentrations2utility/RM_Concentrations2Utility/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_create/RM_Create/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_createmapping/RM_CreateMapping/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_decodeerror/RM_DecodeError/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_destroy/RM_Destroy/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_dumpmodule/RM_DumpModule/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_errormessage/RM_ErrorMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_findcomponents/RM_FindComponents/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getbackwardmapping/RM_GetBackwardMapping/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getchemistrycellcount/RM_GetChemistryCellCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getcomponentcount/RM_GetComponentCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getcomponent/RM_GetComponent/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getconcentrations/RM_GetConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getcurrentselectedoutputusernumber/RM_GetCurrentSelectedOutputUserNumber/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getdensitycalculated/RM_GetDensityCalculated/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getendcell/RM_GetEndCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_geterrorstringlength/RM_GetErrorStringLength/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_geterrorstring/RM_GetErrorString/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getfileprefix/RM_GetFilePrefix/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgascompmoles/RM_GetGasCompMoles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgascompphi/RM_GetGasCompPhi/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgascomppressures/RM_GetGasCompPressures/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgasphasevolume/RM_GetGasPhaseVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgfw/RM_GetGfw/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgridcellcount/RM_GetGridCellCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getiphreeqcid/RM_GetIPhreeqcId/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getmpimyself/RM_GetMpiMyself/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getmpitasks/RM_GetMpiTasks/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getnthselectedoutputusernumber/RM_GetNthSelectedOutputUserNumber/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getporosity/RM_GetPorosity/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getpressure/RM_GetPressure/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsaturationcalculated/RM_GetSaturationCalculated/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getselectedoutputcolumncount/RM_GetSelectedOutputColumnCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getselectedoutputcount/RM_GetSelectedOutputCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getselectedoutputheading/RM_GetSelectedOutputHeading/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getselectedoutputrowcount/RM_GetSelectedOutputRowCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getselectedoutput/RM_GetSelectedOutput/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsolutionvolume/RM_GetSolutionVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciesconcentrations/RM_GetSpeciesConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciescount/RM_GetSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciesd25/RM_GetSpeciesD25/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspecieslog10molalities/RM_GetSpeciesLog10Molalities/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciesname/RM_GetSpeciesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciessaveon/RM_GetSpeciesSaveOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspeciesz/RM_GetSpeciesZ/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getstartcell/RM_GetStartCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getthreadcount/RM_GetThreadCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_gettime/RM_GetTime/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_gettimeconversion/RM_GetTimeConversion/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_gettimestep/RM_GetTimeStep/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_initializeyaml/RM_InitializeYAML/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_initialphreeqc2concentrations/RM_InitialPhreeqc2Concentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_initialphreeqc2module/RM_InitialPhreeqc2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_initialphreeqccell2module/RM_InitialPhreeqcCell2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_initialphreeqc2speciesconcentrations/RM_InitialPhreeqc2SpeciesConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_loaddatabase/RM_LoadDatabase/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_logmessage/RM_LogMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_mpiworkerbreak/RM_MpiWorkerBreak/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_mpiworker/RM_MpiWorker/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_openfiles/RM_OpenFiles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_outputmessage/RM_OutputMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_runcells/RM_RunCells/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_runfile/RM_RunFile/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_runstring/RM_RunString/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_screenmessage/RM_ScreenMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setcomponenth2o/RM_SetComponentH2O/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setconcentrations/RM_SetConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setcurrentselectedoutputusernumber/RM_SetCurrentSelectedOutputUserNumber/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setdensityuser/RM_SetDensityUser/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setdumpfilename/RM_SetDumpFileName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_seterrorhandlermode/RM_SetErrorHandlerMode/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_seterroron/RM_SetErrorOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setfileprefix/RM_SetFilePrefix/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setgascompmoles/RM_SetGasCompMoles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setgasphasevolume/RM_SetGasPhaseVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setmpiworkercallback/RM_SetMpiWorkerCallback/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setnthselectedoutput/RM_SetNthSelectedOutput/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setpartitionuzsolids/RM_SetPartitionUZSolids/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setporosity/RM_SetPorosity/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setprintchemistrymask/RM_SetPrintChemistryMask/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setprintchemistryon/RM_SetPrintChemistryOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setpressure/RM_SetPressure/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setrebalancefraction/RM_SetRebalanceFraction/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setrebalancebycell/RM_SetRebalanceByCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setrepresentativevolume/RM_SetRepresentativeVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setsaturationuser/RM_SetSaturationUser/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setscreenon/RM_SetScreenOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setselectedoutputon/RM_SetSelectedOutputOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setspeciessaveon/RM_SetSpeciesSaveOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_settemperature/RM_SetTemperature/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_settimeconversion/RM_SetTimeConversion/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_settimestep/RM_SetTimeStep/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_settime/RM_SetTime/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitsexchange/RM_SetUnitsExchange/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitsgasphase/RM_SetUnitsGasPhase/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitskinetics/RM_SetUnitsKinetics/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitsppassemblage/RM_SetUnitsPPassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitssolution/RM_SetUnitsSolution/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_setunitssurface/RM_SetUnitsSurface/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_speciesconcentrations2module/RM_SpeciesConcentrations2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_usesolutiondensityvolume/RM_UseSolutionDensityVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_warningmessage/RM_WarningMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getspecieslog10gammas/RM_GetSpeciesLog10Gammas/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getexchangespeciescount/RM_GetExchangeSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getexchangespeciesname/RM_GetExchangeSpeciesName/g"   %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getexchangename/RM_GetExchangeName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsurfacespeciescount/RM_GetSurfaceSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsurfacespeciesname/RM_GetSurfaceSpeciesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsurfacetype/RM_GetSurfaceType/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsurfacename/RM_GetSurfaceName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getequilibriumphasescount/RM_GetEquilibriumPhasesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getequilibriumphasesname/RM_GetEquilibriumPhasesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgascomponentscount/RM_GetGasComponentsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getgascomponentsname/RM_GetGasComponentsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getkineticreactionscount/RM_GetKineticReactionsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getkineticreactionsname/RM_GetKineticReactionsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsolidsolutioncomponentscount/RM_GetSolidSolutionComponentsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsolidsolutioncomponentsname/RM_GetSolidSolutionComponentsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsolidsolutionname/RM_GetSolidSolutionName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsicount/RM_GetSICount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_getsiname/RM_GetSIName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_statesave/RM_StateSave/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_stateapply/RM_StateApply/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/rm_statedelete/RM_StateDelete/g" %%f 
  @del sed*
)

#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_finalize\b/RM_BMI_Finalize/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getcomponentname\b/RM_BMI_GetComponentName/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getcurrenttime\b/RM_BMI_GetCurrentTime/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getendtime\b/RM_BMI_GetEndTime/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getinputitemcount\b/RM_BMI_GetInputItemCount/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getinputvarnames\b/RM_BMI_GetInputVarNames/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getoutputitemcount\b/RM_BMI_GetOutputItemCount/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getoutputvarnames\b/RM_BMI_GetOutputVarNames/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_gettimestep\b/RM_BMI_GetTimeStep/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_gettimeunits\b/RM_BMI_GetTimeUnits/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getvalue\b/RM_BMI_GetValue/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getvaritemsize\b/RM_BMI_GetVarItemsize/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getvarnbytes\b/RM_BMI_GetVarNbytes/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getvartype\b/RM_BMI_GetVarType/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_getvarunits\b/RM_BMI_GetVarUnits/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_initialize\b/RM_BMI_Initialize/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_setvalue\b/RM_BMI_SetValue/g" %%f
#  @del sed*
#)
#@FOR %%f IN (.\html\*) DO @(
#@sed -i "s/\brm_bmi_update\b/RM_BMI_Update/g" %%f
#  @del sed*
#)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getcomponents\b/RM_GetComponents/g" %%f
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getgridcellcountyaml\b/RM_GetGridCellCountYAML/g" %%f
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getselectedoutputheadings\b/RM_GetSelectedOutputHeadings/g" %%f
  @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_gettemperature\b/RM_GetTemperature/g" %%f
  @del sed*
)

@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getviscosity\b/RM_GetViscosity/g" %%f
  @del sed*
)

@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getithconcentration\b/RM_GetIthConcentration/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_getithspeciesconcentration\b/RM_GetIthSpeciesConcentration/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_setithconcentration\b/RM_SetIthConcentration/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_setithspeciesconcentration\b/RM_SetIthSpeciesConcentration/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialsolutions2module\b/RM_InitialSolutions2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialequilibriumphases2module\b/RM_InitialEquilibriumPhases2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialexchanges2module\b/RM_InitialExchanges2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialsurfaces2module\b/RM_InitialSurfaces2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialgasphases2module\b/RM_InitialGasPhases2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialsolidsolutions2module\b/RM_InitialSolidSolutions2Module/g" %%f
	@del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\brm_initialkinetics2module\b/RM_InitialKinetics2Module/g" %%f
	@del sed*
)



@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\bdestroyyamlphreeqcrm\b/DestroyYAMLPhreeqcRM/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\bcreateyamlphreeqcrm\b/CreateYAMLPhreeqcRM/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\bwriteyamldoc\b/WriteYAMLDoc/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamladdoutputvars\b/YAMLAddOutputVars/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlclear\b/YAMLClear/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlclosefiles\b/YAMLCloseFiles/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlcreatemapping\b/YAMLCreateMapping/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamldumpmodule\b/YAMLDumpModule/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlfindcomponents\b/YAMLFindComponents/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialsolutions2module\b/YAMLInitialSolutions2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialequilibriumphases2module\b/YAMLInitialEquilibriumPhases2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialexchanges2module\b/YAMLInitialExchanges2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialsurfaces2module\b/YAMLInitialSurfaces2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialgasphases2module\b/YAMLInitialGasPhases2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialsolidsolutions2module\b/YAMLInitialSolidSolutions2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialkinetics2module\b/YAMLInitialKinetics2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialphreeqc2module\b/YAMLInitialPhreeqc2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialphreeqc2module_mix\b/YAMLInitialPhreeqc2Module_mix/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlinitialphreeqccell2module\b/YAMLInitialPhreeqcCell2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlloaddatabase\b/YAMLLoadDatabase/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamllogmessage\b/YAMLLogMessage/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlopenfiles\b/YAMLOpenFiles/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamloutputmessage\b/YAMLOutputMessage/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlruncells\b/YAMLRunCells/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlrunfile\b/YAMLRunFile/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlrunstring\b/YAMLRunString/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlscreenmessage\b/YAMLScreenMessage/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetcomponenth2o\b/YAMLSetComponentH2O/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetconcentrations\b/YAMLSetConcentrations/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetcurrentselectedoutputusernumber\b/YAMLSetCurrentSelectedOutputUserNumber/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetdensityuser\b/YAMLSetDensityUser/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetdumpfilename\b/YAMLSetDumpFileName/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlseterrorhandlermode\b/YAMLSetErrorHandlerMode/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlseterroron\b/YAMLSetErrorOn/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetfileprefix\b/YAMLSetFilePrefix/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetgascompmoles\b/YAMLSetGasCompMoles/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetgasphasevolume\b/YAMLSetGasPhaseVolume/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetgridcellcount\b/YAMLSetGridCellCount/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetnthselectedoutput\b/YAMLSetNthSelectedOutput/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetpartitionuzsolids\b/YAMLSetPartitionUZSolids/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetporosity\b/YAMLSetPorosity/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetpressure\b/YAMLSetPressure/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetprintchemistrymask\b/YAMLSetPrintChemistryMask/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetprintchemistryon\b/YAMLSetPrintChemistryOn/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetrebalancebycell\b/YAMLSetRebalanceByCell/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetrebalancefraction\b/YAMLSetRebalanceFraction/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetrepresentativevolume\b/YAMLSetRepresentativeVolume/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetsaturationuser\b/YAMLSetSaturationuser/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetscreenon\b/YAMLSetScreenOn/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetselectedoutputon\b/YAMLSetSelectedOutputOn/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetspeciessaveon\b/YAMLSetSpeciesSaveOn/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsettemperature\b/YAMLSetTemperature/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsettimeconversion\b/YAMLSetTimeConversion/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsettimestep\b/YAMLSetTimeStep/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsettime\b/YAMLSetTime/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitsexchange\b/YAMLSetUnitsExchange/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitsgasphase\b/YAMLSetUnitsGasPhase/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitskinetics\b/YAMLSetUnitsKinetics/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitsppassemblage\b/YAMLSetUnitsPPassemblage/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitssolution\b/YAMLSetUnitsSolution/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitsssassemblage\b/YAMLSetUnitsSSassemblage/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlsetunitssurface\b/YAMLSetUnitsSurface/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlspeciesconcentrations2module\b/YAMLSpeciesConcentrations2Module/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlstatesave\b/YAMLStateSave/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlstateapply\b/YAMLStateApply/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlstatedelete\b/YAMLStateDelete/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlthreadcount\b/YAMLThreadCount/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlusesolutiondensityvolume\b/YAMLUseSolutionDensityVolume/g" %%f
     @del sed*
)
@FOR %%f IN (.\html\*) DO @(
@sed -i "s/\byamlwarningmessage\b/YAMLWarningMessage/g" %%f
     @del sed*
)



@REM Rename PhreeqcRM to mixed case 
@FOR %%f IN (.\html\*) DO @(
  @sed -i "s/\bphreeqcrm\b/PhreeqcRM/g" %%f
  @del sed*
)

@REM Change Titles
@sed -i "s/PhreeqcRM Module Reference/PhreeqcRM Fortran Module Reference/g" ./html/classphreeqcrm.html
@sed -i "s/PhreeqcRM Module Reference/PhreeqcRM Fortran Module Reference/g" ./namespacephreeqcrm.html
@sed -i "s/RM_interface_C.h File Reference/PhreeqcRM C Reference/g" ./html/_r_m__interface___c_8h.html
@sed -i "s/PhreeqcRM Class Reference/PhreeqcRM C++ Class Reference/g" ./html/class_phreeqc_r_m.html

echo finished
