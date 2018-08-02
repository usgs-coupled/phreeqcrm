@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_abort/RM_Abort/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_closefiles/RM_CloseFiles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_concentrations2utility/RM_Concentrations2Utility/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_create/RM_Create/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_createmapping/RM_CreateMapping/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_decodeerror/RM_DecodeError/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_destroy/RM_Destroy/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_dumpmodule/RM_DumpModule/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_errormessage/RM_ErrorMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_findcomponents/RM_FindComponents/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getbackwardmapping/RM_GetBackwardMapping/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getchemistrycellcount/RM_GetChemistryCellCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getcomponent/RM_GetComponent/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getcomponentcount/RM_GetComponentCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getconcentrations/RM_GetConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getdensity/RM_GetDensity/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getendcell/RM_GetEndCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_geterrorstring/RM_GetErrorString/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_geterrorstringlength/RM_GetErrorStringLength/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getfileprefix/RM_GetFilePrefix/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getgfw/RM_GetGfw/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getgridcellcount/RM_GetGridCellCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getiphreeqcid/RM_GetIPhreeqcId/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getmpimyself/RM_GetMpiMyself/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getmpitasks/RM_GetMpiTasks/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getnthselectedoutputusernumber/RM_GetNthSelectedOutputUserNumber/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsaturation/RM_GetSaturation/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getselectedoutput/RM_GetSelectedOutput/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getselectedoutputcolumncount/RM_GetSelectedOutputColumnCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getselectedoutputcount/RM_GetSelectedOutputCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getselectedoutputheading/RM_GetSelectedOutputHeading/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getselectedoutputrowcount/RM_GetSelectedOutputRowCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsolutionvolume/RM_GetSolutionVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciesconcentrations/RM_GetSpeciesConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciescount/RM_GetSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciesd25/RM_GetSpeciesD25/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciesname/RM_GetSpeciesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciessaveon/RM_GetSpeciesSaveOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspeciesz/RM_GetSpeciesZ/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getstartcell/RM_GetStartCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getthreadcount/RM_GetThreadCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_gettime/RM_GetTime/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_gettimeconversion/RM_GetTimeConversion/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_gettimestep/RM_GetTimeStep/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_initialphreeqc2concentrations/RM_InitialPhreeqc2Concentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_initialphreeqc2module/RM_InitialPhreeqc2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_initialphreeqccell2module/RM_InitialPhreeqcCell2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_initialphreeqc2speciesconcentrations/RM_InitialPhreeqc2SpeciesConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_loaddatabase/RM_LoadDatabase/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_logmessage/RM_LogMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_mpiworker/RM_MpiWorker/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_mpiworkerbreak/RM_MpiWorkerBreak/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_openfiles/RM_OpenFiles/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_outputmessage/RM_OutputMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_runcells/RM_RunCells/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_runfile/RM_RunFile/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_runstring/RM_RunString/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_screenmessage/RM_ScreenMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setcomponenth2o/RM_SetComponentH2O/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setconcentrations/RM_SetConcentrations/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setcurrentselectedoutputusernumber/RM_SetCurrentSelectedOutputUserNumber/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setdensity/RM_SetDensity/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setdumpfilename/RM_SetDumpFileName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_seterrorhandlermode/RM_SetErrorHandlerMode/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setfileprefix/RM_SetFilePrefix/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setmpiworkercallback/RM_SetMpiWorkerCallback/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setpartitionuzsolids/RM_SetPartitionUZSolids/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setporosity/RM_SetPorosity/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setprintchemistrymask/RM_SetPrintChemistryMask/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setprintchemistryon/RM_SetPrintChemistryOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setpressure/RM_SetPressure/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setrebalancefraction/RM_SetRebalanceFraction/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setrebalancebycell/RM_SetRebalanceByCell/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setrepresentativevolume/RM_SetRepresentativeVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setsaturation/RM_SetSaturation/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setscreenon/RM_SetScreenOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setselectedoutputon/RM_SetSelectedOutputOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setspeciessaveon/RM_SetSpeciesSaveOn/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_settemperature/RM_SetTemperature/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_settime/RM_SetTime/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_settimeconversion/RM_SetTimeConversion/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_settimestep/RM_SetTimeStep/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitsexchange/RM_SetUnitsExchange/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitsgasphase/RM_SetUnitsGasPhase/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitskinetics/RM_SetUnitsKinetics/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitsppassemblage/RM_SetUnitsPPassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitssolution/RM_SetUnitsSolution/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitsssassemblage/RM_SetUnitsSSassemblage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_setunitssurface/RM_SetUnitsSurface/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_speciesconcentrations2module/RM_SpeciesConcentrations2Module/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_usesolutiondensityvolume/RM_UseSolutionDensityVolume/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_warningmessage/RM_WarningMessage/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getspecieslog10gammas/RM_GetSpeciesLog10Gammas/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getexchangespeciescount/RM_GetExchangeSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getexchangespeciesname/RM_GetExchangeSpeciesName/g"   %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getexchangename/RM_GetExchangeName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsurfacespeciescount/RM_GetSurfaceSpeciesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsurfacespeciesname/RM_GetSurfaceSpeciesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsurfacetype/RM_GetSurfaceType/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsurfacename/RM_GetSurfaceName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getequilibriumphasescount/RM_GetEquilibriumPhasesCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getequilibriumphasesname/RM_GetEquilibriumPhasesName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getgascomponentscount/RM_GetGasComponentsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getgascomponentsname/RM_GetGasComponentsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getkineticreactionscount/RM_GetKineticReactionsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getkineticreactionsname/RM_GetKineticReactionsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsolidsolutioncomponentscount/RM_GetSolidSolutionComponentsCount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsolidsolutioncomponentsname/RM_GetSolidSolutionComponentsName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsolidsolutionname/RM_GetSolidSolutionName/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsicount/RM_GetSICount/g" %%f 
  @del sed*
)
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/rm_getsiname/RM_GetSIName/g" %%f 
  @del sed*
)

@REM Rename PhreeqcRM to mixed case 
@FOR %%f IN (.\html\*.html) DO @(
  @sed -i "s/\bphreeqcrm\b/PhreeqcRM/g" %%f
  @del sed*
)

@REM Change Titles
@sed -i "s/PhreeqcRM Module Reference/PhreeqcRM Fortran Module Reference/g" ./html/classphreeqcrm.html
@sed -i "s/RM_interface_C.h File Reference/PhreeqcRM C Reference/g" ./html/_r_m__interface___c_8h.html
@sed -i "s/PhreeqcRM Class Reference/PhreeqcRM C++ Class Reference/g" ./html/class_phreeqc_r_m.html

echo finished
