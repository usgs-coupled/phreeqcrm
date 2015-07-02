set RM_INTERFACE_F90=..\src\RM_interface.F90
REM sed -i -e "s/PhreeqcRM/phreeqcrm/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_Abort/rm_abort/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_Abort/rm_abort/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_CloseFiles/rm_closefiles/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_Concentrations2Utility/rm_concentrations2utility/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_CreateMapping/rm_createmapping/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_Create/rm_create/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_DecodeError/rm_decodeerror/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_Destroy/rm_destroy/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_DumpModule/rm_dumpmodule/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_ErrorMessage/rm_errormessage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_FindComponents/rm_findcomponents/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetBackwardMapping/rm_getbackwardmapping/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetChemistryCellCount/rm_getchemistrycellcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetComponentCount/rm_getcomponentcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetComponent/rm_getcomponent/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetConcentrations/rm_getconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetDensity/rm_getdensity/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetErrorString/rm_geterrorstring/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetErrorStringLength/rm_geterrorstringlength/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetFilePrefix/rm_getfileprefix/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetGfw/rm_getgfw/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetGridCellCount/rm_getgridcellcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetIPhreeqcId/rm_getiphreeqcid/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetMpiMyself/rm_getmpimyself/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetMpiTasks/rm_getmpitasks/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetNthSelectedOutputUserNumber/rm_getnthselectedoutputusernumber/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSaturation/rm_getsaturation/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSelectedOutputColumnCount/rm_getselectedoutputcolumncount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSelectedOutputCount/rm_getselectedoutputcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSelectedOutputHeading/rm_getselectedoutputheading/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSelectedOutputRowCount/rm_getselectedoutputrowcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSelectedOutput/rm_getselectedoutput/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSolutionVolume/rm_getsolutionvolume/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesConcentrations/rm_getspeciesconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesCount/rm_getspeciescount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesD25/rm_getspeciesd25/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesLogGammas/rm_getspeciesloggammas/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesName/rm_getspeciesname/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesSaveOn/rm_getspeciessaveon/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSpeciesZ/rm_getspeciesz/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSurfaceDiffuseLayerArea/rm_getsurfacediffuselayerarea/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSurfaceDiffuseLayerConcentrations/rm_getsurfacediffuselayerconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSurfaceDiffuseLayerCount/rm_getsurfacediffuselayercount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSurfaceDiffuseLayerName/rm_getsurfacediffuselayername/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetSurfaceDiffuseLayerThickness/rm_getsurfacediffuselayerthickness/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetThreadCount/rm_getthreadcount/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetTimeConversion/rm_gettimeconversion/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetTimeStep/rm_gettimestep/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_GetTime/rm_gettime/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_InitialPhreeqc2Concentrations/rm_initialphreeqc2concentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_InitialPhreeqc2Module/rm_initialphreeqc2module/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_InitialPhreeqcCell2Module/rm_initialphreeqccell2module/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_InitialPhreeqc2SpeciesConcentrations/rm_initialphreeqc2speciesconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_LoadDatabase/rm_loaddatabase/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_LogMessage/rm_logmessage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_MpiWorkerBreak/rm_mpiworkerbreak/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_MpiWorker/rm_mpiworker/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_OpenFiles/rm_openfiles/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_OutputMessage/rm_outputmessage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_RunCells/rm_runcells/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_RunFile/rm_runfile/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_RunString/rm_runstring/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_ScreenMessage/rm_screenmessage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetComponentH2O/rm_setcomponenth2o/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetConcentrations/rm_setconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetCurrentSelectedOutputUserNumber/rm_setcurrentselectedoutputusernumber/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetDensity/rm_setdensity/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetDumpFileName/rm_setdumpfilename/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetErrorHandlerMode/rm_seterrorhandlermode/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetFilePrefix/rm_setfileprefix/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetMpiWorkerCallback/rm_setmpiworkercallback/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetPartitionUZSolids/rm_setpartitionuzsolids/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetPorosity/rm_setporosity/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetPrintChemistryMask/rm_setprintchemistrymask/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetPrintChemistryOn/rm_setprintchemistryon/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetPressure/rm_setpressure/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetRebalanceFraction/rm_setrebalancefraction/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetRebalanceByCell/rm_setrebalancebycell/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetRepresentativeVolume/rm_setrepresentativevolume/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetSaturation/rm_setsaturation/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetSelectedOutputOn/rm_setselectedoutputon/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetSpeciesSaveOn/rm_setspeciessaveon/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetSurfaceDiffuseLayerConcentrations/rm_setsurfacediffuselayerconcentrations/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetTemperature/rm_settemperature/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetTimeConversion/rm_settimeconversion/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetTimeStep/rm_settimestep/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetTime/rm_settime/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsExchange/rm_setunitsexchange/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsGasPhase/rm_setunitsgasphase/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsKinetics/rm_setunitskinetics/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsPPassemblage/rm_setunitsppassemblage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsSolution/rm_setunitssolution/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsSSassemblage/rm_setunitsssassemblage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsSSassemblage/rm_setunitsssassemblage/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SetUnitsSurface/rm_setunitssurface/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_SpeciesConcentrations2Module/rm_speciesconcentrations2module/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_UseSolutionDensityVolume/rm_usesolutiondensityvolume/g" %RM_INTERFACE_F90%
sed -i -e "s/RM_WarningMessage/rm_warningmessage/g" %RM_INTERFACE_F90%
del sed*