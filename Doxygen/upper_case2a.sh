#!/bin/sh
sed -i "s/getgridcellcountyaml/GetGridCellCountYAML/g" $(egrep -l getgridcellcountyaml ./html/*)
sed -i "s/setgridcellcountyaml/SetGridCellCountYAML/g" $(egrep -l setgridcellcountyaml ./html/*)
#!/bin/sh
for f in ./html/BMI__interface_8F90.html ./html/BMI__interface_8F90.js \
./html/functions_func_g.html ./html/functions_g.html ./html/namespacebmiphreeqcrm.html \
./html/namespacebmiphreeqcrm.js ./html/namespacemembers_func_g.html \
./html/namespacemembers_g.html ./html/structbmiphreeqcrm_1_1bmi-members.html \
./html/structbmiphreeqcrm_1_1bmi.html ./html/structbmiphreeqcrm_1_1bmi.js; do
sed -i "s/abort/Abort/g" $f
sed -i "s/closefiles/CloseFiles/g" $f
sed -i "s/concentrations2utility/Concentrations2Utility/g" $f
sed -i "s/create/Create/g" $f
sed -i "s/createmapping/CreateMapping/g" $f
sed -i "s/decodeerror/DecodeError/g" $f
sed -i "s/destroy/Destroy/g" $f
sed -i "s/dumpmodule/DumpModule/g" $f
sed -i "s/errormessage/ErrorMessage/g" $f
sed -i "s/findcomponents/FindComponents/g" $f
sed -i "s/getbackwardmapping/GetBackwardMapping/g" $f
sed -i "s/getchemistrycellcount/GetChemistryCellCount/g" $f
sed -i "s/getcomponentcount/GetComponentCount/g" $f
sed -i "s/getcomponent/GetComponent/g" $f
sed -i "s/getconcentrations/GetConcentrations/g" $f
sed -i "s/getcurrentselectedoutputusernumber/GetCurrentSelectedOutputUserNumber/g" $f
sed -i "s/getdensitycalculated/GetDensityCalculated/g" $f
sed -i "s/getdensity/GetDensity/g" $f
sed -i "s/getendcell/GetEndCell/g" $f
sed -i "s/geterrorstringlength/GetErrorStringLength/g" $f
sed -i "s/geterrorstring/GetErrorString/g" $f
sed -i "s/getfileprefix/GetFilePrefix/g" $f
sed -i "s/getgascompmoles/GetGasCompMoles/g" $f
sed -i "s/getgascomppressures/GetGasCompPressures/g" $f
sed -i "s/getgascompphi/GetGasCompPhi/g" $f
sed -i "s/getgasphasevolume/GetGasPhaseVolume/g" $f
sed -i "s/setgascompmoles/SetGasCompMoles/g" $f
sed -i "s/setgasphasevolume/SetGasPhaseVolume/g" $f
sed -i "s/getspecieslog10molalities/GetSpeciesLog10Molalities/g" $f
sed -i "s/getgfw/GetGfw/g" $f
sed -i "s/getgridcellcount/GetGridCellCount/g" $f
sed -i "s/getiphreeqcid/GetIPhreeqcId/g" $f
sed -i "s/getmpimyself/GetMpiMyself/g" $f
sed -i "s/getmpitasks/GetMpiTasks/g" $f
sed -i "s/getnthselectedoutputusernumber/GetNthSelectedOutputUserNumber/g" $f
sed -i "s/getporosity/GetPorosity/g" $f
sed -i "s/getpressure/GetPressure/g" $f
sed -i "s/getsaturationcalculated/GetSaturationCalculated/g" $f
sed -i "s/getsaturation/GetSaturation/g" $f
sed -i "s/getselectedoutputcolumncount/GetSelectedOutputColumnCount/g" $f
sed -i "s/getselectedoutputcount/GetSelectedOutputCount/g" $f
sed -i "s/getselectedoutputheading/GetSelectedOutputHeading/g" $f
sed -i "s/getselectedoutputrowcount/GetSelectedOutputRowCount/g" $f
sed -i "s/getselectedoutput/GetSelectedOutput/g" $f
sed -i "s/getsolutionvolume/GetSolutionVolume/g" $f
sed -i "s/getspeciesconcentrations/GetSpeciesConcentrations/g" $f
sed -i "s/getspeciescount/GetSpeciesCount/g" $f
sed -i "s/getspeciesd25/GetSpeciesD25/g" $f
sed -i "s/getspeciesname/GetSpeciesName/g" $f
sed -i "s/getspeciessaveon/GetSpeciesSaveOn/g" $f
sed -i "s/getspeciesz/GetSpeciesZ/g" $f
sed -i "s/getstartcell/GetStartCell/g" $f
sed -i "s/getthreadcount/GetThreadCount/g" $f
sed -i "s/gettime/GetTime/g" $f
sed -i "s/gettimeconversion/GetTimeConversion/g" $f
sed -i "s/gettimestep/GetTimeStep/g" $f
sed -i "s/initializeyaml/InitializeYAML/g" $f
sed -i "s/initialphreeqc2concentrations/InitialPhreeqc2Concentrations/g" $f
sed -i "s/initialphreeqc2module/InitialPhreeqc2Module/g" $f
sed -i "s/initialphreeqccell2module/InitialPhreeqcCell2Module/g" $f
sed -i "s/initialphreeqc2speciesconcentrations/InitialPhreeqc2SpeciesConcentrations/g" $f
sed -i "s/loaddatabase/LoadDatabase/g" $f
sed -i "s/logmessage/LogMessage/g" $f
sed -i "s/mpiworkerbreak/MpiWorkerBreak/g" $f
sed -i "s/mpiworker/MpiWorker/g" $f
sed -i "s/openfiles/OpenFiles/g" $f
sed -i "s/outputmessage/OutputMessage/g" $f
sed -i "s/runcells/RunCells/g" $f
sed -i "s/runfile/RunFile/g" $f
sed -i "s/runstring/RunString/g" $f
sed -i "s/screenmessage/ScreenMessage/g" $f
sed -i "s/setcomponenth2o/SetComponentH2O/g" $f
sed -i "s/setconcentrations/SetConcentrations/g" $f
sed -i "s/setcurrentselectedoutputusernumber/SetCurrentSelectedOutputUserNumber/g" $f
sed -i "s/setdensityuser/SetDensityUser/g" $f
sed -i "s/setdensity/SetDensity/g" $f
sed -i "s/setdumpfilename/SetDumpFileName/g" $f
sed -i "s/seterrorhandlermode/SetErrorHandlerMode/g" $f
sed -i "s/seterroron/SetErrorOn/g" $f
sed -i "s/setfileprefix/SetFilePrefix/g" $f
sed -i "s/setmpiworkercallback/SetMpiWorkerCallback/g" $f
sed -i "s/setnthselectedoutput/SetNthSelectedOutput/g" $f
sed -i "s/setpartitionuzsolids/SetPartitionUZSolids/g" $f
sed -i "s/setporosity/SetPorosity/g" $f
sed -i "s/setprintchemistrymask/SetPrintChemistryMask/g" $f
sed -i "s/setprintchemistryon/SetPrintChemistryOn/g" $f
sed -i "s/setpressure/SetPressure/g" $f
sed -i "s/setrebalancefraction/SetRebalanceFraction/g" $f
sed -i "s/setrebalancebycell/SetRebalanceByCell/g" $f
sed -i "s/setrepresentativevolume/SetRepresentativeVolume/g" $f
sed -i "s/setsaturationuser/SetSaturationUser/g" $f
sed -i "s/setsaturation/SetSaturation/g" $f
sed -i "s/setscreenon/SetScreenOn/g" $f
sed -i "s/setselectedoutputon/SetSelectedOutputOn/g" $f
sed -i "s/setspeciessaveon/SetSpeciesSaveOn/g" $f
sed -i "s/settemperature/SetTemperature/g" $f
sed -i "s/settimeconversion/SetTimeConversion/g" $f
sed -i "s/settimestep/SetTimeStep/g" $f
sed -i "s/settime/SetTime/g" $f
sed -i "s/setunitsexchange/SetUnitsExchange/g" $f
sed -i "s/setunitsgasphase/SetUnitsGasPhase/g" $f
sed -i "s/setunitskinetics/SetUnitsKinetics/g" $f
sed -i "s/setunitsppassemblage/SetUnitsPPassemblage/g" $f
sed -i "s/setunitssolution/SetUnitsSolution/g" $f
sed -i "s/setunitsssassemblage/SetUnitsSSassemblage/g" $f
sed -i "s/setunitsssassemblage/SetUnitsSSassemblage/g" $f
sed -i "s/setunitssurface/SetUnitsSurface/g" $f
sed -i "s/speciesconcentrations2module/SpeciesConcentrations2Module/g" $f
sed -i "s/usesolutiondensityvolume/UseSolutionDensityVolume/g" $f
sed -i "s/warningmessage/WarningMessage/g" $f
sed -i "s/getspecieslog10gammas/GetSpeciesLog10Gammas/g"  $f
sed -i "s/getexchangespeciescount/GetExchangeSpeciesCount/g"  $f
sed -i "s/getexchangespeciesname/GetExchangeSpeciesName/g"  $f
sed -i "s/getexchangename/GetExchangeName/g"  $f
sed -i "s/getsurfacespeciescount/GetSurfaceSpeciesCount/g"  $f
sed -i "s/getsurfacespeciesname/GetSurfaceSpeciesName/g"  $f
sed -i "s/getsurfacetype/GetSurfaceType/g"  $f
sed -i "s/getsurfacename/GetSurfaceName/g"  $f
sed -i "s/getequilibriumphasescount/GetEquilibriumPhasesCount/g"  $f
sed -i "s/getequilibriumphasesname/GetEquilibriumPhasesName/g"  $f
sed -i "s/getgascomponentscount/GetGasComponentsCount/g"  $f
sed -i "s/getgascomponentsname/GetGasComponentsName/g"  $f
sed -i "s/getkineticreactionscount/GetKineticReactionsCount/g"  $f
sed -i "s/getkineticreactionsname/GetKineticReactionsName/g"  $f
sed -i "s/getsolidsolutioncomponentscount/GetSolidSolutionComponentsCount/g"  $f
sed -i "s/getsolidsolutioncomponentsname/GetSolidSolutionComponentsName/g"  $f
sed -i "s/getsolidsolutionname/GetSolidSolutionName/g"  $f
sed -i "s/getsicount/GetSICount/g"  $f
sed -i "s/getsiname/GetSIName/g"  $f
sed -i "s/statesave/StateSave/g"  $f
sed -i "s/stateapply/StateApply/g"  $f
sed -i "s/statedelete/StateDelete/g"  $f
sed -i "s/getcomponents/GetComponents/g" $f
sed -i "s/getgridcellcountyaml/GetGridCellCountYAML/g" $f
sed -i "s/getselectedoutputheadings/GetSelectedOutputHeadings/g" $f
sed -i "s/gettemperature/GetTemperature/g" $f
sed -i "s/getviscosity/GetViscosity/g" $f
sed -i "s/getithconcentration/GetIthConcentration/g" $f
sed -i "s/getithspeciesconcentration/GetIthSpeciesConcentration/g" $f
sed -i "s/setithconcentration/SetIthConcentration/g" $f
sed -i "s/setithspeciesconcentration/SetIthSpeciesConcentration/g" $f
sed -i "s/initialsolutions2module/InitialSolutions2Module/g" $f
sed -i "s/initialequilibriumphases2module/InitialEquilibriumPhases2Module/g" $f
sed -i "s/initialexchanges2module/InitialExchanges2Module/g" $f
sed -i "s/initialsurfaces2module/InitialSurfaces2Module/g" $f
sed -i "s/initialgasphases2module/InitialGasPhases2Module/g" $f
sed -i "s/initialsolidsolutions2module/InitialSolidSolutions2Module/g" $f
sed -i "s/initialkinetics2module/InitialKinetics2Module/g" $f
done
echo finished