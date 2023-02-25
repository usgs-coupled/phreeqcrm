import phreeqcrm

def units_tester():
    nxyz = 3
    nthreads = 3
    phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)

    # Set properties
    status = phreeqc_rm.SetErrorOn(True)
    print(status)
    status = phreeqc_rm.SetErrorHandlerMode(1)
    print(status)
    status = phreeqc_rm.SetFilePrefix("Units_InitialPhreeqc_1")
    print(status)

    if (phreeqc_rm.GetMpiMyself() == 0):
        phreeqc_rm.OpenFiles()

    status = phreeqc_rm.SetUnitsSolution(1)      # 1, mg/L; 2, mol/L; 3, kg/kgs
    print(status)
    status = phreeqc_rm.SetUnitsPPassemblage(2)  # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)
    status = phreeqc_rm.SetUnitsExchange(1)      # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)
    status = phreeqc_rm.SetUnitsSurface(1)       # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)
    status = phreeqc_rm.SetUnitsGasPhase(1)      # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)
    status = phreeqc_rm.SetUnitsSSassemblage(1)  # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)
    status = phreeqc_rm.SetUnitsKinetics(1)      # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
    print(status)


    # Set representative volume
    rv = phreeqcrm.vectord(nxyz, 1.0)
    status = phreeqc_rm.SetRepresentativeVolume(rv)
    print(status)

    # Set current porosity
    por = phreeqcrm.vectord(nxyz, 0.2)
    status = phreeqc_rm.SetPorosity(por)
    print(status)

    # Set saturation
    sat = phreeqcrm.vectord(nxyz, 1.0)
    status = phreeqc_rm.SetSaturation(sat)
    print(status)

    # Set printing of chemistry file
    status = phreeqc_rm.SetPrintChemistryOn(False, True, False) # workers, initial_phreeqc, utility
    print(status)

    # --------------------------------------------------------------------------
    # Set initial conditions
    # --------------------------------------------------------------------------

    # Load database
    status = phreeqc_rm.LoadDatabase("phreeqc.dat")
    # Run file to define solutions and reactants for initial conditions, selected output
    workers = True
    initial_phreeqc = True
    utility = False
    status = phreeqc_rm.RunFile(workers, initial_phreeqc, utility, "units.pqi")
    print(status)

    input = "DELETE; -all"
    status = phreeqc_rm.RunString(True, False, True, input)
    print(status)

    status = phreeqc_rm.SetFilePrefix("Units_InitialPhreeqc_2")
    if (phreeqc_rm.GetMpiMyself() == 0):
        phreeqc_rm.OpenFiles()

    # Set reference to components
    ncomps = phreeqc_rm.FindComponents()
    components = phreeqc_rm.GetComponents()
    print(components)

    # Set initial conditions
    cell_numbers = phreeqcrm.vectori()
    cell_numbers.push_back(0)
   
    status = phreeqc_rm.InitialPhreeqcCell2Module(1, cell_numbers)
    print(status)
    cell_numbers[0] = 1
    status = phreeqc_rm.InitialPhreeqcCell2Module(2, cell_numbers)
    print(status)
    cell_numbers[0] = 2
    status = phreeqc_rm.InitialPhreeqcCell2Module(3, cell_numbers)
    print(status)
    # Retrieve concentrations
    c = phreeqcrm.vectord()
    status = phreeqc_rm.SetFilePrefix("Units_Worker");
    print(status)
    if (phreeqc_rm.GetMpiMyself() == 0):
        phreeqc_rm.OpenFiles()

    print_mask = phreeqcrm.vectori(3, 1)
    phreeqc_rm.SetPrintChemistryMask(print_mask)
    phreeqc_rm.SetPrintChemistryOn(True, True, True)
    status = phreeqc_rm.RunCells()
    print(status)

    status = phreeqc_rm.GetConcentrations(c)
    print(status)
    so = phreeqcrm.vectord()
    status = phreeqc_rm.GetSelectedOutput(so)
    print(status)
    print("Cell {}".format("TODO"))
    for i in range(nxyz):
        print("{}   {}".format(i, so[i]))

    # # Use utility instance of PhreeqcRM
    # # std::vector<double> tc, p_atm;
    # #tc.resize(nxyz, 25.0);
    # #p_atm.resize(nxyz, 1.0);
    # # IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c, tc, p_atm);
    # # std::string input;
    # input = "RUN_CELLS; -cells 0-2"
    # # Output goes to new file
    # int iphreeqc_result;
    # util_ptr->SetOutputFileName("Units_utility.out");
    # util_ptr->SetOutputFileOn(true);
    # iphreeqc_result = util_ptr->RunString(input.c_str());

    status = phreeqc_rm.MpiWorkerBreak()
    return 0

# def main():
#     mpi_tasks = 1
# 	mpi_myself = 0

#     rm = phreeqcrm.PhreeqcRM()

if __name__ == "__main__":
    units_tester()