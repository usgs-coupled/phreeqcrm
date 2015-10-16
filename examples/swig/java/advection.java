public class advection extends BasicCallback
{
	// data members
	public PhreeqcRM phreeqcrm;
	public int nxyz;
	public int threads;
	public DoubleVector hydraulic_K;


	// load the dll/so
	static {
		try {
			System.loadLibrary("phreeqcrm_java");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. See the chapter on Dynamic Linking Problems in the SWIG Java documentation for help.\n" + e);
			System.exit(1);
		}
	}

	public advection() {
		phreeqcrm   = null;
		nxyz        = 40;
		threads     = 3;
		hydraulic_K = null;
	}

    public double Callback(double x1, double x2, String str) {
		/*
		Use of a callback is optional.

		The callback provides a way to obtain data from a Basic program
		through the variables x1, x2, and str1, and send data to a
		Basic program through the return value of the callback.

		The callback function is called whenever CALLBACK(x1, x2, str$)
		is used in a Basic program (usually USER_PUNCH). See file "ic".
		*/
		int rm_cell_number = (int) x1;
		if (rm_cell_number >= 0 && rm_cell_number < phreeqcrm.GetChemistryCellCount()) {
			IntVectorVector back = phreeqcrm.GetBackwardMapping();
			if (str.compareTo("HYDRAULIC_K") == 0) {
				return hydraulic_K.get(back.get(rm_cell_number).get(0));
			}
		}
		return -999.9;
	}

	public void Exec() {

		// Based on PHREEQC Example 11

		// --------------------------------------------------------------------------
		// Create PhreeqcRM
		// --------------------------------------------------------------------------

		hydraulic_K = new DoubleVector();
		for (int i = 0; i < nxyz; ++i) {
			hydraulic_K.add(i*2.0);
		}

		// Create module
		phreeqcrm = new PhreeqcRM(nxyz, threads);

		IRM_RESULT status;

		// Set properties
		status = phreeqcrm.SetErrorHandlerMode(0);
		status = phreeqcrm.SetComponentH2O(false);
		status = phreeqcrm.SetRebalanceFraction(0.5);
		status = phreeqcrm.SetRebalanceByCell(true);
		phreeqcrm.UseSolutionDensityVolume(false);
		phreeqcrm.SetPartitionUZSolids(false);

		// Open files
		status = phreeqcrm.SetFilePrefix("Advect_cpp");
		phreeqcrm.OpenFiles();

		// Set concentration units
		status = phreeqcrm.SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqcrm.SetUnitsPPassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqcrm.SetUnitsExchange(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqcrm.SetUnitsSurface(1);            // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqcrm.SetUnitsGasPhase(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqcrm.SetUnitsSSassemblage(1);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqcrm.SetUnitsKinetics(1);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		// Set conversion from seconds to user units (days)
		double time_conversion = 1.0 / 86400;
		status = phreeqcrm.SetTimeConversion(time_conversion);

		// Set representative volume
		DoubleVector rv = new DoubleVector(nxyz);
		for (int i=0; i < nxyz; ++i) {
			rv.set(i, 1.0);
		}
		status = phreeqcrm.SetRepresentativeVolume(rv);

		// Set initial porosity
		DoubleVector por = new DoubleVector(nxyz);
		for (int i=0; i < nxyz; ++i) {
			por.set(i, 0.2);
		}
		status = phreeqcrm.SetPorosity(por);

		// Set initial saturation
		DoubleVector sat = new DoubleVector(nxyz);
		for (int i=0; i < nxyz; ++i) {
			sat.set(i, 1.0);
		}
		status = phreeqcrm.SetSaturation(sat);

		// Set cells to print chemistry when print chemistry is turned on
		IntVector print_chemistry_mask = new IntVector(nxyz);
		for (int i = 0; i < nxyz; i++) {
			print_chemistry_mask.set(i, 0);
		}
		for (int i = 0; i < nxyz/2; i++) {
			print_chemistry_mask.set(i, 1);
		}
		status = phreeqcrm.SetPrintChemistryMask(print_chemistry_mask);

		// test getters
		IntVector print_chemistry_mask1 = phreeqcrm.GetPrintChemistryMask();
		BoolVector print_on             = phreeqcrm.GetPrintChemistryOn();
		boolean rebalance               = phreeqcrm.GetRebalanceByCell();
		double f_rebalance              = phreeqcrm.GetRebalanceFraction();
		boolean so_on                   = phreeqcrm.GetSelectedOutputOn();
		int units_exchange              = phreeqcrm.GetUnitsExchange();
		int units_gas_phase             = phreeqcrm.GetUnitsGasPhase();
		int units_kinetics              = phreeqcrm.GetUnitsKinetics();
		int units_pp_assemblage         = phreeqcrm.GetUnitsPPassemblage();
		int units_solution              = phreeqcrm.GetUnitsSolution();
		int units_ss_exchange           = phreeqcrm.GetUnitsSSassemblage();
		int units_surface               = phreeqcrm.GetUnitsSurface();

		// Demonstation of mapping, two equivalent rows by symmetry
		IntVector grid2chem = new IntVector(nxyz);
		for (int i = 0; i < nxyz; ++i) {
			grid2chem.set(i, -1);
		}

		for (int i = 0; i < nxyz/2; i++) {
			grid2chem.set(i, i);
			grid2chem.set(i + nxyz/2, i);
		}
		status = phreeqcrm.CreateMapping(grid2chem);
		if (status != IRM_RESULT.IRM_OK) phreeqcrm.DecodeError(status.swigValue());
		int nchem = phreeqcrm.GetChemistryCellCount();

		// --------------------------------------------------------------------------
		// Set initial conditions
		// --------------------------------------------------------------------------

		// Set printing of chemistry file
		status = phreeqcrm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility

		// Load database
		phreeqcrm.LoadDatabase("phreeqc.dat");

		// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
		register_basic_callback();

		// Demonstration of error handling if ErrorHandlerMode is 0
		if (status != IRM_RESULT.IRM_OK) {
			System.err.println(phreeqcrm.GetErrorString());
			System.exit(1);
		}

		// Run file to define solutions and reactants for initial conditions, selected output
		boolean workers         = true;     // Worker instances do the reaction calculations for transport
		boolean initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		boolean utility         = true;     // Utility instance is available for processing
		status = phreeqcrm.RunFile(workers, initial_phreeqc, utility, "advect.pqi");

		// Clear contents of workers and utility
		initial_phreeqc = false;
		String input = "DELETE; -all";
		status = phreeqcrm.RunString(workers, initial_phreeqc, utility, input);

		// Determine number of components to transport
		int ncomps = phreeqcrm.FindComponents();

		// Print some of the reaction module information
		{
			StringBuffer sb = new StringBuffer();
			sb.append("Database:                                         " + phreeqcrm.GetDatabaseFileName() + "\n");
			sb.append("Number of threads:                                " + phreeqcrm.GetThreadCount() + "\n");
			sb.append("Number of MPI processes:                          " + phreeqcrm.GetMpiTasks() + "\n");
			sb.append("MPI task number:                                  " + phreeqcrm.GetMpiMyself() + "\n");
			sb.append("File prefix:                                      " + phreeqcrm.GetFilePrefix() + "\n");
			sb.append("Number of grid cells in the user's model:         " + phreeqcrm.GetGridCellCount() + "\n");
			sb.append("Number of chemistry cells in the reaction module: " + phreeqcrm.GetChemistryCellCount() + "\n");
			sb.append("Number of components for transport:               " + phreeqcrm.GetComponentCount() + "\n");
			sb.append("Partioning of UZ solids:                          " + phreeqcrm.GetPartitionUZSolids() + "\n");
			sb.append("Error handler mode:                               " + phreeqcrm.GetErrorHandlerMode() + "\n");
			phreeqcrm.OutputMessage(sb.toString());
		}

		IntVector f_map = phreeqcrm.GetForwardMapping();
		// Get component information
		StringVector components = phreeqcrm.GetComponents();
		DoubleVector gfw = phreeqcrm.GetGfw();
		for (int i = 0; i < ncomps; i++) {
			phreeqcrm.OutputMessage(String.format("%10s    %10f \n", components.get(i), gfw.get(i)));
		}
		phreeqcrm.OutputMessage("\n");

		// Set array of initial conditions
		IntVector ic1 = new IntVector(nxyz*7);
		IntVector ic2 = new IntVector(nxyz*7);
		DoubleVector f1 = new DoubleVector(nxyz*7);
		for (int i = 0; i < nxyz*7; ++i) {
			ic1.set(i, -1);
			ic2.set(i, -1);
			f1.set(i, 1.0);
		}
		for (int i = 0; i < nxyz; ++i) {
			ic1.set(i,           1);    // Solution 1
			ic1.set(nxyz + i,   -1);    // Equilibrium phases none
			ic1.set(2*nxyz + i,  1);    // Exchange 1
			ic1.set(3*nxyz + i, -1);    // Surface none
			ic1.set(4*nxyz + i, -1);    // Gas phase none
			ic1.set(5*nxyz + i, -1);    // Solid solutions none
			ic1.set(6*nxyz + i, -1);    // Kinetics none
		}
		status = phreeqcrm.InitialPhreeqc2Module(ic1, ic2, f1);

		// No mixing is defined, so the following is equivalent
		// status = phreeqcrm.InitialPhreeqc2Module(ic1.data());

		// alternative for setting initial conditions
		// cell number in first argument (-1 indicates last solution, 40 in this case)
		// in advect.pqi and any reactants with the same number--
		// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
		// will be written to cells 18 and 19 (0 based)
		IntVector module_cells = new IntVector();
		module_cells.add(18);
		module_cells.add(19);
		status = phreeqcrm.InitialPhreeqcCell2Module(-1, module_cells);

		// Get temperatures
		DoubleVector tempc = phreeqcrm.GetTemperature();

		// get current saturation
		DoubleVector current_sat = new DoubleVector();
		status = phreeqcrm.GetSaturation(current_sat);

		// Initial equilibration of cells
		double time = 0.0;
		double time_step = 0.0;
		DoubleVector c = new DoubleVector(nxyz * components.size());
		status = phreeqcrm.SetTime(time);
		status = phreeqcrm.SetTimeStep(time_step);
		status = phreeqcrm.RunCells();
		status = phreeqcrm.GetConcentrations(c);

		// --------------------------------------------------------------------------
		// Set boundary condition
		// --------------------------------------------------------------------------

		int nbound = 1;
		DoubleVector bc_conc = new DoubleVector(nbound);
		DoubleVector bc_f1 = new DoubleVector(nbound);
		IntVector bc1 = new IntVector(nbound);
		IntVector bc2 = new IntVector(nbound);

		for (int i = 0; i > 0; ++i) {
			bc1.set(i, 0);                          // solution 0 from Initial IPhreeqc instance
			bc2.set(i, -1);                         // no bc2 solution for mixing
			bc_f1.set(i, 1.0);                      // mixing fraction for bc1
		}
		status = phreeqcrm.InitialPhreeqc2Concentrations(bc_conc, bc1, bc2, bc_f1);

		// --------------------------------------------------------------------------
		// Transient loop
		// --------------------------------------------------------------------------

		int nsteps = 10;
		DoubleVector initial_density = new DoubleVector(nxyz);
		DoubleVector temperature = new DoubleVector(nxyz);
		DoubleVector pressure = new DoubleVector(nxyz);
		for (int i=0; i < nxyz; ++i) {
			initial_density.set(i, 1.0);
			temperature.set(i, 20.0);
			pressure.set(i, 2.0);
		}
		phreeqcrm.SetDensity(initial_density);
		phreeqcrm.SetTemperature(temperature);
		phreeqcrm.SetPressure(pressure);
		time_step = 86400.;
		status = phreeqcrm.SetTimeStep(time_step);
		for (int steps = 0; steps < nsteps; ++steps) {

			// Transport calculation here
			{
				StringBuffer strm = new StringBuffer();
				strm.append("Beginning transport calculation             " + phreeqcrm.GetTime() * phreeqcrm.GetTimeConversion() + " days\n");
				strm.append("          Time step                         " + phreeqcrm.GetTimeStep() * phreeqcrm.GetTimeConversion() + " days\n");
				phreeqcrm.LogMessage(strm.toString());
				phreeqcrm.SetScreenOn(true);
				phreeqcrm.ScreenMessage(strm.toString());
			}
			AdvectCpp(c, bc_conc, ncomps, nxyz, nbound);

			// Transfer data to PhreeqcRM for reactions
			boolean print_selected_output_on = (steps == nsteps - 1) ? true : false;
			boolean print_chemistry_on = (steps == nsteps - 1) ? true : false;
			status = phreeqcrm.SetSelectedOutputOn(print_selected_output_on);
			status = phreeqcrm.SetPrintChemistryOn(print_chemistry_on, false, false); // workers, initial_phreeqc, utility
			status = phreeqcrm.SetPorosity(por);             // If pororosity changes due to compressibility
			status = phreeqcrm.SetSaturation(sat);           // If saturation changes
			status = phreeqcrm.SetTemperature(temperature);  // If temperature changes
			status = phreeqcrm.SetPressure(pressure);        // If pressure changes
			status = phreeqcrm.SetConcentrations(c);         // Transported concentrations
			status = phreeqcrm.SetTimeStep(time_step);		 // Time step for kinetic reactions
			time += time_step;
			status = phreeqcrm.SetTime(time);

			// Run cells with transported conditions
			{
				StringBuffer strm = new StringBuffer();
				strm.append("Beginning reaction calculation              " + time * phreeqcrm.GetTimeConversion() + " days\n");
				phreeqcrm.LogMessage(strm.toString());
				phreeqcrm.ScreenMessage(strm.toString());
			}
			status = phreeqcrm.RunCells();

			// Transfer data from PhreeqcRM for transport
			status = phreeqcrm.GetConcentrations(c);
			DoubleVector density = new DoubleVector();
			status = phreeqcrm.GetDensity(density);
			DoubleVector volume = phreeqcrm.GetSolutionVolume();

			// Print results at last time step
			if (print_chemistry_on)
			{
				{
					StringBuffer oss = new StringBuffer();
					oss.append("Current distribution of cells for workers\n");
					oss.append("Worker      First cell        Last Cell\n");
					int n;
					n = phreeqcrm.GetThreadCount() * phreeqcrm.GetMpiTasks();
					for (int i = 0; i < n; i++)
					{
						oss.append(i + "           " + phreeqcrm.GetStartCell().get(i) + "                 "
							+ phreeqcrm.GetEndCell().get(i) + "\n");
					}
					phreeqcrm.OutputMessage(oss.toString());
				}
				for (int isel = 0; isel < phreeqcrm.GetSelectedOutputCount(); isel++) {
					// Loop through possible multiple selected output definitions
					int n_user = phreeqcrm.GetNthSelectedOutputUserNumber(isel);
					status = phreeqcrm.SetCurrentSelectedOutputUserNumber(n_user);
					System.err.println("Selected output sequence number: " + isel);
					System.err.println("Selected output user number:     " + n_user);

					// Get double array of selected output values
					DoubleVector so = new DoubleVector();
					int col = phreeqcrm.GetSelectedOutputColumnCount();
					status = phreeqcrm.GetSelectedOutput(so);

					// Print results
					for (int i = 0; i < phreeqcrm.GetSelectedOutputRowCount()/2; i++) {
						System.err.println("Cell number " + i);
						System.err.println("     Density: " + density.get(i));
						System.err.println("     Volume:  " + volume.get(i));
						System.err.println("     Components: ");
						for (int j = 0; j < ncomps; j++) {
							System.err.println("          " + j + " " + components.get(j) + ": " + c.get(j*nxyz + i));
						}
						StringVector headings = new StringVector(col);
						System.err.println("     Selected output: ");
						for (int j = 0; j < col; ++j) {
							headings.set(j, phreeqcrm.GetSelectedOutputHeading(j));
							System.err.format("           %d %s: %f\n", j, headings.get(j), so.get(j*nxyz + i));
						}
					}
				}
			}
		}

		// --------------------------------------------------------------------------
		// Additional features and finalize
		// --------------------------------------------------------------------------

 		// Use utility instance of PhreeqcRM to calculate pH of a mixture
		DoubleVector c_well = new DoubleVector(ncomps);
		for (int i = 0; i < ncomps; ++i) {
			c_well.set(i, 0.5 * c.get(0 + nxyz*i) + 0.5 * c.get(9 + nxyz*i));
		}
		DoubleVector tc = new DoubleVector(1);
		DoubleVector p_atm = new DoubleVector(1);
		tc.set(0, 15.0);
		p_atm.set(0, 3.0);

		IPhreeqc util_ptr = phreeqcrm.Concentrations2Utility(c_well, tc, p_atm);
		input = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1";
		VRESULT iphreeqc_result;
		util_ptr.SetOutputFileName("utility_cpp.txt");
		util_ptr.SetOutputFileOn(true);
		int ir = util_ptr.RunString(input);

		// Alternatively, utility pointer is worker nthreads + 1
		IPhreeqc util_ptr1 = phreeqcrm.GetIPhreeqcPointer(phreeqcrm.GetThreadCount() + 1);
		if (ir != 0) {
			phreeqcrm.ErrorHandler(IRM_RESULT.IRM_FAIL.swigValue(), "IPhreeqc RunString failed");
		}
		util_ptr.SetCurrentSelectedOutputUserNumber(5);
		VAR v = new VAR();
		iphreeqc_result = util_ptr.GetSelectedOutputValue(1, 0, v);

		// Dump results
		boolean dump_on = true;
		boolean append = false;
		status = phreeqcrm.SetDumpFileName("advection_cpp.dmp");
		status = phreeqcrm.DumpModule(dump_on, append);    // gz disabled unless compiled with #define USE_GZ

		// Get pointer to worker
		IPhreeqcPhastVector w = phreeqcrm.GetWorkers();
		w.get(0).AccumulateLine("Delete; -all");
		ir = w.get(0).RunAccumulated();

		// Clean up
		status = phreeqcrm.CloseFiles();
		status = phreeqcrm.MpiWorkerBreak();
	}

	public void	AdvectCpp(DoubleVector c, DoubleVector bc_conc, int ncomps, int nxyz, int dim) {
		for (int i = nxyz/2 - 1 ; i > 0; --i) {
			for (int j = 0; j < ncomps; ++j) {
				c.set(j * nxyz + i, c.get(j * nxyz + i - 1));                    // component j
			}
		}
		// Cell zero gets boundary condition
		for (int j = 0; j < ncomps; ++j) {
			c.set(j * nxyz, bc_conc.get(j * dim));                               // component j
		}
	}

	public void register_basic_callback() {
		IPhreeqcPhastVector w = phreeqcrm.GetWorkers();
		for (int i = 0; i < (int) w.size(); i++) {
			w.get(i).SetBasicCallback(this);
		}
	}

	public static void main(String argv[]) {
		advection a = new advection();
		a.Exec();
	}
}
