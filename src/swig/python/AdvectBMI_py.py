import phreeqcrm
import numpy as np
#ifdef USE_YAML
    #module mydata
    #  double precision, dimension(:), pointer :: K_ptr
    #  integer                                 :: rm_id
    #end module mydata

def AdvectBMI_py():

	# Based on PHREEQC Example 11

	# --------------------------------------------------------------------------
	# Create PhreeqcRM
	# --------------------------------------------------------------------------
	yaml_file = "AdvectBMI_py.yaml"
	
	# phreeqc_rm.GetGridCellCountYAML must be called BEFORE
	# the PhreeqcRM instance is created. The
	# return value can be used to create the
	# PhreeqcRM instance.
	#
	# If the YAML file does not contain
	# a node "SetGridCellCount:" (usually written
	# using the YAMLPhreeqcRM class and the method
	# YAMLSetGridCellCount), the return
	# value is zero.
	###nxyz = GetGridCellCountYAML(yaml_file)
	#nxyz = 40
	# Bogus conductivity field for Basic callback demonstration
	#hydraulic_K = [i*2.0 for i in range(nxyz)] 
	#nthreads = 3
	
	bmi = phreeqcrm.BMIPhreeqcRM(40, 3)
	# Initialize with YAML file
	status = bmi.Initialize(yaml_file)

    # Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
    #TODO CALL register_basic_callback_fortran()
#ifdef USE_MPI
    # Optional callback for MPI
    #TODO status = do_something()   # only root is calling do_something here
#endif
	#nxyz = 0
	components = bmi.GetValue("Components")
	ncomps = bmi.GetValue("ComponentCount")
	nxyz = bmi.GetValue("GridCellCount")
	nxyz1 = bmi.GetGridCellCount()
	print(f"nxyz1={nxyz1}")
	time = bmi.GetValue("Time")
	print(f"time={time}")

	hydraulic_K = [0.0] * nxyz
	for i in range(nxyz):
		hydraulic_K[i] = i*2.0
	# Print some of the reaction module information
	nthreads = bmi.GetThreadCount()
	#print(f"Number of threads:                                {nthreads}")
	string1 = f"Number of threads:                                {nthreads}"
	print(string1)
	string1 = f"MPI task number:                                  {bmi.GetMpiMyself()}"
	print(string1)
	string1 = f"File prefix:                                      {bmi.GetValue('FilePrefix')}"
	print(string1)
	string1 = f"Number of grid cells in the user's model:         {nxyz}"
	print(string1)
	nchem = bmi.GetChemistryCellCount()
	string1 = f"Number of chemistry cells in the reaction module: {nchem}"
	print(string1)
	string1 = f"Number of components for transport:               {ncomps}"
	print(string1)
	# Get component information)
	gfw = bmi.GetValue("Gfw")
	for i in range(ncomps):
		print(f"{components[i].rjust(10,' ')}  {gfw[i]}")
	print()
	# Get initial temperatures
	temperature = bmi.GetValue("Temperature")
	# Get initial temperature
	sat = bmi.GetValue("Saturation")
	# Get initial porosity
	por = bmi.GetValue("Porosity")
	# Get initial temperature
	volume = bmi.GetValue("SolutionVolume")
	# Get initial concentrations
	c = bmi.GetValue("Concentrations")
	#c_dbl_vect = phreeqcrm.DoubleVector(nxyz * len(components), 0.0)
	#for i in range(nxyz):
	#	c_dbl_vect[i] = c[i]
	# Set density, pressure, and temperature (previously allocated)
	density = [1.0] * nxyz
	bmi.SetValue("Density", density)
	pressure = [2.0] * nxyz
	bmi.SetValue("Pressure", pressure)  
	temperature = [20.0] * nxyz
	bmi.SetValue("Temperature", temperature)  
    # --------------------------------------------------------------------------
    # Set boundary condition
    # --------------------------------------------------------------------------
	nbound = 1
	bc1 = [0]           # solution 0 from Initial IPhreeqc instance
	bc2 = [-1]          # no bc2 solution for mixing
	bc_f1 = [1.0]       # mixing fraction for bc1
	#bc_conc = phreeqcrm.DoubleVector()
	#phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, nbound, bc1, bc2, bc_f1)
	#bc_conc = bmi.InitialPhreeqc2Concentrations(bc1, bc2, bc_f1)
	out = bmi.InitialPhreeqc2Concentrations(bc1)
	status = out[0]
	bc_conc = list(out[1])
	
	#for i in range(nxyz*ncomps):
	#	c_dbl_vect[i] = c[i]
    # --------------------------------------------------------------------------
    # Transient loop
    # --------------------------------------------------------------------------
	nsteps = 10
	time = 0.0
	bmi.SetValue("Time", time)
	time_step = 86400.0
	bmi.SetValue("TimeStep", time_step)  
	bmi.SetScreenOn(True)	
	for step in range(nsteps):
		print(f"Beginning transport calculation {time*bmi.GetTimeConversion()} days")
		print(f"          Time step             {time_step*bmi.GetTimeConversion()} days")
		# Transport calculation here, changes c
		advectionbmi_py(c, bc_conc, ncomps, nxyz, nbound)

		# Print selected output and chemistry on last step
		if (step == (nsteps - 1)):
			print_chemistry_on = True
		else:
			print_chemistry_on = False
		bmi.SetValue("SelectedOutputOn", print_chemistry_on)
		bmi.SetPrintChemistryOn(print_chemistry_on, False, False)  # workers, initial_phreeqc, utility
		time += time_step
		status = bmi.SetTime(time)
		# Transfer data to PhreeqcRM after transport      
		bmi.SetValue("Concentrations", c)   # Transported concentrations
		# Optionally, if values changed during transport
		bmi.SetValue("Porosity", por)              
		bmi.SetValue("Saturation", sat)            
		bmi.SetValue("Temperature", temperature) 
		bmi.SetValue("Pressure", pressure)          
		bmi.SetValue("TimeStep", time_step) 
		
		# Run cells with transported conditions
		print(f"Beginning reaction calculation  {time*bmi.GetTimeConversion()} days")
		bmi.Update()

		# Get new data calculated by PhreeqcRM for transport
		c = bmi.GetValue("Concentrations")
		density = bmi.GetValue("Density")
		volume = bmi.GetValue("SolutionVolume")   
		# Print results at last time step
		if (step == (nsteps - 1)):
			#print("Current distribution of cells for workers")
			#print("Worker      First cell        Last Cell")
			#n = phreeqc_rm.GetThreadCount() * phreeqc_rm.GetMpiTasks()
			#sc = phreeqc_rm.IntVector()
			#ec = phreeqc_rm.IntVector()
			#phreeqc_rm.GetStartCell(sc)
			#phreeqc_rm.GetEndCell(ec)
			#for i in range(n):
			#	print(i,"           ", sc(i),"                 ",ec(i))
			
			# Loop through possible multiple selected output definitions
			n = bmi.GetValue("SelectedOutputCount")
			for isel in range(n): 
				i = isel
				bmi.SetValue("NthSelectedOutput", i)
				n_user = bmi.GetValue("CurrentSelectedOutputUserNumber")
				print(f"Selected output sequence number: {isel}")
				print(f"Selected output user number:     {n_user}")
				# Get 2D array of selected output values
				col = bmi.GetValue("SelectedOutputColumnCount")
				rows = bmi.GetValue("SelectedOutputRowCount")
				#selected_out = phreeqc_rm.DoubleVector()
				# Get headings
				headings = bmi.GetValue("SelectedOutputHeadings")
				# Get selected output
				selected_out = bmi.GetValue("SelectedOutput")
				# Print results
				for i in range(rows//2):
					print("Cell number ", i)
					print(f"     Density:    f{density[i]}")
					print(f"     Volume:     {volume[i]}")
					print(f"     Components: ")
					for j in range(ncomps):
						print(f"{j}, {components[j]}: {c[j * nxyz + i]}")
					print(f"     Selected output: ")
					for j in range(col):
						print(f"{j}, {headings[j]}, {selected_out[j * nxyz + i]}")
	# Clean up
	bmi.Finalize()						

def advectionbmi_py(c, bc_conc, ncomps, nxyz, dim):
    # Advect
    for i in range(nxyz - 1, 0, -1):
        for j in range(ncomps):
            c[j * nxyz + i] = c[j * nxyz + i - 1]              # component j
    
    # Cell zero gets boundary condition
    for j in range(ncomps):
        c[j * nxyz] = bc_conc[j * dim];                        # component j


if __name__ == '__main__':
	AdvectBMI_py()

#endif # USE_YAML
