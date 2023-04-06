import phreeqcrm
import numpy as np
#ifdef USE_YAML
    #module mydata
    #  double precision, dimension(:), pointer :: K_ptr
    #  integer                                 :: rm_id
    #end module mydata

def AdvectBMI_f90():

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
	nxyz = 40
	# Bogus conductivity field for Basic callback demonstration
	#hydraulic_K = [i*2.0 for i in range(nxyz)] 
	hydraulic_K = [0.0] * nxyz
	for i in range(nxyz):
		hydraulic_K[i] = i*2.0
	nthreads = 3
	phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, nthreads)
	# Initialize with YAML file
	status = phreeqc_rm.InitializeYAML(yaml_file)
	###phreeqc_rm.BMI_Initialize(yaml_file)

    # Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
    #TODO CALL register_basic_callback_fortran()
#ifdef USE_MPI
    # Optional callback for MPI
    #TODO status = do_something()   # only root is calling do_something here
#endif

	phreeqc_rm.BMI_GetValue("ComponentCount", ncomps)
	# Print some of the reaction module information
	print("Number of threads:                                ", phreeqc_rm.GetThreadCount())
	phreeqc_rm.OutputMessage(string1)
	#write(string1, "(A,I10)") "Number of MPI processes:                          ", phreeqc_rm.GetMpiTasks()
	print("MPI task number:                                  ", phreeqc_rm.GetMpiMyself())
	print("File prefix:                                      ", phreeqc_rm.BMI_GetValue("FilePrefix", prefix))
	print("Number of grid cells in the user's model:         ", nxyz)
	print("Number of chemistry cells in the reaction module: ", nchem)
	print("Number of components for transport:               ", ncomps)
	# Get component information
	phreeqc_rm.BMI_GetValue("Components", components)
	phreeqc_rm.BMI_GetValue("Gfw", gfw)
	for i in range(ncomps):
		print(components(i)), gfw(i)
	print()
	# Get initial temperatures
	phreeqc_rm.BMI_GetValue("Temperature", temperature)
	# Get initial temperature
	phreeqc_rm.BMI_GetValue("Saturation", sat)
	# Get initial porosity
	phreeqc_rm.BMI_GetValue("Porosity", por)
	# Get initial temperature
	phreeqc_rm.BMI_GetValue("SolutionVolume", volume)
	# Get initial concentrations
	c = phreeqcrm.DoubleVector()
	phreeqc_rm.BMI_GetValue("Concentrations", c)
	# Set density, pressure, and temperature (previously allocated)
	density = [1.0] * nxyz
	phreeqc_rm.BMI_SetValue("Density", density)
	pressure = [2.0] * nxyz
	phreeqc_rm.BMI_SetValue("Pressure", pressure)  
	temperature = [20.0] * nxyz
	phreeqc_rm.BMI_SetValue("Temperature", temperature)  
    # --------------------------------------------------------------------------
    # Set boundary condition
    # --------------------------------------------------------------------------
	nbound = 1
	bc1 = [0]           # solution 0 from Initial IPhreeqc instance
	bc2 = [-1]          # no bc2 solution for mixing
	bc_f1 = [1.0]       # mixing fraction for bc1
	bc_conc = phreeqcrm.DoubleVector()
	phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, nbound, bc1, bc2, bc_f1)
    # --------------------------------------------------------------------------
    # Transient loop
    # --------------------------------------------------------------------------
	nsteps = 10
	time = 0.0
	phreeqc_rm.BMI_SetValue("Time", time)
	time_step = 86400.0
	phreeqc_rm.BMI_SetValue("TimeStep", time_step)   
	for i in range(nsteps):
		print("Beginning transport calculation ", time/86400., " days")
		phreeqc_rm.SetScreenOn(1)
		print("          Time step             ", time_step / 86400., " days")
		# Transport calculation here, changes c
		advectionbmi_py(c, bc_conc, ncomps, nxyz)

		# Transfer data to PhreeqcRM for reactions
		print_selected_output_on = (steps == nsteps - 1)
		print_chemistry_on = (steps == nsteps - 1)
		phreeqc_rm.SetSelectedOutputOn(print_selected_output_on)
		phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, False, False)  # workers, initial_phreeqc, utility
		phreeqc_rm.SetConcentrations(c)                                   # Transported concentrations
		time += time_step
		status = phreeqc_rm.SetTime(time)
		# Transfer data to PhreeqcRM after transport      
		phreeqc_rm.BMI_SetValue("Concentrations", c)  # Transported concentrations
		# Optionally, if values changed during transport
		phreeqc_rm.BMI_SetValue("Porosity", por)              
		phreeqc_rm.BMI_SetValue("Saturation", sat)            
		phreeqc_rm.BMI_SetValue("Temperature", temperature) 
		phreeqc_rm.BMI_SetValue("Pressure", pressure)          
		phreeqc_rm.BMI_SetValue("TimeStep", time_step) 
		# Set new time
		time = time + time_step              
		phreeqc_rm.BMI_SetValue("Time", time)  # Current time
		# Run cells with transported conditions
		print("Beginning reaction calculation  ", time / 86400., " days")
		# Demonstration of state 
		phreeqc_rm.StateSave(1)
		phreeqc_rm.StateApply(1)
		phreeqc_rm.StateDelete(1)
		# Run chemistry
		phreeqc_rm.BMI_Update()
		# Get new data calculated by PhreeqcRM for transport
		phreeqc_rm.BMI_GetValue("Concentrations", c)
		density = phreeqc_rm.DoubleVector()
		phreeqc_rm.BMI_GetValue("Density", density)
		volume = phreeqc_rm.DoubleVector()
		phreeqc_rm.BMI_GetValue("SolutionVolume", volume)   
		# Print results at last time step
		if (isteps == nsteps):
			print("Current distribution of cells for workers")
			print("Worker      First cell        Last Cell")
			n = phreeqc_rm.GetThreadCount() * phreeqc_rm.GetMpiTasks()
			sc = phreeqc_rm.IntVector()
			ec = phreeqc_rm.IntVector()
			phreeqc_rm.GetStartCell(sc)
			phreeqc_rm.GetEndCell(ec)
			for i in range(n):
				print(i,"           ", sc(i),"                 ",ec(i))
			
			# Loop through possible multiple selected output definitions
			phreeqc_rm.BMI_GetValue("SelectedOutputCount", n)
			
			for isel in range(n): 
				i = isel
				phreeqc_rm.BMI_SetValue("NthSelectedOutput", i)
				phreeqc_rm.BMI_GetValue("CurrentSelectedOutputUserNumber", n_user)
				print("Selected output sequence number: ", isel)
				print("Selected output user number:     ", n_user)
				# Get 2D array of selected output values
				phreeqc_rm.BMI_GetValue("SelectedOutputColumnCount", col)
				phreeqc_rm.BMI_GetValue("SelectedOutputRowCount", rows)
				selected_out = phreeqc_rm.DoubleVector()
				# Get headings
				phreeqc_rm.BMI_GetValue("SelectedOutputHeadings", headings)
				# Get selected output
				phreeqc_rm.BMI_GetValue("SelectedOutput", selected_out)
				# Print results
				for i in range(rows//2):
					print("Cell number ", i)
					print("     Density: ", density[i])
					print("     Volume:  ", volume[i])
					print("     Components: ")
					for j in range(ncomps):
						print(j, " ",components[j], ": ", c[j * nxyz + i])
					print("     Selected output: ")
					for j in range(col):
						print(j, " ", headings[j],": ", selected_out[j * nxyz + i])
						

	# Clean up
	phreeqc_rm.CloseFiles()
	phreeqc_rm.MpiWorkerBreak()
	phreeqc_rm.BMI_Finalize()


def advectionbmi_py(c, bc_conc, ncomps, nxyz):
    # Advect
    for i in range(nxyz - 1, 0, -1):
        for j in range(ncomps):
            c[j * nxyz + i] = c[j * nxyz + i - 1]              # component j
    
    # Cell zero gets boundary condition
    for j in range(ncomps):
        c[j * nxyz] = bc_conc[j * dim];                        # component j


if __name__ == '__main__':
    AdvectBMI_f90()


#endif # USE_YAML
