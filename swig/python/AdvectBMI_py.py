import phreeqcrm
import numpy as np

from constants import FilePaths

	#module mydata
	#  double precision, dimension(:), pointer :: K_ptr
	#  integer                                 :: rm_id
	#end module mydata

class AdvectBMI(phreeqcrm.BMIPhreeqcRM):

	def __init__(self, yaml=""):
		phreeqcrm.BMIPhreeqcRM.__init__(self)
		self.initialize(yaml)

	def display_results(self):

		nxyz = self.get_value_ptr("GridCellCount")[0]
		components = self.get_value_ptr("Components")
		density = self.get_value_ptr("DensityCalculated")
		volume = self.get_value_ptr("SolutionVolume")
		c = self.get_value_ptr("Concentrations")

		ncomps = len(components)

		n = self.get_scalar("SelectedOutputCount")
		for isel in range(n):
			self.set_scalar("NthSelectedOutput", isel)
			n_user = self.get_scalar("CurrentSelectedOutputUserNumber")

			print(f"Selected output sequence number: {isel}")
			print(f"Selected output user number:     {n_user}")

			# Get 2D array of selected output values
			cols = self.get_scalar("SelectedOutputColumnCount")
			rows = self.get_scalar("SelectedOutputRowCount")

			# Get headings
			headings = self.get_selected_output_headings()

			# Get selected output
			selected_out = self.get_selected_output()
			
			# Print results
			for i in range(rows//2):
				print("Cell number ", i)
				print(f"     Density:    {density[i]}")
				print(f"     Volume:     {volume[i]}")
				print(f"     Components: ")
				for j in range(ncomps):
					print(f"{j}, {components[j]}: {c[j * nxyz + i]}")
				print(f"     Selected output: ")
				for j in range(cols):
					print(f"{j}, {headings[j]}, {selected_out[j * nxyz + i]}")

	def get_scalar(self, var_name):
		itemsize = self.get_var_itemsize(var_name)
		nbytes = self.get_var_nbytes(var_name)
		dim = nbytes // itemsize
	
		if dim != 1:
			raise ValueError(f"{var_name} is not a scalar")
	
		vtype = self.get_var_type(var_name)
		dest = np.empty(1, dtype=vtype)
		x = self.get_value(var_name, dest)
		return x[0]
	
	def set_scalar(self, var_name, value):
		itemsize = self.get_var_itemsize(var_name)
		nbytes = self.get_var_nbytes(var_name)
		dim = nbytes // itemsize
	
		if dim != 1:
			raise ValueError(f"{var_name} is not a scalar")
		
		vtype = self.get_var_type(var_name)
		dest = np.empty(1, dtype=vtype)
		dest[0] = value
		x = self.set_value(var_name, dest)

	def get_selected_output_headings(self):
		# @todo use bmi methods
		return self.GetSelectedOutputHeadings()
	
	def get_selected_output(self):
		# @todo use bmi methods
		return self.GetSelectedOutput()



def AdvectBMI_py():

	# Based on PHREEQC Example 11

	# --------------------------------------------------------------------------
	# Create PhreeqcRM
	# --------------------------------------------------------------------------
	yaml_file = FilePaths.YAML
	
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
	
	##bmi = phreeqcrm.BMIPhreeqcRM()
	bmi = AdvectBMI(yaml_file)

	# Initialize with YAML file
	##status = bmi.initialize(yaml_file)

    # Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
    #TODO CALL register_basic_callback_fortran()
#ifdef USE_MPI
    # Optional callback for MPI
    #TODO status = do_something()   # only root is calling do_something here
#endif

	# print(bmi.get_pointable_var_names())

	components = bmi.get_value_ptr("Components")
	ncomps = bmi.get_value_ptr("ComponentCount")[0]
	nxyz = bmi.get_value_ptr("GridCellCount")[0]

	# print(components)
	# print(ncomps)
	# print(nxyz)

	time = bmi.get_value_ptr("Time")
	time_step = bmi.get_value_ptr("TimeStep")


	hydraulic_K = [0.0] * nxyz
	for i in range(nxyz):
		hydraulic_K[i] = i*2.0

	nthreads = bmi.GetThreadCount()

	print(f"Number of threads:                                {nthreads}")

	nbytes = bmi.get_var_nbytes("FilePrefix")
	dest = np.full(1, " " * nbytes)
	prefix = bmi.get_value('FilePrefix', dest)[0]

	nchem = bmi.GetChemistryCellCount()
	print(nchem)

	print(f"Number of components for transport:               {ncomps}")

	# Get component information)
	gfw = bmi.get_value_ptr("Gfw")
	for i in range(ncomps):
		print(f"{components[i].rjust(10,' ')}  {gfw[i]}")
	print()

	# Get initial temperatures
	temperature = bmi.get_value_ptr("Temperature")
	# Get initial saturation
	sat = bmi.get_value_ptr("SaturationCalculated")
	# Get initial porosity
	por = bmi.get_value_ptr("Porosity")

	# print(temperature)
	# print(sat)
	# print(por)

	# Get initial volume
	volume = bmi.get_value_ptr("SolutionVolume")
	# Get initial concentrations
	c = bmi.get_value_ptr("Concentrations")

	# print("volume")
	# print(volume)
	# print("c")
	# print(c)

	# Set density, pressure, and temperature
	density = [1.0] * nxyz
	bmi.set_value("DensityUser", density)
	pressure = [2.0] * nxyz
	bmi.set_value("Pressure", pressure)
	temperature = [20.0] * nxyz
	bmi.set_value("Temperature", temperature)

    # --------------------------------------------------------------------------
    # Set boundary condition
    # --------------------------------------------------------------------------
	nbound = 1
	bc1 = [0]           # solution 0 from Initial IPhreeqc instance
	bc2 = [-1]          # no bc2 solution for mixing
	bc_f1 = [1.0]       # mixing fraction for bc1

	bc_conc = bmi.InitialPhreeqc2Concentrations(bc1)

	bmi.SetScreenOn(True)
	time[0] = 0.0
	time_step[0] = 86400.0

	# --------------------------------------------------------------------------
	# Transient loop
	# --------------------------------------------------------------------------
	dummy_int_scalar = np.empty((1,), dtype=int)
	selected_output = bmi.get_value_ptr("SelectedOutputOn")
	nsteps = 10
	for step in range(nsteps):
		print(f"Beginning transport calculation {time*bmi.GetTimeConversion()} days")
		print(f"          Time step             {time_step*bmi.GetTimeConversion()} days")

		advectionbmi_py(c, bc_conc, ncomps, nxyz, nbound)

		if (step == (nsteps - 1)):
			print_chemistry_on = 1
		else:
			print_chemistry_on = 0
		selected_output[0] = print_chemistry_on

		bmi.SetPrintChemistryOn(print_chemistry_on==1, False, False)  # workers, initial_phreeqc, utility
		time[0] += time_step
		# bmi.set_value("Time", time)

		# Optionally, if values changed during transport
		##bmi.set_value("Porosity", por)
		bmi.set_value("SaturationUser", sat)
		##bmi.set_value("Temperature", temperature)
		##bmi.set_value("Pressure", pressure)
		##bmi.set_value("TimeStep", time_step)
		
		# Run cells with transported conditions
		print(f"Beginning reaction calculation  {time*bmi.GetTimeConversion()} days")
		bmi.update()

		# Print results at last time step
		if (step == (nsteps - 1)):
			bmi.display_results()

	# Clean up
	bmi.finalize()
	print("Done.")
	
def advectionbmi_py(c: np.ndarray, bc_conc: np.ndarray, ncomps: int, nxyz: int, dim: int):
	# Advect
	for i in range(nxyz - 1, 0, -1):
		for j in range(ncomps):
			c[j * nxyz + i] = c[j * nxyz + i - 1]              # component j
	
	# Cell zero gets boundary condition
	for j in range(ncomps):
		c[j * nxyz] = bc_conc[j * dim];                        # component j

if __name__ == '__main__':
	AdvectBMI_py()
