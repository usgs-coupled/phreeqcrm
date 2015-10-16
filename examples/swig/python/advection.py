import sys
import phreeqcrm

class advection(phreeqcrm.BasicCallback):
	def __init__(self):
		# must call base class ctor
		phreeqcrm.BasicCallback.__init__(self)
		
		# Create module
		self.nxyz        = 40
		self.threads     = 3
		self.hydraulic_K = phreeqcrm.DoubleVector()
		self.rm          = phreeqcrm.PhreeqcRM(self.nxyz, self.threads)

	def Callback(self, x1, x2, s1):
		rm_cell_number = int(x1);
		if (rm_cell_number >= 0) and (rm_cell_number < self.rm.GetChemistryCellCount()):
			back = self.rm.GetBackwardMapping();
			if (s1 == "HYDRAULIC_K"):
				return self.hydraulic_K.get(back.get(rm_cell_number).get(0))		
		return -999.9
		

	def Exec(self):
		# Based on PHREEQC Example 11

		# Set properties
		status = self.rm.SetErrorHandlerMode(0)
		status = self.rm.SetComponentH2O(False)
		status = self.rm.SetRebalanceFraction(0.5)
		status = self.rm.SetRebalanceByCell(True)
		self.rm.UseSolutionDensityVolume(False)
		self.rm.SetPartitionUZSolids(False)

		# Open files
		status = self.rm.SetFilePrefix("Advect_cpp")
		self.rm.OpenFiles()

		# Set concentration units
		status = self.rm.SetUnitsSolution(2)           # 1, mg/L; 2, mol/L; 3, kg/kgs
		status = self.rm.SetUnitsPPassemblage(1)       # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = self.rm.SetUnitsExchange(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = self.rm.SetUnitsSurface(1)            # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = self.rm.SetUnitsGasPhase(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = self.rm.SetUnitsSSassemblage(1)       # 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = self.rm.SetUnitsKinetics(1)           # 0, mol/L cell; 1, mol/L water; 2 mol/L rock

		# Set conversion from seconds to user units (days)
		time_conversion = 1.0 / 86400
		status = self.rm.SetTimeConversion(time_conversion)

		# Set representative volume
		rv = phreeqcrm.DoubleVector(self.nxyz)
		for i in range(self.nxyz):
			rv[i] = 1.0
		status = self.rm.SetRepresentativeVolume(rv)

		# Set initial porosity
		por = phreeqcrm.DoubleVector(self.nxyz)
		for i in range(self.nxyz):
			por[i] = 0.2
		status = self.rm.SetPorosity(por)

		# Set initial saturation
		sat = phreeqcrm.DoubleVector(self.nxyz)
		for i in range(self.nxyz):
			sat[i] = 1.0
		status = self.rm.SetSaturation(sat)	
	
		# Set cells to print chemistry when print chemistry is turned on
		print_chemistry_mask = phreeqcrm.IntVector(self.nxyz)
		for i in range(self.nxyz):
			print_chemistry_mask[i] = 0			
		for i in range(self.nxyz/2):
			print_chemistry_mask[i] = 1
		status = self.rm.SetPrintChemistryMask(print_chemistry_mask);	
	
		# test getters
		print_chemistry_mask1 = self.rm.GetPrintChemistryMask();
		print_on              = self.rm.GetPrintChemistryOn();
		rebalance             = self.rm.GetRebalanceByCell();
		f_rebalance           = self.rm.GetRebalanceFraction();
		so_on                 = self.rm.GetSelectedOutputOn();
		units_exchange        = self.rm.GetUnitsExchange();
		units_gas_phase       = self.rm.GetUnitsGasPhase();
		units_kinetics        = self.rm.GetUnitsKinetics();
		units_pp_assemblage   = self.rm.GetUnitsPPassemblage();
		units_solution        = self.rm.GetUnitsSolution();
		units_ss_exchange     = self.rm.GetUnitsSSassemblage();
		units_surface         = self.rm.GetUnitsSurface();	
	
		# Demonstation of mapping, two equivalent rows by symmetry
		grid2chem = phreeqcrm.IntVector(self.nxyz)
		for i in range(self.nxyz):
			grid2chem[i] = -1

		for i in range(self.nxyz/2):
			grid2chem[i] = i
			grid2chem[i + self.nxyz/2] = i

		status = self.rm.CreateMapping(grid2chem);
		if (status != phreeqcrm.IRM_OK):
			self.rm.DecodeError(status)
		nchem = self.rm.GetChemistryCellCount()

		# --------------------------------------------------------------------------
		# Set initial conditions
		# --------------------------------------------------------------------------	
	
		# Set printing of chemistry file
		status = self.rm.SetPrintChemistryOn(False, True, False) # workers, initial_phreeqc, utility
	
		# Load database
		self.rm.LoadDatabase('phreeqc.dat')



if __name__ == '__main__':
	a = advection()
	a.Exec()