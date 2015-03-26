#include "PhreeqcRM.h"
#include "RM_interface_F.h"
#include "IPhreeqcPhastLib.h"
#include "Phreeqc.h"
#include "PHRQ_io.h"
#include <string>
#include <map>
#ifdef USE_MPI
#include "mpi.h"
#endif

static void
rmpadfstring(char *dest, const char *src, unsigned int len)
{
    size_t sofar;

    for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
        *dest++ = *src++;

    while (sofar++ < len)
        *dest++ = ' ';
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_Abort(int *id, int *result, const char * str)
/* ---------------------------------------------------------------------- */
{
	// decodes error
	// writes any error messages
	// exits 
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->DecodeError(*result);
		Reaction_module_ptr->ErrorMessage(str);
		Reaction_module_ptr->MpiAbort();
		Reaction_module_ptr->DestroyReactionModule(*id);
	    exit(4);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_CloseFiles(int *id)
/* ---------------------------------------------------------------------- */
{
	// closes output and log file
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->CloseFiles();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_Concentrations2Utility(int *id, double *c, int *n, double *tc, double *p_atm)
/* ---------------------------------------------------------------------- */
{
	// set of concentrations c is imported to SOLUTIONs in the Utility IPhreeqc
	// n is the number of sets of concentrations
	// c is of size n times count_components, equivalent to Fortran definition (n, ncomps)
	// tc is array of dimension n of temperatures 
	// p_atm is array of dimension n pressure
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> c_vector, tc_vector, p_atm_vector;
		size_t ncomps = Reaction_module_ptr->GetComponents().size();
		c_vector.resize(*n * ncomps, 0.0);

		for (size_t i = 0; i < (size_t) *n; i++)
		{
			for (size_t j = 0; j < ncomps; j++)
			{
				c_vector[j * (*n) + i] = c[j * (*n) + i];
			}
			tc_vector.push_back(tc[i]);
			p_atm_vector.push_back(p_atm[i]);
		}
		IPhreeqc * util_ptr = Reaction_module_ptr->Concentrations2Utility(c_vector, tc_vector, p_atm_vector);
		if (util_ptr != NULL)
		{
			return util_ptr->GetId();
		}
		return IRM_FAIL;
	}
	return IRM_BADINSTANCE;
}
#ifdef USE_MPI
/* ---------------------------------------------------------------------- */
int
RMF_Create(int *nxyz, int *nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
  return PhreeqcRM::CreateReactionModule(*nxyz, MPI_Comm_f2c(*nthreads));
}
#else
/* ---------------------------------------------------------------------- */
int
RMF_Create(int *nxyz, int *nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return PhreeqcRM::CreateReactionModule(*nxyz, *nthreads);
}
#endif
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_CreateMapping(int *id, int *grid2chem)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates mapping from all grid cells to only cells for chemistry
	// Excludes inactive cells (negative values) and cells that are redundant by symmetry (many-to-one mapping)
	// (1D or 2D chemistry)
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<int> grid2chem_vector;
		grid2chem_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&grid2chem_vector.front(), grid2chem, (size_t) (Reaction_module_ptr->GetGridCellCount() * sizeof(int)));
		return Reaction_module_ptr->CreateMapping(grid2chem_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_DecodeError(int *id, int *e)
/* ---------------------------------------------------------------------- */
{
	// Prints the error message for IRM_RESULT e
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->DecodeError(*e);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_Destroy(int *id)
/* ---------------------------------------------------------------------- */
{
	// Destroys reaction module
	return PhreeqcRM::DestroyReactionModule(*id);
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_DumpModule(int *id, int *dump_on, int *append)
/* ---------------------------------------------------------------------- */
{	
	// Dumps raw format of all chemistry cells to file defined by
	// SetDumpFileName
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->DumpModule(*dump_on != 0, *append != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_ErrorMessage(int *id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	// writes an error message
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string e_string(err_str);
		trim_right(e_string);
		e_string.append("\n");
		Reaction_module_ptr->ErrorMessage(e_string);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_FindComponents(int *id)
/* ---------------------------------------------------------------------- */
{
	// Accumulates a list of components from the definitions in InitialPhreeqc
	// Multiple calls will add new components to the list
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return (Reaction_module_ptr->FindComponents());
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetChemistryCellCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of chemistry cells <= number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetChemistryCellCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetComponent(int * id, int * num, char *chem_name, int * l1)
/* ---------------------------------------------------------------------- */
{
	// Retrieves the component name in position num to chem_name
	// num <= result of FindComponent
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (chem_name != NULL)
		{
			if (l1 > 0 && *num > 0 && *num <= Reaction_module_ptr->GetComponentCount())
			{
				rmpadfstring(chem_name, Reaction_module_ptr->GetComponents()[*num - 1].c_str(), (unsigned int) *l1);
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetComponentCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of components 
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetComponentCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetConcentrations(int *id, double * c)
/* ---------------------------------------------------------------------- */
{
	// Retrieves concentrations for all grid nodes to c
	// size of c must be the number of grid nodes time the number of components
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> c_vector;
		c_vector.resize(Reaction_module_ptr->GetGridCellCount() * Reaction_module_ptr->GetComponentCount());
		IRM_RESULT return_value = Reaction_module_ptr->GetConcentrations(c_vector);
		if (return_value == IRM_OK)
		{
			memcpy(c, &c_vector.front(), c_vector.size() * sizeof(double));
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetDensity(int *id, double * d)
/* ---------------------------------------------------------------------- */
{
	// Retrieves density for all grid nodes in d
	// size of d must be the number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		IRM_RESULT return_value = IRM_OK;
		std::vector <double> density_vector;
		Reaction_module_ptr->GetDensity(density_vector);
		if ((int) density_vector.size() == Reaction_module_ptr->GetGridCellCount())
		{
			memcpy(d, &density_vector.front(), (size_t) (Reaction_module_ptr->GetGridCellCount()*sizeof(double)));
		}
		else
		{
			for (int i = 0; i < Reaction_module_ptr->GetGridCellCount(); i++)
			{
				d[i] = INACTIVE_CELL_VALUE;
			}
			return_value = IRM_FAIL;
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetErrorString(int * id, char *errstr, int * l)
/* ---------------------------------------------------------------------- */
{
	// Retrieves file prefix in prefix
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		rmpadfstring(errstr, Reaction_module_ptr->GetErrorString().c_str(), (unsigned int) *l);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_GetErrorStringLength(int * id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves file prefix in prefix
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return (int) (Reaction_module_ptr->GetErrorString().size());
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetFilePrefix(int * id, char *prefix, int *l)
/* ---------------------------------------------------------------------- */
{
	// Retrieves file prefix in prefix
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		rmpadfstring(prefix, Reaction_module_ptr->GetFilePrefix().c_str(), (unsigned int) *l);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetGfw(int *id, double * gfw)
/* ---------------------------------------------------------------------- */
{
	// Retrieves gram formula weights
	// size of d must be the number of components
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		size_t ncomps = Reaction_module_ptr->GetComponents().size();
		if (ncomps > 0)
		{
			memcpy(gfw, &Reaction_module_ptr->GetGfw().front(), ncomps * sizeof(double));
			return IRM_OK;
		}
		return IRM_FAIL;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_GetGridCellCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetGridCellCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_GetIPhreeqcId(int * id, int * i)
	/* ---------------------------------------------------------------------- */
{
	// Returns an integer for an IPhreeqc instance
	// If using threads, there are nthreads + 2 IPhreeqc instances
	//      0 through nthreads - 1 are the worker IPhreeqcs
	//      nthreads is InitialPhreeqc IPhreeqc
	//      nthreads + 1 is Utility IPhreeqc
	// If using MPI,
	//      0 is the worker IPhreeqc
	//      1 is InitialIPhreeqc
	//      2 is Utility IPhreeqc
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		IPhreeqc * iptr = Reaction_module_ptr->GetIPhreeqcPointer(*i);
		if (iptr != NULL)
		{
			return iptr->GetId();
		}
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetMpiMyself(int * id)
	/* ---------------------------------------------------------------------- */
{
	// Returns 0 for threaded version and
	// process number for MPI version, 0 for root, > 0 for workers
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetMpiMyself();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetMpiTasks(int * id)
	/* ---------------------------------------------------------------------- */
{
	// Returns 1 for threaded and
	// the number of processes for MPI
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetMpiTasks();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetNthSelectedOutputUserNumber(int * id, int * i)
	/* ---------------------------------------------------------------------- */
{
	// This method returns the user number for the nth selected output
	// RMF_GetSelectedOutputCount is number of selected outputs
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetNthSelectedOutputUserNumber(*i - 1);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSaturation(int *id, double * sat)
/* ---------------------------------------------------------------------- */
{
	// Retrieves saturation for all grid nodes in sat
	// size of sat must be the number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		IRM_RESULT return_value = IRM_OK;
		std::vector <double> sat_vector;
		Reaction_module_ptr->GetSaturation(sat_vector);
		if ((int) sat_vector.size() == Reaction_module_ptr->GetGridCellCount())
		{
			memcpy(sat, &sat_vector.front(), (size_t) (Reaction_module_ptr->GetGridCellCount()*sizeof(double)));
		}
		else
		{
			for (int i = 0; i < Reaction_module_ptr->GetGridCellCount(); i++)
			{
				sat[i] = INACTIVE_CELL_VALUE;
			}
			return_value = IRM_FAIL;
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSelectedOutput(int * id, double * so)
/* ---------------------------------------------------------------------- */
{
	// Retrieves selected output for the currently selected selected output
	// as an array of doubles
	// The size of the array so must be nxyz * RMF_GetSelectedOutputColumnCount
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> so_vector;
		so_vector.resize(Reaction_module_ptr->GetSelectedOutputColumnCount() * 
			Reaction_module_ptr->GetSelectedOutputRowCount());
		IRM_RESULT return_value = Reaction_module_ptr->GetSelectedOutput(so_vector);
		if (return_value == IRM_OK)
		{
			memcpy(so, &so_vector.front(), so_vector.size() * sizeof(double));
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetSelectedOutputColumnCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of columns for the currently selected selected output
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputColumnCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetSelectedOutputCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of selected output definitions (different user numbers) that have been defined
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSelectedOutputHeading(int * id, int *icol, char *heading, int *length)
/* ---------------------------------------------------------------------- */
{
	// Returns the heading at position icol (numbered from 0)
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string head;
		IRM_RESULT return_value = Reaction_module_ptr->GetSelectedOutputHeading(*icol - 1, head);
		if (return_value == IRM_OK)
		{
			rmpadfstring(heading, head.c_str(), (unsigned int) *length);
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetSelectedOutputRowCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of rows for the currently selected selected output = number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputRowCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSolutionVolume(int *id, double * v)
/* ---------------------------------------------------------------------- */
{
	// Retrieves solution volumes for each grid cell
	// size of v is number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		IRM_RESULT return_value = IRM_OK;
		const std::vector<double> &v_vector = Reaction_module_ptr->GetSolutionVolume();
		if ((int) v_vector.size() == Reaction_module_ptr->GetGridCellCount())
		{
			memcpy(v, &v_vector.front(), v_vector.size() * sizeof(double));
		}
		else
		{
			for (int i = 0; i < Reaction_module_ptr->GetGridCellCount(); i++)
			{
				v[i] = INACTIVE_CELL_VALUE;
			}
			return_value = IRM_FAIL;
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSpeciesConcentrations(int *id, double * species_conc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		IRM_RESULT return_value = IRM_OK;
		std::vector<double> species_conc_vector;
		return_value = Reaction_module_ptr->GetSpeciesConcentrations(species_conc_vector);
		if (return_value == IRM_OK)
		{
			memcpy(species_conc, &species_conc_vector.front(), species_conc_vector.size()*sizeof(double));
		}
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_GetSpeciesCount(int *id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSpeciesCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSpeciesD25(int *id, double * diffc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		const std::vector<double> & diffc_vector = Reaction_module_ptr->GetSpeciesD25();
		memcpy(diffc, &diffc_vector.front(), diffc_vector.size()*sizeof(double));
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSpeciesName(int *id, int *i_in, char *name, int *length)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		int i = *i_in - 1;
		const std::vector<std::string> & names = Reaction_module_ptr->GetSpeciesNames();
		if (i >= 0 && i < (int) names.size())
		{
			rmpadfstring(name, names[i].c_str(), (unsigned int) *length);  
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RMF_GetSpeciesSaveOn(int *id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return (Reaction_module_ptr->GetSpeciesSaveOn() ? 1 : 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_GetSpeciesZ(int *id, double * z)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		const std::vector<double> & z_vector = Reaction_module_ptr->GetSpeciesZ();
		memcpy(z, &z_vector.front(), z_vector.size()*sizeof(double));
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RMF_GetThreadCount(int * id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of threads for threaded version
	// 1 for MPI

	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetThreadCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double
RMF_GetTime(int * id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current simulation time, in seconds
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTime();
	}
	return (double) IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
double
RMF_GetTimeConversion(int * id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves conversion factor such that seconds times conversion factor is user time
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTimeConversion();
	}
	return (double) IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
double RMF_GetTimeStep(int * id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves current time step, in seconds
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTimeStep();
	}
	return (double) IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2Concentrations(
			int *id,
			double *boundary_c,
			int *n_boundary,
			int *boundary_solution1)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine takes a list of solution numbers and returns a set of
 *   concentrations
 *   Input: n_boundary - number of boundary conditions in list
 *          boundary_solution1 - list of first solution numbers to be mixed
 *          Optionally, boundary_solution2 - list of second solution numbers to be mixed
 *          Optionally, fraction1 - list of mixing fractions of solution 1
 *
 *          fraction1 - fraction of first solution 0 <= f <= 1
 *          boundary_solution2 and fraction1 may be omitted if no mixing
 *
 *   Output: boundary_c - concentrations for boundary conditions
 *                      - dimension must be of size n_boundary x n_comp
 *
 */
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
		std::vector < double > destination_c, fraction1_vector;
		boundary_solution1_vector.resize(*n_boundary);
		memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (*n_boundary * sizeof(int)));
		IRM_RESULT return_value = Reaction_module_ptr->InitialPhreeqc2Concentrations(
			destination_c,
			boundary_solution1_vector,
			boundary_solution2_vector,
			fraction1_vector);
		if (return_value == 0)
		{
			memcpy(boundary_c, &destination_c.front(), destination_c.size() * sizeof(double));
		}       
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2Concentrations2(
			int *id,
			double *boundary_c,
			int *n_boundary,
			int *boundary_solution1,  
			int *boundary_solution2, 
			double *fraction1)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine takes a list of solution numbers and returns a set of
 *   concentrations
 *   Input: n_boundary - number of boundary conditions in list
 *          boundary_solution1 - list of first solution numbers to be mixed
 *          Optionally, boundary_solution2 - list of second solution numbers to be mixed
 *          Optionally, fraction1 - list of mixing fractions of solution 1
 *
 *          fraction1 - fraction of first solution 0 <= f <= 1
 *          boundary_solution2 and fraction1 may be omitted if no mixing
 *
 *   Output: boundary_c - concentrations for boundary conditions
 *                      - dimension must be of size n_boundary x n_comp
 *
 */
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
		std::vector < double > destination_c, fraction1_vector;
		boundary_solution1_vector.resize(*n_boundary);
		memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (*n_boundary * sizeof(int)));
		if (boundary_solution2 != NULL)
		{
			boundary_solution2_vector.resize(*n_boundary);
			memcpy(&boundary_solution2_vector.front(), boundary_solution2, (size_t) (*n_boundary * sizeof(int)));
		}
		if (fraction1 != NULL)
		{
			fraction1_vector.resize(*n_boundary);
			memcpy(&fraction1_vector.front(), fraction1, (size_t) (*n_boundary * sizeof(double)));
		}
		IRM_RESULT return_value = Reaction_module_ptr->InitialPhreeqc2Concentrations(
			destination_c,
			boundary_solution1_vector,
			boundary_solution2_vector,
			fraction1_vector);
		if (return_value == 0)
		{
			memcpy(boundary_c, &destination_c.front(), destination_c.size() * sizeof(double));
		}       
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2Module(int *id,
						  int *initial_conditions1)		// 7 x nxyz end-member 1	
/* ---------------------------------------------------------------------- */
{
	// 7 sets of indices for initial conditions
	// 0 solution
	// 1 ppassemblage
	// 2 exchange
	// 3 surface
	// 4 gas phase
	// 5 ss_assemblage
	// 6 kinetics
	// array is equivalent to Fortran ic(nxyz, 7), where nxyz is number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<int> i1_vector, i2_vector;
		std::vector<double> f1_vector;
		int nxyz = Reaction_module_ptr->GetGridCellCount();
		i1_vector.resize(nxyz * 7);
		i2_vector.resize(nxyz * 7, -1);
		f1_vector.resize(nxyz * 7, 1.0);
		memcpy(&i1_vector.front(), initial_conditions1, (size_t) (nxyz * 7 * sizeof(int)));
		return Reaction_module_ptr->InitialPhreeqc2Module(
			i1_vector,
			i2_vector,
			f1_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2Module2(int *id,
							  int *initial_conditions1,		// 7 x nxyz end-member 1
							  int *initial_conditions2,		// 7 x nxyz end-member 2
							  double *fraction1)			// 7 x nxyz fraction of end-member 1
/* ---------------------------------------------------------------------- */
{
	// 7 sets of indices for initial conditions
	// 0 solution
	// 1 ppassemblage
	// 2 exchange
	// 3 surface
	// 4 gas phase
	// 5 ss_assemblage
	// 6 kinetics
	// array is equivalent to Fortran ic(nxyz, 7), where nxyz is number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<int> i1_vector, i2_vector;
		std::vector<double> f1_vector;
		int nxyz = Reaction_module_ptr->GetGridCellCount();
		i1_vector.resize(nxyz * 7);
		i2_vector.resize(nxyz * 7, -1);
		f1_vector.resize(nxyz * 7, 1.0);
		memcpy(&i1_vector.front(), initial_conditions1, (size_t) (nxyz * 7 * sizeof(int)));
		if (initial_conditions2 != NULL)
		{
			memcpy(&i2_vector.front(), initial_conditions2, (size_t) (nxyz * 7 * sizeof(int)));
		}
		if (fraction1 != NULL)
		{
			memcpy(&f1_vector.front(), fraction1, (size_t) (nxyz * 7 * sizeof(double)));
		}
		return Reaction_module_ptr->InitialPhreeqc2Module(
			i1_vector,
			i2_vector,
			f1_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqcCell2Module(int *id,
                int *n,		                            // InitialPhreeqc cell number
                int *module_numbers,		            // Module cell numbers
                int *dim_module_numbers)			    // Number of module cell numbers
/* ---------------------------------------------------------------------- */
{
	// Copies a cell number n from InitialPhreeqc to the module
	// A cell includes any of the following types reactants with the number n
	// solution, exchange, equilibrium_phases, gas_phase, solid_solution, surface, and kinetics
	// If n is negative, n is set to the largest numbered SOLUTION or MIX in InitialPhreeqc
	// MIX definitions are converted to solutions and copied
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector <int> module_numbers_vector;
		module_numbers_vector.resize(*dim_module_numbers);
		memcpy(&module_numbers_vector.front(), module_numbers, (size_t) (*dim_module_numbers) * sizeof(int));
		return Reaction_module_ptr->InitialPhreeqcCell2Module(
			*n,
			module_numbers_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2SpeciesConcentrations(
			int *id,
			double *species_c,
			int *n_boundary,
			int *boundary_solution1)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine takes a list of solution numbers and returns a set of
 *   aqueous species concentrations
 *   Input: n_boundary - number of boundary conditions in list
 *          boundary_solution1 - list of first solution numbers to be mixed
 *          boundary_solution2 - list of second solution numbers to be mixed
 *          fraction1 - list of mixing fractions of solution 1
 *
 *   Output: species_c - aqueous species concentrations for boundary conditions
 *                     - dimensions must be n_boundary x n_species
 *
 */
	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
		std::vector < double > destination_c, fraction1_vector;
		boundary_solution1_vector.resize(*n_boundary);
		memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (*n_boundary * sizeof(int)));
		IRM_RESULT return_value = Reaction_module_ptr->InitialPhreeqc2SpeciesConcentrations(
			destination_c,
			boundary_solution1_vector,
			boundary_solution2_vector,
			fraction1_vector);		
		if (return_value == 0)
		{
			memcpy(species_c, &destination_c.front(), destination_c.size() * sizeof(double));
		}       
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_InitialPhreeqc2SpeciesConcentrations2(
			int *id,
			double *species_c,
			int *n_boundary,
			int *boundary_solution1,  
			int *boundary_solution2, 
			double *fraction1)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine takes a list of solution numbers and returns a set of
 *   aqueous species concentrations
 *   Input: n_boundary - number of boundary conditions in list
 *          boundary_solution1 - list of first solution numbers to be mixed
 *          boundary_solution2 - list of second solution numbers to be mixed
 *          fraction1 - list of mixing fractions of solution 1
 *
 *   Output: species_c - aqueous species concentrations for boundary conditions
 *                     - dimensions must be n_boundary x n_species
 *
 */
	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
		std::vector < double > destination_c, fraction1_vector;
		boundary_solution1_vector.resize(*n_boundary);
		memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (*n_boundary * sizeof(int)));
		if (boundary_solution2 != NULL)
		{
			boundary_solution2_vector.resize(*n_boundary);
			memcpy(&boundary_solution2_vector.front(), boundary_solution2, (size_t) (*n_boundary * sizeof(int)));
		}
		if (fraction1 != NULL)
		{
			fraction1_vector.resize(*n_boundary);
			memcpy(&fraction1_vector.front(), fraction1, (size_t) (*n_boundary * sizeof(double)));
		}
		IRM_RESULT return_value = Reaction_module_ptr->InitialPhreeqc2SpeciesConcentrations(
			destination_c,
			boundary_solution1_vector,
			boundary_solution2_vector,
			fraction1_vector);		
		if (return_value == 0)
		{
			memcpy(species_c, &destination_c.front(), destination_c.size() * sizeof(double));
		}       
		return return_value;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_LoadDatabase(int * id, const char *db_name)
/* ---------------------------------------------------------------------- */
{
	// Loads a database, must be done before any simulations
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string db(db_name);
		size_t strEnd = db.find_last_not_of(" \t\n");
		db = db.substr(0, strEnd + 1);
		return Reaction_module_ptr->LoadDatabase(db.c_str());
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_LogMessage(int *id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	// write a message to the log file
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string e_string(err_str);
		trim_right(e_string);
		e_string.append("\n");
		Reaction_module_ptr->LogMessage(e_string);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_MpiWorker(int *id)
/* ---------------------------------------------------------------------- */
{
	// Starts a listener for a worker MPI process
	// The listener will receive a message from root and perform the
	// necessary calculations and data transfers
	// Return from the method occurs on a RMF_MpiWorkerBreak() message
	// from root
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->MpiWorker();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_MpiWorkerBreak(int *id)
/* ---------------------------------------------------------------------- */
{
	// Sends a message from root to all workers to return from the MpiWorker method
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->MpiWorkerBreak();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_OpenFiles(int *id)
/* ---------------------------------------------------------------------- */
{
	// Opens output and log files on root
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->OpenFiles();
	}
	return IRM_BADINSTANCE;
}
/* 
--------------------------------------------------------------------- */
IRM_RESULT
RMF_OutputMessage(int *id, const char *str)
/* ---------------------------------------------------------------------- */
{
	// writes a message to the output file
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string e_string(str);
		e_string.append("\n");
		Reaction_module_ptr->OutputMessage(e_string);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;

}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_RunCells(int *id)
/* ---------------------------------------------------------------------- */
{
	// Runs reactions for each of the chemistry cells
	// Rebalances the load among threads or MPI processes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->RunCells(); 
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_RunFile(int *id, int *workers, int *initial_phreeqc, int *utility, const char *chem_name)
/* ---------------------------------------------------------------------- */
{
	// Runs a PHREEQC input file
	// for all workers (thread or process) if workers != 0
	// for InitialPhreeqc if initial_phreeqc != 0
	// for Utility if utility != 0
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str(chem_name);
		size_t strEnd = str.find_last_not_of(" \t\n");
		str = str.substr(0, strEnd + 1);
		return Reaction_module_ptr->RunFile((*workers != 0), (*initial_phreeqc != 0), (*utility != 0), str.c_str());
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_RunString(int *id, int *workers, int *initial_phreeqc, int *utility, const char *input_string)
/* ---------------------------------------------------------------------- */
{	
	// Runs a PHREEQC input string
	// for all workers (thread or process) if workers != 0
	// for InitialPhreeqc if initial_phreeqc != 0
	// for Utility if utility != 0
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{	
		std::string str(input_string);
		return Reaction_module_ptr->RunString((*workers != 0), (*initial_phreeqc != 0), (*utility != 0), input_string);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_ScreenMessage(int *id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	// writes a message to the screen
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string e_string(err_str);
		Reaction_module_ptr->ScreenMessage(e_string);
		Reaction_module_ptr->ScreenMessage("\n");
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetComponentH2O(int *id, int *tf)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetComponentH2O(*tf != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetConcentrations(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the concentrations for the solutions in the cells
	// This method would be called to transfer transported concentrations
	// to the reaction module
	// size of t is number of grid nodes times number of components
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> c_vector;
		c_vector.resize(Reaction_module_ptr->GetGridCellCount() * Reaction_module_ptr->GetComponentCount());
		memcpy(&c_vector.front(), t, c_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetConcentrations(c_vector);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetCurrentSelectedOutputUserNumber(int * id, int * i)
/* ---------------------------------------------------------------------- */
{
	// i is the user number of the selected output that is to be 
	// processed
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetCurrentSelectedOutputUserNumber(*i);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetDensity(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the density that may be used to convert concentrations from
	// transport to the number of moles in chemistry cell solutions
	// size of t is the number of grid nodes.
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
			std::vector<double> d_vector;
			d_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&d_vector.front(), t, d_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetDensity(d_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetDumpFileName(int *id, const char *name)
/* ---------------------------------------------------------------------- */
{
	// Sets the name of the dump file
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str(name);		
		size_t strEnd = str.find_last_not_of(" \t\n");
		str = str.substr(0, strEnd + 1);
		return Reaction_module_ptr->SetDumpFileName(str.c_str());
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetErrorHandlerMode(int *id, int *mode)
/* ---------------------------------------------------------------------- */
{
	// Sets the action on encountering an error
	// mode = 0, return an error return code from each method
	// mode = 1, throw an exception,
	// mode = 2, exit
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetErrorHandlerMode(*mode);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetFilePrefix(int *id, const char *name)
/* ---------------------------------------------------------------------- */
{
	// Sets the file prefix for the output and log files that are opened with
	// RMF_OpenFiles on the root process
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string str(name);
		size_t strEnd = str.find_last_not_of(" \t\n");
		str = str.substr(0, strEnd + 1);
		return Reaction_module_ptr->SetFilePrefix(str.c_str());
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetMpiWorkerCallback(int *id, int (*fcn)(int *x1))
/* ---------------------------------------------------------------------- */
{
	// Registers a callback, effectively extending the set of messages
	// that are responded to by RMF_MpiWorker()
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetMpiWorkerCallbackFortran(fcn);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetPartitionUZSolids(int *id, int *t)
/* ---------------------------------------------------------------------- */
{
	// Sets a flag that determines whether the ratio of water to solid
	// changes in unsaturated reactive transport simulations
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		bool tf = false;
		if (t != NULL)
		{
			tf = (*t != 0);
		}
		return Reaction_module_ptr->SetPartitionUZSolids(tf);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetPorosity(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current porosity of cells, which may change due to 
	// compressibility of the water and media.
	// size of t is the number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> v_vector;
		v_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&v_vector.front(), t, v_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetPorosity(v_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetPressure(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current pressure for each cell
	// size of t is the number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> p_vector;
		p_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&p_vector.front(), t, p_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetPressure(p_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetPrintChemistryOn(int *id,	 int *worker, int *ip, int *utility)
/* ---------------------------------------------------------------------- */
{
	// Sets flag for whether to print PHREEQC output  
	// worker != 0 turns on printing for workers 
	// ip != 0 turns on printing for InitialPhreeqc
	// ip != 0 turns on printing for Utility
	// Warning: the output could be huge for the workers if all cells are selected
	// RMF_SetPrintChemistryMask can be used to limit printing from workers to
	// selected cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetPrintChemistryOn(*worker != 0, *ip != 0, *utility != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetPrintChemistryMask(int *id, int *t)
/* ---------------------------------------------------------------------- */
{
	// Sets flags to print or not print chemistry for each cell
	// size of t is number of grid cells
	// if a value for t(i) is 0, no printing occurs for the associated chemistry cell
	// if a value for t(i) is != 0, printing occurs for the associated chemistry cell
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{		
		std::vector<int> m_vector;
		m_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&m_vector.front(), t, m_vector.size() * sizeof(int));
		return Reaction_module_ptr->SetPrintChemistryMask(m_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetRebalanceFraction(int *id, double *f)
/* ---------------------------------------------------------------------- */
{
	// Scalar value 0.0-1.0, for rebalanceing cells among threads or processes
	// A value of zero turns of rebalancing, otherwise a fraction of the
	// calculated optimum redistribution of cells is tranferred among
	// threads or processes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetRebalanceFraction(*f);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetRebalanceByCell(int *id, int *method)
/* ---------------------------------------------------------------------- */
{
	// alternative rebalancing method is used if method != 0
	// It is not clear if one method is any better than the other
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetRebalanceByCell(*method != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetRepresentativeVolume(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the pore volume of a cell, can be an absolute volume
	// or a relative volume
	// size of t is number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> v_vector;
		v_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&v_vector.front(), t, v_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetRepresentativeVolume(v_vector);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetSaturation(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current saturation for each cell
	// size of t is the number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> s_vector;
		s_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&s_vector.front(), t, s_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetSaturation(s_vector);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetSelectedOutputOn(int *id, int *selected_output_on)
/* ---------------------------------------------------------------------- */
{
	// Specifies whether selected output will be retrieved for this time step
	// selected_output_on != 0 indicates selected output will be retrieved
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetSelectedOutputOn(*selected_output_on != 0);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetSpeciesSaveOn(int *id, int *save_on)
/* ---------------------------------------------------------------------- */
{
	// pass pointers from Fortran to the Reaction module
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetSpeciesSaveOn(*save_on != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetTemperature(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current temperature for each cell
	// size of t is the number of grid cells
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::vector<double> t_vector;
		t_vector.resize(Reaction_module_ptr->GetGridCellCount());
		memcpy(&t_vector.front(), t, t_vector.size() * sizeof(double));
		return Reaction_module_ptr->SetTemperature(t_vector);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetTime(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current time in seconds
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTime(*t);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetTimeConversion(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	//
	// multiply seconds to convert to user time units
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTimeConversion(*t);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetTimeStep(int *id, double *t)
/* ---------------------------------------------------------------------- */
{
	// Sets the current time step, seconds
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTimeStep(*t);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsExchange (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsExchange(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsGasPhase (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsGasPhase(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsKinetics (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsKinetics(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsPPassemblage (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsPPassemblage(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsSolution (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units of transport 1 mg/L, 2 mmol/L, 3 kg/kgs
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSolution(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsSSassemblage (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSSassemblage(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SetUnitsSurface (int *id, int *u)
/* ---------------------------------------------------------------------- */
{	
	// units for solid, 1 is per liter water, 2 is per liter rock
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSurface(*u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RMF_SpeciesConcentrations2Module(int *id, double * species_conc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		if (species_conc)
		{
			IRM_RESULT return_value = IRM_OK;
			std::vector<double> species_conc_vector;
			species_conc_vector.resize(Reaction_module_ptr->GetGridCellCount() * Reaction_module_ptr->GetSpeciesCount());
			memcpy(&species_conc_vector.front(), species_conc, species_conc_vector.size()*sizeof(double));
			return_value = Reaction_module_ptr->SpeciesConcentrations2Module(species_conc_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* --------------------------------------------------------------------- */
IRM_RESULT
RMF_UseSolutionDensityVolume(int *id, int *tf)
/* ---------------------------------------------------------------------- */
{
	// writes a warning message to screen, log, and output files
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->UseSolutionDensityVolume(*tf != 0);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;

}
/* --------------------------------------------------------------------- */
IRM_RESULT
RMF_WarningMessage(int *id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	// writes a warning message to screen, log, and output files
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(*id);
	if (Reaction_module_ptr)
	{
		std::string e_string(err_str);
		trim_right(e_string);
		e_string.append("\n");
		Reaction_module_ptr->WarningMessage(e_string);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;

}
