#ifdef USE_MPI
#include "mpi.h"
#endif
#include "PhreeqcRM.h"
#include "RM_interface_C.h"
#include "IPhreeqcPhastLib.h"
#include "Phreeqc.h"
#include "PHRQ_io.h"
#include <string>
#include <map>

/* ---------------------------------------------------------------------- */
IRM_RESULT RM_Abort(int id, int result, const char * str)
/* ---------------------------------------------------------------------- */
{
	// decodes error
	// writes any error messages
	// exits 
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->DecodeError(result);
		Reaction_module_ptr->ErrorMessage(str);
		Reaction_module_ptr->MpiAbort();
		Reaction_module_ptr->DestroyReactionModule(id);
		exit(4);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_CloseFiles(int id)
/* ---------------------------------------------------------------------- */
{	
	// closes output and log file
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->CloseFiles();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RM_Concentrations2Utility(int id, double *c, int n, double *tc, double *p_atm)
/* ---------------------------------------------------------------------- */
{
	// set of concentrations c is imported to SOLUTIONs in the Utility IPhreeqc
	// n is the number of sets of concentrations
	// c is of size n times count_components, equivalent to Fortran definition (n, ncomps)
	// tc is array of dimension n of temperatures 
	// p_atm is array of dimension n pressure
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (c != NULL && tc != NULL && p_atm != NULL)
		{
			std::vector<double> c_vector, tc_vector, p_atm_vector;
			size_t ncomps = Reaction_module_ptr->GetComponents().size();
			c_vector.resize(n * ncomps, 0.0);

			for (size_t i = 0; i < (size_t) n; i++)
			{
				for (size_t j = 0; j < ncomps; j++)
				{
					c_vector[j * n + i] = c[j * n + i];
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
#ifdef USE_MPI
/* ---------------------------------------------------------------------- */
int RM_Create(int nxyz, MPI_Comm comm)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return PhreeqcRM::CreateReactionModule(nxyz, comm);
}
#else
/* ---------------------------------------------------------------------- */
int RM_Create(int nxyz, int nthreads)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates reaction module, called by root and MPI workers
	//
	return PhreeqcRM::CreateReactionModule(nxyz, nthreads);
}
#endif
/* ---------------------------------------------------------------------- */
IRM_RESULT RM_CreateMapping(int id, int *grid2chem)
/* ---------------------------------------------------------------------- */
{
	//
	// Creates mapping from all grid cells to only cells for chemistry
	// Excludes inactive cells (negative values) and cells that are redundant by symmetry (many-to-one mapping)
	// (1D or 2D chemistry)
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (grid2chem != NULL)
		{
			std::vector<int> grid2chem_vector;
			grid2chem_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&grid2chem_vector.front(), grid2chem, (size_t) (Reaction_module_ptr->GetGridCellCount() * sizeof(int)));
			return Reaction_module_ptr->CreateMapping(grid2chem_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT RM_DecodeError(int id, int e)
/* ---------------------------------------------------------------------- */
{
	// Prints the error message for IRM_RESULT e
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->DecodeError(e);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT RM_Destroy(int id)
/* ---------------------------------------------------------------------- */
{
	return PhreeqcRM::DestroyReactionModule(id);
}

/* ---------------------------------------------------------------------- */
IRM_RESULT RM_DumpModule(int id, int dump_on, int append)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->DumpModule((dump_on != 0), (append != 0));
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_ErrorMessage(int id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (err_str)
		{
			std::string e_string(err_str);
			trim_right(e_string);
			Reaction_module_ptr->ErrorMessage(e_string);
			return IRM_OK;
		}
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int
RM_FindComponents(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return (Reaction_module_ptr->FindComponents());
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_DLL_EXPORT IRM_RESULT RM_GetBackwardMapping(int id, int n, int *list, int *size)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (n >= 0 && n < Reaction_module_ptr->GetChemistryCellCount() && list != NULL)
		{
			const std::vector < std::vector<int> > & back = Reaction_module_ptr->GetBackwardMapping();
			if (*size >= (int) back[n].size())
			{
				*size = (int) back[n].size();
				for (int i = 0; i < (int) back[n].size(); i++)
				{
					list[i] = back[n][i];
				}
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int RM_GetChemistryCellCount(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetChemistryCellCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetComponent(int id, int num, char *chem_name, int l1)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (chem_name != NULL)
		{
			//if (l1 >= 0)
			if (l1 > 0 && num >= 0 && num < Reaction_module_ptr->GetComponentCount())
			{
				strncpy(chem_name, Reaction_module_ptr->GetComponents()[num].c_str(), l1);
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int RM_GetComponentCount(int id)
/* ---------------------------------------------------------------------- */
{
	// Returns the number of components 
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetComponentCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_GetConcentrations(int id, double * c)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (c != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetDensity(int id, double * d)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (d != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
IRM_RESULT 
RM_GetEndCell(int id, int *ec)
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		const std::vector <int> & l = Reaction_module_ptr->GetEndCell();
		memcpy(ec, &l.front(), l.size() * sizeof(int));
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_GetErrorString(int id, char *errstr, int l)
/* ---------------------------------------------------------------------- */
{
	// Retrieves file prefix in prefix
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		strncpy(errstr, Reaction_module_ptr->GetErrorString().c_str(), l);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int 
RM_GetErrorStringLength(int id)
/* ---------------------------------------------------------------------- */
{
	// Retrieves file prefix in prefix
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return (int) (Reaction_module_ptr->GetErrorString().size());
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_GetFilePrefix(int id, char *prefix, int l)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (prefix)
		{
			strncpy(prefix, Reaction_module_ptr->GetFilePrefix().c_str(), l);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetGfw(int id, double * gfw)
/* ---------------------------------------------------------------------- */
{
	// Retrieves gram formula weights
	// size of d must be the number of components
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (gfw != NULL)
		{
			size_t ncomps = Reaction_module_ptr->GetComponents().size();
			if (ncomps > 0)
			{
				memcpy(gfw, &Reaction_module_ptr->GetGfw().front(), ncomps * sizeof(double));
				return IRM_OK;
			}
			return IRM_FAIL;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int RM_GetGridCellCount(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetGridCellCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int 
RM_GetIPhreeqcId(int id, int i)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		IPhreeqc * iptr = Reaction_module_ptr->GetIPhreeqcPointer(i);
		if (iptr != NULL)
		{
			return iptr->GetId();
		}
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int 
RM_GetMpiMyself(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetMpiMyself();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int 
RM_GetMpiTasks(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetMpiTasks();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int 
RM_GetNthSelectedOutputUserNumber(int id, int i)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetNthSelectedOutputUserNumber(i);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetSaturation(int id, double * sat)
/* ---------------------------------------------------------------------- */
{
	// Retrieves saturation for all grid nodes in sat
	// size of sat must be the number of grid nodes
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
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
RM_GetSelectedOutput(int id, double * so)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (so != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int 
RM_GetSelectedOutputColumnCount(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputColumnCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int RM_GetSelectedOutputCount(int id)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_GetSelectedOutputHeading(int id, int icol, char *heading, int length)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (heading != NULL)
		{
			std::string head;
			IRM_RESULT return_value = Reaction_module_ptr->GetSelectedOutputHeading(icol, head);
			if (return_value >= 0)
			{
				strncpy(heading, head.c_str(), length);
			}
			return return_value;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
int RM_GetSelectedOutputRowCount(int id)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSelectedOutputRowCount();
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetSolutionVolume(int id, double * v)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (v != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetSpeciesConcentrations(int id, double * species_conc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (species_conc != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int
RM_GetSpeciesCount(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetSpeciesCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetSpeciesD25(int id, double * diffc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (diffc != NULL)
		{
			const std::vector<double> & diffc_vector = Reaction_module_ptr->GetSpeciesD25();
			memcpy(diffc, &diffc_vector.front(), diffc_vector.size()*sizeof(double));
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_GetSpeciesName(int id, int i, char *name, int length)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (name != NULL)
		{
			const std::vector<std::string> & names = Reaction_module_ptr->GetSpeciesNames();
			if (i >= 0 && i < (int) names.size())
			{
				strncpy(name, names[i].c_str(), length);
				return IRM_OK;
			}
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int RM_GetSpeciesSaveOn(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return (Reaction_module_ptr->GetSpeciesSaveOn() ? 1 : 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_GetSpeciesZ(int id, double * z)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (z != NULL)
		{
			const std::vector<double> & z_vector = Reaction_module_ptr->GetSpeciesZ();
			memcpy(z, &z_vector.front(), z_vector.size()*sizeof(double));
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
IRM_RESULT 
RM_GetStartCell(int id, int *sc)
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		const std::vector <int> & l = Reaction_module_ptr->GetStartCell();
		memcpy(sc, &l.front(), l.size() * sizeof(int));
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
int 
RM_GetThreadCount(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetThreadCount();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
double RM_GetTime(int id)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTime();
	}
	return (double) IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
double RM_GetTimeConversion(int id)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTimeConversion();
	}
	return (double) IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
double RM_GetTimeStep(int id)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->GetTimeStep();
	}
	return (double) IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_InitialPhreeqc2Concentrations(
			int id,
			double *boundary_c,
			int n_boundary,
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
	*          boundary_solution2 - list of second solution numbers to be mixed
	*          fraction1 - list of mixing fractions of solution 1
	*
	*   Output: boundary_c - concentrations for boundary conditions
	*                      - dimensions must be >= n_boundary x n_comp
	*
	*/

	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (boundary_c != NULL && boundary_solution1 != NULL)
		{
			std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
			std::vector < double > destination_c, fraction1_vector;
			boundary_solution1_vector.resize(n_boundary);
			memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (n_boundary * sizeof(int)));
			if (boundary_solution2 != NULL)
			{
				boundary_solution2_vector.resize(n_boundary);
				memcpy(&boundary_solution2_vector.front(), boundary_solution2, (size_t) (n_boundary * sizeof(int)));
			}
			if (fraction1 != NULL)
			{
				fraction1_vector.resize(n_boundary);
				memcpy(&fraction1_vector.front(), fraction1, (size_t) (n_boundary * sizeof(double)));
			}
			IRM_RESULT return_value = Reaction_module_ptr->InitialPhreeqc2Concentrations(
				destination_c,
				boundary_solution1_vector,
				boundary_solution2_vector,
				fraction1_vector);		
			if (return_value == IRM_OK)
			{
				memcpy(boundary_c, &destination_c.front(), destination_c.size() * sizeof(double));
			}       
			return return_value;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_InitialPhreeqc2Module(int id,
							  int *initial_conditions1,		// 7 x nxyz end-member 1
							  int *initial_conditions2,		// 7 x nxyz end-member 2
							  double *fraction1)			// 7 x nxyz fraction of end-member 1
/* ---------------------------------------------------------------------- */
{
	// 7 indices for initial conditions
	// 0 solution
	// 1 ppassemblage
	// 2 exchange
	// 3 surface
	// 4 gas phase
	// 5 ss_assemblage
	// 6 kinetics
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (initial_conditions1 != NULL)
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_InitialPhreeqcCell2Module(int id,
                int n,		                            // InitialPhreeqc cell number
                int *module_numbers,		            // Module cell numbers
                int dim_module_numbers)			    // Number of module cell numbers
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (module_numbers != NULL)
		{
			std::vector <int> module_numbers_vector;
			module_numbers_vector.resize(dim_module_numbers);
			memcpy(&module_numbers_vector.front(), module_numbers, (size_t) (dim_module_numbers) * sizeof(int));
			return Reaction_module_ptr->InitialPhreeqcCell2Module(
				n,
				module_numbers_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_InitialPhreeqc2SpeciesConcentrations(
			int id,
			double *species_c,
			int n_boundary,
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

	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (species_c && boundary_solution1)
		{
			std::vector < int > boundary_solution1_vector, boundary_solution2_vector;
			std::vector < double > destination_c, fraction1_vector;
			boundary_solution1_vector.resize(n_boundary);
			memcpy(&boundary_solution1_vector.front(), boundary_solution1, (size_t) (n_boundary * sizeof(int)));
			if (boundary_solution2 != NULL)
			{
				boundary_solution2_vector.resize(n_boundary);
				memcpy(&boundary_solution2_vector.front(), boundary_solution2, (size_t) (n_boundary * sizeof(int)));
			}
			if (fraction1 != NULL)
			{
				fraction1_vector.resize(n_boundary);
				memcpy(&fraction1_vector.front(), fraction1, (size_t) (n_boundary * sizeof(double)));
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
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_LoadDatabase(int id, const char *db_name)
	/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (db_name != NULL)
		{
			std::string db = PhreeqcRM::Char2TrimString(db_name);
			return Reaction_module_ptr->LoadDatabase(db.c_str());
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_LogMessage(int id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (err_str)
		{
			std::string e_string(err_str);
			Reaction_module_ptr->LogMessage(e_string);
			return IRM_OK;
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_MpiWorker(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->MpiWorker();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_MpiWorkerBreak(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->MpiWorkerBreak();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_OpenFiles(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->OpenFiles();
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_OutputMessage(int id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (err_str)
		{
			std::string e_string(err_str);
			Reaction_module_ptr->OutputMessage(e_string);
		}
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_RunCells(int id)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->RunCells(); 
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_RunFile(int id, int workers, int initial_phreeqc, int utility, const char *chem_name)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (chem_name != NULL)
		{
			std::string str = PhreeqcRM::Char2TrimString(chem_name);
			return Reaction_module_ptr->RunFile((workers != 0), (initial_phreeqc != 0), (utility != 0), str.c_str());
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_RunString(int id, int workers, int initial_phreeqc, int utility, const char *input_string)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		std::string str = PhreeqcRM::Char2TrimString(input_string);
		return Reaction_module_ptr->RunString((workers != 0), (initial_phreeqc != 0), (utility != 0), input_string);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_ScreenMessage(int id, const char *err_str)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (err_str)
		{
			std::string e_string(err_str);
			Reaction_module_ptr->ScreenMessage(e_string);
		}
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetComponentH2O(int id, int tf)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetComponentH2O(tf != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetConcentrations(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> c_vector;
			c_vector.resize(Reaction_module_ptr->GetGridCellCount() * Reaction_module_ptr->GetComponentCount());
			memcpy(&c_vector.front(), t, c_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetConcentrations(c_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetCurrentSelectedOutputUserNumber(int id, int i)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetCurrentSelectedOutputUserNumber(i);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetDensity(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> d_vector;
			d_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&d_vector.front(), t, d_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetDensity(d_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetDumpFileName(int id, const char *name)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (name != NULL)
		{
			std::string str = PhreeqcRM::Char2TrimString(name);
			return Reaction_module_ptr->SetDumpFileName(str.c_str());
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetErrorHandlerMode(int id, int mode)
/* ---------------------------------------------------------------------- */
{
	// pass pointers from Fortran to the Reaction module
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetErrorHandlerMode(mode);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetFilePrefix(int id, const char *name)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (name != NULL)
		{
			std::string str = PhreeqcRM::Char2TrimString(name);
			return Reaction_module_ptr->SetFilePrefix(str.c_str());
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetMpiWorkerCallback(int id, int (*fcn)(int *x1, void *cookie))
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetMpiWorkerCallbackC(fcn);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetMpiWorkerCallbackCookie(int id, void *cookie)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetMpiWorkerCallbackCookie(cookie);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetPartitionUZSolids(int id, int t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetPartitionUZSolids(t != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetPorosity(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> v_vector;
			v_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&v_vector.front(), t, v_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetPorosity(v_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetPressure(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> p_vector;
			p_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&p_vector.front(), t, p_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetPressure(p_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetPrintChemistryOn(int id,	 int worker, int ip, int utility)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetPrintChemistryOn(worker != 0, ip != 0, utility != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetPrintChemistryMask(int id, int *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<int> m_vector;
			m_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&m_vector.front(), t, m_vector.size() * sizeof(int));
			return Reaction_module_ptr->SetPrintChemistryMask(m_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetRebalanceFraction(int id, double f)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetRebalanceFraction(f);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetRebalanceByCell(int id, int method)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetRebalanceByCell(method != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetRepresentativeVolume(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> v_vector;
			v_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&v_vector.front(), t, v_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetRepresentativeVolume(v_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetSaturation(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> s_vector;
			s_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&s_vector.front(), t, s_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetSaturation(s_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetScreenOn(int id, int tf)
/* ---------------------------------------------------------------------- */
{
	// pass pointers from Fortran to the Reaction module
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetSelectedOutputOn(tf != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetSelectedOutputOn(int id, int selected_output_on)
/* ---------------------------------------------------------------------- */
{
	// pass pointers from Fortran to the Reaction module
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetSelectedOutputOn(selected_output_on != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SetSpeciesSaveOn(int id, int save_on)
/* ---------------------------------------------------------------------- */
{
	// pass pointers from Fortran to the Reaction module
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetSpeciesSaveOn(save_on != 0);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetTemperature(int id, double *t)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (t != NULL)
		{
			std::vector<double> t_vector;
			t_vector.resize(Reaction_module_ptr->GetGridCellCount());
			memcpy(&t_vector.front(), t, t_vector.size() * sizeof(double));
			return Reaction_module_ptr->SetTemperature(t_vector);
		}
		return IRM_INVALIDARG;
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetTime(int id, double t)
/* ---------------------------------------------------------------------- */
{
	//
	// multiply seconds to convert to user time units
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTime(t);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetTimeConversion(int id, double t)
/* ---------------------------------------------------------------------- */
{
	//
	// multiply seconds to convert to user time units
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTimeConversion(t);
	}
	return IRM_BADINSTANCE;
}

/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetTimeStep(int id, double t)
/* ---------------------------------------------------------------------- */
{
	//
	// multiply seconds to convert to user time units
	//
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetTimeStep(t);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsExchange (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsExchange(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsGasPhase (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsGasPhase(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsKinetics (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsKinetics(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsPPassemblage (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsPPassemblage(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsSolution (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSolution(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsSSassemblage (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSSassemblage(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT 
RM_SetUnitsSurface (int id, int u)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		return Reaction_module_ptr->SetUnitsSurface(u);
	}
	return IRM_BADINSTANCE;
}
/* ---------------------------------------------------------------------- */
IRM_RESULT
RM_SpeciesConcentrations2Module(int id, double * species_conc)
/* ---------------------------------------------------------------------- */
{
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
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
RM_UseSolutionDensityVolume(int id, int tf)
/* ---------------------------------------------------------------------- */
{
	// writes a warning message to screen, log, and output files
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		Reaction_module_ptr->UseSolutionDensityVolume(tf != 0);
		return IRM_OK;
	}
	return IRM_BADINSTANCE;

}
/* --------------------------------------------------------------------- */
IRM_RESULT
RM_WarningMessage(int id, const char *err_str)
/* ---------------------------------------------------------------------- */
{	
	PhreeqcRM * Reaction_module_ptr = PhreeqcRM::GetInstance(id);
	if (Reaction_module_ptr)
	{
		if (err_str)
		{
			std::string e_string(err_str);
			trim_right(e_string);
			Reaction_module_ptr->WarningMessage(e_string);
		}
		return IRM_OK;
	}
	return IRM_BADINSTANCE;
}

