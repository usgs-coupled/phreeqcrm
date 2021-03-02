#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>


#if defined(__cplusplus)
extern "C" {
#endif
	
extern void advection_f90(void);
extern void advection_c(void);
extern void species_f90(void);
extern void species_c(void);
extern void gas_c(void);
extern void gas_f90(void);

#if defined(__cplusplus)
}
#endif

// C++ function
extern int advection_cpp();
extern int species_cpp();
extern int units_tester();
extern int gas_cpp();

int main(int argc, char* argv[])
{
	int mpi_tasks;
	int mpi_myself;

#if defined(USE_MPI)
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		return EXIT_FAILURE;
	}

	if (MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks) != MPI_SUCCESS)
	{
		return EXIT_FAILURE;
	}

	if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
	{
		return EXIT_FAILURE;
	}
#else
	mpi_tasks = 1;
	mpi_myself = 0;
#endif
	std::cerr << mpi_myself << std::endl;
	bool root = (mpi_myself == 0);
	units_tester();
	if (root) std::cerr << "Done units_tester.===================================" << std::endl;
	advection_cpp();
	if (root) std::cerr << "Done advection_cpp.==================================" << std::endl;
	advection_c();
	if (root) std::cerr << "Done advection_c.====================================" << std::endl;
	species_cpp();
	if (root) std::cerr << "Done species_cpp.====================================" << std::endl;
	species_c();
	if (root) std::cerr << "Done species_c.======================================" << std::endl;
	gas_cpp();
	if (root) std::cerr << "Done gas_ccp.========================================" << std::endl;
	gas_c();
	if (root) std::cerr << "Done gas_c.==========================================" << std::endl;
#if defined(TEST_FORTRAN)
	advection_f90();
	if (root) std::cerr << "Done advection_f90.==================================" << std::endl;
	species_f90();
	if (root) std::cerr << "Done species_f90.====================================" << std::endl;
	gas_f90();
	if (root) std::cerr << "Done gas_f90.========================================" << std::endl;
#endif
#if defined(USE_MPI)
	MPI_Finalize();
#endif
	if (mpi_myself == 0)
	{
		std::cerr << "Done with tests." << std::endl;
	}
	return EXIT_SUCCESS;
}
