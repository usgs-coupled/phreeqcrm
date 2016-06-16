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

#if defined(__cplusplus)
}
#endif

// C++ function
extern int advection_cpp();
extern int species_cpp();
extern int units_tester();

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

	units_tester();
	advection_cpp();
#if defined(TEST_FORTRAN)
	advection_f90();
#endif
	advection_c();
	species_cpp();
#if defined(TEST_FORTRAN)
	species_f90();
#endif
	species_c();

#if defined(USE_MPI)
	MPI_Finalize();
#endif
	if (mpi_myself == 0)
	{
		std::cerr << "Done with tests." << std::endl;
	}
	return EXIT_SUCCESS;
}
