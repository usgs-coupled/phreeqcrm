#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>


#if defined(__cplusplus)
extern "C" {
#endif

extern void Advect_c(void);
extern void Advect_f90(void);
extern void AdvectBMI_f90(void);
extern void AdvectBMI_f90_test(void);
extern void Gas_c(void);
extern void Gas_f90(void);
extern void SimpleAdvect_c(void);
extern void SimpleAdvect_f90(void);
extern void Species_c(void);
extern void Species_f90(void);
extern void TestAllMethods_c(void);
extern void TestAllMethods_f90(void);
extern void WriteYAMLFile_f90(void);
extern void WriteYAMLFile_f90_test(void);

#if defined(__cplusplus)
}
#endif

// C++ function
extern int SimpleAdvect_cpp();
extern int Advect_cpp();
extern int AdvectBMI_cpp();
extern int AdvectBMI_cpp_test();
extern int Species_cpp();
extern int units_tester();
extern int Gas_cpp();
extern void TestAllMethods_cpp();
extern void WriteYAMLFile_cpp();
extern void WriteYAMLFile_cpp_test();

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
	bool root = (mpi_myself == 0);
	units_tester();
	if (root) std::cerr << "Done units_tester.===================================" << std::endl;
	SimpleAdvect_cpp();
	if (root) std::cerr << "Done SimpleAdvection_cpp.============================" << std::endl;
	Advect_cpp();
	if (root) std::cerr << "Done Advect_cpp.=====================================" << std::endl;
#ifdef USE_YAML
	if (root) WriteYAMLFile_cpp();
#if defined(USE_MPI)
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) exit(3); //  make sure yaml is finished being written
#endif
	if (root) std::cerr << "Done WriteYAMLFile_cpp.==============================" << std::endl;
	AdvectBMI_cpp();
	if (root) std::cerr << "Done AdvectBMI_cpp.==================================" << std::endl;
	if (root) WriteYAMLFile_cpp_test();
	AdvectBMI_cpp_test();
	if (root) std::cerr << "Done AdvectBMI_cpp_test.=============================" << std::endl;
	TestAllMethods_cpp();
	if (root) std::cerr << "Done TestAllMethods_cpp.=============================" << std::endl;
	TestAllMethods_c();
	if (root) std::cerr << "Done TestAllMethods_c.=============================" << std::endl;
#endif
	SimpleAdvect_c();
	if (root) std::cerr << "Done SimpleAdvect_c.=================================" << std::endl;
	Advect_c();
	if (root) std::cerr << "Done Advect_c.=======================================" << std::endl;
	Species_cpp();
	if (root) std::cerr << "Done Species_cpp.====================================" << std::endl;
	Species_c();
	if (root) std::cerr << "Done Species_c.======================================" << std::endl;
	Gas_cpp();
	if (root) std::cerr << "Done Gas_ccp.========================================" << std::endl;
	Gas_c();
	if (root) std::cerr << "Done Gas_c.==========================================" << std::endl;
#if defined(TEST_FORTRAN)
	SimpleAdvect_f90();
	if (root) std::cerr << "Done SimpleAdvect_f90.===============================" << std::endl;
	Advect_f90();
	if (root) std::cerr << "Done Advect_f90.=====================================" << std::endl;
#ifdef USE_YAML
    if (root) WriteYAMLFile_f90();
#if defined(USE_MPI)
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) exit(3); //  make sure yaml is finished being written
#endif
	if (root) std::cerr << "Done WriteYAMLFile_f90.==============================" << std::endl;
	AdvectBMI_f90();
	if (root) std::cerr << "Done AdvectBMI_f90.==================================" << std::endl;
	if (root) WriteYAMLFile_f90_test();
	AdvectBMI_f90_test();
	if (root) std::cerr << "Done AdvectBMI_f90_test.=============================" << std::endl;
	TestAllMethods_f90();
	if (root) std::cerr << "Done TestAllMethods_f90.=============================" << std::endl;
#endif
	Species_f90();
	if (root) std::cerr << "Done Species_f90.====================================" << std::endl;
	Gas_f90();
	if (root) std::cerr << "Done Gas_f90.========================================" << std::endl;
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
