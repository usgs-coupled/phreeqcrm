#include <cstdlib>
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include "PhreeqcRM.h"

int main(int argc, char* argv[])
{
#if defined(USE_MPI)
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        return EXIT_FAILURE;
    }
    PhreeqcRM* rm = new PhreeqcRM(10, MPI_COMM_WORLD);
#else
    PhreeqcRM* rm = new PhreeqcRM(10, 2);
#endif
    int idx = rm->GetIndex();
    assert(rm == PhreeqcRM::GetInstance(idx));
    delete rm;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
#if defined(USE_MPI)
    if (MPI_Finalize() != MPI_SUCCESS)
    {
        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}
