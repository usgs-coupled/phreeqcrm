#include <cstdlib>
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include "PhreeqcRM.h"

int main(int argc, char* argv[])
{
#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif
    PhreeqcRM* rm = new PhreeqcRM(10, 2);
    int idx = rm->GetIndex();
    assert(rm == PhreeqcRM::GetInstance(idx));
    delete rm;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
    return EXIT_SUCCESS;
}
