#include <cstdlib>
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include "BMIPhreeqcRM.h"


// The template method PhreeqcRM::GetInstance<T>()) doesn't work
// when PhreeqcRM is compiled as a dll on windows
int main(int argc, char* argv[])
{
#if defined(USE_MPI)
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        return EXIT_FAILURE;
    }
#endif
    BMIPhreeqcRM* bmi = new BMIPhreeqcRM();
    int idx = bmi->GetIndex();
    assert(bmi == PhreeqcRM::GetInstance(idx));
    // assert(bmi == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    assert(bmi == BMIPhreeqcRM::GetInstance(idx));
    delete bmi;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
    // assert(nullptr == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    assert(nullptr == BMIPhreeqcRM::GetInstance(idx));
#if defined(USE_MPI)
    if (MPI_Finalize() != MPI_SUCCESS)
    {
        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}
