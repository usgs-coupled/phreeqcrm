#include <cstdlib>
#include "BMIPhreeqcRM.h"

// The template method PhreeqcRM::GetInstance<T>()) doesn't work
// when PhreeqcRM is compiled as a dll on windows
int main(void)
{
    BMIPhreeqcRM* bmi = new BMIPhreeqcRM;
    int idx = bmi->GetIndex();
    assert(bmi == PhreeqcRM::GetInstance(idx));
    // assert(bmi == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    assert(bmi == BMIPhreeqcRM::GetInstance(idx));
    delete bmi;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
    // assert(nullptr == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    assert(nullptr == BMIPhreeqcRM::GetInstance(idx));
    return EXIT_SUCCESS;
}
