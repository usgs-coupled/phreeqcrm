#include <cstdlib>
#include "BMIPhreeqcRM.h"

int main(void)
{
    BMIPhreeqcRM* bmi = new BMIPhreeqcRM;
    int idx = bmi->GetIndex();
    assert(bmi == PhreeqcRM::GetInstance(idx));
    assert(bmi == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    delete bmi;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
    assert(nullptr == PhreeqcRM::GetInstance<BMIPhreeqcRM>(idx));
    return EXIT_SUCCESS;
}
