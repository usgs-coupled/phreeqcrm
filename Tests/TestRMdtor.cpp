#include <cstdlib>
#include "PhreeqcRM.h"

int main(void)
{
    PhreeqcRM* rm = new PhreeqcRM(10, 2);
    int idx = rm->GetIndex();
    assert(rm == PhreeqcRM::GetInstance(idx));
    delete rm;
    assert(nullptr == PhreeqcRM::GetInstance(idx));
    return EXIT_SUCCESS;
}
