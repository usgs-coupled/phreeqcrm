#include <cstdlib>
#include "PhreeqcRM.h"

int main(void)
{
    PhreeqcRM* rm = new PhreeqcRM(10, 2);
    delete rm;
    return EXIT_SUCCESS;
}
