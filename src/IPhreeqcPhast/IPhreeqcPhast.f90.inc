
      include "IPhreeqc.f90.inc"

       INTERFACE
        FUNCTION CreateIPhreeqcPhast()
         INTEGER(KIND=4) :: CreateIPhreeqcPhast
        END FUNCTION CreateIPhreeqcPhast
       END INTERFACE


       INTERFACE
        FUNCTION DestroyIPhreeqcPhast(ID)
         INTEGER(KIND=4), INTENT(IN) :: ID
         INTEGER(KIND=4)             :: DestroyIPhreeqcPhast
        END FUNCTION DestroyIPhreeqcPhast
       END INTERFACE