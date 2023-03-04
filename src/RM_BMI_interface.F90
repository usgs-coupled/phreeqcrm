!*MODULE BMI_PhreeqcRM PhreeqcRM BMI interface
!> @brief Fortran Documentation for the geochemical reaction module PhreeqcRM. 
!> @par "" 
!> "USE PhreeqcRM" is included in Fortran source code to define the PhreeqcRM functions.
!> For Windows, define the module by including the file RM_interface.F90 in your project.
!> For Linux, configure, compile, and install the PhreeqcRM library and module file. 
!> You will need installed include directory (-I) added to the project) to reference the module file.
!> You will need to link to the library to produce the executable for your code.
!>
    MODULE BMI_PhreeqcRM
    IMPLICIT NONE
    
    CONTAINS
    
INTEGER FUNCTION RM_BMI_GetComponentName(id, component_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetComponentName(id, component_name, l) &
            BIND(C, NAME='RMF_BMI_GetComponentName')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: component_name(*)
        END FUNCTION RMF_BMI_GetComponentName 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: component_name
    RM_BMI_GetComponentName = RMF_BMI_GetComponentName(id, component_name, len(component_name))
    return
END FUNCTION RM_BMI_GetComponentName 

DOUBLE PRECISION FUNCTION RM_BMI_GetCurrentTime(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetCurrentTime(id) &
            BIND(C, NAME='RMF_BMI_GetCurrentTime')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetCurrentTime
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_GetCurrentTime = RMF_BMI_GetCurrentTime(id)
END FUNCTION RM_BMI_GetCurrentTime

INTEGER FUNCTION RM_BMI_GetInputItemCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetInputItemCount(id) &
            BIND(C, NAME='RMF_BMI_GetInputItemCount')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetInputItemCount
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_GetInputItemCount = RMF_BMI_GetInputItemCount(id)
END FUNCTION RM_BMI_GetInputItemCount

INTEGER FUNCTION RM_BMI_GetInputVarNames(id, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetInputVarNames(id, var_names, l) &
            BIND(C, NAME='RMF_BMI_GetInputVarNames')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: var_names(*)
        END FUNCTION RMF_BMI_GetInputVarNames 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: var_names
    RM_BMI_GetInputVarNames = RMF_BMI_GetInputVarNames(id, var_names, len(var_names))
    return
END FUNCTION RM_BMI_GetInputVarNames 

INTEGER FUNCTION RM_BMI_GetOutputItemCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetOutputItemCount(id) &
            BIND(C, NAME='RMF_BMI_GetOutputItemCount')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetOutputItemCount
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_GetOutputItemCount = RMF_BMI_GetOutputItemCount(id)
END FUNCTION RM_BMI_GetOutputItemCount

INTEGER FUNCTION RM_BMI_GetOutputVarNames(id, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetOutputVarNames(id, var_names, l) &
            BIND(C, NAME='RMF_BMI_GetOutputVarNames')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: var_names(*)
        END FUNCTION RMF_BMI_GetOutputVarNames 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: var_names
    RM_BMI_GetOutputVarNames = RMF_BMI_GetOutputVarNames(id, var_names, len(var_names))
    return
END FUNCTION RM_BMI_GetOutputVarNames 

DOUBLE PRECISION FUNCTION RM_BMI_GetTimeStep(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetTimeStep(id) &
            BIND(C, NAME='RMF_BMI_GetTimeStep')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetTimeStep
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_GetTimeStep = RMF_BMI_GetTimeStep(id)
END FUNCTION RM_BMI_GetTimeStep

INTEGER FUNCTION RM_BMI_GetTimeUnits(id, time_units)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetTimeUnits(id, time_units, l) &
            BIND(C, NAME='RMF_BMI_GetTimeUnits')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: time_units(*)
        END FUNCTION RMF_BMI_GetTimeUnits 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: time_units
    RM_BMI_GetTimeUnits = RMF_BMI_GetTimeUnits(id, time_units, len(time_units))
    return
END FUNCTION RM_BMI_GetTimeUnits 


INTEGER FUNCTION RM_BMI_GetValue(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            CHARACTER(KIND=C_CHAR), INTENT(out) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(out) :: dest
    RM_BMI_GetValue = RMF_BMI_GetValue(id, var, dest)
    return
END FUNCTION RM_BMI_GetValue 

INTEGER FUNCTION RM_BMI_GetVarItemsize(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetVarItemsize(id, var) &
            BIND(C, NAME='RMF_BMI_GetVarItemsize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        END FUNCTION RMF_BMI_GetVarItemsize
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    RM_BMI_GetVarItemsize = RMF_BMI_GetVarItemsize(id, var)
END FUNCTION RM_BMI_GetVarItemsize

INTEGER FUNCTION RM_BMI_GetVarNbytes(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetVarNbytes(id, var) &
            BIND(C, NAME='RMF_BMI_GetVarNbytes')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        END FUNCTION RMF_BMI_GetVarNbytes
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    RM_BMI_GetVarNbytes = RMF_BMI_GetVarNbytes(id, var)
END FUNCTION RM_BMI_GetVarNbytes




END MODULE BMI_PhreeqcRM


    
