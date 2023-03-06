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
INTERFACE RM_BMI_GetValue
    procedure RM_BMI_GetValue_c
    procedure RM_BMI_GetValue_c1
    procedure RM_BMI_GetValue_d
    procedure RM_BMI_GetValue_d1
    procedure RM_BMI_GetValue_d2
    procedure RM_BMI_GetValue_i
    procedure RM_BMI_GetValue_i1
    procedure RM_BMI_GetValue_i2
END INTERFACE RM_BMI_GetValue

INTERFACE RM_BMI_SetValue
  procedure RM_BMI_SetValue_c
  procedure RM_BMI_SetValue_d
  procedure RM_BMI_SetValue_d1
  procedure RM_BMI_SetValue_d2
  procedure RM_BMI_SetValue_i
  procedure RM_BMI_SetValue_i1
  procedure RM_BMI_SetValue_i2
END INTERFACE RM_BMI_SetValue 

    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarItemsize(id, var) &
            BIND(C, NAME='RMF_BMI_GetVarItemsize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        END FUNCTION RMF_BMI_GetVarItemsize

        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarNbytes(id, var) &
            BIND(C, NAME='RMF_BMI_GetVarNbytes')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        END FUNCTION RMF_BMI_GetVarNbytes
    END INTERFACE
    CONTAINS
    
INTEGER FUNCTION RM_BMI_Finalize(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RM_Destroy(id) &
            BIND(C, NAME='RM_Destroy')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RM_Destroy 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_Finalize = RM_Destroy(id)
    return
END FUNCTION RM_BMI_Finalize 

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
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetInputItemCount(id) &
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
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetOutputItemCount(id) &
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


INTEGER FUNCTION RM_BMI_GetValue_c(id, var, dest)
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
    RM_BMI_GetValue_c = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_GetValue_c
INTEGER FUNCTION RM_BMI_GetValue_c1(id, var, dest)
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
    CHARACTER(len=*), INTENT(out) :: dest(:)
    if ((RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) .ne. SIZE(dest)) &
        .or. (len(dest(1)) .ne. RMF_BMI_GetVarItemsize(id, var))) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_GetValue_c1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
END FUNCTION RM_BMI_GetValue_c1

INTEGER FUNCTION RM_BMI_GetValue_d(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            REAL(KIND=C_DOUBLE), INTENT(out) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(out) :: dest
    RM_BMI_GetValue_d = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_GetValue_d

INTEGER FUNCTION RM_BMI_GetValue_d1(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
        BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(out) :: dest(:)
    if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_GetValue_d1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
END FUNCTION RM_BMI_GetValue_d1

INTEGER FUNCTION RM_BMI_GetValue_d2(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
        BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
   
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(out) :: dest(:,:)
    if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest, 1)*SIZE(dest, 2) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_GetValue_d2 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
END FUNCTION RM_BMI_GetValue_d2

INTEGER FUNCTION RM_BMI_GetValue_i(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest
    RM_BMI_GetValue_i = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_GetValue_i

INTEGER FUNCTION RM_BMI_GetValue_i1(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(out) :: dest(:)
    if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_GetValue_i1 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
END FUNCTION RM_BMI_GetValue_i1

INTEGER FUNCTION RM_BMI_GetValue_i2(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
    INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
        BIND(C, NAME='RMF_BMI_GetValue')
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(KIND=C_INT), INTENT(in) :: id
        CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
        INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
   
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(out) :: dest(:,:)
    if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest, 1)*SIZE(dest, 2) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_GetValue_i2 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
END FUNCTION RM_BMI_GetValue_i2

INTEGER FUNCTION RM_BMI_GetVarItemsize(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    RM_BMI_GetVarItemsize = RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR)
END FUNCTION RM_BMI_GetVarItemsize

INTEGER FUNCTION RM_BMI_GetVarNbytes(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE

    INTEGER, INTENT(in) :: id
    CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
    RM_BMI_GetVarNbytes = RMF_BMI_GetVarNbytes(id, var)
END FUNCTION RM_BMI_GetVarNbytes

INTEGER FUNCTION RM_BMI_GetVarType(id, var, vtype)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarType(id, var, vtype, l) &
            BIND(C, NAME='RMF_BMI_GetVarType')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            CHARACTER(KIND=C_CHAR), INTENT(out) :: vtype(*)
        END FUNCTION RMF_BMI_GetVarType 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(out) :: vtype
    RM_BMI_GetVarType = RMF_BMI_GetVarType(id, trim(var)//C_NULL_CHAR, vtype, len(vtype))
    return
END FUNCTION RM_BMI_GetVarType 

INTEGER FUNCTION RM_BMI_GetVarUnits(id, var, units)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetVarUnits(id, var, units, l) &
            BIND(C, NAME='RMF_BMI_GetVarUnits')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            CHARACTER(KIND=C_CHAR), INTENT(out) :: units(*)
        END FUNCTION RMF_BMI_GetVarUnits 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(out) :: units
    RM_BMI_GetVarUnits = RMF_BMI_GetVarUnits(id, trim(var)//C_NULL_CHAR, units, len(units))
    return
END FUNCTION RM_BMI_GetVarUnits 

INTEGER FUNCTION RM_BMI_Initialize(id, config_file)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_Initialize(id, config_file) &
            BIND(C, NAME='RMF_BMI_Initialize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: config_file(*)
        END FUNCTION RMF_BMI_Initialize 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: config_file
    RM_BMI_Initialize = RMF_BMI_Initialize(id, trim(config_file//C_NULL_CHAR))
    return
END FUNCTION RM_BMI_Initialize

INTEGER FUNCTION RM_BMI_SetValue_c(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            CHARACTER(KIND=C_CHAR), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(out) :: dest
    RM_BMI_SetValue_c = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_SetValue_c

INTEGER FUNCTION RM_BMI_SetValue_i(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(out) :: dest
    RM_BMI_SetValue_i = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_SetValue_i

INTEGER FUNCTION RM_BMI_SetValue_i1(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(out) :: dest(:)
	if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_SetValue_i1 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
END FUNCTION RM_BMI_SetValue_i1

INTEGER FUNCTION RM_BMI_SetValue_i2(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            INTEGER(KIND=C_INT), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(out) :: dest(:,:)
	if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest,1)*SIZE(dest,2) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_SetValue_i2 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
END FUNCTION RM_BMI_SetValue_i2

INTEGER FUNCTION RM_BMI_SetValue_d(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            REAL(KIND=C_DOUBLE), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(inout) :: dest
    RM_BMI_SetValue_d = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_SetValue_d

INTEGER FUNCTION RM_BMI_SetValue_d1(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            REAL(KIND=C_DOUBLE), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(out) :: dest(:)
	if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_SetValue_d1 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest(1))
    return
END FUNCTION RM_BMI_SetValue_d1

INTEGER FUNCTION RM_BMI_SetValue_d2(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            REAL(KIND=C_DOUBLE), INTENT(out) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE


    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(out) :: dest(:,:)
	if (RMF_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var) &
        .ne. SIZE(dest,1)*SIZE(dest,2) ) then
        stop "Dimension error in RMF_BMI_GetValue"
    endif
    RM_BMI_SetValue_d2 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
END FUNCTION RM_BMI_SetValue_d2

INTEGER FUNCTION RM_BMI_Update(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_Update(id) &
            BIND(C, NAME='RMF_BMI_Update')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_Update
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_Update = RMF_BMI_Update(id)
    return
END FUNCTION RM_BMI_Update

INTEGER FUNCTION RM_GetGridCellCountYAML(config_file)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetGridCellCountYAML(config_file) &
            BIND(C, NAME='RMF_GetGridCellCountYAML')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(in) :: config_file(*)
        END FUNCTION RMF_GetGridCellCountYAML 
    END INTERFACE
    CHARACTER(len=*), INTENT(in) :: config_file
    RM_GetGridCellCountYAML = RMF_GetGridCellCountYAML(trim(config_file)//C_NULL_CHAR)
    return
END FUNCTION RM_GetGridCellCountYAML 



END MODULE BMI_PhreeqcRM




    
