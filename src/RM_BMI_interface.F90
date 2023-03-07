
    MODULE BMI_PhreeqcRM
    IMPLICIT NONE
INTERFACE RM_BMI_GetValue
    procedure RM_BMI_GetValue_b
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
  procedure RM_BMI_SetValue_b
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

!> Basic Model Interface method that returns the component name--PhreeqcRM. The BMI interface to PhreeqcRM is
!> only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the BMI_
!> prefix) provide a complete interface, and it is expected that the native methods will be used in preference to the BMI_
!> methods.
!> 
!> @param id            The instance @a id returned from @ref RM_Create.
!> @param component_name is filled with "PhreeqcRM", the name of the component.
!> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_BMI_GetComponentName(id, component_name)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.
INTEGER FUNCTION RM_BMI_GetComponentName(id, component_name)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetComponentName(id, component_name, l) &
            BIND(C, NAME='RMF_BMI_GetComponentName')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: component_name(*)
        END FUNCTION RMF_BMI_GetComponentName 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: component_name
    RM_BMI_GetComponentName = RMF_BMI_GetComponentName(id, component_name, len(component_name))
    return
END FUNCTION RM_BMI_GetComponentName 


!> Basic Model Interface method that returns the current simulation time, in seconds. (Same as @ref GetTime.)
!> The reaction module does not change the time value, so the
!> returned value is equal to the default (0.0) or the last time set by
!> @ref BMI_SetValue("Time", time) or @ref RM_SetTime.
!> @param id            The instance @a id returned from @ref RM_Create.
!> @retval                 The current simulation time, in seconds.
!> @see
!> @ref RM_BMI_GetEndTime,
!> @ref RM_BMI_GetTimeStep,
!> @ref RM_BMI_SetValue,
!> @ref RM_GetTime,
!> @ref RM_GetTimeStep,
!> @ref RM_SetTime,
!> @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> time = RM_BMI_GetCurrentTime(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

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

!> Basic Model Interface method that returns @ref RM_BMI_GetCurrentTime plus 
!> @ref RM_BMI_GetTimeStep, in seconds.
!> @param id            The instance @a id returned from @ref RM_Create.
!> @retval                 The end of the time step, in seconds.
!> @see
!> @ref RM_BMI_GetCurrentTime,
!> @ref RM_BMI_GetTimeStep,
!> @ref RM_BMI_SetValue,
!> @ref RM_GetTime,
!> @ref RM_GetTimeStep,
!> @ref RM_SetTime,
!> @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> time = RM_BMI_GetEndTime(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

DOUBLE PRECISION FUNCTION RM_BMI_GetEndTime(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_BMI_GetEndTime(id) &
            BIND(C, NAME='RMF_BMI_GetEndTime')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_BMI_GetEndTime
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_BMI_GetEndTime = RMF_BMI_GetEndTime(id)
END FUNCTION RM_BMI_GetEndTime

!> Basic Model Interface method that returns count of input variables that 
!> can be set with @ref RM_BMI_SetValue.
!> @retval  Count of input variables that can be set with @ref RM_BMI_SetValue.
!> @param id            The instance @a id returned from @ref RM_Create.
!> 
!> @see
!> @ref RM_BMI_GetInputVarNames,
!> @ref RM_BMI_GetVarItemsize,
!> @ref RM_BMI_GetVarNbytes,
!> @ref RM_BMI_GetVarType,
!> @ref RM_BMI_GetVarUnits,
!> @ref RM_BMI_SetValue.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> 		std::vector<std::string> InputVarNames = phreeqc_rm.BMI_GetInputVarNames();
!> 		int count = phreeqc_rm.BMI_GetInputItemCount();
!> 		oss << "BMI_SetValue variables:\n";
!> 		for (size_t i = 0; i < count; i++)
!> 		{
!> 			oss << "  " << i << "  " << InputVarNames[i] << "\n";
!> 			oss << "     Type:        " << phreeqc_rm.BMI_GetVarType(InputVarNames[i]) << "\n";
!> 			oss << "     Units:       " << phreeqc_rm.BMI_GetVarUnits(InputVarNames[i]) << "\n";
!> 			oss << "     Total bytes: " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) << "\n";
!> 			oss << "     Item bytes:  " << phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
!> 			oss << "     Dim:         " << phreeqc_rm.BMI_GetVarNbytes(InputVarNames[i]) /
!> 										   phreeqc_rm.BMI_GetVarItemsize(InputVarNames[i]) << "\n";
!> 		}
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.
!> 
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: var_names(*)
        END FUNCTION RMF_BMI_GetInputVarNames 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: var_names
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: var_names(*)
        END FUNCTION RMF_BMI_GetOutputVarNames 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: var_names
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: time_units(*)
        END FUNCTION RMF_BMI_GetTimeUnits 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: time_units
    RM_BMI_GetTimeUnits = RMF_BMI_GetTimeUnits(id, time_units, len(time_units))
    return
END FUNCTION RM_BMI_GetTimeUnits 

INTEGER FUNCTION RM_BMI_GetValue_b(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_GetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_GetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            LOGICAL(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL(KIND=4), INTENT(inout) :: dest
   RM_BMI_GetValue_b = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_GetValue_b

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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: dest
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: dest(:)
    if ((RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) .ne. SIZE(dest)) &
        .or. (len(dest(1)) .ne. RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR))) then
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
            REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(inout) :: dest
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
    double precision, INTENT(inout) :: dest(:)
    if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
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
    double precision, INTENT(inout) :: dest(:,:)
    if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
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
            INTEGER(KIND=C_INT), INTENT(inout) :: dest
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
            INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:)
    if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
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
        INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_GetValue 
    END INTERFACE
   
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:,:)
    if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
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
    CHARACTER(len=*), INTENT(in) :: var
    RM_BMI_GetVarNbytes = RMF_BMI_GetVarNbytes(id, trim(var)//C_NULL_CHAR)
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: vtype(*)
        END FUNCTION RMF_BMI_GetVarType 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: vtype
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
            CHARACTER(KIND=C_CHAR), INTENT(inout) :: units(*)
        END FUNCTION RMF_BMI_GetVarUnits 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(inout) :: units
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

INTEGER FUNCTION RM_BMI_SetValue_b(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_BMI_SetValue(id, var, dest) &
            BIND(C, NAME='RMF_BMI_SetValue')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: var(*)
            LOGICAL(KIND=C_BOOL), INTENT(in) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL, INTENT(in) :: dest
    RM_BMI_SetValue_b = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest)
    return
END FUNCTION RM_BMI_SetValue_b

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
            CHARACTER(KIND=C_CHAR), INTENT(in) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    CHARACTER(len=*), INTENT(in) :: dest
    RM_BMI_SetValue_c = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, trim(dest)//C_NULL_CHAR)
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
            INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest
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
            INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:)
	if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_SetValue"
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
            INTEGER(KIND=C_INT), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE
    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:,:)
	if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
        .ne. SIZE(dest,1)*SIZE(dest,2) ) then
        stop "Dimension error in RMF_BMI_SetValue"
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
            REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
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
            REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(inout) :: dest(:)
	if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
        .ne. SIZE(dest) ) then
        stop "Dimension error in RMF_BMI_SetValue"
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
            REAL(KIND=C_DOUBLE), INTENT(inout) :: dest
        END FUNCTION RMF_BMI_SetValue 
    END INTERFACE

    integer, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    double precision, INTENT(inout) :: dest(:,:)
	if (RMF_BMI_GetVarNBytes(id, trim(var)//C_NULL_CHAR) / &
        RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR) &
        .ne. SIZE(dest,1)*SIZE(dest,2) ) then
        stop "Dimension error in RMF_BMI_SetValue"
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




    
