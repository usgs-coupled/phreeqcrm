
        !> Destroys a reaction module, same as @ref RM_Destroy.
    !> @param id               The instance @a id returned from @ref RM_Create.
    !> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
    !> @see
    !> @ref RM_Create.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = RM_BMI_Finalize(id)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root and workers.
    
    INTEGER FUNCTION RM_BMI_Finalize(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
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
    !> integer inputvarcount
    !> inputvarcount = RM_BMI_GetInputItemCount(id);
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

    !> Basic Model Interface method that returns a list of the variable names that can be set
    !> with @ref RM_BMI_SetValue.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var_names     Deferred length, allocatable, 1D character vector.
    !> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
    !>
    !> @see
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetValue,
    !> @ref RM_BMI_GetVarItemsize,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_GetVarType,
    !> @ref RM_BMI_GetVarUnits,
    !> @ref RM_BMI_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), dimension(:), allocatable          :: inputvars
    !> status = RM_BMI_GetInputVarNames(id, inputvars)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

    INTEGER FUNCTION RM_BMI_GetInputVarNames(id, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), allocatable, INTENT(inout) :: var_names(:)
    RM_BMI_GetInputVarNames = RM_BMI_GetValue(id, "inputvarnames", var_names)
    return
    END FUNCTION RM_BMI_GetInputVarNames

    !> Basic Model Interface method that returns count of output variables that can be
    !> retrieved with @ref RM_BMI_GetValue.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @retval  Count of output variables that can be retrieved with @ref RM_BMI_GetValue.
    !>
    !> @see
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetValue,
    !> @ref RM_BMI_GetVarItemsize,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_GetVarType,
    !> @ref RM_BMI_GetVarUnits.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer outputvarcount
    !> outputvarcount = RM_BMI_GetOutputItemCount(id);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.


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

    !> Basic Model Interface method that returns a list of the variable names that can be
    !> retrieved with @ref RM_BMI_GetValue.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var_names     Deferred length, allocatable, 1D character vector.
    !> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
    !>
    !> @see
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetValue,
    !> @ref RM_BMI_GetVarItemsize,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_GetVarType,
    !> @ref RM_BMI_GetVarUnits.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(len=:), allocatable, dimension(:) :: var_names
    !> var_names = RM_BMI_GetOutputVarNames(id);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

    INTEGER FUNCTION RM_BMI_GetOutputVarNames(id, var_names)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), allocatable,  INTENT(inout) :: var_names(:)
    RM_BMI_GetOutputVarNames = RM_BMI_GetValue(id, "outputvarnames", var_names)
    return
    END FUNCTION RM_BMI_GetOutputVarNames

    !> Basic Model Interface method that returns the current simulation time step,
    !> in seconds. (Same as @ref RM_GetTimeStep.)
    !> The reaction module does not change the time-step value, so the
    !> returned value is equal to the last time step set by
    !> @ref RM_BMI_SetValue("TimeStep", time_step) or @ref SetTimeStep.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @retval                 The current simulation time step, in seconds.
    !> @see
    !> @ref RM_BMI_GetCurrentTime,
    !> @ref RM_BMI_GetEndTime,
    !> @ref RM_BMI_SetValue,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> time_step = RM_BMI_GetTimeStep(id)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

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

    !> Basic Model Interface method that returns the time units of PhreeqcRM.
    !> All time units are seconds for PhreeqcRM.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param time_units    The instance @a id returned from @ref RM_Create.
    !> @retval                 Returns the string "seconds".
    !> @see
    !> @ref RM_BMI_GetCurrentTime,
    !> @ref RM_BMI_GetEndTime,
    !> @ref RM_BMI_GetTimeStep,
    !> @ref RM_BMI_SetValue,
    !> @ref RM_GetTime,
    !> @ref RM_GetTimeStep,
    !> @ref RM_SetTime,
    !> @ref RM_SetTimeStep,
    !> @ref RM_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> character(20) time_units
    !> status = RM_BMI_GetTimeUnits(id, time_units) << ".\n";
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

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

    !> Basic Model Interface method that retrieves the units of a
    !> variable that can be set with
    !> @ref RM_BMI_SetValue or retrieved with @ref RM_BMI_GetValue.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetInputVarNames can be set.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve total bytes.
    !> @param vtype Units of the variable.
    !> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
    !> @see
    !> @ref RM_BMI_GetInputVarNames,
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetVarType.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do i = 1, size(inputvars)
    !>     write(*,"(1x, I4, A40)") i, trim(inputvars(i))
    !>     status = RM_RM_BMI_GetVarUnits(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = RM_RM_BMI_GetVarType(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarItemsize(id, inputvars(i))
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarNbytes(id, inputvars(i))
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.


    !> Basic Model Interface method that retrieves model variables. Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved. The BMI interface to PhreeqcRM is
    !> only partial, and provides only the most basic functions. The native PhreeqcRM methods (those without the the RM_BMI_
    !> prefix) provide a complete interface.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve.
    !> @param dest Variable in which to place results.
    !> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
    !>
    !> Variable names for the second argument (@a name) and variable type of the
    !> third argument (@a dest).
    !> @n "ComponentCount", @a dest: integer;
    !> @n "Components", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "Concentrations", @a dest: double precision, allocatable, dimension(:,:);
    !> @n "CurrentSelectedOutputUserNumber", @a dest: integer;
    !> @n "Density", @a dest: double precision, allocatable, dimension(:);
    !> @n "ErrorString", @a dest: character;
    !> @n "FilePrefix", @a dest: character;
    !> @n "Gfw", @a dest: double precision, allocatable, dimension(:);
    !> @n "GridCellCount", @a dest: integer;
    !> @n "InputVarNames", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "OutputVarNames", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "Porosity", @a dest: double precision, allocatable, dimension(:);
    !> @n "Pressure", @a dest: double precision, allocatable, dimension(:);
    !> @n "Saturation", @a dest: double precision, allocatable, dimension(:);
    !> @n "SelectedOutput", @a dest: double precision, allocatable, dimension(:,:);
    !> @n "SelectedOutputColumnCount", @a dest: integer;
    !> @n "SelectedOutputCount", @a dest: integer;
    !> @n "SelectedOutputHeadings", @a dest: character(len=:), allocatable, dimension(:);
    !> @n "SelectedOutputOn", @a dest: logical;
    !> @n "SelectedOutputRowCount", @a dest: integer;
    !> @n "SolutionVolume", @a dest: double precision, allocatable, dimension(:);
    !> @n "Temperature", @a dest: double precision, allocatable, dimension(:);
    !> @n "Time",	@a dest: double precision;
    !> @n "TimeStep",	@a dest: double precision.
    !>
    !> @see
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetVarItemsize,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_GetVarType,
    !> @ref RM_BMI_GetVarUnits,
    !> @ref RM_BMI_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> double precision, allocatable, dimension(:) :: bmi_density
    !> character(len=:), allocatable, dimension(:) :: bmi_comps
    !> status = RM_BMI_GetValue(id, "Density", bmi_density)
    !> status = RM_BMI_GetValue("Components", bmi_comps)
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref MpiWorker.

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
    character(20) :: vartype
    integer :: status
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "logical") then
        stop "Variable type error."
    endif
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
    CHARACTER(len=:), allocatable, INTENT(inout) :: dest
    character(20) :: vartype
    integer :: bytes, status
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "character") then
        stop "Variable type error."
    endif
    bytes = RM_BMI_GetVarItemsize(id, var)
    if (len(var) < bytes) then
        if (allocated(dest)) then
            deallocate(dest)
        endif
        allocate(character(len=bytes) :: dest)
    endif
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
    CHARACTER(len=:), allocatable, dimension(:), INTENT(inout) :: dest
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "character,1d") then
        stop "Variable type error."
    endif
    itemsize = RM_BMI_GetVarItemsize(id, var)
    dim = RM_BMI_GetVarNBytes(id, var) / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
        if (dim2 > 0) then
            dim1 = len(dest(1))
        endif
    endif
    if (dim1 .ne. itemsize .or. dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(character(len=itemsize) :: dest(dim))
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double") then
        stop "Variable type error."
    endif
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
    double precision, allocatable, dimension(:), INTENT(inout) :: dest
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double,1d") then
        stop "Variable type error."
    endif
    itemsize = RM_BMI_GetVarItemsize(id, var)
    dim = RM_BMI_GetVarNBytes(id, var) / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
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
    double precision, allocatable, INTENT(inout) :: dest(:,:)
    character(20) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double,2d") then
        stop "Variable type error."
    endif
    varname = Lower(var)
    need_alloc = .true.
    if (varname .eq. "concentrations") then
        dim1 = RM_GetGridCellCount(id)
        dim2 = RM_GetComponentCount(id)
    else if (varname .eq. "selectedoutput") then
        dim1 = RM_GetSelectedOutputRowCount(id)
        dim2 = RM_GetSelectedOutputColumnCount(id)
    else
        stop "Unknown 2d variable"
    endif
    if (allocated(dest)) then
        if ((size(dest,1) .eq. dim1) .and. &
            (size(dest,2) .eq. dim2)) then
            need_alloc = .false.
        else
            deallocate(dest)
        endif
    endif
    if (need_alloc) then
        allocate(dest(dim1, dim2))
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "integer") then
        stop "Variable type error."
    endif
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
    integer, allocatable, INTENT(inout) :: dest(:)
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim, itemsize
    integer :: dim1, dim2
    dim1 = 0
    dim2 = 0
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "integer,1d") then
        stop "Variable type error."
    endif
    itemsize = RM_BMI_GetVarItemsize(id, var)
    dim = RM_BMI_GetVarNBytes(id, var) / itemsize
    if (allocated(dest)) then
        dim2 = size(dest)
    endif
    if (dim2 .ne. dim) then
        if(allocated(dest)) deallocate(dest)
        allocate(dest(dim))
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
    integer, allocatable, INTENT(inout) :: dest(:,:)
    character(20) :: vartype
    character(40) :: varname
    integer :: status
    integer :: dim1, dim2
    logical :: need_alloc
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double,2d") then
        stop "Variable type error."
    endif
    varname = Lower(varname)
    need_alloc = .true.
    stop "Unknown 2d variable"
    RM_BMI_GetValue_i2 = RMF_BMI_GetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION RM_BMI_GetValue_i2

    !> Basic Model Interface method that retrieves the size of an
    !> individual item that can be set or retrived.
    !> Sizes may be sizeof(int), sizeof(double),
    !> or a character length for string variables. Only variables in the list
    !> provided by @ref RM_BMI_GetInputVarNames can be set.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve.
    !> @retval   Size of one element of the variable.
    !>
    !> @see
    !> @ref RM_BMI_GetInputVarNames,
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetValue,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer nbytes, item_size, dim
    !> double precision, allocatable, dimension(:) :: bmi_temperature
    !> nbytes = RM_BMI_GetVarNbytes(id, "Temperature")
    !> item_size = RM_BMI_GetVarItemSize(id, "Temperature");
    !> int dim = nbytes/item_size;
    !> allocate(bmi_temperature(dim))
    !> bmi_temperature = 25.0
    !> status = RM_BMI_SetValue("Temperature", bmi_temperature);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION RM_BMI_GetVarItemsize(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    RM_BMI_GetVarItemsize = RMF_BMI_GetVarItemsize(id, trim(var)//C_NULL_CHAR)
    END FUNCTION RM_BMI_GetVarItemsize
    !> Basic Model Interface method that retrieves the size of an
    !> individual item that can be set or retrived.
    !> Sizes may be sizeof(int), sizeof(double),
    !> or a character length for string variables. Only variables in the list
    !> provided by @ref RM_BMI_GetInputVarNames can be set.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve size.
    !> @retval Size of one element of the variable.
    !>
    !> @see
    !> @ref RM_BMI_GetInputVarNames,
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetValue,
    !> @ref RM_BMI_GetVarNbytes,
    !> @ref RM_BMI_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> integer nbytes, item_size, dim
    !>  double precision, allocatable, dimension(:) :: bmi_temperature
    !> nbytes = RM_BMI_GetVarNbytes(id, "Temperature")
    !> item_size = RM_BMI_GetVarItemSize(id, "Temperature");
    !> int dim = nbytes/item_size;
    !> allocate(bmi_temperature(dim))
    !>  bmi_temperature = 25.0
    !> status = RM_BMI_SetValue("Temperature", bmi_temperature);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
    INTEGER FUNCTION RM_BMI_GetVarNbytes(id, var)
    USE ISO_C_BINDING
    IMPLICIT NONE

    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    RM_BMI_GetVarNbytes = RMF_BMI_GetVarNbytes(id, trim(var)//C_NULL_CHAR)
    END FUNCTION RM_BMI_GetVarNbytes

    !> Basic Model Interface method that retrieves the type of a variable that can be set with
    !> @ref RM_BMI_SetValue or retrieved with @ref RM_BMI_GetValue. Types are "character","double",
    !> "int", or "logical",
    !> followed by the dimension if the variable is a array.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetInputVarNames can be set.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve total bytes.
    !> @param vtype Type of the variable.
    !> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
    !> @see
    !> @ref RM_BMI_GetInputVarNames,
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetVarUnits.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do i = 1, size(inputvars)
    !>     write(*,"(1x, I4, A40)") i, trim(inputvars(i))
    !>     status = RM_RM_BMI_GetVarUnits(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = RM_RM_BMI_GetVarType(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarItemsize(id, inputvars(i))
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarNbytes(id, inputvars(i))
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.

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
    !> Basic Model Interface method that retrieves the units of a
    !> variable that can be set with
    !> @ref RM_BMI_SetValue or retrieved with @ref RM_BMI_GetValue.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetInputVarNames can be set.
    !> Only variables in the list
    !> provided by @ref RM_BMI_GetOutputVarNames can be retrieved.
    !> @param id            The instance @a id returned from @ref RM_Create.
    !> @param var Name of the variable to retrieve total bytes.
    !> @param vtype Units of the variable.
    !> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
    !> @see
    !> @ref RM_BMI_GetInputVarNames,
    !> @ref RM_BMI_GetInputItemCount,
    !> @ref RM_BMI_GetOutputVarNames,
    !> @ref RM_BMI_GetOutputItemCount,
    !> @ref RM_BMI_GetVarType.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> do i = 1, size(inputvars)
    !>     write(*,"(1x, I4, A40)") i, trim(inputvars(i))
    !>     status = RM_RM_BMI_GetVarUnits(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     status = RM_RM_BMI_GetVarType(id, inputvars(i), string)
    !>     write(*,"(5x, A15)") trim(string)
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarItemsize(id, inputvars(i))
    !>     write(*, "(5x, I15)") RM_RM_BMI_GetVarNbytes(id, inputvars(i))
    !> enddo
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root.
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
	!> A YAML file can be used to initialize an instance of PhreeqcRM. Same as
	!> @ref RM_InitializeYAML.
    !> @param id               The instance @a id returned from @ref RM_Create.
    !> @param yaml_name         String containing the YAML file name.
    !> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
    !> @par
    !> The file contains a YAML map of PhreeqcRM methods
    !> and the arguments corresponding to the methods.
    !> Note that the PhreeqcRM methods do not have the "RM_" prefix
    !> and the id argument is not included.
    !> For example,
    !> @par
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> LoadDatabase: phreeqc.dat
    !> RunFile:
    !> workers: true
    !> initial_phreeqc: true
    !> utility: true
    !> chemistry_name: advect.pqi
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par
    !> @ref RM_BMI_Initialize will read the YAML file and execute the specified methods with
    !> the specified arguments. Using YAML
    !> terminology, the argument(s) for a method may be a scalar, a sequence, or a map,
    !> depending if the argument is
    !> a single item, a single vector, or there are multiple arguments.
    !> In the case of a map, the name associated
    !> with each argument (for example "chemistry_name" above) is arbitrary.
    !> The names of the map keys for map
    !> arguments are not used in parsing the YAML file; only the order of
    !> the arguments is important.
    !> @par
    !> The following list gives the PhreeqcRM methods that can be specified in a YAML file
    !> and the arguments that are required. The arguments are described with C++ formats, which
    !> are sufficient to identify which arguments are YAML scalars (single bool/logical,
    !> int, double, string/character argument),
    !> sequences (single vector argument), or maps (multiple arguments).
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> @n CloseFiles(void);
    !> @n CreateMapping(std::vector< int >& grid2chem);
    !> @n DumpModule();
    !> @n FindComponents();
    !> @n InitialPhreeqc2Module(std::vector< int > initial_conditions1);
    !> @n InitialPhreeqc2Module(std::vector< int > initial_conditions1, std::vector< int > initial_conditions2, std::vector< double > fraction1);
    !> @n InitialPhreeqcCell2Module(int n, std::vector< int > cell_numbers);
    !> @n LoadDatabase(std::string database);
    !> @n OpenFiles(void);
    !> @n OutputMessage(std::string str);
    !> @n RunCells(void);
    !> @n RunFile(bool workers, bool initial_phreeqc, bool utility, std::string chemistry_name);
    !> @n RunString(bool workers, bool initial_phreeqc, bool utility, std::string input_string);
    !> @n ScreenMessage(std::string str);
    !> @n SetComponentH2O(bool tf);
    !> @n SetConcentrations(std::vector< double > c);
    !> @n SetCurrentSelectedOutputUserNumber(int n_user);
    !> @n SetDensity(std::vector< double > density);
    !> @n SetDumpFileName(std::string dump_name);
    !> @n SetErrorHandlerMode(int mode);
    !> @n SetErrorOn(bool tf);
    !> @n SetFilePrefix(std::string prefix);
    !> @n SetGasCompMoles(std::vector< double > gas_moles);
    !> @n SetGasPhaseVolume(std::vector< double > gas_volume);
    !> @n SetPartitionUZSolids(bool tf);
    !> @n SetPorosity(std::vector< double > por);
    !> @n SetPressure(std::vector< double > p);
    !> @n SetPrintChemistryMask(std::vector< int > cell_mask);
    !> @n SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
    !> @n SetRebalanceByCell(bool tf);
    !> @n SetRebalanceFraction(double f);
    !> @n SetRepresentativeVolume(std::vector< double > rv);
    !> @n SetSaturation(std::vector< double > sat);
    !> @n SetScreenOn(bool tf);
    !> @n SetSelectedOutputOn(bool tf);
    !> @n SetSpeciesSaveOn(bool save_on);
    !> @n SetTemperature(std::vector< double > t);
    !> @n SetTime(double time);
    !> @n SetTimeConversion(double conv_factor);
    !> @n SetTimeStep(double time_step);
    !> @n SetUnitsExchange(int option);
    !> @n SetUnitsGasPhase(int option);
    !> @n SetUnitsKinetics(int option);
    !> @n SetUnitsPPassemblage(int option);
    !> @n SetUnitsSolution(int option);
    !> @n SetUnitsSSassemblage(int option);
    !> @n SetUnitsSurface(int option);
    !> @n SpeciesConcentrations2Module(std::vector< double > species_conc);
    !> @n StateSave(int istate);
    !> @n StateApply(int istate);
    !> @n StateDelete(int istate);
    !> @n UseSolutionDensityVolume(bool tf);
    !> @n WarningMessage(std::string warnstr);
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !>
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> id = RM_Create(nxyz, MPI_COMM_WORLD)
    !> status = RM_InitializeYAML(id, "myfile.yaml")
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
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
    LOGICAL(KIND=C_INT), INTENT(in) :: dest
    END FUNCTION RMF_BMI_SetValue
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    LOGICAL(kind=4), INTENT(in) :: dest
    CHARACTER(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "logical") then
        stop "Variable type error."
    endif
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "character") then
        stop "Variable type error."
    endif
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "integer") then
        stop "Variable type error."
    endif
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "integer,1d") then
        stop "Variable type error."
    endif
    dim = RM_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var)
    if (dim .ne. size(dest)) then
        stop "Variable dimension error"
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "integer,2d") then
        stop "Variable type error."
    endif
    dim = RM_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var)
    if (dim .ne. size(dest,1)*size(dest,2)) then
        stop "Variable dimension error"
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double") then
        stop "Variable type error."
    endif
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double,1d") then
        stop "Variable type error."
    endif
    dim = RM_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var)
    if (dim .ne. size(dest)) then
        stop "Variable dimension error"
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
    character(20) :: vartype
    integer :: bytes, nbytes, status, dim
    status = RM_BMI_GetVarType(id, var, vartype)
    if (vartype .ne. "double,2d") then
        stop "Variable type error."
    endif
    dim = RM_BMI_GetVarNBytes(id, var) / RMF_BMI_GetVarItemsize(id, var)
    if (dim .ne. size(dest,1)*size(dest,2)) then
        stop "Variable dimension error"
    endif
    RM_BMI_SetValue_d2 = RMF_BMI_SetValue(id, trim(var)//C_NULL_CHAR, dest(1,1))
    return
    END FUNCTION RM_BMI_SetValue_d2
	!> Runs a reaction step for all of the cells in the reaction module.
	!> Same as @ref RunCells.
    !> Normally, tranport concentrations are transferred to the reaction cells 
	!> (@ref RM_BMI_SetValue "Concentrations" before
    !> reaction calculations are run. The length of time over which kinetic 
	!> reactions are integrated is set
    !> by @ref RM_BMI_SetTimeStep. Other properties that may need to be updated 
	!> as a result of the transport
    !> calculations include porosity (@ref RM_BMI_SetValue "Porosity"), 
	!> pressure (@ref RM_BMI_SetValue "Pressure"),
	!> saturation (@ref RM_BMI_SetValue "Saturation"),
    !> temperature (@ref RM_BMI_SetValue "Temperature").
    !> @param id               The instance @a id returned from @ref RM_BMI_Create.
    !> @retval IRM_BMI_RESULT      0 is success, negative is failure (See @ref RM_BMI_DecodeError).
    !> @see
    !> @ref RM_BMI_SetValue.
    !> @par Fortran Example:
    !> @htmlonly
    !> <CODE>
    !> <PRE>
    !> status = RM_BMI_SetValue(id, "Porosity", por)                ! If pore volume changes
    !> status = RM_BMI_SetValue(id, "Saturation", sat)              ! If saturation changes
    !> status = RM_BMI_SetValue(id, "Temperature", temperature)     ! If temperature changes
    !> status = RM_BMI_SetValue(id, "Pressure", pressure)           ! If pressure changes
    !> status = RM_BMI_SetValue(id, "Concentrations", c)            ! Transported concentrations
    !> status = RM_BMI_SetValue(id, "TimeStep", time_step)          ! Time step for kinetic reactions
    !> status = RM_BMI_Update(id)
    !> status = RM_BMI_GetValue(id, "Concentrations", c)            ! Concentrations after reaction
    !> status = RM_BMI_GetValue(id, "Density", density)             ! Density after reaction
    !> status = RM_BMI_GetValue(id, "SolutionVolume", volume)       ! Solution volume after reaction
    !> </PRE>
    !> </CODE>
    !> @endhtmlonly
    !> @par MPI:
    !> Called by root, workers must be in the loop of @ref RM_MpiWorker.
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

    FUNCTION Lower(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i

    DO i = 1,LEN(s1)
        ch = s1(i:i)
        IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
        s2(i:i) = ch
    END DO
    END FUNCTION Lower
    END MODULE PhreeqcRM