
	! ====================================================
    ! Functions not implemented
	! ====================================================
    integer FUNCTION bmif_get_var_location(id, var, location)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
	CHARACTER(len=*), INTENT(in) :: var
	CHARACTER(len=*), INTENT(inout) :: location
    bmif_get_var_location = BMI_FAILURE
    END FUNCTION bmif_get_var_location

    !> \overload
    INTEGER FUNCTION bmif_get_value_float(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=4), INTENT(inout) :: dest
    bmif_get_value_float = BMI_FAILURE
    return
    END FUNCTION bmif_get_value_float   
    
    !> \overload
    INTEGER FUNCTION bmif_get_value_ptr_float(id, var, dest)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=4), pointer, INTENT(inout) :: dest
    bmif_get_value_ptr_float = BMI_FAILURE
    return 
    END FUNCTION bmif_get_value_ptr_float
    
    !> \overload
    INTEGER FUNCTION get_value_at_indices_double(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT(inout) :: dest(:)
    integer, INTENT(inout) :: inds(:)
    get_value_at_indices_double = BMI_FAILURE
    return 
    END FUNCTION get_value_at_indices_double
    
    !> \overload
    INTEGER FUNCTION get_value_at_indices_float(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=4), INTENT(inout) :: dest(:)
    integer, INTENT(in) :: inds(:)
    get_value_at_indices_float = BMI_FAILURE
    return 
    END FUNCTION get_value_at_indices_float
    
    !> \overload
    INTEGER FUNCTION get_value_at_indices_int(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:)
    integer, INTENT(in) :: inds(:)
    get_value_at_indices_int = BMI_FAILURE
    return 
    END FUNCTION get_value_at_indices_int
    
    !> \overload
    INTEGER FUNCTION bmif_set_value_float(id, var, src)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=4), INTENT(inout) :: src
    bmif_set_value_float = BMI_FAILURE
    return
    END FUNCTION bmif_set_value_float
    
    !> \overload
    INTEGER FUNCTION set_value_at_indices_double(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=8), INTENT(inout) :: dest(:)
    integer, INTENT(inout) :: inds(:)
    set_value_at_indices_double = BMI_FAILURE
    return 
    END FUNCTION set_value_at_indices_double
    
    !> \overload
    INTEGER FUNCTION set_value_at_indices_float(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    real(kind=4), INTENT(inout) :: dest(:)
    integer, INTENT(in) :: inds(:)
    set_value_at_indices_float = BMI_FAILURE
    return 
    END FUNCTION set_value_at_indices_float
    
    !> \overload
    INTEGER FUNCTION set_value_at_indices_int(id, var, dest, inds)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: var
    integer, INTENT(inout) :: dest(:)
    integer, INTENT(in) :: inds(:)
    set_value_at_indices_int = BMI_FAILURE
    return 
    END FUNCTION set_value_at_indices_int    
    
    
    ! Get the dimensions of the computational grid.
    integer function bmif_get_grid_shape(grid, shape) result(bmi_status)
    integer, intent(in) :: grid
    integer, dimension(:), intent(inout) :: shape
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_shape

    ! Get distance between nodes of the computational grid.
    integer function bmif_get_grid_spacing(grid, spacing) result(bmi_status)
    integer, intent(in) :: grid
    double precision, dimension(:), intent(inout) :: spacing
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_spacing

    ! Get coordinates of the origin of the computational grid.
    integer function bmif_get_grid_origin(grid, origin) result(bmi_status)
    integer, intent(in) :: grid
    double precision, dimension(:), intent(inout) :: origin
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_origin

    ! Get the x-coordinates of the nodes of a computational grid.
    integer function bmif_get_grid_x(grid, x) result(bmi_status)
    integer, intent(in) :: grid
    double precision, dimension(:), intent(inout) :: x
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_x

    ! Get the y-coordinates of the nodes of a computational grid.
    integer function bmif_get_grid_y(grid, y) result(bmi_status)
    integer, intent(in) :: grid
    double precision, dimension(:), intent(inout) :: y
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_y

    ! Get the z-coordinates of the nodes of a computational grid.
    integer function bmif_get_grid_z(grid, z) result(bmi_status)
    integer, intent(in) :: grid
    double precision, dimension(:), intent(inout) :: z
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_z

    ! Get the number of nodes in an unstructured grid.
    integer function bmif_get_grid_node_count(grid, count) result(bmi_status)
    integer, intent(in) :: grid
    integer, intent(inout) :: count
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_node_count

    ! Get the number of edges in an unstructured grid.
    integer function bmif_get_grid_edge_count(grid, count) result(bmi_status)
    integer, intent(in) :: grid
    integer, intent(inout) :: count
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_edge_count

    ! Get the number of faces in an unstructured grid.
    integer function bmif_get_grid_face_count(grid, count) result(bmi_status)
    integer, intent(in) :: grid
    integer, intent(inout) :: count
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_face_count

    ! Get the edge-node connectivity.
    integer function bmif_get_grid_edge_nodes(grid, edge_nodes) &
        result(bmi_status)
    integer, intent(in) :: grid
    integer, dimension(:), intent(inout) :: edge_nodes
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_edge_nodes

    ! Get the face-edge connectivity.
    integer function bmif_get_grid_face_edges(grid, face_edges) &
        result(bmi_status)
    integer, intent(in) :: grid
    integer, dimension(:), intent(inout) :: face_edges
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_face_edges

    ! Get the face-node connectivity.
    integer function bmif_get_grid_face_nodes(grid, face_nodes) &
        result(bmi_status)
    integer, intent(in) :: grid
    integer, dimension(:), intent(inout) :: face_nodes
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_face_nodes

    ! Get the number of nodes for each face.
    integer function bmif_get_grid_nodes_per_face(grid, nodes_per_face) &
        result(bmi_status)
    integer, intent(in) :: grid
    integer, dimension(:), intent(inout) :: nodes_per_face
    bmi_status = BMI_FAILURE
    return
    end function bmif_get_grid_nodes_per_face