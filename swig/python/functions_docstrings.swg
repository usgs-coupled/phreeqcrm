%ignore basic_callback_shim;
%ignore mpi_worker_callback_shim;


%feature("docstring") has_mpi
"""
Determine whether phreeqcrm was compiled with MPI support.

This function reports whether distributed-memory parallel
execution via MPI is available in the current phreeqcrm build.

Returns:
    bool
        True if MPI support is enabled, False otherwise.
"""


%feature("docstring") has_openmp
"""
Determine whether phreeqcrm was compiled with OpenMP support.

This function reports whether parallel execution via OpenMP
is available in the current phreeqcrm build.

Returns:
    bool
        True if OpenMP support is enabled, False otherwise.
"""
