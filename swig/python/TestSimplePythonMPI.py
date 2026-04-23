import phreeqcrm
import sys
import numpy as np
# installed msmpisetup.exe 10.1.3 (10.1.12498.52)
# sudo apt install openmpi-bin openmpi-common libopenmpi-dev  # version 4.1.6-5ubuntu1

from mpi4py import MPI
import pytest

def test_comm():
    # assumes mpiexec -n 3 python TestSimplePythonMPI.py

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()
    print(f"Hello, MPI! I am process {rank} of {size} on {name}.", flush=True)
    assert size == 3

    # Prepare array to gather all ranks
    all_ranks = np.empty(size, dtype=int)
    
    # Gather ranks from all processes
    MPI.COMM_WORLD.Allgather(np.array([rank], dtype=int), all_ranks)
    
    # Assert that ranks are unique and match [0, 1, 2]
    assert sorted(all_ranks) == [0, 1, 2]

    # Create PhreeqcRM instance and test MPI functions
    nxyz = 20
    phreeqc_rm = phreeqcrm.PhreeqcRM(nxyz, MPI.COMM_WORLD)
    tasks = phreeqc_rm.GetMpiTasks()
    mpi_myself = phreeqc_rm.GetMpiMyself()
    print(f"Hello, PhreeqcRM! I am process {mpi_myself} of {tasks}.", flush=True)   
    assert tasks == size
    assert mpi_myself == rank

    # Prepare array to gather all phreeqcrm ranks
    all_phreeqc_rm_ranks = np.empty(tasks, dtype=int)
    
    # Gather mpi_myselfs from all processes
    MPI.COMM_WORLD.Allgather(np.array([mpi_myself], dtype=int), all_phreeqc_rm_ranks)
    
    # Assert that ranks are unique and match [0, 1, 2]
    assert sorted(all_phreeqc_rm_ranks) == [0, 1, 2]

if __name__ == "__main__":
    print("run using: mpiexec -n 3 python TestSimplePythonMPI.py")
    test_comm()
