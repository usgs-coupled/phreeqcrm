# #!/usr/bin/env python3

import pytest
import os
import sys

# Add the project directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def phreeqcrm_module():
    """Import the callback module."""
    try:
        import phreeqcrm
        return phreeqcrm
    except ImportError as e:
        pytest.skip(f"Cannot import phreeqcrm module: {e}")

class TestModuleImport:
    """Tests for module import functionality."""
    
    def test_import_module(self, phreeqcrm_module):
        """Test that the module can be imported."""
        assert phreeqcrm_module is not None

    def test_has_functions(self, phreeqcrm_module):
        """Test that the expected methods exist."""
        # Check has_mpi
        assert hasattr(phreeqcrm_module, 'has_mpi')
        assert callable(phreeqcrm_module.has_mpi)
        assert isinstance(phreeqcrm_module.has_mpi(), bool)  # Should return a boolean
        assert phreeqcrm_module.has_mpi.__doc__ is not None  # Should have a docstring
        # Check has_openmp
        assert hasattr(phreeqcrm_module, 'has_openmp')
        assert callable(phreeqcrm_module.has_openmp)
        assert isinstance(phreeqcrm_module.has_openmp(), bool)  # Should return a boolean
        assert phreeqcrm_module.has_openmp.__doc__ is not None  # Should have a docstring

    def test_shims_are_ignored(self, phreeqcrm_module):
        """Test that shim functions are not present in the module."""
        # Check basic_callback_shim
        assert not hasattr(phreeqcrm_module, 'basic_callback_shim')
        # Check mpi_worker_callback_shim
        assert not hasattr(phreeqcrm_module, 'mpi_worker_callback_shim')
