# #!/usr/bin/env python3
# """
# Comprehensive test suite for the callback module using pytest.

# Tests cover:
# - Module import
# - Callback registration
# - Callback execution with various argument types
# - Return value handling
# - Multiple callback scenarios
# - Edge cases
# """

import sys
import os
import pytest
import numpy as np

from constants import FilePaths

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
    
    def test_import_callback_module(self, phreeqcrm_module):
        """Test that the callback module can be imported."""
        assert phreeqcrm_module is not None

    def test_classes_exist(self, phreeqcrm_module):
        """Test that the expected classes are present in the module."""
        assert hasattr(phreeqcrm_module, 'PhreeqcRM')
        assert hasattr(phreeqcrm_module, 'BMIPhreeqcRM')

    def test_classes_have_set_basic_callback(self, phreeqcrm_module):
        """Test that the expected classes have set_basic_callback methods."""
        # Check PhreeqcRM class
        assert hasattr(phreeqcrm_module.PhreeqcRM, 'set_basic_callback')
        assert callable(phreeqcrm_module.PhreeqcRM.set_basic_callback)
        assert phreeqcrm_module.PhreeqcRM.set_basic_callback.__doc__ is not None  # Should have a docstring

        # Check BMIPhreeqcRM class
        assert hasattr(phreeqcrm_module.BMIPhreeqcRM, 'set_basic_callback')
        assert callable(phreeqcrm_module.BMIPhreeqcRM.set_basic_callback)
        assert phreeqcrm_module.BMIPhreeqcRM.set_basic_callback.__doc__ is not None  # Should have a docstring

    def test_classes_have_set_mpi_worker_callback(self, phreeqcrm_module):
        """Test that the expected classes have set_mpi_worker_callback methods."""
        # Check PhreeqcRM class
        assert hasattr(phreeqcrm_module.PhreeqcRM, 'set_mpi_worker_callback')
        assert callable(phreeqcrm_module.PhreeqcRM.set_mpi_worker_callback)
        assert phreeqcrm_module.PhreeqcRM.set_mpi_worker_callback.__doc__ is not None  # Should have a docstring

        # Check BMIPhreeqcRM class
        assert hasattr(phreeqcrm_module.BMIPhreeqcRM, 'set_mpi_worker_callback')
        assert callable(phreeqcrm_module.BMIPhreeqcRM.set_mpi_worker_callback)
        assert phreeqcrm_module.BMIPhreeqcRM.set_mpi_worker_callback.__doc__ is not None  # Should have a docstring


class TestBasicCallbackRegistration:
    """Tests for callback registration."""

    def test_set_callback_throws_before_initialize(self, phreeqcrm_module):
        """Test setting a callback before initialization raises an error."""
        def dummy_callback(val1, val2, message, cookie):
            return 0.0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        # with pytest.raises(RuntimeError, match="must call initialize first"):
        #     model.set_basic_callback(dummy_callback)

        # Note: For compatibility with MPI (which can't call initialize for the
        #  mpi workers), set_basic_callback can be called before initialize. 
        #  The callback will be set on workers towards the end of Construct().
        model.set_basic_callback(dummy_callback)

    def test_set_callback_with_valid_function(self, phreeqcrm_module):
        """Test setting a callback with a valid function."""
        def dummy_callback(val1, val2, message, cookie):
            return 0.0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_basic_callback(dummy_callback)

    def test_set_callback_with_valid_function_and_cookie(self, phreeqcrm_module):
        """Test setting a callback with a valid function and cookie."""
        def dummy_callback(val1, val2, message, cookie):
            return 0.0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_basic_callback(dummy_callback, model)  # Using the model instance as a cookie

    def test_set_callback_with_lambda(self, phreeqcrm_module):
        """Test setting a callback with a lambda function."""

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_basic_callback(lambda v1, v2, m, cookie: 0.0)

    def test_set_callback_with_lambda_and_cookie(self, phreeqcrm_module):
        """Test setting a callback with a lambda function and cookie."""

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_basic_callback(lambda v1, v2, m, cookie: 0.0, model)

    def test_set_callback_with_non_callable_raises_error(self, phreeqcrm_module):
        """Test that setting a non-callable raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_basic_callback("not a function")

    def test_set_callback_with_none_raises_error(self, phreeqcrm_module):
        """Test that setting None raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_basic_callback(None)

    def test_set_callback_with_integer_raises_error(self, phreeqcrm_module):
        """Test that setting an integer raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_basic_callback(42)


class TestMpiWorkerCallbackRegistration:
    """Tests for callback registration."""

    def test_set_callback_throws_before_initialize(self, phreeqcrm_module):
        """Test setting a callback before initialization raises an error."""
        def dummy_callback(val, cookie):
            return 0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        # with pytest.raises(RuntimeError, match="must call initialize first"):
        #     model.set_basic_callback(dummy_callback)

        # Note: For compatibility with MPI (which can't call initialize for the
        #  mpi workers), set_basic_callback can be called before initialize. 
        #  The callback will be set on workers towards the end of Construct().
        model.set_mpi_worker_callback(dummy_callback)

    def test_set_mpi_worker_callback_with_valid_function(self, phreeqcrm_module):
        """Test setting a callback with a valid function."""
        def dummy_callback(val, cookie):
            return 0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_mpi_worker_callback(dummy_callback)

    def test_set_callback_with_valid_function_and_cookie(self, phreeqcrm_module):
        """Test setting a callback with a valid function and cookie."""
        def dummy_callback(val, cookie):
            return 0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_mpi_worker_callback(dummy_callback, model)  # Using the model instance as a cookie

    def test_set_callback_with_lambda(self, phreeqcrm_module):
        """Test setting a callback with a lambda function."""

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_mpi_worker_callback(lambda v, cookie: 0)

    def test_set_callback_with_lambda_and_cookie(self, phreeqcrm_module):
        """Test setting a callback with a lambda function and cookie."""

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()

        # Should not raise an exception
        model.set_mpi_worker_callback(lambda v, cookie: 0, model)  # Using the model instance as a cookie

    def test_set_callback_with_non_callable_raises_error(self, phreeqcrm_module):
        """Test that setting a non-callable raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_mpi_worker_callback("not a function")

    def test_set_callback_with_none_raises_error(self, phreeqcrm_module):
        """Test that setting None raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_mpi_worker_callback(None)

    def test_set_callback_with_integer_raises_error(self, phreeqcrm_module):
        """Test that setting an integer raises TypeError."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        with pytest.raises(TypeError, match="Callback must be callable"):
            model.set_mpi_worker_callback(42)


class TestBasicCallbackExecution:
    """Tests for callback execution."""

    def test_execute_with_simple_return(self, phreeqcrm_module):
        """Test executing a callback that returns a simple double."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: 3.14)
        result = model._execute_basic_callback(1.5, 2.5, "test")
        assert result == pytest.approx(3.14)

    def test_execute_with_sum_of_arguments(self, phreeqcrm_module):
        """Test executing a callback that sums the double arguments."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(10.0, 20.0, "sum")
        assert result == pytest.approx(30.0)
    
    def test_execute_with_multiplication(self, phreeqcrm_module):
        """Test executing a callback that multiplies arguments."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 * v2)
        result = model._execute_basic_callback(5.0, 4.0, "multiply")
        assert result == pytest.approx(20.0)

    def test_execute_with_message_based_computation(self, phreeqcrm_module):
        """Test executing a callback that uses message length in computation."""
        def callback(v1, v2, msg, cookie):
            return v1 + v2 + len(msg)
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(callback)
        result = model._execute_basic_callback(1.0, 2.0, "hello")
        # 1.0 + 2.0 + len("hello") = 1.0 + 2.0 + 5 = 8.0
        assert result == pytest.approx(8.0)

    def test_execute_with_zero_values(self, phreeqcrm_module):
        """Test executing a callback with zero values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(0.0, 0.0, "zeros")
        assert result == pytest.approx(0.0)

    def test_execute_with_negative_values(self, phreeqcrm_module):
        """Test executing a callback with negative values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(-5.0, 3.0, "negative")
        assert result == pytest.approx(-2.0)
    
    def test_execute_with_large_values(self, phreeqcrm_module):
        """Test executing a callback with large double values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(1e10, 2e10, "large")
        assert result == pytest.approx(3e10)
    
    def test_execute_with_small_values(self, phreeqcrm_module):
        """Test executing a callback with very small double values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(1e-10, 2e-10, "small")
        assert result == pytest.approx(3e-10, abs=1e-20)

    def test_execute_with_empty_message(self, phreeqcrm_module):
        """Test executing a callback with an empty message string."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(5.0, 3.0, "")
        assert result == pytest.approx(8.0)
    
    def test_execute_with_long_message(self, phreeqcrm_module):
        """Test executing a callback with a long message string."""
        long_msg = "a" * 1000
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2 + len(m))
        result = model._execute_basic_callback(1.0, 2.0, long_msg)
        assert result == pytest.approx(1003.0)
    
    def test_execute_with_special_characters_in_message(self, phreeqcrm_module):
        """Test executing a callback with special characters in message."""
        special_msg = "!@#$%^&*()"
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(1.0, 2.0, special_msg)
        assert result == pytest.approx(3.0)


class TestMpiWorkerCallbackExecution:
    """Tests for callback execution."""

    def test_set_callback_with_dummy_callback(self, phreeqcrm_module):
        """Test setting a callback before initialization raises an error."""
        def dummy_callback(val, cookie):
            return 0
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.set_mpi_worker_callback(dummy_callback)
        result = model._execute_mpi_worker_callback(301)
        assert result == 0

    def test_execute_with_simple_return(self, phreeqcrm_module):
        """Test executing a callback that returns a simple int."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: 303)
        result = model._execute_mpi_worker_callback(301)
        assert result == 303

    def test_execute_with_multiplication(self, phreeqcrm_module):
        """Test executing a callback that multiplies arguments."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: 3 * v)
        result = model._execute_mpi_worker_callback(101)
        assert result == 303

    def test_execute_with_zero_values(self, phreeqcrm_module):
        """Test executing a callback with zero values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: 555)
        result = model._execute_mpi_worker_callback(0)
        assert result == 555

    def test_execute_with_negative_values(self, phreeqcrm_module):
        """Test executing a callback with negative values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: v)
        result = model._execute_mpi_worker_callback(-5)
        assert result == -5
    
    def test_execute_with_large_values(self, phreeqcrm_module):
        """Test executing a callback with large integer value."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: v)
        result = model._execute_mpi_worker_callback(2147483647)
        assert result == 2147483647

class TestBasicCallbackStateManagement:
    """Tests for managing callback state."""
    
    def test_multiple_set_basic_callback_calls(self, phreeqcrm_module):
        """Test that multiple set_basic_callback calls override the previous callback."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: 1.0)
        result1 = model._execute_basic_callback(1.0, 1.0, "first")
        assert result1 == pytest.approx(1.0)
        
        model.set_basic_callback(lambda v1, v2, m, cookie: 2.0)
        result2 = model._execute_basic_callback(1.0, 1.0, "second")
        assert result2 == pytest.approx(2.0)

    def test_callback_with_state_via_closure(self, phreeqcrm_module):
        """Test that callbacks can maintain state via closures."""
        call_count = [0]  # Use a list to allow modification in nested function
        
        def stateful_callback(v1, v2, m, cookie):
            call_count[0] += 1
            return v1 + v2 + call_count[0]
        
        if phreeqcrm_module.has_mpi():
            from mpi4py import MPI
            model = phreeqcrm_module.BMIPhreeqcRM(10, MPI.COMM_WORLD)
            model.initialize()
            assert MPI.COMM_WORLD.Get_size() == 1
        elif phreeqcrm_module.has_openmp():
            model = phreeqcrm_module.BMIPhreeqcRM(10, 1)
            model.initialize()
            assert model.GetThreadCount() == 1
        else:
            # serial case - default to 1 thread
            model = phreeqcrm_module.BMIPhreeqcRM(10, 1)
            model.initialize()
            assert model.GetThreadCount() == 1

        assert model is not None
        model.initialize()
        model.set_basic_callback(stateful_callback)
        
        result1 = model._execute_basic_callback(1.0, 2.0, "first")
        assert result1 == pytest.approx(6.0)  # 1 + 2 + threads+2
        
        result2 = model._execute_basic_callback(1.0, 2.0, "second")
        assert result2 == pytest.approx(9.0)  # 1 + 2 + 2*(threads+2)

    def test_callback_still_gets_called_after_initialize(self, phreeqcrm_module):
        """Test that callbacks set before initialization are still called."""

        def basic_callback(v1, v2, m, cookie):
            print("basic_callback called with v1:", v1, "v2:", v2, "m:", m, "cookie:", cookie)
            return v1 + v2

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.set_basic_callback(basic_callback)  # Set before initialize
        model.initialize(FilePaths.MINIMUM_YAML)  # Initialize after setting callback

        result1 = model._execute_basic_callback(1.0, 2.0, "first")
        assert result1 == pytest.approx(3.0)  # 1 + 2

    def test_callback_still_gets_called_after_initialize_alt(self, phreeqcrm_module):
        """Test that callbacks set before initialization are still called."""

        def basic_callback(v1, v2, m, cookie):
            print("basic_callback called with v1:", v1, "v2:", v2, "m:", m, "cookie:", cookie)
            return v1 + v2

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.set_basic_callback(basic_callback)  # Set before initialize
        model.initialize()  # Initialize after setting callback
        model.LoadDatabase("minimum.dat")

        result1 = model._execute_basic_callback(1.0, 2.0, "first")
        assert result1 == pytest.approx(3.0)  # 1 + 2

class TestMpiWorkerCallbackStateManagement:
    """Tests for managing callback state."""
    
    def test_multiple_set_basic_callback_calls(self, phreeqcrm_module):
        """Test that multiple set_mpi_worker_callback calls override the previous callback."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(lambda v, cookie: 1)
        result1 = model._execute_mpi_worker_callback(1)
        assert result1 == 1
        
        model.set_mpi_worker_callback(lambda v, cookie: 2)
        result2 = model._execute_mpi_worker_callback(1)
        assert result2 == 2

    def test_callback_with_state_via_closure(self, phreeqcrm_module):
        """Test that callbacks can maintain state via closures."""
        call_count = [0]  # Use a list to allow modification in nested function

        def stateful_callback(v, cookie):
            call_count[0] += 1
            return v + call_count[0]

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_mpi_worker_callback(stateful_callback)

        result1 = model._execute_mpi_worker_callback(1)
        assert result1 == 2  # 1 + 1

        result2 = model._execute_mpi_worker_callback(1)
        assert result2 == 3  # 1 + 2

    # def test_callback_still_gets_called_after_initialize(self, phreeqcrm_module):
    #     """Test that callbacks set before initialization are still called."""
    #     call_count = [0]

    #     def stateful_callback(v, cookie):
    #         call_count[0] += 1
    #         return v + call_count[0]

    #     model = phreeqcrm_module.BMIPhreeqcRM()
    #     model.set_mpi_worker_callback(stateful_callback)  # Set before initialize
    #     model.initialize()  # Initialize after setting callback

    #     result1 = model._execute_mpi_worker_callback(1)
    #     assert result1 == 2  # 1 + 1

    #     result2 = model._execute_mpi_worker_callback(1)
    #     assert result2 == 3  # 1 + 2

    def test_callback_still_gets_called_after_initialize(self, phreeqcrm_module):
        """Test that callbacks set before initialization are still called."""

        def mpi_worker_callback(v, cookie):
            print("basic_callback called with v:", v, "cookie:", cookie)
            return v

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.set_mpi_worker_callback(mpi_worker_callback)  # Set before initialize
        model.initialize(FilePaths.MINIMUM_YAML)  # Initialize after setting callback

        result1 = model._execute_mpi_worker_callback(301)
        assert result1 == 301

    def test_callback_still_gets_called_after_initialize_alt(self, phreeqcrm_module):
        """Test that callbacks set before initialization are still called."""

        def mpi_worker_callback(v, cookie):
            print("basic_callback called with v:", v, "cookie:", cookie)
            return v

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.set_mpi_worker_callback(mpi_worker_callback)  # Set before initialize
        model.initialize()  # Initialize after setting callback
        model.LoadDatabase("minimum.dat")

        result1 = model._execute_mpi_worker_callback(707)
        assert result1 == 707

class TestReturnTypes:
    """Tests for return type handling."""
    
    def test_callback_returns_float_value(self, phreeqcrm_module):
        """Test that callback return value is properly converted to float."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: 3.14159)
        result = model._execute_basic_callback(0.0, 0.0, "pi")
        assert isinstance(result, float)
        assert result == pytest.approx(3.14159)
    
    def test_callback_returns_int_converted_to_float(self, phreeqcrm_module):
        """Test that integer return values are converted to float."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: 42)
        result = model._execute_basic_callback(0.0, 0.0, "int")
        assert isinstance(result, (float, int))
        assert result == pytest.approx(42.0)
    
    def test_callback_returns_zero(self, phreeqcrm_module):
        """Test that callback can return zero."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: 0)
        result = model._execute_basic_callback(1.0, 2.0, "zero")
        assert result == pytest.approx(0.0)
    
    def test_callback_returns_negative(self, phreeqcrm_module):
        """Test that callback can return negative values."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: -99.99)
        result = model._execute_basic_callback(0.0, 0.0, "negative")
        assert result == pytest.approx(-99.99)


class TestComplexCallbacks:
    """Tests for more complex callback scenarios."""
    
    def test_callback_with_conditional_logic(self, phreeqcrm_module):
        """Test callback with conditional logic."""
        def conditional_callback(v1, v2, msg, cookie):
            if v1 > v2:
                return v1 - v2
            else:
                return v2 - v1
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(conditional_callback)
        
        result1 = model._execute_basic_callback(10.0, 5.0, "greater")
        assert result1 == pytest.approx(5.0)
        
        result2 = model._execute_basic_callback(3.0, 8.0, "less")
        assert result2 == pytest.approx(5.0)
    
    def test_callback_with_string_checking(self, phreeqcrm_module):
        """Test callback that checks message string."""
        def string_aware_callback(v1, v2, msg, cookie):
            if msg == "double":
                return (v1 + v2) * 2
            else:
                return v1 + v2

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(string_aware_callback)

        result1 = model._execute_basic_callback(5.0, 3.0, "double")
        assert result1 == pytest.approx(16.0)

        result2 = model._execute_basic_callback(5.0, 3.0, "normal")
        assert result2 == pytest.approx(8.0)

    def test_callback_with_exception_handling(self, phreeqcrm_module):
        """Test callback that handles edge cases gracefully."""
        def safe_callback(v1, v2, msg, cookie):
            try:
                if len(msg) == 0:
                    return v1 + v2
                else:
                    return (v1 + v2) / len(msg)
            except ZeroDivisionError:
                return 0.0

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(safe_callback)

        result1 = model._execute_basic_callback(10.0, 5.0, "")
        assert result1 == pytest.approx(15.0)

        result2 = model._execute_basic_callback(10.0, 5.0, "ab")
        assert result2 == pytest.approx(7.5)

class TestUnicodeHandling:
    """Tests for unicode string handling in callbacks."""
    
    def test_callback_with_basic_unicode(self, phreeqcrm_module):
        """Test callback with basic unicode characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        result = model._execute_basic_callback(1.0, 2.0, "café")
        assert result == pytest.approx(3.0)
    
    def test_callback_with_emoji(self, phreeqcrm_module):
        """Test callback with emoji characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "😀🎉🚀")
        assert result == pytest.approx(3.0)
    
    def test_callback_with_chinese_characters(self, phreeqcrm_module):
        """Test callback with Chinese characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "你好世界")
        assert result == pytest.approx(4.0)
    
    def test_callback_with_japanese_characters(self, phreeqcrm_module):
        """Test callback with Japanese characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "こんにちは")
        assert result == pytest.approx(5.0)
    
    def test_callback_with_arabic_characters(self, phreeqcrm_module):
        """Test callback with Arabic characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "مرحبا")
        assert result == pytest.approx(5.0)
    
    def test_callback_with_greek_characters(self, phreeqcrm_module):
        """Test callback with Greek characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "Ελληνικά")
        assert result == pytest.approx(8.0)
    
    def test_callback_with_cyrillic_characters(self, phreeqcrm_module):
        """Test callback with Cyrillic (Russian) characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "Привет")
        assert result == pytest.approx(6.0)
    
    def test_callback_with_mixed_unicode_and_ascii(self, phreeqcrm_module):
        """Test callback with mixed ASCII and unicode characters."""
        def mixed_callback(v1, v2, msg, cookie):
            ascii_count = sum(1 for c in msg if ord(c) < 128)
            return ascii_count
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(mixed_callback)
        result = model._execute_basic_callback(0.0, 0.0, "Hello世界!")
        # "Hello!" = 6 ASCII characters
        assert result == pytest.approx(6.0)
    
    def test_callback_with_unicode_diacritics(self, phreeqcrm_module):
        """Test callback with accented characters and diacritics."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        result = model._execute_basic_callback(0.0, 0.0, "Français Español Português")
        assert result == pytest.approx(26.0)
    
    def test_callback_with_combining_characters(self, phreeqcrm_module):
        """Test callback with combining diacritical marks."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        # Using a combining character sequence
        combining_string = "é"  # e with combining acute accent
        result = model._execute_basic_callback(0.0, 0.0, combining_string)
        assert result >= 1.0  # Should have at least the base character
    
    def test_callback_with_zero_width_characters(self, phreeqcrm_module):
        """Test callback with zero-width characters."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        # Zero-width space (U+200B)
        zero_width_string = "test\u200bstring"
        result = model._execute_basic_callback(0.0, 0.0, zero_width_string)
        assert result == pytest.approx(11.0)
    
    def test_callback_with_bidi_text(self, phreeqcrm_module):
        """Test callback with bidirectional text (LTR and RTL)."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        bidi_string = "Hello שלום"  # English and Hebrew
        result = model._execute_basic_callback(0.0, 0.0, bidi_string)
        # "Hello " (6) + "שלום" (4) = 10 characters
        assert result == pytest.approx(10.0)
    
    def test_callback_with_unicode_math_symbols(self, phreeqcrm_module):
        """Test callback with unicode mathematical symbols."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: len(m))
        math_string = "∑∫∂∇"
        result = model._execute_basic_callback(0.0, 0.0, math_string)
        assert result == pytest.approx(4.0)
    
    def test_callback_with_unicode_in_conditional_logic(self, phreeqcrm_module):
        """Test callback using unicode string in conditional logic."""
        def unicode_aware_callback(v1, v2, msg, cookie):
            if msg == "你好":
                return 100.0
            elif msg == "مرحبا":
                return 200.0
            else:
                return 0.0

        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(unicode_aware_callback)
        
        result1 = model._execute_basic_callback(0.0, 0.0, "你好")
        assert result1 == pytest.approx(100.0)
        
        result2 = model._execute_basic_callback(0.0, 0.0, "مرحبا")
        assert result2 == pytest.approx(200.0)


class TestIntegration:
    """Integration tests combining multiple operations."""
    
    def test_workflow_set_and_execute(self, phreeqcrm_module):
        """Test a complete workflow of setting and executing callbacks."""
        def workflow_callback(v1, v2, msg, cookie):
            return v1 * v2 + len(msg)
        
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(workflow_callback)
        result = model._execute_basic_callback(3.0, 4.0, "test")
        # 3 * 4 + len("test") = 12 + 4 = 16
        assert result == pytest.approx(16.0)
    
    def test_sequence_of_executions(self, phreeqcrm_module):
        """Test executing the same callback multiple times."""
        model = phreeqcrm_module.BMIPhreeqcRM()
        model.initialize()
        model.set_basic_callback(lambda v1, v2, m, cookie: v1 + v2)
        
        results = [
            model._execute_basic_callback(1.0, 2.0, "a"),
            model._execute_basic_callback(2.0, 3.0, "b"),
            model._execute_basic_callback(3.0, 4.0, "c"),
        ]
        
        expected = [3.0, 5.0, 7.0]
        assert all(r == pytest.approx(e) for r, e in zip(results, expected))


class TestCallbacks:
    def test_set_python_callback_and_invoke(self, phreeqcrm_module):
        """Test setting a Python callback and invoking it."""
        called = []

        def pycb(x, y, msg, cookie):
            called.append((x, y, msg, cookie))
            return 123.5
        
        if phreeqcrm_module.has_mpi():
            from mpi4py import MPI
            model = phreeqcrm_module.BMIPhreeqcRM(10, MPI.COMM_WORLD)
            model.initialize()
            assert MPI.COMM_WORLD.Get_size() == 1
        elif phreeqcrm_module.has_openmp():
            model = phreeqcrm_module.BMIPhreeqcRM(10, 1)
            model.initialize()
            assert model.GetThreadCount() == 1
        else:
            model = phreeqcrm_module.BMIPhreeqcRM(10, 1)
            model.initialize()
            assert model.GetThreadCount() == 1

        assert model is not None
        model.initialize()
        model.set_basic_callback(pycb, None)

        value = model._execute_basic_callback(3.0, 4.0, "test")
        assert value == pytest.approx(123.5)
        # nthreads + 2
        assert called == [(3.0, 4.0, "test", None), (3.0, 4.0, "test", None), (3.0, 4.0, "test", None)]

    def test_default_cookie_is_none(self, phreeqcrm_module):
        """Test setting callback without cookie sets cookie to None."""
        model = phreeqcrm_module.BMIPhreeqcRM()

        def pycb(x, y, msg, cookie):
            print(f"Callback called with x={x}, y={y}, msg='{msg}', cookie={cookie}")
            assert cookie is None
            return x + y

        model.initialize()
        model.set_basic_callback(pycb)
        result = model._execute_basic_callback(3.0, 4.0, "No cookie")
        assert result == pytest.approx(7.0)

    def test_callback_with_cookie(self, phreeqcrm_module):
        """Test setting callback with a cookie value."""
        model = phreeqcrm_module.BMIPhreeqcRM()

        def pycb(x, y, msg, cookie):
            print(f"Callback called with x={x}, y={y}, msg='{msg}', cookie={cookie}")
            assert cookie == "test_cookie"
            return x + y + cookie.count("t")

        model.initialize()
        model.set_basic_callback(pycb, "test_cookie")
        result = model._execute_basic_callback(3.0, 4.0, "With cookie")
        assert result == pytest.approx(9.0)  # 3 + 4 + 2 (count of 't' in "test_cookie")

    def test_callback_with_model_cookie(self, phreeqcrm_module):
        """Test setting callback with a model as cookie value."""
        model = phreeqcrm_module.BMIPhreeqcRM()

        def pycb(x, y, msg, cookie):
            print(f"Callback called with x={x}, y={y}, msg='{msg}', cookie={cookie}")
            assert cookie is model
            return 0.0

        model.initialize()
        model.set_basic_callback(pycb, model)  # Using the model instance as cookie
        result = model._execute_basic_callback(3.0, 4.0, "With model as cookie")
        assert result == pytest.approx(0.0)



if __name__ == "__main__":
    pytest.main([__file__, "-v"])
