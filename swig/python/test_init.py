import os
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_less
import pytest

from phreeqcrm import BMIPhreeqcRM, IRM_OK, State

from constants import FilePaths

ERROR_MUST_INITIALIZE             = "must call initialize first"
ERROR_NOT_IMPLEMENTED             = "Not Implemented"

###############################################################################
def test_initialized():
    model = BMIPhreeqcRM()
    assert model._state == State.UNINITIALIZED

def test_finalize_not_initialized():
    model = BMIPhreeqcRM()

    assert model._state == State.UNINITIALIZED
    model.finalize()
    assert model._state == State.UNINITIALIZED

def test_initialize_AdvectBMI():
    model = BMIPhreeqcRM()

    assert model._state == State.UNINITIALIZED
    model.initialize(FilePaths.YAML)
    assert model._state == State.INITIALIZED
###############################################################################

#    def finalize(self) -> None:
def test_finalize_no_throw():
    model = BMIPhreeqcRM()
    model.finalize()

#    def get_component_name(self) -> str:
def test_get_component_name_no_throw():
    model = BMIPhreeqcRM()
    component_name = model.get_component_name()

#    def get_current_time(self) -> float:
def test_get_current_time_no_throw():
    model = BMIPhreeqcRM()
    current_time = model.get_current_time()

#    def get_end_time(self) -> float:
def test_get_end_time_raises_uninitialized():
    model = BMIPhreeqcRM()
    end_time = model.get_end_time()

#    def get_grid_edge_count(self, grid: int) -> int:
def test_get_grid_edge_count_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_edge_count(0)

#    def get_grid_edge_nodes(self, grid: int, edge_nodes: np.ndarray) -> np.ndarray:
def test_get_grid_edge_nodes_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_edge_nodes(0, None)

#    def get_grid_face_count(self, grid: int) -> int:
def test_get_grid_face_count_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_face_count(0)

#    def get_grid_face_edges(self, grid: int, face_edges: np.ndarray) -> np.ndarray:
def test_get_grid_face_edges_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_face_edges(0, None)

#    def get_grid_face_nodes(self, grid: int, face_nodes: np.ndarray) -> np.ndarray:
def test_get_grid_face_nodes_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_face_nodes(0, None)

#    def get_grid_node_count(self, grid: int) -> int:
def test_get_grid_node_count_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_node_count(0)

#    def get_grid_nodes_per_face(self, grid: int, nodes_per_face: np.ndarray) -> np.ndarray:
def test_get_grid_nodes_per_face_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_nodes_per_face(0, None)

#    def get_grid_origin(self, grid: int, origin: np.ndarray) -> np.ndarray:
def test_get_grid_origin_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_origin(0, None)

#    def get_grid_rank(self, grid: int) -> int:
def test_get_grid_rank_no_throw():
    model = BMIPhreeqcRM()
    assert model.get_grid_rank(0) == 1

#    def get_grid_shape(self, grid: int, shape: np.ndarray) -> np.ndarray:
def test_get_grid_shape_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_shape(0, None)

#    def get_grid_size(self, grid: int) -> int:
def test_get_grid_size_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_grid_size(0)

#    def get_grid_spacing(self, grid: int, spacing: np.ndarray) -> np.ndarray:
def test_get_grid_spacing_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_spacing(0, None)

#    def get_grid_type(self, grid: int) -> str:
def test_get_grid_type_no_throw():
    model = BMIPhreeqcRM()
    assert model.get_grid_type(0) == "points"

#    def get_grid_x(self, grid: int, x: np.ndarray) -> np.ndarray:
def test_get_grid_x_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_x(0, None)

#    def get_grid_y(self, grid: int, y: np.ndarray) -> np.ndarray:
def test_get_grid_y_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_y(0, None)

#    def get_grid_z(self, grid: int, z: np.ndarray) -> np.ndarray:
def test_get_grid_z_raises_notimplemented():
    model = BMIPhreeqcRM()
    with pytest.raises(NotImplementedError, match=ERROR_NOT_IMPLEMENTED):
        model.get_grid_z(0, None)

#    def get_input_item_count(self) -> int:
def test_get_input_item_count_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_input_item_count()

#    def get_input_var_names(self) -> Tuple[str]:
def test_get_input_var_names_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_input_var_names()

#    def get_output_item_count(self) -> int:
def test_get_output_item_count_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_output_item_count()

#    def get_output_var_names(self) -> Tuple[str]:
def test_get_output_var_names_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_output_var_names()

#    def get_start_time(self) -> float:
def test_get_start_time_no_throw():
    model = BMIPhreeqcRM()
    start_time = model.get_start_time()

#    def get_time_step(self) -> float:
def test_get_time_step_no_throw():
    model = BMIPhreeqcRM()
    time_step = model.get_time_step()

#    def get_time_units(self) -> str:
def test_get_time_units_no_throw():
    model = BMIPhreeqcRM()
    assert model.get_time_units() == "seconds"

#    def get_value(self, name: str, dest: np.ndarray) -> np.ndarray:
def test_get_value_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_value("Temperature", None)

#    def get_value_at_indices(self, name: str, dest: np.ndarray, inds: np.ndarray) -> np.ndarray:
def test_get_value_at_indices_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_value_at_indices("Temperature", None, None)

#    def get_value_ptr(self, name: str) -> np.ndarray:
def test_get_value_ptr_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_value_ptr("Temperature")

#    def get_var_grid(self, name: str) -> int:
def test_get_var_grid_no_throw():
    model = BMIPhreeqcRM()
    assert model.get_var_grid("Temperature") == 0

#    def get_var_itemsize(self, name: str) -> int:
def test_get_var_itemsize_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_var_itemsize("Temperature")

#    def get_var_location(self, name: str) -> str:
def test_get_var_location_no_throw():
    model = BMIPhreeqcRM()
    assert model.get_var_location("Temperature") == "Unknown"

#    def get_var_nbytes(self, name: str) -> int:
def test_get_var_nbytes_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_var_nbytes("Temperature")

#    def get_var_type(self, name: str) -> str:
def test_get_var_type_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_var_type("Temperature")

#    def get_var_units(self, name: str) -> str:
def test_get_var_units_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.get_var_units("Temperature")

#    def initialize(self, config_file: str) -> None:
def test_initialize_no_throw():
    model = BMIPhreeqcRM()
    model.initialize()

#    def set_value(self, name: str, src: np.ndarray) -> None:
def test_set_value_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.set_value("Temperature", None)

#    def set_value_at_indices(self, name: str, inds: np.ndarray, src: np.ndarray) -> None:
def test_set_value_at_indices_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.set_value_at_indices("Temperature", None, None)

#    def update(self) -> None:
def test_update_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.update()

#    def update_until(self, time: float) -> None:
def test_update_until_raises_uninitialized():
    model = BMIPhreeqcRM()
    with pytest.raises(RuntimeError, match=ERROR_MUST_INITIALIZE):
        model.update_until(1.0)
