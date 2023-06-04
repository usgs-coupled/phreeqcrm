#include <iostream>
#include <memory>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "BMIPhreeqcRM.h"
#include "VarManager.h"

#include "docstrings.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

#define PYBIND11_DETAILED_ERROR_MESSAGES

/*
* Get a copy of values of the given variable
*/
py::array BMIPhreeqcRM::get_value(std::string name, py::array dest)
{
    if (!this->_initialized) throw NotIntialized();

    assert(this->var_man);
    RMVARS v_enum = this->var_man->GetEnum(name);
    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");

    BMIVariant& bv = this->var_man->VariantMap[v_enum];

    assert(this->language == "Py");

    if (bv.GetHasPtr())
    {
        // this call is REQUIRED for the BMIVariant::GetVoidPtr() to be valid
        this->GetValuePtr(name);
    }

    assert(this->GetVarType(name) == bv.GetPType());

    // dest dtype
    py::ssize_t dst_ndim = dest.ndim();
    py::dtype dst_dt = dest.dtype();
    std::string dst_dt_str = std::string(dst_dt.str());

    // check for scalars
    int dim = bv.GetDim();
    if (!bv.GetHasPtr() && dim == 1)
    {
        // scalar
        py::array a;
        if (bv.GetPType() == "float64")
        {
            assert(bv.GetItemsize() == sizeof(double));
            assert(bv.GetNbytes() == sizeof(double));

            void* dst_ptr = dest.mutable_data();
            const void* src_ptr = bv.GetDVarPtr();
            assert(src_ptr);

            std::memcpy(dst_ptr, src_ptr, bv.GetNbytes());

            return dest;
        }
        else if (bv.GetPType() == "int32")
        {
            assert(bv.GetItemsize() == sizeof(int));
            assert(bv.GetNbytes() == sizeof(int));

            void* dst_ptr = dest.mutable_data();
            const void* src_ptr = bv.GetIVarPtr();
            assert(src_ptr);

            std::memcpy(dst_ptr, src_ptr, bv.GetNbytes());

            return dest;
        }
        else if (bv.GetPType().rfind("<U", 0) == 0)
        {
            // This is a hack to convert the string to unicode
            std::vector<std::string> strings;
            strings.push_back(bv.GetStringVar());
            py::array src = py::cast(strings);

            ssize_t nsrc = src.request().itemsize;
            ssize_t ndst = dest.request().itemsize;

            if (src.request().itemsize > dest.request().itemsize)
            {
                throw std::runtime_error("buffer too small");
            }

            std::memset(dest.mutable_data(), 0, dest.request().itemsize);
            std::memcpy(dest.mutable_data(), src.request().ptr, src.request().itemsize);

            return dest;
        }
        else
        {
            assert(false);
        }
    }

    // src array
    auto src = get_value_ptr(name);

    // src dtype
    py::ssize_t src_ndim = src.ndim();
    py::dtype src_dt = src.dtype();
    std::string src_dt_str = std::string(src_dt.str());

    if (dst_dt != src_dt)
    {
        throw std::runtime_error("bad array dtype");
    }

    auto src_buf = src.request();  // obtain buffer information for source array
    auto dst_buf = dest.request();  // obtain buffer information for destination array
    void* src_ptr = src_buf.ptr;
    void* dst_ptr = dest.mutable_data();

    std::memcpy(dst_ptr, src_ptr, src_buf.size * src_buf.itemsize);

    return dest;
}

py::array BMIPhreeqcRM::get_value_at_indices(std::string name, py::array dest, py::array indices)
{
    // name and initilized validated in get_value_ptr
    py::array src = get_value_ptr(name);

    // validate dest and indices
    py::buffer_info dest_info = dest.request();
    py::buffer_info indices_info = indices.request();

    if (dest_info.ndim != 1 || indices_info.ndim != 1)
    {
        throw std::runtime_error("Both dest and indices must be 1-dimensional arrays.");
    }

    if (dest_info.shape[0] != indices_info.shape[0])
    {
        throw std::runtime_error("The size of dest and indices must be the same.");
    }

    // src dtype
    py::dtype src_dt = src.dtype();

    // dest dtype
    py::dtype dst_dt = dest.dtype();

    // indices dtype
    py::dtype indices_dt = indices.dtype();

    // validate dtypes
    if (dst_dt != src_dt || std::string(indices_dt.str()) != "int32")
    {
        throw std::runtime_error("bad array dtype");
    }


    int* indices_data = static_cast<int*>(indices_info.ptr);

    std::string src_dt_str = std::string(src_dt.str());
    if (src_dt_str == "float64")
    {
        double* dest_data = static_cast<double*>(dest_info.ptr);

        py::array_t<double> double_array(src);
        auto src_data = double_array.unchecked<1>();

        for (py::ssize_t i = 0; i < dest_info.shape[0]; ++i)
        {
            dest_data[i] = src_data(static_cast<size_t>(indices_data[i]));
        }
    }
    else if (src_dt_str == "int32")
    {
        int* dest_data = static_cast<int*>(dest_info.ptr);

        py::array_t<int> int_array(src);
        auto src_data = int_array.unchecked<1>();

        for (py::ssize_t i = 0; i < dest_info.shape[0]; ++i)
        {
            dest_data[i] = src_data(static_cast<size_t>(indices_data[i]));
        }
    }
    else if (src_dt_str.rfind("<U", 0) == 0)
    {
        //py::array_t<std::string> str_array(src);
        //auto src_data = int_array.unchecked<1>();

        //for (py::ssize_t i = 0; i < dest_info.shape[0]; ++i)
        //{
        //    dest_data[i] = src_data(static_cast<size_t>(indices_data[i]));
        //}
        assert(false);
        throw std::runtime_error("str are unsupported in get_value_at_indices");
    }
    else
    {
        assert(false);
        throw std::runtime_error("bad array dtype");
    }


    return dest;
}

py::array BMIPhreeqcRM::get_value_ptr(std::string name)
{
    if (!this->_initialized) throw NotIntialized();

    assert(this->var_man);
    RMVARS v_enum = this->var_man->GetEnum(name);
    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");

    BMIVariant& bv = this->var_man->VariantMap[v_enum];

    if (bv.HasPyArray())
    {
        return bv.GetPyArray();
    }

    // this call is REQUIRED for the BMIVariant::GetVoidPtr() to be valid
    // and BMIVariant::GetHasPtr()
    void* ptr = this->GetValuePtr(name);

    if (!bv.GetHasPtr() && !(bv.GetPType().rfind("<U", 0) == 0))
    {
        throw std::runtime_error("This variable does not support get_value_ptr.");
    }

    assert(this->language == "Py");

    //// this call is REQUIRED for the BMIVariant::GetVoidPtr() to be valid
    //void* ptr = this->GetValuePtr(name);

    assert(this->GetVarType(name) == bv.GetPType());

    bool read_only = !(bv.GetHasSetter());

    py::array a = py::none();

    if (bv.GetPType() == "float64")
    {
        assert(bv.GetItemsize() == sizeof(double));
        assert(bv.GetVoidPtr() == this->GetValuePtr(name));

        py::array_t<double> arr({ bv.GetDim() }, (double*)bv.GetVoidPtr(), py::cast(bv.GetVoidPtr()));
        a = arr;

        // Check that the array shares its data
        assert(!arr.owndata());
        assert(arr.data() == bv.GetVoidPtr());

        py::dtype dt = arr.dtype();
        assert(std::string(dt.str()) == bv.GetPType());
    }
    else if (bv.GetPType() == "int32")
    {
        assert(bv.GetItemsize() == sizeof(int));
        assert(bv.GetVoidPtr() == this->GetValuePtr(name));

        py::array_t<int> arr({ bv.GetDim() }, (int*)bv.GetVoidPtr(), py::cast(bv.GetVoidPtr()));
        a = arr;

        // Check that the array shares its data
        assert(!arr.owndata());
        assert(arr.data() == bv.GetVoidPtr());

        py::dtype dt = arr.dtype();
        assert(std::string(dt.str()) == bv.GetPType());
    }
    else if (bv.GetPType().rfind("<U", 0) == 0)
    {
        std::vector<std::string> &strings = bv.GetStringVectorRef();
        assert(strings.size() > 0);

        a = py::cast(strings);
        // string arrays are always readonly
        read_only = true;

        py::dtype dt = a.dtype();
        assert(std::string(dt.str()) == bv.GetPType());
    }
    else
    {
        assert(false);
    }

    if (read_only)
    {
        assert(a.writeable());
        py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
        assert(!a.writeable());
    }
    return bv.SetPyArray(a);
}

void BMIPhreeqcRM::set_value(std::string name, py::array arr)
{
    if (!this->_initialized) throw NotIntialized();

    RMVARS v_enum = this->var_man->GetEnum(name);
    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");

    BMIVariant& bv = this->var_man->VariantMap[v_enum];

    if (!bv.GetInitialized())
    {
        this->var_man->task = VarManager::VAR_TASKS::Info;
        ((*this->var_man).*bv.GetFn())();
    }

    py::ssize_t ndim = arr.ndim();
    py::dtype dt = arr.dtype();
    std::string dt_str = std::string(dt.str());

    std::string bv_str = bv.GetPType();
    int dim = bv.GetDim();

    const py::ssize_t* shape = arr.shape();
    py::ssize_t sz = arr.size();

    // Check dimension
    if (!(ndim > 0 && shape[0] == dim))
    {
        std::ostringstream oss;
        oss << "dimension error in set_value: " << name;
        this->ErrorMessage(oss.str());
        throw std::runtime_error(oss.str());
    }

    if (dt_str == "float64")
    {
        //std::vector<double>& ref = this->var_man->VarExchange.GetDoubleVectorRef();
        //ref.resize(arr.size());        
        //std::memcpy(ref.data(), arr.data(), arr.size() * sizeof(double));
        assert(bv.GetCType() == "double" && dim > 1);
        this->SetValue(name, (void*)arr.data());
    }
    else if (dt_str == "int32")
    {
        //std::vector<int>& ref = this->var_man->VarExchange.GetIntVectorRef();
        //ref.resize(arr.size());
        //std::memcpy(ref.data(), arr.data(), arr.size() * sizeof(int));
        assert(bv.GetCType() == "int" && dim > 1);
        this->SetValue(name, (void*)arr.data());
    }
    else if (dt_str.rfind("<U", 0) == 0)
    {
        assert(dt.kind() == 'U');

        std::vector<std::string> result;
        for (auto item : arr)
        {
            result.push_back(py::cast<std::string>(item));
        }

        size_t n = result.size();
        assert(n > 0);
        assert(shape[0] == n);
        assert(ndim == 1);
        assert(sz == n);

        if (n == 1)
        {
            this->SetValue(name, result[0]);

            std::string value;
            this->GetValue(name, value);
            assert(value == "prefix");
        }
        else
        {
            // this should have been caught in check dims above
            assert(false);
        }
    }
    else if (dt_str == "bool")
    {
        std::vector<bool> result;
        for (auto item : arr)
        {
            bool b = py::cast<bool>(item);
            result.push_back(py::cast<bool>(item));
        }

        size_t n = result.size();
        assert(n > 0);
        assert(shape[0] == n);
        assert(ndim == 1);
        assert(sz == n);

        if (n == 1)
        {
            this->SetValue(name, result[0]);
        }
        else
        {
            // this should have been caught in check dims above
            assert(false);
        }
    }
    return;
}

void BMIPhreeqcRM::set_value_at_indices(std::string name, py::array indices, py::array src)
{
    // name and initilized validated in get_value_ptr
    py::array dest = get_value_ptr(name);
    //{{
    py::buffer_info dest_info = dest.request();
    //}}

    // validate src and indices
    py::buffer_info src_info = src.request();
    py::buffer_info indices_info = indices.request();

    if (src_info.ndim != 1 || indices_info.ndim != 1)
    {
        throw std::runtime_error("Both src and indices must be 1-dimensional arrays.");
    }

    if (src_info.shape[0] != indices_info.shape[0])
    {
        throw std::runtime_error("The size of src and indices must be the same.");
    }

    // src dtype
    py::dtype src_dt = src.dtype();

    // dest dtype
    py::dtype dst_dt = dest.dtype();

    // indices dtype
    py::dtype indices_dt = indices.dtype();

    // validate dtypes
    if (dst_dt != src_dt || std::string(indices_dt.str()) != "int32")
    {
        throw std::runtime_error("bad array dtype");
    }

    int* indices_data = static_cast<int*>(indices_info.ptr);

    std::string src_dt_str = std::string(src_dt.str());
    if (src_dt_str == "float64")
    {
        double* dest_data = static_cast<double*>(dest_info.ptr);

        py::array_t<double> double_array(src);
        auto src_data = double_array.unchecked<1>();

        for (py::ssize_t i = 0; i < src_info.shape[0]; ++i)
        {
            dest_data[indices_data[i]] = src_data(static_cast<size_t>(i));
        }
    }
    else if (src_dt_str == "int32")
    {
        int* dest_data = static_cast<int*>(dest_info.ptr);

        py::array_t<int> int_array(src);
        auto src_data = int_array.unchecked<1>();

        for (py::ssize_t i = 0; i < src_info.shape[0]; ++i)
        {
            dest_data[indices_data[i]] = src_data(static_cast<size_t>(i));
        }
    }
    else if (src_dt_str.rfind("<U", 0) == 0)
    {
        //py::array_t<std::string> str_array(src);
        //auto src_data = int_array.unchecked<1>();

        //for (py::ssize_t i = 0; i < dest_info.shape[0]; ++i)
        //{
        //    dest_data[i] = src_data(static_cast<size_t>(indices_data[i]));
        //}
        assert(false);
        throw std::runtime_error("str are unsupported in get_value_at_indices");
    }
    else
    {
        assert(false);
        throw std::runtime_error("bad array dtype");
    }
}

PYBIND11_MODULE(phreeqcrm, m) {

    py::class_<BMIPhreeqcRM>(m, "bmi_phreeqcrm" , py::dynamic_attr())

        // Constructor
        .def(py::init())

        //
        // Model control functions.
        //

        // virtual void Initialize(std::string config_file) = 0;
        // def initialize(self, config_file: str) -> None:
        .def("initialize",
            &BMIPhreeqcRM::Initialize,
            initialize_docstring.c_str(),
            py::arg("config_file") = py::str("")
        )

        // virtual void Update() = 0;
        // def update(self) -> None:
        .def("update",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                self.Update();
            },
            update_docstring.c_str()
        )

        // virtual void UpdateUntil(double time) = 0;
        // def update_until(self, time: float) -> None: 
        .def("update_until",
            [](BMIPhreeqcRM& self, double time) {
                if (!self._initialized) throw NotIntialized();
                self.UpdateUntil(time);
            },
            update_until_docstring.c_str(),
            py::arg("time")
        )

        // virtual void Finalize() = 0;
        // def finalize(self) -> None:
        .def("finalize",
            [](BMIPhreeqcRM& self) {
                self.Finalize();
            }
        )

        //
        // Model information functions.
        // 

        // virtual std::string GetComponentName() = 0;
        // def get_component_name(self) -> str:
        .def("get_component_name",
            [](BMIPhreeqcRM& self) {
                return self.GetComponentName();
            },
            get_component_name_docstring.c_str()
        )

        // virtual int GetInputItemCount() = 0;
        // def get_input_item_count(self) -> int:
        .def("get_input_item_count",
            &BMIPhreeqcRM::GetInputItemCount
        )

        // virtual int GetOutputItemCount() = 0;
        // def get_input_item_count(self) -> int:
        .def("get_output_item_count",
            &BMIPhreeqcRM::GetOutputItemCount
        )

        // virtual std::vector<std::string> GetInputVarNames() = 0;
        // def get_input_var_names(self) -> Tuple[str]:
        .def("get_input_var_names",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                auto output = self.GetInputVarNames();
                py::tuple out = py::cast(output);
                return out;
            }
        )

        // virtual std::vector<std::string> GetOutputVarNames() = 0;
        // def get_output_var_names(self) -> Tuple[str]:
        .def("get_output_var_names",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                auto output = self.GetOutputVarNames();
                py::tuple out = py::cast(output);
                return out;
            }
        )

        //
        // Variable information functions
        //

        // virtual int GetVarGrid(std::string name) = 0;
        // def get_var_grid(self, name: str) -> int:
        .def("get_var_grid",
            &BMIPhreeqcRM::GetVarGrid,
            py::arg("name")
        )

        // virtual std::string GetVarType(std::string name) = 0;
        // def get_var_type(self, name: str) -> str:
        .def("get_var_type",
            &BMIPhreeqcRM::GetVarType,
            py::arg("name")
        )

        // virtual std::string GetVarUnits(std::string name) = 0;
        // def get_var_units(self, name: str) -> str:
        .def("get_var_units",
            &BMIPhreeqcRM::GetVarUnits,
            py::arg("name")
        )

        // virtual int GetVarItemsize(std::string name) = 0;
        // def get_var_itemsize(self, name: str) -> int:
        .def("get_var_itemsize",
            &BMIPhreeqcRM::GetVarItemsize,
            py::arg("name")
        )

        
        // virtual int GetVarNbytes(std::string name) = 0;
        // def get_var_nbytes(self, name: str) -> int:
        .def("get_var_nbytes",
            &BMIPhreeqcRM::GetVarNbytes,
            py::arg("name")
        )

        // virtual std::string GetVarLocation(std::string name) = 0;
        // def get_var_location(self, name: str) -> str:
        .def("get_var_location",
            [](BMIPhreeqcRM& self, std::string name) {
                if (!self._initialized) throw NotIntialized();
                return self.GetVarLocation(name);
            },
            py::arg("name")
        )

        // virtual double GetCurrentTime() = 0;
        // def get_current_time(self) -> float:
        .def("get_current_time",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                return self.GetCurrentTime();
            }
        )

        // virtual double GetStartTime() = 0;
        // def get_start_time(self) -> float:
        .def("get_start_time",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                return self.GetStartTime();
            }
        )

        // virtual double GetEndTime() = 0;
        // def get_end_time(self) -> float:
        .def("get_end_time",
            [](BMIPhreeqcRM& self) {
                if (!self._initialized) throw NotIntialized();
                return self.GetEndTime();
            }
        )

        // virtual std::string GetTimeUnits() = 0;
        // def get_time_units(self) -> str:
        .def("get_time_units",
            [](BMIPhreeqcRM& self) {
                return self.GetTimeUnits();
            }
        )

        // virtual double GetTimeStep() = 0;
        // def get_time_step(self) -> float:
        .def("get_time_step",
            [](BMIPhreeqcRM& self) {
                return self.GetTimeStep();
            }
        )

        //
        // Variable getters
        //

        // virtual void GetValue(std::string name, void *dest) = 0;
        // def get_value(self, name: str, dest: np.ndarray) -> np.ndarray:
        .def("get_value",
            &BMIPhreeqcRM::get_value,
            py::arg("name"),
            py::arg("dest")
        )

        // virtual void *GetValuePtr(std::string name) = 0;
        // def get_value_ptr(self, name: str) -> np.ndarray:
        .def("get_value_ptr",
            &BMIPhreeqcRM::get_value_ptr,
            py::arg("name")
        )

        // virtual void GetValueAtIndices(std::string name, void *dest, int *inds, int count) = 0;
        // def get_value_at_indices(
        //     self, name: str, dest : np.ndarray, inds : np.ndarray
        // )->np.ndarray:
        .def("get_value_at_indices",
            &BMIPhreeqcRM::get_value_at_indices,
            py::arg("name"),
            py::arg("dest"),
            py::arg("inds")
        )

        //
        // Variable setters
        //

        // virtual void SetValue(std::string name, void *src) = 0;
        // def set_value(self, name: str, src: np.ndarray) -> None:
        .def("set_value",
            &BMIPhreeqcRM::set_value,
            py::arg("name"),
            py::arg("src")
        )

        // virtual void SetValueAtIndices(std::string name, int *inds, int count, void *src) = 0;
        // def set_value_at_indices(
        //     self, name: str, inds : np.ndarray, src : np.ndarray
        // )->None:
        .def("set_value_at_indices",
            &BMIPhreeqcRM::set_value_at_indices,
            py::arg("name"),
            py::arg("inds"),
            py::arg("src")
        )

        //
        // Grid information functions
        // 

        // virtual int GetGridRank(const int grid) = 0;
        // def get_grid_rank(self, grid: int) -> int:
        .def("get_grid_rank",
            &BMIPhreeqcRM::GetGridRank,
            py::arg("grid")
        )

        // virtual int GetGridSize(const int grid) = 0;
        // def get_grid_size(self, grid: int) -> int:
        .def("get_grid_size",
            [](BMIPhreeqcRM& self, int grid) {
                if (!self._initialized) throw NotIntialized();
                return self.GetGridSize(grid);
            },
            get_grid_size_docstring.c_str(),
            py::arg("grid")
        )

        // virtual std::string GetGridType(const int grid) = 0;
        // def get_grid_type(self, grid: int) -> str:
        .def("get_grid_type",
            &BMIPhreeqcRM::GetGridType,
            py::arg("grid")
        )

        // void GetGridShape(const int grid, int* shape) override
        // def get_grid_shape(self, grid: int, shape: np.ndarray) -> np.ndarray:
        .def("get_grid_shape",
            [](BMIPhreeqcRM& self, int grid, py::array_t<int> shape) {
                throw NotImplemented();
                if (!self._initialized) throw NotIntialized();

                int* data = shape.mutable_data();
                py::array::ShapeContainer shapeContainer{1};
                shape.reshape(shapeContainer);
                data[0] = self.GetGridCellCount();
                return shape;
            }
            ,
            //get_grid_shape_docstring.c_str(),
            py::arg("grid"),
            py::arg("shape")
        )

        // virtual void GetGridSpacing(const int grid, double *spacing) = 0;
        // def get_grid_spacing(self, grid: int, spacing: np.ndarray) -> np.ndarray:
        .def("get_grid_spacing",
            [](BMIPhreeqcRM& self, int grid, py::array_t<double> spacing) {
                throw NotImplemented();
                if (!self._initialized) throw NotIntialized();
                return spacing;
            }
            ,
            //get_grid_spacing_docstring.c_str(),
            py::arg("grid"),
            py::arg("spacing")
        )

        // virtual void GetGridOrigin(const int grid, double *origin) = 0;
        // def get_grid_origin(self, grid: int, origin: np.ndarray) -> np.ndarray:
        .def("get_grid_origin",
            [](BMIPhreeqcRM& self, int grid, py::array_t<double> origin) {
                throw NotImplemented();
                if (!self._initialized) throw NotIntialized();
                return origin;
            }
            ,
            //get_grid_origin_docstring.c_str(),
            py::arg("grid"),
            py::arg("origin")
        )

        // virtual void GetGridX(const int grid, double* x) = 0;
        // def get_grid_x(self, grid: int, x: np.ndarray) -> np.ndarray:
        .def("get_grid_x",
            [](BMIPhreeqcRM& self, int grid, py::array x) {
                throw NotImplemented();
                return x;
            },
            //get_grid_x_docstring.c_str(),
            py::arg("grid"),
            py::arg("x")
        )

        // virtual void GetGridY(const int grid, double* y) = 0;
        // def get_grid_y(self, grid: int, y: np.ndarray) -> np.ndarray:
        .def("get_grid_y",
            [](BMIPhreeqcRM& self, int grid, py::array y) {
                throw NotImplemented();
                return y;
            },
            //get_grid_y_docstring.c_str(),
            py::arg("grid"),
            py::arg("y")
        )

        // virtual void GetGridZ(const int grid, double* z) = 0;
        // def get_grid_y(self, grid: int, z: np.ndarray) -> np.ndarray:
        .def("get_grid_z",
            [](BMIPhreeqcRM& self, int grid, py::array z) {
                throw NotImplemented();
                return z;
            },
            //get_grid_z_docstring.c_str(),
            py::arg("grid"),
            py::arg("z")
        )

        // virtual int GetGridNodeCount(const int grid) = 0;
        // def get_grid_node_count(self, grid: int) -> int:
        .def("get_grid_node_count",
            [](BMIPhreeqcRM& self, int grid) {
                //throw NotImplemented();
                return self.GetGridNodeCount(grid);
            },
            // get_grid_node_count_docstring.c_str(),
            py::arg("grid")
        )

        // virtual int GetGridEdgeCount(const int grid) = 0;
        // def get_grid_edge_count(self, grid: int) -> int:
        .def("get_grid_edge_count",
            [](BMIPhreeqcRM& self, int grid) {
                // throw NotImplemented();
                return self.GetGridEdgeCount(grid);
            },
            // get_grid_edge_count_docstring.c_str(),
            py::arg("grid")
        )

        // virtual int GetGridFaceCount(const int grid) = 0;
        // def get_grid_face_count(self, grid: int) -> int:
        .def("get_grid_face_count",
            [](BMIPhreeqcRM& self, int grid) {
                // throw NotImplemented();
                return self.GetGridFaceCount(grid);
            },
            // get_grid_face_count_docstring.c_str(),
            py::arg("grid")
        )

        // virtual void GetGridEdgeNodes(const int grid, int* edge_nodes) = 0;
        // def get_grid_edge_nodes(self, grid: int, edge_nodes: np.ndarray) -> np.ndarray:
        .def("get_grid_edge_nodes",
            [](BMIPhreeqcRM& self, int grid, py::array_t<int> edge_nodes) {
                throw NotImplemented();
                return edge_nodes;
            },
            // get_grid_edge_nodes_docstring.c_str(),
            py::arg("grid"),
            py::arg("edge_nodes")
        )

        // virtual void GetGridFaceEdges(const int grid, int* face_edges) = 0;
        // def get_grid_face_edges(self, grid: int, face_edges: np.ndarray) -> np.ndarray:
        .def("get_grid_face_edges",
            [](BMIPhreeqcRM& self, int grid, py::array_t<int> face_edges) {
                throw NotImplemented();
                return face_edges;
            },
            // get_grid_face_edges_docstring.c_str(),
            py::arg("grid"),
            py::arg("face_edges")
        )

        // virtual void GetGridFaceNodes(const int grid, int* face_nodes) = 0;
        // def get_grid_face_nodes(self, grid: int, face_nodes: np.ndarray) -> np.ndarray:
        .def("get_grid_face_nodes",
            [](BMIPhreeqcRM& self, int grid, py::array_t<int> face_nodes) {
                throw NotImplemented();
                return face_nodes;
            },
            // get_grid_face_nodes_docstring.c_str(),
            py::arg("grid"),
            py::arg("face_nodes")
        )

        // virtual void GetGridNodesPerFace(const int grid, int* nodes_per_face) = 0;
        // def get_grid_nodes_per_face(
        //     self, grid: int, nodes_per_face: np.ndarray
        // ) -> np.ndarray:
        .def("get_grid_nodes_per_face",
            [](BMIPhreeqcRM& self, int grid, py::array_t<int> nodes_per_face) {
                throw NotImplemented();
                return nodes_per_face;
            },
            // get_grid_nodes_per_face_docstring.c_str(),
            py::arg("grid"),
            py::arg("nodes_per_face")
        )

        //
        // Extras
        //


        .def("get_components",
            [](BMIPhreeqcRM& self) {
                // initialization not necessary
                auto output = self.GetComponents();
                py::tuple out = py::cast(output);
                return out;
            },
            get_components_docstring.c_str()
        )

        // int GetGridCellCount(void)
        // def GetGridCellCount(self: phreeqcrm.bmi_phreeqcrm) -> int:
        .def("GetGridCellCount",
            &BMIPhreeqcRM::GetGridCellCount
        )

        // int GetMpiMyself(void) const
        // def GetMpiMyself(self: phreeqcrm.bmi_phreeqcrm) -> int:
        .def("GetMpiMyself",
            &BMIPhreeqcRM::GetMpiMyself
        )

        // int GetThreadCount(void)
        // def GetThreadCount(self: phreeqcrm.bmi_phreeqcrm) -> int:
        .def("GetThreadCount",
            &BMIPhreeqcRM::GetThreadCount
        )


        // IRM_RESULT SetComponentH2O(bool tf);
        // SetComponentH2O(self: phreeqcrm.bmi_phreeqcrm, arg0: bool) -> IRM_RESULT
        .def("SetComponentH2O",
            &BMIPhreeqcRM::SetComponentH2O,
            py::arg("tf")
            )

        .def_readonly("_initialized", &BMIPhreeqcRM::_initialized)
        ;

        py::enum_<IRM_RESULT>(m, "IRM_RESULT")
            .value("IRM_OK",          IRM_OK)
            .value("IRM_OUTOFMEMORY", IRM_OUTOFMEMORY)
            .value("IRM_BADVARTYPE",  IRM_BADVARTYPE)
            .value("IRM_INVALIDARG",  IRM_INVALIDARG)
            .value("IRM_INVALIDROW",  IRM_INVALIDROW)
            .value("IRM_INVALIDCOL",  IRM_INVALIDCOL)
            .value("IRM_BADINSTANCE", IRM_BADINSTANCE)
            .value("IRM_FAIL",        IRM_FAIL)
            .export_values();
}
