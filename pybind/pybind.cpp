#include <vector>
#include <memory>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
// #include <chrono>
// #include <thread>

#include "BMIPhreeqcRM.h"
#include "VarManager.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;



//{{
#define PYBIND11_DETAILED_ERROR_MESSAGES
//}}

//#include <pybind11/numpy.h>
//#include <pybind11/pybind11.h>
//#include <vector>
//#include <memory>
//{{
//#include <iostream>
//#include <pybind11/stl.h>
//#include <pybind11/complex.h>
//#include <pybind11/functional.h>
//#include <pybind11/chrono.h>
//}}

namespace py = pybind11;

// Function that returns a shared_ptr to a std::vector<double> with some values
std::shared_ptr<std::vector<double>> get_shared_vector() {
    auto v = std::make_shared<std::vector<double>>(std::initializer_list<double>{1.0, 2.0, 3.0, 4.0, 5.0});
    return v;
}

// Function that wraps a numpy array around a shared_ptr to a std::vector<double>
py::array_t<double> wrap_shared_vector(const std::shared_ptr<std::vector<double>>& v) {
    auto size = v->size();
    auto data = v->data();
    py::capsule free_when_done(data, [](void* f) {
        double* data = reinterpret_cast<double*>(f);
        std::shared_ptr<std::vector<double>>* v = reinterpret_cast<std::shared_ptr<std::vector<double>>*>(data);
        delete v;
        });
    return py::array(size, data, free_when_done);
}

//// Module definition
//PYBIND11_MODULE(phreeqcrm, m) {
//    //m.def("get_shared_vector", &get_shared_vector, "Get a shared_ptr to a std::vector<double>");
//    //m.def("wrap_shared_vector", &wrap_shared_vector, "Wrap a numpy array around a shared_ptr to a std::vector<double>");
//    m.def("echo3", [](py::str s) {
//        if (s.is_none()) std::cout << "empty";
//        else std::cout << s;
//        }, "A function that outputs a string.", py::arg("s").none(true) = py::none());
//}

//template <typename T>
//py::array_t<T> get_value(std::string name, py::array_t<T, py::array::c_style | py::array::forcecast> array)
//{
//    auto r = array.mutable_unchecked<1>();
//    for (py::ssize_t i = 0; i < r.shape(0); i++) {
//        r(i) = 2 * r(i);
//    }
//    return array;
//}


py::array_t<double> use_mutable_unchecked(std::string name, py::array_t<double, py::array::c_style | py::array::forcecast> array)
{
    auto r = array.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < r.shape(0); i++) {
        r(i) = 2 * r(i);
    }
    return array;
}

py::array get_value(std::string name, py::array array)
{
    py::buffer_info buffer = array.request();
    std::cout << "itemsize = " << buffer.itemsize << std::endl;
    std::cout << "size = " << buffer.size << std::endl;
    std::cout << "format = " << buffer.format << std::endl;
    std::cout << "ndim = " << buffer.ndim << std::endl;
    std::cout << "shape.size() = " << buffer.shape.size() << std::endl;

    std::cout << "shape: (";
    for (auto dim : buffer.shape) {
        std::cout << dim << ", ";
    }
    std::cout << ")" << std::endl;


    std::cout << "strides.size() = " << buffer.strides.size() << std::endl;
    std::cout << "readonly = " << buffer.readonly << std::endl;

    if (buffer.ndim != 1) {
        throw std::domain_error("array has incorrect number of dimensions: "
            + std::to_string(buffer.ndim) + "; expected "
            + std::to_string(1));
    }

    return array;
}

py::array get_value_ptr(std::string name)
{
    // numpy.float64
    static std::vector<double> v = { 1.2, 2.2, 3.2, 4.2, 5.2 };

    py::buffer_info buf(
        v.data(),
        sizeof(double),
        py::format_descriptor<double>::format(),
        1,
        { v.size() },
        { sizeof(double) }
    );

    ///py::detail::get_numpy_internals().get_type_info<double>(true)->format_str;

    ///py::dtype numpy_dtype_double(py::format_descriptor<double>::format());
    // numpy_dtype_double.type().str()

    std::cout << py::format_descriptor<double>::format() << " = " << py::detail::get_numpy_internals().get_type_info<double>(true)->format_str << std::endl;
    std::cout << py::format_descriptor<int>::format() << " = " << py::detail::get_numpy_internals().get_type_info<int>(true)->format_str << std::endl;

    return py::array_t<double>(buf);
}

//std::string get_dtype_string(py::array_t<double>& arr) {
//    py::dtype dtype = arr.dtype();  // Get the dtype of the input array
//    py::dtype::type type = dtype.type();  // Get the numpy dtype object for the dtype
//    return type.str();  // Get the numpy dtype string for the numpy dtype object
//}

// this requires #include <numpy/arrayobject.h>
//std::string double_to_numpy_dtype() {
//    std::string c_format = py::format_descriptor<double>::format();  // Get the C type format string for double
//    PyObject* numpy_dtype = PyArray_DescrFromType(NPY_NOTYPE);  // Create a NumPy dtype object
//    PyObject* resolved_dtype = PyArray_DescrFromTypestr(c_format.c_str());  // Resolve the NumPy dtype from the C type format string
//    PyObject* dtype_repr = PyObject_Repr(resolved_dtype);  // Get the string representation of the resolved dtype
//    std::string numpy_dtype_str = PyUnicode_AsUTF8(dtype_repr);  // Convert the dtype representation to a C++ string
//    Py_XDECREF(numpy_dtype);
//    Py_XDECREF(resolved_dtype);
//    Py_XDECREF(dtype_repr);
//    return numpy_dtype_str;
//}

//py::array BMIPhreeqcRM::get_value_ptr(std::string name)
//{
//    // numpy.float64
//    //static std::vector<double> v = { 1.2, 2.2, 3.2, 4.2, 5.2 };
//
//    //assert(name == "Temperature");
//    //std::vector<double> v;
//    //this->GetValue(name, v);
//
//
//    RMVARS v_enum = this->var_man->GetEnum(name);
//    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");
//
//    BMIVariant& bv = this->var_man->VariantMap[v_enum];
//
//    // int dim = bv.GetDim();
//    // std::cout << "0 dim = " << dim << std::endl;
//
//    bool read_only = !bv.GetHasSetter();
//
//    //assert(name == "ComponentCount");
//    //std::vector<int> v;
//    //// this->GetValue(name, v);
//    //int n;
//    //this->GetValue(name, n);
//    //v.push_back(n);
//
//    void* ptr = this->GetValuePtr(name);            // GetValuePtr takes a variable name and returns a pointer to a current copy of the variable values.Unlike the buffer returned from @ref GetValue, the reference always points to the current values of the variable, even if the model's state has changed.
//    int nbytes = this->GetVarNbytes(name);          // retrieves the total number of bytes that are set for a variable with @ref SetValue
//    int itemsize = this->GetVarItemsize(name);      // GetVarItemsize retrieves size of an individual item that can be set or retrived.Sizes may be sizeof(int), sizeof(double), or a character length for string variables.
//    std::string type = this->GetVarType(name);      // std::string GetVarType(const std::string name)
//
//    // betting that GetDim call needs to be after first call to GetValuePtr()
//    int dim = bv.GetDim();
//    std::cout << "0 dim = " << dim << std::endl;
//
//
//    //ssize_t itemsize{ sizeof(double) };
//    //std::string format(py::format_descriptor<double>::format());
//    std::string format;
//    std::vector<ssize_t> strides;
//    if (type == "float64")
//    {
//        assert(itemsize == sizeof(double));
//        format = py::format_descriptor<double>::format();
//        strides.push_back(sizeof(double));
//    }
//    else if (type == "int32")
//    {
//        assert(itemsize == sizeof(int));
//        format = py::format_descriptor<int>::format();
//        strides.push_back(sizeof(int));
//    }
//    else
//    {
//        assert(false);
//        throw std::runtime_error("Unknown dtype");
//    }
//
//    int count = nbytes / itemsize;
//    //assert(count == v.size());
//    std::cout << "count = " << count << std::endl;
//    std::cout << "dim = " << dim << std::endl;
//    assert(count == dim);
//
//    //buffer_info_format;
//
//    int nn = *((int*)ptr);
//
//    //py::buffer_info buf(
//    //    ptr,
//    //    itemsize, // sizeof(double),
//    //    format, // py::format_descriptor<double>::format(),
//    //    1,
//    //    { count },
//    //    strides,  //{ sizeof(int) } // { sizeof(double) }
//    //    read_only
//    //);
//
//    py::buffer_info buf(
//        ptr,
//        itemsize, // sizeof(double),
//        format, // py::format_descriptor<double>::format(),
//        1,
//        { count },
//        strides  //{ sizeof(int) } // { sizeof(double) }
//    );
//
//
//    ///py::detail::get_numpy_internals().get_type_info<double>(true)->format_str;
//
//    ///py::dtype numpy_dtype_double(py::format_descriptor<double>::format());
//    // numpy_dtype_double.type().str()
//
//    //std::cout << py::format_descriptor<double>::format() << " = " << py::detail::get_numpy_internals().get_type_info<double>(true)->format_str << std::endl;
//    //std::cout << py::format_descriptor<int>::format() << " = " << py::detail::get_numpy_internals().get_type_info<int>(true)->format_str << std::endl;
//
//    if (type == "int")
//    {
//        //if (read_only) return py::array_t<int>(buf, py::capsule(ptr, [](void* ptr) { /* no-op */}));
//
//        //if (read_only) {
//        //    py::array_t<int> a = py::array_t<int>(buf);
//
//        //    bool writeable = a.writeable();
//
//        //    //reinterpret_cast<pybind11::detail::PyArray_Proxy*>(a.ptr())->flags &= ~pybind11::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//        //    reinterpret_cast<py::detail::PyArray_Proxy*>(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//
//        //    writeable = a.writeable();
//
//        //    return a;
//        //}
//
//
//        if (read_only) {
//            py::array_t<int> a = py::array_t<int>(buf);
//            assert(a.writeable());
//            py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//            assert(!a.writeable());
//            //reinterpret_cast<py::detail::PyArray_Proxy*>(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//            return a;
//        }
//
//        return py::array_t<int>(buf);
//    }
//
//    if (read_only) {
//        py::array_t<double> a = py::array_t<double>(buf);
//        reinterpret_cast<py::detail::PyArray_Proxy*>(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//        return a;
//    }
//    return py::array_t<double>(buf);
//}

//py::array BMIPhreeqcRM::get_value_ptr(std::string name)
//{
//    RMVARS v_enum = this->var_man->GetEnum(name);
//    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");
//
//    BMIVariant& bv = this->var_man->VariantMap[v_enum];
//    // if (!bv.GetHasPtr()) throw std::runtime_error("This variable does not support get_value_ptr.");   // this throws on "ComponentCount"
//
//
//    void* ptr = this->GetValuePtr(name);            // GetValuePtr takes a variable name and returns a pointer to a current copy of the variable values.Unlike the buffer returned from @ref GetValue, the reference always points to the current values of the variable, even if the model's state has changed.
//    int nbytes = this->GetVarNbytes(name);          // retrieves the total number of bytes that are set for a variable with @ref SetValue
//    int itemsize = this->GetVarItemsize(name);      // GetVarItemsize retrieves size of an individual item that can be set or retrived.Sizes may be sizeof(int), sizeof(double), or a character length for string variables.
//    std::string type = this->GetVarType(name);      // std::string GetVarType(const std::string name)
//
//    // must call GetDim after GetValuePtr
//    int dim = bv.GetDim();
//    
//    bool read_only = !bv.GetHasSetter();
//
//    std::string format;
//    std::vector<ssize_t> strides;
//    if (type == "float64")
//    {
//        assert(itemsize == sizeof(double));
//        format = py::format_descriptor<double>::format();
//        strides.push_back(sizeof(double));
//    }
//    else if (type == "int32")
//    {
//        assert(itemsize == sizeof(int));
//        format = py::format_descriptor<int>::format();
//        strides.push_back(sizeof(int));
//    }
//    else
//    {
//        assert(false);
//        throw std::runtime_error("Unknown dtype");
//    }
//
//    int count = nbytes / itemsize;
//    assert(count == dim);
//
//    py::buffer_info buf(
//        ptr,
//        itemsize,
//        format,
//        1,
//        { count },
//        strides
//    );
//
//
//    if (type == "int32")
//    {
//        if (read_only) {
//            py::array_t<int> a = py::array_t<int>(buf);
//            assert(a.writeable());
//            py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//            assert(!a.writeable());
//            return a;
//        }
//        return py::array_t<int>(buf);
//    }
//
//    if (read_only) {
//        py::array_t<double> a = py::array_t<double>(buf);
//        assert(a.writeable());
//        py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//        assert(!a.writeable());
//        return a;
//    }
//
//
//    //return py::array_t<double>(buf);
//    py::array_t<double> a = py::array_t<double>(buf);
//
//    //auto r = a.mutable_unchecked<1>();
//    ////r.data(0);
//    //assert(r.data(0) == ptr);
//    //for (py::ssize_t i = 0; i < r.shape(0); i++) {
//    //    r(i) = 2 * r(i);
//    //}
//
//    
//    //if (true) {
//    //    py::array_t<double> a = py::array_t<double>(buf);
//    //    assert(a.owndata());
//    //    py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_OWNDATA_;
//    //    assert(!a.owndata());
//    //    return a;
//    //    // auto view = py::array_t<double>::ensure(arr.ptr(), nullptr);
//    //}
//
//    //if (true) {
//    //    // py::array_t<double>::ensure()
//    //    py::array_t<double> arr({ count }, (const double*)ptr, );
//    //    assert(arr.data() == ptr);
//    //    return arr;
//    //}
//
//    if (true) {
//        py::array_t<double> arr({ 2, 3 });
//        auto view = py::array_t<double>::ensure(arr.ptr());
//        // assert(!view.owndata());
//        assert(arr.data() == view.data());
//    }
//
//    // std::vector<double> vec = { 1.2, 1.3, 1.4 };
//    if (type == "float64") {
//        //// Get a pointer to the vector data.
//        //const double* data_ptr = vec.data();
//
//        //// Get the size of the vector.
//        ////const size_t size = vec.size();
//        //const ssize_t size = vec.size();
//
//        //// Create a NumPy array that shares its data with the vector.
//        //py::array_t<double> arr({ size }, data_ptr);
//
//        // Check that the array shares its data with the vector.
//        //assert(arr.data() == data_ptr);
//
//        std::vector<double> &vec = bv.GetDoubleVectorRef();
//
//        ///py::cast(vec, pybind11::return_value_policy::reference);
//
//        //py::array_t<double> arr( { static_cast<py::ssize_t>(vec.size()) }, vec.data(), py::cast(vec, pybind11::return_value_policy::reference));
//
//        //py::array_t<double> arr({ static_cast<py::ssize_t>(vec.size()) }, vec.data(), py::cast(ptr, pybind11::return_value_policy::reference));
//
//        assert(ptr == bv.GetVoidPtr());
//
//        //py::array_t<double> arr({ static_cast<py::ssize_t>(vec.size()) }, vec.data(), py::cast(ptr, pybind11::return_value_policy::reference));
//        //py::array_t<double> arr({ bv.GetDim() }, vec.data(), py::cast(ptr, pybind11::return_value_policy::reference));
//        ///py::array_t<double> arr({ bv.GetDim() }, (double*)bv.GetVoidPtr(), py::cast(bv.GetVoidPtr(), pybind11::return_value_policy::reference));
//        py::array_t<double> arr({ bv.GetDim() }, (double*)bv.GetVoidPtr(), py::cast(bv.GetVoidPtr()));
//
//        assert(arr.data() == vec.data());
//        assert(!arr.owndata());
//
//        py::dtype dt = arr.dtype();
//        
//        std::cout << "Array data type : " << dt.str() << std::endl;
//
//
//        if (read_only) {
//            assert(arr.writeable());
//            py::detail::array_proxy(arr.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
//            assert(!arr.writeable());
//            return arr;
//        }
//
//        return arr;
//    }
//
//    return a;
//}


//void copy_array(const py::array_t<double>& src, const py::array_t<double>& dst) {
//    auto src_buf = src.request();  // obtain buffer information for source array
//    auto dst_buf = dst.request();  // obtain buffer information for destination array
//    double* src_ptr = static_cast<double*>(src_buf.ptr);  // get pointer to source data
//    double* dst_ptr = static_cast<double*>(dst.mutable_data());  // get pointer to mutable destination data
//    std::memcpy(dst_ptr, src_ptr, src_buf.size * sizeof(double));  // copy data from source to destination
//}

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
    ///if (!bv.GetHasPtr()) throw std::runtime_error("This variable does not support get_value_ptr.");   // this throws on "ComponentCount"

    assert(this->language == "Py");

    // this call is REQUIRED for the BMIVariant::GetVoidPtr() to be valid
    this->GetValuePtr(name);

    assert(this->GetVarType(name) == bv.GetPType());

    // dest dtype
    py::ssize_t dst_ndim = dest.ndim();
    py::dtype dst_dt = dest.dtype();
    std::string dst_dt_str = std::string(dst_dt.str());

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

py::array BMIPhreeqcRM::get_value_ptr(std::string name)
{
    if (!this->_initialized) throw NotIntialized();

    assert(this->var_man);
    RMVARS v_enum = this->var_man->GetEnum(name);
    if (v_enum == RMVARS::NotFound) throw std::runtime_error("Unknown variable name");

    BMIVariant& bv = this->var_man->VariantMap[v_enum];

    //{{
    if (bv.HasPyArray())
    {
        return bv.GetPyArray();
    }
    //}}

    ///if (!bv.GetHasPtr()) throw std::runtime_error("This variable does not support get_value_ptr.");   // this throws on "ComponentCount"

    assert(this->language == "Py");

    // this call is REQUIRED for the BMIVariant::GetVoidPtr() to be valid
    void* ptr = this->GetValuePtr(name);

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


// def get_grid_shape(self, grid: int, shape : np.ndarray)->np.ndarray:


// see https://stackoverflow.com/questions/3898572/what-are-the-most-common-python-docstring-formats
// see https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html#example-google
// see https://numpydoc.readthedocs.io/en/latest/example.html


// def foo(var1, var2, *args, long_var_name = "hi", only_seldom_used_keyword = 0, **kwargs) :
const std::string numpydoc_example_docstring =
R"pbdoc(Summarize the function in one line.

Several sentences providing an extended description.Refer to
variables using back - ticks, e.g. `var`.

Parameters
----------
var1 : array_like
    Array_like means all those objects -- lists, nested lists, etc. --
    that can be converted to an array.We can also refer to
    variables like `var1`.
var2 : int
    The type above can either refer to an actual Python type
    (e.g. ``int``), or describe the type of the variable in more
        detail, e.g. ``(N, ) ndarray`` or ``array_like``.
* args : iterable
    Other arguments.
long_var_name : {'hi', 'ho'}, optional
    Choices in brackets, default first when optional.

Returns
-------
type
    Explanation of anonymous return value of type ``type``.
describe : type
    Explanation of return value named `describe`.
out : type
    Explanation of `out`.
type_without_description

Other Parameters
----------------
only_seldom_used_keyword : int, optional
    Infrequently used parameters can be described under this optional
    section to prevent cluttering the Parameters section.
** kwargs : dict
    Other infrequently used keyword arguments.Note that all keyword
    arguments appearing after the first parameter specified under the
    Other Parameters section, should also be described under this
    section.

Raises
------
BadException
    Because you shouldn't have done that.

See Also
--------
numpy.array : Relationship(optional).
numpy.ndarray : Relationship(optional), which could be fairly long, in
    which case the line wraps here.
numpy.dot, numpy.linalg.norm, numpy.eye

Notes
-----
Notes about the implementation algorithm(if needed).

This can have multiple paragraphs.

You may include some math :

..math::X(e^ { j\omega }) = x(n)e^ { -j\omega n }

And even use a Greek symbol like : math:`\omega` inline.

References
----------
Cite the relevant literature, e.g.[1]_.You may also cite these
references in the notes section above.

.. [1] O.McNoleg, "The integration of GIS, remote sensing,
    expert systems and adaptive co - kriging for environmental habitat
    modelling of the Highland Haggis using object - oriented, fuzzy - logic
    and neural - network techniques, " Computers & Geosciences, vol. 22,
    pp. 585 - 588, 1996.

Examples
--------
These are written in doctest format, and should illustrate how to
use the function.

>>> a = [1, 2, 3]
>>> print([x + 3 for x in a])
[4, 5, 6]
>>> print("a\nb")
a
b
)pbdoc";


const std::string Numpydoc_docstring =
R"pbdoc(My numpydoc description of a kind
of very exhautive numpydoc format docstring.

Parameters
----------
first : array_like
    the 1st param name `first`
second :
    the 2nd param
    third : {'value', 'other'}, optional
    the 3rd param, by default 'value'

Returns
-------
string
    a value in a string

Raises
------
KeyError
    when a key error
OtherError
    when an other error
)pbdoc";


// Finalize closes any files open in the BMIPhreeqcRM instance.

const std::string finalize_docstring =
R"pbdoc(Closes any files open in the bmi_phreeqcrm instance.
)pbdoc";


const std::string initialize_docstring =
R"pbdoc(Initializes the bmi_phreeqcrm instance.

A YAML file used for initialization contains a YAML map of
bmi_phreeqcrm methods and the arguments corresponding to
each method.

Parameters
----------
config_file : str, optional
    The path to the model configuration file.
)pbdoc";


const std::string get_components_docstring =
R"pbdoc(Component list that was generated by calls to @FindComponents@.

Returns
-------
tuple(str)
    Each str is a component name.

Examples
--------

>>> import phreeqcrm as rm
>>> bmi = rm.bmi_phreeqcrm()
>>> bmi.initialize("AdvectBMI_py.yaml")
>>> components = bmi.get_components()
>>> print(components)
('H', 'O', 'Charge', 'Ca', 'Cl', 'K', 'N', 'Na')
)pbdoc";


const std::string get_component_name_docstring =
R"pbdoc(Name of the component

get_component_name returns the component name--"BMI PhreeqcRM".
BMI PhreeqcRM is a partial interface to PhreeqcRM, and provides 
the methods to implement chemical reactions in a multicomponent
transport model. All of the native PhreeqcRM methods (non BMI
methods) are also available, which provides a complete
interface to PhreeqcRM.

Returns
-------
str
    The string "BMI PhreeqcRM".

Examples
--------

>>> import phreeqcrm as rm
>>> bmi = rm.bmi_phreeqcrm()
>>> print(bmi.get_component_name())
BMI PhreeqcRM
)pbdoc";

const std::string get_grid_size_docstring =
R"pbdoc(Get the total number of cells defined.

get_grid_size returns the number of cells specified
by the YAML configuration file. If not specified a
default value of 10 is used.

Parameters
----------
grid : int
    Grid number, only grid 0 is considered.

Returns
-------
int
    Same value as GetGridCellCount is returned for grid 0;
    0 for all other values of `grid`.
)pbdoc";

const std::string update_docstring =
R"pbdoc(Advance model state by one time step.

Runs PhreeqcRM for one time step. This method is equivalent to
@ref RunCells. PhreeqcRM will equilibrate the solutions with all equilibrium 
reactants (EQUILIBRIUM_PHASES, EXCHANGE, GAS_PHASE, SOLID_SOLUTIONS, and SURFACE) 
and integrate KINETICS reactions for the specified time step 
(@ref SetValue "TimeStep" or @ref SetTimeStep).
)pbdoc";

const std::string update_until_docstring =
R"pbdoc(Advance model state until the given time.

Parameters
----------
end_time : float
    Time at the end of the time step. 

update_until is the same as @ref update, except the time step is calculated
from the argument @a end_time. The time step is calculated to be @a end_time minus 
the current time (@ref GetCurrentTime).
)pbdoc";


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Define a templated function that accepts an array-like variable of ints or doubles
template<typename T>
void my_function(py::array_t<T> input_array) {
    // Get a pointer to the data buffer
    auto ptr = static_cast<T*>(input_array.request().ptr);

    // Get the shape of the array
    auto shape = input_array.shape();

    // Get the number of dimensions
    auto ndim = input_array.ndim();

    // Print some information about the array
    std::cout << "Array shape: (";
    for (int i = 0; i < ndim; i++) {
        std::cout << shape[i];
        if (i < ndim - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;

    // Access the elements of the array
    for (int i = 0; i < shape[0]; i++) {
        for (int j = 0; j < shape[1]; j++) {
            auto value = *(ptr + i * shape[1] + j);
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}

//// Define the Python module
//PYBIND11_MODULE(example, m) {
//    m.def("my_function", &my_function<double>, "A function that accepts an array-like variable of doubles");
//    m.def("my_function", &my_function<int>, "A function that accepts an array-like variable of ints");
//}

//void my_func(py::array_t<py::str> arr) {
void my_func(py::list strings)
{
    std::vector<std::string> result;
    for (auto item : strings)
    {
        std::string s = py::cast<std::string>(item);
        result.push_back(py::cast<std::string>(item));
    }
    //return result;
}

py::sequence BMIPhreeqcRM::process_sequence(py::sequence seq)
{
    // Get the length of the sequence
    ssize_t size = py::len(seq);

    py::handle type = seq.get_type();
    std::string string = std::string(type.str());

    py::sequence return_val = seq;

    py::dtype int_dt = py::dtype::of<int>();
    py::dtype double_dt = py::dtype::of<double>();
    //py::dtype string_dt = py::dtype::of<std::string>();


    if (py::isinstance<py::list>(seq))
    {
        py::list list(size);
        //return_val = list;

        for (ssize_t i = 0; i < size; ++i)
        {
            // Access the elements of the sequence
            py::object element = seq[i];

            int modified_element = py::cast<int>(element) * 2;

            //list.append(modified_element);
            list[i] = modified_element;
        }
        return list;
    }
    else if (py::isinstance<py::tuple>(seq))
    {
        //py::tuple tuple(size);
        py::tuple tuple(0);
        return_val = tuple;

        for (ssize_t i = 0; i < size; ++i)
        {
            // Access the elements of the sequence
            py::object element = seq[i];

            int modified_element = py::cast<int>(element) * 2;

            //tuple[i] = modified_element;
            tuple = tuple + py::make_tuple(modified_element);
        }
        return tuple;
    }
    else if (py::isinstance<py::array>(seq))
    {
        //py::array a = seq.cast(py::array)();
        py::array arr = seq.cast<py::array>();
        py::dtype dt = arr.dtype();
        std::string dt_str = std::string(dt.str());

        py::object item = seq[0];

        py::handle type = item.get_type();
        std::string string = std::string(type.str());

        int int_element = py::cast<int>(item);

        double double_element = py::cast<double>(item);

        //if (py::isinstance<py::float_>(item))
        if (dt_str == "float64")
        {
            py::array_t<double> d_array(size);
            return_val = d_array;
        }
        //else if (py::isinstance<py::int_>(item))
        else if (dt_str == "int32")
        {
            py::array_t<int> i_array(size);
            return_val = i_array;
        }
    }
    else
    {
        // unsupported sequence
        return py::none();
    }

    //py::list mult_by_2;

    // Iterate over the sequence
    //for (ssize_t i = 0; i < size; ++i)
    //{
    //    // Access the elements of the sequence
    //    py::object element = seq[i];

    //    int modified_element = py::cast<int>(element) * 2;

    //    mult_by_2.append(modified_element);

    //    // Process the element
    //    // ...
    //    py::handle type = element.get_type();
    //    std::string string = std::string(type.str());
    //}
    //return mult_by_2;
    return return_val;
}

PYBIND11_MODULE(phreeqcrm, m) {

    m.def("my_function", &my_function<double>, "A function that accepts an array-like variable of doubles");
    m.def("my_function", &my_function<int>, "A function that accepts an array-like variable of ints");

    m.def("my_func",
        &my_func,
        "A function that takes a NumPy array of strings."
    );

    py::class_<BMIPhreeqcRM>(m, "bmi_phreeqcrm" , py::dynamic_attr())
        .def(py::init())

        //
        // testing
        //

        .def("process_sequence",
            &BMIPhreeqcRM::process_sequence,
            py::arg("seq")
        )

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

        //.def("get_value_at_indices",
        //    [](BMIPhreeqcRM& self, std::string name, py::array dest) {
        //        auto resultobj = self.get_value_ptr(name).attr("take");
        //        return resultobj();
        //    }
        //)
        .def("get_value_at_indices",
            &BMIPhreeqcRM::get_value_at_indices,
            py::arg("name"),
            py::arg("dest"),
            py::arg("inds")
        )


        //.def("my_sqrt",
        //    [](BMIPhreeqcRM& self, double value) {
        //        auto math = py::module::import("math");
        //        auto resultobj = math.attr("sqrt")(value);
        //        return resultobj.cast<double>();
        //    },
        //    py::arg("value")
        //)


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


        .def_readonly("_initialized", &BMIPhreeqcRM::_initialized)

        ;

    //m.def("get_value", &get_value, "Get a shared_ptr to a std::vector<double>");
    //m.def("get_value", &get_value<int>, "std::vector<int>");
    //m.def("get_value", &get_value<double>, "std::vector<double>");
    //m.def("get_value", &get_value, "");
    //m.def("get_value_ptr", &get_value_ptr, "");
    //m.def("use_mutable_unchecked", &use_mutable_unchecked, "");
}





// PYBIND11_MODULE(phreeqcrm, m) {
//   m.doc() = "I'm a docstring hehe";
//   m.def("some_fn_python_name", &some_fn);
//   py::class_<BMIPhreeqcRM>(
// 			m, "bmi_phreeqcrm"
// 			)
//       .def(py::init())
//       .def("initialize", &BMIPhreeqcRM::Initialize)
//       ;
// }

int add(int x, int y=0) {
    return x + y;
}

const std::string bmi_initialize_docstring = 
R"pbdoc(A YAML file can be used to initialize an instance of bmi_phreeqcrm.

Parameters
----------
filename : str, optional
    Path to name of input file.
)pbdoc";


const std::string bmi_set_value_docstring = 
R"pbdoc("Use the vector of concentrations (c) to set the moles of components in each 
reaction cell. The volume of water in a cell is the product of porosity 
(:meth:`SetPorosity`), saturation (:meth:`SetSaturation`), and reference volume 
(:meth:`SetRepresentativeVolume`). The moles of each component are determined 
by the volume of water and per liter concentrations. If concentration units 
(:meth:`SetUnitsSolution`) are mass fraction, the density (as specified by 
:meth:`SetDensity`) is used to convert from mass fraction to per mass per liter.

:param c: a list or numpy array of doubles.

Args:
	c (double list, numpy.ndarray, or tuple): Component concentrations. Size of 
		vector is ncomps times nxyz, where ncomps is the number of components as 
		determined by :meth:`FindComponents` or :meth:`GetComponentCount` and nxyz 
		is the number of grid cells in the user model (:meth:`GetGridCellCount`).
)pbdoc";


#include <stdexcept>

// https://github.com/Deltares/xmipy/blob/f3954fd888e42d24920adc58543185b52f80a6c7/xmipy/xmiwrapper.py#L313

//py::array BMIPhreeqcRM::get_value_test(std::string name, py::array dest/* = py::none()*/)
////py::array BMIPhreeqcRM::get_value_test(std::string name, py::array_t<double> dest = py::none())
//{
//    //{{
//    // make sure that optional array is of correct layout:
//    if (!dest.is_none())
//    {
//        if (dest.flags() & py::array::f_style)
//        {
//            throw std::runtime_error("Array should have C layout");
//        }
//    }
//
//    RMVARS v_enum = this->var_man->GetEnum(name);
//    return py::none();
//    //return dest;
//
//
//    //if (v_enum != RMVARS::NotFound)
//    //{
//    //    // first deal with scalars
//    //    BMIVariant& bv = this->var_man->VariantMap[v_enum];
//
//    //    int dims = this->var_man->VarExchange.GetDim();
//
//    //    if (dims > 0)
//    //    {
//    //        py::array out = py::cast(this->var_man->VarExchange.GetDoubleVectorRef());
//    //        return out;
//    //    }
//    //    else
//    //    {
//    //        dest = this->var_man->VarExchange.GetDVar();
//    //        double output = 3.14;
//    //        return py::cast();
//    //    }
//
//    //}
//
//
//    ////}}
//
//
//    //RMVARS v_enum = this->var_man->GetEnum(name);
//    //if (v_enum != RMVARS::NotFound)
//    //{
//    //    BMIVariant& bv = this->var_man->VariantMap[v_enum];
//    //    //VarManager::VarFunction fn = this->var_man->GetFn(v_enum);
//    //    if (!bv.GetInitialized())
//    //    {
//    //        this->var_man->task = VarManager::VAR_TASKS::Info;
//    //        ((*this->var_man).*bv.GetFn())();
//    //    }
//    //    this->var_man->task = VarManager::VAR_TASKS::GetVar;
//    //    ((*this->var_man).*bv.GetFn())();
//    //    
//    //    int dims = this->var_man->VarExchange.GetDim();
//
//    //    if (this->var_man->VarExchange.GetCType() == "double")
//    //    {
//    //        if (dims > 0)
//    //        {
//    //            py::array out = py::cast(this->var_man->VarExchange.GetDoubleVectorRef());
//    //            return out;
//    //        }
//    //        else
//    //        {
//    //            dest = this->var_man->VarExchange.GetDVar();
//    //            double output = 3.14;
//    //            return py::cast();
//    //        }
//    //    }
//
//    //    //dest = this->var_man->VarExchange.GetDoubleVectorRef();
//    //    //return;
//    //}
//    //return py::none();
//
//
//    ////if (name.compare("Pi") == 0)
//    ////{
//    ////    double output = 3.14;
//    ////    return py::cast(output);
//    ////}
//    ////else
//    ////{
//    ////    std::vector<double> output(20, 3.14);
//    ////    py::array out = py::cast(output);
//    ////    return out;
//    ////}
//}

//#include <optional>
//
//void echo(std::string s)
//{
//    std::cout << s << std::endl;
//}
//
////int my_len(py::array_t<double> a)
////{
////    if (a.is_none()) return 0;
////
////    py::buffer_info buf = a.request();
////    return buf.shape[0];
////}
//
//int my_len(py::array a)
//{
//    if (a.is_none()) return 0;
//
//    py::buffer_info buf = a.request();
//    return buf.shape[0];
//}
//
//
////using array_i = pybind11::array_t<int, pybind11::array::c_style | pybind11::array::forcecast>;
//using array_i = pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>;
//
//struct myObj_array {
//    array_i ids;
//};
//




//PYBIND11_MODULE(phreeqcrm, m) {
//
//    py::class_<myObj_array>(m, "myObj_array_int")
//        .def(py::init([](array_i py_ids) {
//        std::cout << "ndim " << py_ids.ndim() << std::endl;
//        if (py_ids.ndim())
//            return myObj_array{ py_ids };
//        else
//            return myObj_array{ array_i(0) };
//            }), py::arg("py_ids") = (array_i*) nullptr
//                );
//
//
//
//    //m.def("add", &add, "A function that adds two integers.", py::arg("x"), py::arg("y")=0);
//
//    m.def("echo", &echo, "A function that outputs a string.", py::arg("s") = std::string("Hello"));
//
//    m.def("echo2", [](py::str s) {
//        if (s.is_none()) std::cout << "empty";
//        else std::cout << s;
//        }, "A function that outputs a string.", py::arg("s") = static_cast<py::str>(nullptr);
//
//
//    m.def("echo3", [](py::object s) {
//        if (s.is_none()) std::cout << "empty";
//        else std::cout << s;
//        }, "A function that outputs a string.", py::arg("s") = py::none());
//
//
//    m.def("my_len", &my_len, "A function that outputs a string.", py::arg("a") = py::none());
//
//    m.def("my_len2",
//        [](std::optional<py::array> a) {
//            if (!a.has_value()) {
//                return 0;
//            }
//            return 10;
//        },
//        py::arg("a") = py::none()
//            );
//
//
//
//    m.doc() = R"pbdoc(
//        PhreeqcRM plugin
//        -----------------------
//
//        .. currentmodule:: phreeqcrm
//
//        .. autosummary::
//           :toctree: _generate
//
//    )pbdoc";
//
//    py::class_<BMIPhreeqcRM>(
//			m, "bmi_phreeqcrm" /* , py::dynamic_attr()*/
//    )
//    .def(py::init())
//
//    // Model control functions.
//    // .def("initialize", &BMIPhreeqcRM::Initialize, R"pbdoc(
//    //     A YAML file can be used to initialize an instance of bmi_phreeqcrm.
//
//    //     Parameters
//    //     ----------
//	//     filename : str, optional
//    //         Path to name of input file.
//    // )pbdoc", py::arg("filename"))
//
//    // .def("initialize2", [](BMIPhreeqcRM &self, py::str filename) {
//    //     if (filename.is_none()) {
//    //         return self.Initialize("");
//    //     }
//    //     return self.Initialize(std::string(py::str(filename)));
//    // }, bmi_initialize_docstring2.c_str(), py::arg("filename") = py::none())
//
//    .def("initialize", [](BMIPhreeqcRM &self, py::str filename) {
//        if (filename.is_none()) {
//            return self.Initialize("");
//        }
//        return self.Initialize(std::string(py::str(filename)));
//    }, bmi_initialize_docstring.c_str(), py::arg("filename") = py::none())
//
//
//
//    // Model information functions.
//    .def("get_component_name", &BMIPhreeqcRM::GetComponentName)
//
//    // .def("get_input_var_names", &BMIPhreeqcRM::GetInputVarNames)
//    .def("get_input_var_names", [](BMIPhreeqcRM &self) {
//        auto output = self.GetInputVarNames();
//        py::tuple out = py::cast(output);
//        return out;
//    })
//
//    .def("get_output_var_names", [](BMIPhreeqcRM &self) {
//        auto output = self.GetOutputVarNames();
//        py::tuple out = py::cast(output);
//        return out;
//    })
//
//
//
//    // Variable getters
//    //   virtual void GetValue(std::string name, void *dest) = 0;
//    .def("get_value", [](BMIPhreeqcRM &self, const std::string &input) {
//        std::vector<double> output;
//        self.GetValue(input, output);
//        py::array out = py::cast(output);
//        return out;
//        })
//        //   virtual void *GetValuePtr(std::string name) = 0;
//        //   virtual void GetValueAtIndices(std::string name, void *dest, int *inds, int count) = 0;
//
//    //.def("get_value_test", &BMIPhreeqcRM::get_value_test)
//    //.def("get_value_test", &BMIPhreeqcRM::get_value_test, py::arg("name"), py::arg("dest") = py::none())
//
//    .def("get_value_test", [](BMIPhreeqcRM& self, const std::string& name, py::array dest) {
//            std::vector<double> output;
//            //self.GetValue(input, output);
//            //py::array out = py::cast(dest);
//            //return out;
//            return dest;
//        }, py::arg("name"), py::arg("dest") = py::none())
//
//
//
//    // Variable setters
//    //   virtual void SetValue(std::string name, void *src) = 0;
//    //.def("set_value", py::overload_cast<std::string, std::vector<double>>(&BMIPhreeqcRM::SetValue), py::arg("name", std::string()), py::arg("src", std::vector<double>()))
//    // .def("set_value", py::overload_cast<std::string, std::vector<double>>(&BMIPhreeqcRM::SetValue), py::arg("name"), py::arg("src").noconvert(), py::keep_alive<1, 2>(),
//    //     bmi_set_value_docstring.c_str()
//    // )
//
//
//    .def("set_value", [](BMIPhreeqcRM &self, std::string name, const py::array_t<double>& src) {
//        // Extract the buffer information from the numpy array.
//        py::buffer_info info = src.request();
//
//        // Construct a std::vector<double> from the buffer data.
//        double* ptr = static_cast<double*>(info.ptr);
//        size_t size = info.shape[0];
//        std::vector<double> src_vec(ptr, ptr + size);
//
//        self.SetValue(name, src_vec);
//    },
//        py::arg("name"),
//        //py::arg("src").noconvert(),       // this will fail if a List is used
//        py::arg("src"),
//        bmi_set_value_docstring.c_str()
//    )
//
//
//    .def("set_value", [](BMIPhreeqcRM &self, std::string name, const py::array_t<int>& src) {
//        // Extract the buffer information from the numpy array.
//        py::buffer_info info = src.request();
//
//        // Construct a std::vector<double> from the buffer data.
//        int* ptr = static_cast<int*>(info.ptr);
//        size_t size = info.shape[0];
//        std::vector<int> src_vec(ptr, ptr + size);
//
//        self.SetValue(name, src_vec);
//    },
//        py::arg("name"),
//        //py::arg("src").noconvert(),       // this will fail if a List is used
//        py::arg("src"),
//        bmi_set_value_docstring.c_str()
//    )
//
//
//
//
//    // .def("set_value", [](BMIPhreeqcRM &self, const std::string &input, std::vector<double> values) {
//    //     self.GetValue(input, output);
//    //     py::array out = py::cast(output);
//    //     return out;
//    // })
//
//
//    //   virtual void SetValueAtIndices(std::string name, int *inds, int count, void *src) = 0;    
//
//
//    // Grid information functions    
//    .def("get_grid_size", &BMIPhreeqcRM::GetGridSize)
//    //   virtual int GetGridRank(const int grid) = 0;
//    //   virtual int GetGridSize(const int grid) = 0;
//    //   virtual std::string GetGridType(const int grid) = 0;
//
//    //   virtual void GetGridShape(const int grid, int *shape) = 0;
//    //   virtual void GetGridSpacing(const int grid, double *spacing) = 0;
//    //   virtual void GetGridOrigin(const int grid, double *origin) = 0;
//
//    //   virtual void GetGridX(const int grid, double *x) = 0;
//    //   virtual void GetGridY(const int grid, double *y) = 0;
//    //   virtual void GetGridZ(const int grid, double *z) = 0;
//
//    //   virtual int GetGridNodeCount(const int grid) = 0;
//    //   virtual int GetGridEdgeCount(const int grid) = 0;
//    //   virtual int GetGridFaceCount(const int grid) = 0;
//
//    //   virtual void GetGridEdgeNodes(const int grid, int *edge_nodes) = 0;
//    //   virtual void GetGridFaceEdges(const int grid, int *face_edges) = 0;
//    //   virtual void GetGridFaceNodes(const int grid, int *face_nodes) = 0;
//    //   virtual void GetGridNodesPerFace(const int grid, int *nodes_per_face) = 0;
//
//
//    // Extras
//    .def("get_components", &BMIPhreeqcRM::GetComponents)
//
//    ;
//
//
//
//#ifdef VERSION_INFO
//    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
//#else
//    m.attr("__version__") = "dev";
//#endif
//}



// #include <pybind11/pybind11.h>

// #define STRINGIFY(x) #x
// #define MACRO_STRINGIFY(x) STRINGIFY(x)

// int add(int i, int j) {
//     return i + j;
// }

// namespace py = pybind11;

// PYBIND11_MODULE(cmake_example, m) {
//     m.doc() = R"pbdoc(
//         Pybind11 example plugin
//         -----------------------

//         .. currentmodule:: cmake_example

//         .. autosummary::
//            :toctree: _generate

//            add
//            subtract
//     )pbdoc";

//     m.def("add", &add, R"pbdoc(
//         Add two numbers

//         Some other explanation about the add function.
//     )pbdoc");

//     m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
//         Subtract two numbers

//         Some other explanation about the subtract function.
//     )pbdoc");

// #ifdef VERSION_INFO
//     m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
// #else
//     m.attr("__version__") = "dev";
// #endif
// }
