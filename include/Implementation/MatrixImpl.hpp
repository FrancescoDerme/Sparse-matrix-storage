#ifndef MATRIX_IMPL_HPP
#define MATRIX_IMPL_HPP

#include <cassert>
#include <memory>

#include "Dimensions.hpp"
#include "Matrix.hpp"

using namespace comparators;
namespace algebra {

/// @brief Constructs a matrix with the given arguments.
/// @tparam Args Variadic template for constructor arguments.
/// @param args The arguments to initialize the matrix.
/// @details This constructor initializes the dynamic storage format with the
/// provided arguments and sets the matrix to uncompressed by default.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(Args&&... args)
    : Compressed<T, S>{},
      Dynamic<T, S>{std::false_type{}, std::forward<Args>(args)...} {
    isCompressed = false;
}

/// @brief Constructs a matrix with specified rows and columns.
/// @tparam Args Variadic template for additional arguments.
/// @param r The number of rows.
/// @param c The number of columns.
/// @param args Additional arguments to initialize the matrix.
/// @details This constructor initializes the dimensions of the matrix and
/// sets up the dynamic storage format with the provided arguments.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(ConvertibleToSizeT auto const r,
                    ConvertibleToSizeT auto const c, Args&&... args)
    : Dimensions(r, c),
      Compressed<T, S>{},
      Dynamic<T, S>{std::true_type{}, std::forward<Args>(args)...} {
    isCompressed = false;
}

/// @brief Constructs a compressed matrix with the given arguments.
/// @tparam Args Variadic template for constructor arguments.
/// @param tag A tag to indicate compressed storage.
/// @param args The arguments to initialize the matrix.
/// @details This constructor initializes the compressed storage format and
/// sets the matrix to compressed.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(UseCompressed, Args&&... args)
    : Compressed<T, S>{std::false_type{}, std::forward<Args>(args)...},
      Dynamic<T, S>{} {
    isCompressed = true;
}

/// @brief Constructs a dynamic matrix with the given arguments.
/// @tparam Args Variadic template for constructor arguments.
/// @param tag A tag to indicate dynamic storage.
/// @param args The arguments to initialize the matrix.
/// @details This constructor initializes the dynamic storage format and
/// sets the matrix to uncompressed.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(UseDynamic, Args&&... args)
    : Compressed<T, S>{},
      Dynamic<T, S>{std::false_type{}, std::forward<Args>(args)...} {
    isCompressed = false;
}

/// @brief Constructs a compressed matrix with specified rows and columns.
/// @tparam Args Variadic template for additional arguments.
/// @param tag A tag to indicate compressed storage.
/// @param r The number of rows.
/// @param c The number of columns.
/// @param args Additional arguments to initialize the matrix.
/// @details This constructor initializes the dimensions of the matrix and
/// sets up the compressed storage format with the provided arguments.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(UseCompressed, ConvertibleToSizeT auto const r,
                    ConvertibleToSizeT auto const c, Args&&... args)
    : Dimensions(r, c),
      Compressed<T, S>{std::true_type{}, std::forward<Args>(args)...},
      Dynamic<T, S>{} {
    isCompressed = true;
}

/// @brief Constructs a dynamic matrix with specified rows and columns.
/// @tparam Args Variadic template for additional arguments.
/// @param tag A tag to indicate dynamic storage.
/// @param r The number of rows.
/// @param c The number of columns.
/// @param args Additional arguments to initialize the matrix.
/// @details This constructor initializes the dimensions of the matrix and
/// sets up the dynamic storage format with the provided arguments.
MATRIX_TEMPLATE
template <typename... Args>
MATRIX_TYPE::Matrix(UseDynamic, ConvertibleToSizeT auto const r,
                    ConvertibleToSizeT auto const c, Args&&... args)
    : Dimensions(r, c),
      Compressed<T, S>{},
      Dynamic<T, S>{std::true_type{}, std::forward<Args>(args)...} {
    isCompressed = false;
}

/// @brief Constructs a matrix by reading data from a file.
/// @param file_name The name of the file to read the matrix data from.
/// @details This constructor initializes the dynamic storage format by
/// parsing the contents of the specified file.
MATRIX_TEMPLATE
MATRIX_TYPE::Matrix(std::string& file_name)
    : Compressed<T, S>{}, Dynamic<T, S>{file_name} {}

/// @brief Accesses the element at the specified position (read-only).
/// @param i The row index.
/// @param j The column index.
/// @return The value at the specified position.
/// @details This function checks whether the matrix is compressed or dynamic
/// and retrieves the value at the specified position accordingly.
MATRIX_TEMPLATE
T MATRIX_TYPE::operator()(std::size_t i, std::size_t j) const {
#ifdef DEBUG
    assert(i < this->rows && j < this->columns && i >= 0 && j >= 0 &&
           "Error in call to operator(): indexes out of bounds.\n");
#endif
    if (!isCompressed) {
        return this->find_dynamic_const(i, j);
    }
    else {
        return this->find_compressed_const(i, j);
    }
}

/// @brief Accesses the element at the specified position (read-write).
/// @param i The row index.
/// @param j The column index.
/// @return A reference to the value at the specified position.
/// @details This function checks whether the matrix is compressed or dynamic
/// and retrieves a reference to the value at the specified position
/// accordingly.
MATRIX_TEMPLATE
T& MATRIX_TYPE::operator()(std::size_t i, std::size_t j) {
#ifdef DEBUG
    assert(i < this->rows && j < this->columns && i >= 0 && j >= 0 &&
           "Error in call to operator(): indexes out of bounds.\n");
#endif
    if (!isCompressed) {
        return this->find_dynamic(i, j);
    }
    else {
        return this->find_compressed(i, j);
    }
}

/// @brief Checks if the matrix is in compressed storage format.
/// @return True if the matrix is compressed, false otherwise.
/// @details This function returns the value of the `isCompressed` flag.
MATRIX_TEMPLATE
bool MATRIX_TYPE::is_compressed() const { return isCompressed; }

/// @brief Gets the number of non-zero elements in the matrix.
/// @return The number of non-zero elements.
/// @details This function checks whether the matrix is compressed or dynamic
/// and retrieves the number of non-zero elements accordingly.
MATRIX_TEMPLATE
size_t MATRIX_TYPE::get_num_elements() const {
    if (!isCompressed) {
        return (this)->get_num_elements_dynamic();
    }
    else {
        return (this)->get_num_elements_compressed();
    }
}

/// @brief Compresses the matrix into a compact storage format.
/// @details This function converts the matrix from dynamic to compressed
/// storage format by iterating through its elements and transferring them
/// into the compressed structure.
MATRIX_TEMPLATE
void MATRIX_TYPE::compress() {
#ifdef DEBUG
    assert(isCompressed == false &&
           "Error in call to compress method: matrix already compressed.\n");
#endif

    size_t num_elements = get_num_elements();

    if constexpr (S == rowMajor) {
        this->initialize_compressed(this->get_rows() + 1, num_elements,
                                    num_elements);
    }
    else {
        this->initialize_compressed(this->get_columns() + 1, num_elements,
                                    num_elements);
    }

    auto tupleptr = std::make_unique<std::tuple<size_t, size_t, T>>();
    for (size_t i = 0; i < num_elements; ++i) {
        (this)->compress_from_dynamic(i == 0, tupleptr);
        (this)->compress_from_triplets(i == num_elements - 1, tupleptr);
    }

    this->release_dynamic();
    isCompressed = true;
}

/// @brief Uncompresses the matrix into a dynamic storage format.
/// @details This function converts the matrix from compressed to dynamic
/// storage format by iterating through its elements and transferring them
/// into the dynamic structure.
MATRIX_TEMPLATE
void MATRIX_TYPE::uncompress() {
#ifdef DEBUG
    assert(
        isCompressed == true &&
        "Error in call to uncompress method: matrix already uncompressed.\n");
#endif

    size_t num_elements = get_num_elements();

    this->initialize_dynamic();

    auto tupleptr = std::make_unique<std::tuple<size_t, size_t, T>>();
    for (size_t i = 0; i < num_elements; ++i) {
        (this)->uncompress_from_compressed(i == 0, tupleptr);
        (this)->uncompress_from_triplets(i == 0, tupleptr);
    }

    this->release_compressed();
    isCompressed = false;
}

/// @brief Computes the norm of the matrix.
/// @tparam N The type of norm to compute (infinity, one, or Frobenius).
/// @return The computed norm value.
/// @details This function checks whether the matrix is compressed or dynamic
/// and computes the specified norm accordingly.
MATRIX_TEMPLATE
template <NormType N>
double MATRIX_TYPE::norm() const {
    if (!isCompressed) {
        return (this)->template norm_dynamic<N>();
    }
    else {
        return (this)->template norm_compressed<N>();
    }
}

/// @brief Removes the element at the specified position.
/// @param i The row index.
/// @param j The column index.
/// @return True if the element was removed, false otherwise.
/// @details This function checks whether the matrix is compressed or dynamic
/// and removes the specified element accordingly.
MATRIX_TEMPLATE
bool MATRIX_TYPE::remove(size_t i, size_t j) {
#ifdef DEBUG
    assert(i < this->rows && j < this->columns && i >= 0 && j >= 0 &&
           "Error in call method remove: indexes out of bounds.\n");
#endif
    if (!isCompressed) {
        return (this)->remove_dynamic(i, j);
    }
    else {
        return (this)->remove_compressed(i, j);
    }
}

/// @brief Prints the matrix to the standard output.
/// @details This function checks whether the matrix is compressed or dynamic
/// and prints its elements accordingly.
MATRIX_TEMPLATE
void MATRIX_TYPE::print() const {
    if (!isCompressed) {
        (this)->print_dynamic();
        return;
    }
    else {
        (this)->print_compressed();
        return;
    }
}

/// @brief Performs matrix-vector product.
/// @param m An object of type Matrix representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
MATRIX_TEMPLATE
std::vector<T> operator*(MATRIX_TYPE const& m, std::vector<T> const& v) {
#ifdef DEBUG
    assert(m.columns == v.size() &&
           "Error in call operator *: non-matching dimensions.\n");
#endif

    if (!m.isCompressed) {
        return by_vector_dynamic(static_cast<const Dynamic<T, S>&>(m), v);
    }
    else {
        return by_vector_compressed(static_cast<const Compressed<T, S>&>(m), v);
    }
}

}  // namespace algebra
#endif
