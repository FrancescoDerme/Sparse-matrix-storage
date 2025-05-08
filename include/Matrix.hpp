#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Comparators.hpp"
#include "Concepts.hpp"

#define MATRIX_TEMPLATE                                           \
    template <NumericOrComplex T,                                 \
              template <typename, StorageOrder> class Compressed, \
              template <typename, StorageOrder> class Dynamic, StorageOrder S>

#define MATRIX_TYPE Matrix<T, Compressed, Dynamic, S>

using namespace comparators;
namespace algebra {

/// @brief Type to flag dynamic building.
struct UseDynamic {};

/// @brief Type to flag compressed building.
struct UseCompressed {};

MATRIX_TEMPLATE
class Matrix;

MATRIX_TEMPLATE
std::vector<T> operator*(MATRIX_TYPE const& m, std::vector<T> const& v);

/// @brief Represents a matrix that supports both compressed and dynamic storage
/// formats.
/// @tparam T The type of the matrix elements (e.g., numeric or complex).
/// @tparam Compressed The class template for compressed storage.
/// @tparam Dynamic The class template for dynamic storage.
/// @tparam S The storage order (row-major or column-major).
MATRIX_TEMPLATE
class Matrix : public Compressed<T, S>, public Dynamic<T, S> {
   public:
    /// @brief Constructs a matrix with the given arguments.
    /// @tparam Args Variadic template for constructor arguments.
    /// @param args The arguments to initialize the matrix.
    template <typename... Args>
    Matrix(Args&&... args);

    /// @brief Constructs a matrix with specified rows and columns.
    /// @tparam Args Variadic template for additional arguments.
    /// @param r The number of rows.
    /// @param c The number of columns.
    /// @param args Additional arguments to initialize the matrix.
    template <typename... Args>
    Matrix(ConvertibleToSizeT auto const r, ConvertibleToSizeT auto const c,
           Args&&... args);

    /// @brief Constructs a compressed matrix with the given arguments.
    /// @tparam Args Variadic template for constructor arguments.
    /// @param tag A tag to indicate compressed storage.
    /// @param args The arguments to initialize the matrix.
    template <typename... Args>
    Matrix(UseCompressed, Args&&... args);

    /// @brief Constructs a dynamic matrix with the given arguments.
    /// @tparam Args Variadic template for constructor arguments.
    /// @param tag A tag to indicate dynamic storage.
    /// @param args The arguments to initialize the matrix.
    template <typename... Args>
    Matrix(UseDynamic, Args&&... args);

    /// @brief Constructs a compressed matrix with specified rows and columns.
    /// @tparam Args Variadic template for additional arguments.
    /// @param tag A tag to indicate compressed storage.
    /// @param r The number of rows.
    /// @param c The number of columns.
    /// @param args Additional arguments to initialize the matrix.
    template <typename... Args>
    Matrix(UseCompressed, ConvertibleToSizeT auto const r,
           ConvertibleToSizeT auto const c, Args&&... args);

    /// @brief Constructs a dynamic matrix with specified rows and columns.
    /// @tparam Args Variadic template for additional arguments.
    /// @param tag A tag to indicate dynamic storage.
    /// @param r The number of rows.
    /// @param c The number of columns.
    /// @param args Additional arguments to initialize the matrix.
    template <typename... Args>
    Matrix(UseDynamic, ConvertibleToSizeT auto const r,
           ConvertibleToSizeT auto const c, Args&&... args);

    /// @brief Constructs a matrix by reading data from a file.
    /// @param file_name The name of the file to read the matrix data from.
    Matrix(std::string& file_name);

    /// @brief Accesses the element at the specified position (read-only).
    /// @param i The row index.
    /// @param j The column index.
    /// @return The value at the specified position.
    T operator()(std::size_t i, std::size_t j) const;

    /// @brief Accesses the element at the specified position (read-write).
    /// @param i The row index.
    /// @param j The column index.
    /// @return A reference to the value at the specified position.
    T& operator()(std::size_t i, std::size_t j);

    /// @brief Checks if the matrix is in compressed storage format.
    /// @return True if the matrix is compressed, false otherwise.
    bool is_compressed() const;

    /// @brief Gets the number of non-zero elements in the matrix.
    /// @return The number of non-zero elements.
    size_t get_num_elements() const;

    /// @brief Computes the norm of the matrix.
    /// @tparam N The type of norm to compute (e.g., Infinity, One, or
    /// Frobenius).
    /// @return The computed norm value.
    template <NormType N>
    double norm() const;

    friend std::vector<T> operator*
        <>(MATRIX_TYPE const&, std::vector<T> const&);
    ///
    /// @brief Removes the element at the specified position.
    /// @param i The row index.
    /// @param j The column index.
    /// @return True if the element was removed, false otherwise.
    bool remove(size_t i, size_t j);

    /// @brief Compresses the matrix into a compact storage format.
    void compress();
    /// @brief Uncompresses the matrix into a dynamic storage format.
    void uncompress();

    /// @brief Prints the matrix to the standard output.
    void print() const;

   private:
    /// @brief Indicates whether the matrix is in compressed storage format.
    bool isCompressed;
};

}  // namespace algebra

#endif
