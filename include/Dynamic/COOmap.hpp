#ifndef COOMAP_HPP
#define COOMAP_HPP

#include <map>
#include <memory>

#include "Comparators.hpp"
#include "Concepts.hpp"
#include "Dimensions.hpp"

using namespace comparators;
namespace algebra {

/// @brief Represents a matrix in Coordinate Map (COOmap) format.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
class COOmap;

/// @brief Performs matrix-vector product.
/// @param m An object of type COOmap representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
std::vector<T> by_vector_dynamic(class COOmap<T, S> const& m,
                                 std::vector<T> const& v);

template <NumericOrComplex T, StorageOrder S>
class COOmap : virtual public Dimensions {
    using valuesmap = std::map<std::pair<size_t, size_t>, T,
                               Comparator<S>>;  ///< Map of index-value pairs.
    using ivtuple = std::tuple<size_t, size_t,
                               T>;  ///< Tuple representing a matrix element.

   protected:
    /// @brief Default constructor for the COOmap class.
    COOmap() = default;

    /// @brief Constructs a COOmap matrix from index-value pairs.
    /// @tparam B Boolean constant to indicate wether the matrix' size was given
    /// as input.
    /// @param indexes The container of index pairs.
    /// @param values The container of matrix values.
    template <bool B>
    COOmap(std::bool_constant<B>, SizetPairContainer auto const& indexes,
           NumericContainer auto const& values);

    /// @brief Constructs a COOmap matrix from a map of index-value pairs.
    /// @tparam B Boolean constant to indicate wether the matrix' size was given
    /// as input.
    /// @param m The map of index-value pairs.
    template <bool B>
    COOmap(std::bool_constant<B>,
           std::map<std::pair<size_t, size_t>, T> const& m);

    /// @brief Constructs a COOmap matrix by reading data from a file.
    /// @param file_name The name of the file to read the matrix data from.
    COOmap(std::string& file_name);

    /// @brief In the process of compressing the matrix sends the next tuple.
    /// @param first_it Flags the first iteration.
    /// @param tupleptr Pointer to a tuple representing a matrix element.
    void compress_from_dynamic(bool first_it,
                               std::unique_ptr<ivtuple> const& tupleptr);

    /// @brief Takes a triplet and inserts it into the data structure.
    /// @param first_it Flags the first iteration.
    /// @param tupleptr Pointer to a tuple representing a matrix element.
    void uncompress_from_triplets(bool first_it,
                                  std::unique_ptr<ivtuple> const& tupleptr);

    /// @brief Gets the number of non-zero elements.
    /// @return The number of non-zero elements.
    size_t get_num_elements_dynamic() const;

    /// @brief Initializes the dynamic storage format.
    void initialize_dynamic();

    /// @brief Releases the dynamic storage format.
    void release_dynamic();

    std::unique_ptr<valuesmap>
        matrixptr;  ///< Pointer to the map of index-value pairs.

   public:
    /// @brief Finds the value at the specified position (read-only).
    /// @param i The row index.
    /// @param j The column index.
    /// @return The value at the specified position.
    T find_dynamic_const(size_t i, size_t j) const;

    /// @brief Finds the value at the specified position (read-write).
    /// @param i The row index.
    /// @param j The column index.
    /// @return A reference to the value at the specified position.
    T& find_dynamic(size_t i, size_t j);

    /// @brief Removes the element at the specified position.
    /// @param i The row index.
    /// @param j The column index.
    /// @return True if the element was removed, false otherwise.
    bool remove_dynamic(size_t i, size_t j);

    /// @brief Prints the matrix in dynamic format to the standard output.
    void print_dynamic() const;

    /// @brief Computes the norm of the matrix.
    /// @tparam N The type of norm to compute (Infinity, One, or Frobenius).
    /// @return The computed norm value.
    template <NormType N>
    double norm_dynamic() const;

    friend std::vector<T> by_vector_dynamic<>(COOmap<T, S> const& m,
                                              std::vector<T> const& v);
};

}  // namespace algebra

#endif
