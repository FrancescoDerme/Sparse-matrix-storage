#ifndef YALE_HPP
#define YALE_HPP

#include <memory>
#include <vector>

#include "Comparators.hpp"
#include "Concepts.hpp"
#include "Dimensions.hpp"

using namespace comparators;
namespace algebra {

/// @brief Represents a matrix in YALE (compressed) format.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
class YALE;

/// @brief Performs matrix-vector product.
/// @param m An object of type YALE representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
std::vector<T> by_vector_compressed(class YALE<T, S> const& m,
                                    std::vector<T> const& v);

template <NumericOrComplex T, StorageOrder S>
class YALE : virtual public Dimensions {
    using indexvec = std::vector<size_t>;  ///< Vector of indices.
    using valuesvec = std::vector<T>;      ///< Vector of matrix values.
    using ivtuple = std::tuple<size_t, size_t,
                               T>;  ///< Tuple representing a matrix element.

   protected:
    /// @brief Default constructor for the YALE class.
    YALE() = default;

    /// @brief Constructs a YALE matrix from compressed format data.
    /// @tparam B Boolean constant to indicate wether the matrix' size was given
    /// as input.
    /// @param out The outer index container.
    /// @param in The inner index container.
    /// @param val The container of matrix values.
    template <bool B>
    YALE(std::bool_constant<B>, SizetContainer auto const& out,
         SizetContainer auto const& in, NumericContainer auto const& val);

    /// @brief In the process of compressing the matrix sends the next tuple.
    /// @param first_it Flags the first iteration.
    /// @param tupleptr Pointer to a tuple representing a matrix element.
    void uncompress_from_compressed(bool first_it,
                                    std::unique_ptr<ivtuple> const& tupleptr);

    /// @brief Takes a triplet and inserts it into the data structure.
    /// @param last_it Flags the last iteration.
    /// @param tupleptr Pointer to a tuple representing a matrix element.
    void compress_from_triplets(bool last_it,
                                std::unique_ptr<ivtuple> const& tupleptr);

    /// @brief Gets the number of non-zero elements.
    /// @return The number of non-zero elements.
    size_t get_num_elements_compressed() const;

    /// @brief Initializes the compressed storage format.
    /// @param insize The size of the inner index vector.
    /// @param outsize The size of the outer index vector.
    /// @param valsize The size of the values vector.
    void initialize_compressed(size_t insize, size_t outsize, size_t valsize);

    /// @brief Releases the compressed storage format.
    void release_compressed();

    /// @brief Computes the inner and outer indexes for a given position.
    /// @param i The row index.
    /// @param j The column index.
    /// @return A pair representing the inner and outer indexes.
    std::pair<size_t, size_t> inner_outer(size_t i, size_t j) const;

    std::unique_ptr<indexvec>
        outerindex_ptr;  ///< Pointer to the outer index vector.
    std::unique_ptr<indexvec>
        innerindex_ptr;  ///< Pointer to the inner index vector.
    std::unique_ptr<valuesvec> values_ptr;  ///< Pointer to the values vector.
    Comparator<S> comparator;  ///< Comparator coherent with the storage order.

   public:
    /// @brief Finds the value at the specified position (read-only).
    /// @param i The row index.
    /// @param j The column index.
    /// @return The value at the specified position.
    T find_compressed_const(size_t i, size_t j) const;

    /// @brief Finds the value at the specified position (read-write).
    /// @param i The row index.
    /// @param j The column index.
    /// @return A reference to the value at the specified position.
    T& find_compressed(size_t i, size_t j);

    /// @brief Removes the element at the specified position.
    /// @param i The row index.
    /// @param j The column index.
    /// @return True if the element was removed, false otherwise.
    bool remove_compressed(size_t i, size_t j);

    /// @brief Prints the matrix in compressed format to the standard output.
    void print_compressed() const;

    /// @brief Computes the norm of the matrix.
    /// @tparam N The type of norm to compute (Infinity, One, or Frobenius).
    /// @return The computed norm value.
    template <NormType N>
    double norm_compressed() const;

    friend std::vector<T> by_vector_compressed<>(YALE<T, S> const& m,
                                                 std::vector<T> const& v);
};

}  // namespace algebra
#endif
