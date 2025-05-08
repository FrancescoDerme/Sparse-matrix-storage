#ifndef COOMAPIMPL_HPP
#define COOMAPIMPL_HPP

#include <cassert>
#include <fstream>
#include <iostream>

#include "COOmap.hpp"
#include "Comparators.hpp"
#include "Concepts.hpp"

using namespace comparators;
namespace algebra {

/// @brief Constructs a COOmap matrix from index-value pairs.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam B Boolean constant to indicate wether the matrix' size was given
/// as input.
/// @param indexes The container of index pairs.
/// @param values The container of matrix values.
/// @details This constructor iterates through the provided index-value pairs,
/// validates them and inserts them into the internal map. If bounds checking is
/// disabled, it also determines the matrix dimensions dynamically.
template <NumericOrComplex T, StorageOrder S>
template <bool B>
COOmap<T, S>::COOmap(std::bool_constant<B>,
                     SizetPairContainer auto const& indexes,
                     NumericContainer auto const& values) {
#ifdef DEBUG
    size_t index_size =
        static_cast<size_t>(std::distance(indexes.begin(), indexes.end()));
    size_t values_size =
        static_cast<size_t>(std::distance(values.begin(), values.end()));
    assert(index_size == values_size &&
           "Error in COOmap constructor: sizes don't match.\n");
#endif

    size_t max_index_r = 0, max_index_c = 0;
    matrixptr = std::make_unique<valuesmap>();

    auto it2 = values.begin();
    for (auto it1 = indexes.begin(); it1 != indexes.end(); ++it1) {
#ifdef DEBUG
        if constexpr (B) {
            assert((*it1).first < this->rows && (*it1).second < this->columns &&
                   "Error in COOmap constructor: indexes out of bounds (too "
                   "big).\n");
        }

        assert(
            (*it1).first >= 0 && (*it1).second >= 0 &&
            "Error in COOmap constructor: indexes out of bounds (negative).\n");

        assert((*matrixptr).find(*it1) == (*matrixptr).end() &&
               "Error in COOmap constructor: redefinition of the same element "
               "(equal indexes).\n");
#endif

        if constexpr (!B) {
            max_index_r = std::max(max_index_r, (*it1).first);
            max_index_c = std::max(max_index_c, (*it1).second);
        }

        (*matrixptr)[*it1] = *it2;
        it2++;
    }

    if constexpr (!B) this->resize(max_index_r + 1, max_index_c + 1);
}

/// @brief Constructs a COOmap matrix from a map of index-value pairs.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam B Boolean constant to indicate wether the matrix' size was given
/// as input.
/// @param m The map of index-value pairs.
/// @details This constructor initializes the internal map with the provided
/// index-value pairs. It validates the indices or it dynamically determines the
/// matrix dimensions if they weren't given.
template <NumericOrComplex T, StorageOrder S>
template <bool B>
COOmap<T, S>::COOmap(std::bool_constant<B>,
                     std::map<std::pair<size_t, size_t>, T> const& m) {
    matrixptr = std::make_unique<valuesmap>(Comparator<S>{});

    size_t max_index_r = 0, max_index_c = 0;
    for (auto const& [key, val] : m) {
#ifdef DEBUG
        if constexpr (B) {
            assert(key.first < this->rows && key.second < this->columns &&
                   "Error in COOmap constructor: indexes out of bounds (too "
                   "big).\n");
        }

        assert(
            key.first >= 0 && key.second >= 0 &&
            "Error in COOmap constructor: indexes out of bounds (negative).\n");
#endif

        (*matrixptr)[key] = val;

        if constexpr (!B) {
            max_index_r = std::max(max_index_r, key.first);
            max_index_c = std::max(max_index_c, key.second);
        }
    }

    if constexpr (!B) this->resize(max_index_r + 1, max_index_c + 1);
}

/// @brief Constructs a COOmap matrix by reading data from a file.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param file_name The name of the file to read the matrix data from.
/// @details This constructor reads a `Matrix Market` file, parses its contents,
/// and initializes the matrix. It validates the file format and ensures that
/// the indices are within bounds.
template <NumericOrComplex T, StorageOrder S>
COOmap<T, S>::COOmap(std::string& file_name) {
    std::ifstream file(file_name);
#ifdef DEBUG
    assert(file && "Cannot open file.");
#endif

    std::string line;
    std::getline(file, line);
#ifdef DEBUG
    assert(line.substr(0, 14) == "%%MatrixMarket" &&
           "Invalid Matrix Market file.");
#endif

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] != '%') break;
    }

    std::istringstream dims(line);
    int nnz;
    dims >> this->rows >> this->columns >> nnz;

    valuesmap m;
    size_t row, col;
    T value;
    for (int i = 0; i < nnz && std::getline(file, line); ++i) {
        std::istringstream dims(line);
        dims >> row >> col >> value;
        m[{row - 1, col - 1}] = value;
    }

    this->matrixptr = std::make_unique<valuesmap>(std::move(m));
}

/// @brief Finds the value at the specified position (read-write).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return A reference to the value at the specified position.
/// @details If the position does not exist, this function inserts a new entry
/// into the map and returns a reference to it.
template <NumericOrComplex T, StorageOrder S>
T& COOmap<T, S>::find_dynamic(size_t i, size_t j) {
    return (*matrixptr)[{i, j}];
}

/// @brief Finds the value at the specified position (read-only).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return The value at the specified position.
/// @details This function checks if the specified position exists in the map.
/// If it does, the value is returned. Otherwise, it returns 0.
template <NumericOrComplex T, StorageOrder S>
T COOmap<T, S>::find_dynamic_const(size_t i, size_t j) const {
    if ((*matrixptr).find({i, j}) != (*matrixptr).end())
        return (*matrixptr)[{i, j}];
    else
        return 0;
}

/// @brief In the process of compressing the matrix sends the next tuple.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param first_it Flags the first iteration.
/// @param tupleptr Pointer to a tuple representing a matrix element.
/// @details This function iterates through the dynamic data structure and
/// sends the tuples to a receiver that uses them to build a compressed
/// representation.
template <NumericOrComplex T, StorageOrder S>
void COOmap<T, S>::compress_from_dynamic(
    bool first_it, std::unique_ptr<ivtuple> const& tupleptr) {
    auto static matrixit = (*matrixptr).begin();

    if (first_it) {
        matrixit = (*matrixptr).begin();
    }

    (*tupleptr) = std::make_tuple((*matrixit).first.first,
                                  (*matrixit).first.second, (*matrixit).second);

    ++matrixit;
}

/// @brief Takes a triplet and inserts it into the data structure.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param first_it Flags the first iteration.
/// @param tupleptr Pointer to a tuple representing a matrix element.
/// @details This function takes a tuple and inserts its contents into the
/// dynamic data structure.
template <NumericOrComplex T, StorageOrder S>
void COOmap<T, S>::uncompress_from_triplets(
    bool first_it, std::unique_ptr<ivtuple> const& tupleptr) {
    (*matrixptr)[std::make_pair(std::get<0>(*tupleptr),
                                std::get<1>(*tupleptr))] =
        std::get<2>(*tupleptr);
}

/// @brief Gets the number of non-zero elements.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @return The number of non-zero elements.
/// @details This function returns the size of the map, which represents the
/// number of non-zero elements in the matrix.
template <NumericOrComplex T, StorageOrder S>
size_t COOmap<T, S>::get_num_elements_dynamic() const {
    return (*matrixptr).size();
}

/// @brief Initializes the dynamic storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function allocates memory for the map.
template <NumericOrComplex T, StorageOrder S>
void COOmap<T, S>::initialize_dynamic() {
    matrixptr = std::make_unique<valuesmap>();
}

/// @brief Releases the dynamic storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function deallocates the memory used by the map.
template <NumericOrComplex T, StorageOrder S>
void COOmap<T, S>::release_dynamic() {
    matrixptr.reset();
}

/// @brief Computes the norm of the matrix.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam N The type of norm to compute (Infinity, One, or Frobenius).
/// @return The computed norm value.
/// @details This function calculates the specified norm by iterating through
/// the map and applying the appropriate formula.
template <NumericOrComplex T, StorageOrder S>
template <NormType N>
double COOmap<T, S>::norm_dynamic() const {
    if constexpr (N == Infinity) {
        std::vector<double> partial_res(this->rows, 0);

        for (auto& el : *matrixptr) {
            partial_res[el.first.first] += std::abs(el.second);
        }
        return *std::max_element(partial_res.begin(), partial_res.end());
    }

    else if constexpr (N == One) {
        std::vector<double> partial_res(this->columns, 0);

        for (auto& el : *matrixptr) {
            partial_res[el.first.second] += std::abs(el.second);
        }
        return *std::max_element(partial_res.begin(), partial_res.end());
    }
    else {
        double sum = 0.0;

        for (auto& el : *matrixptr) {
            sum += std::abs(el.second) * std::abs(el.second);
        }
        return std::sqrt(sum);
    }
}

/// @brief Removes the element at the specified position.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return True if the element was removed, false otherwise.
/// @details This function removes the specified position from the map. If the
/// position does not exist, it returns false.
template <NumericOrComplex T, StorageOrder S>
bool COOmap<T, S>::remove_dynamic(size_t i, size_t j) {
    return (*matrixptr).erase({i, j});
}

/// @brief Prints the matrix in dynamic format to the standard output.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function iterates through the matrix and prints its elements
/// in a human-readable format. It handles both row-major and column-major
/// storage.
template <NumericOrComplex T, StorageOrder S>
void COOmap<T, S>::print_dynamic() const {
    auto it = matrixptr->begin();
    size_t dist = matrixptr->size();
    size_t count = 0;

    if constexpr (S == rowMajor) {
        for (size_t i = 0; i < this->rows; ++i) {
            for (size_t j = 0; j < this->columns; ++j) {
                if ((*it).first.first == i && (*it).first.second == j) {
                    std::cout << (*it).second << " ";
                    if (count < dist) ++it;
                    count++;
                }
                else {
                    std::cout << "0 ";
                }
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout
            << "Printing the transpose matrix (since it is stored column-wise)."
            << std::endl;
        for (size_t i = 0; i < this->columns; ++i) {
            for (size_t j = 0; j < this->rows; ++j) {
                if ((*it).first.first == j && (*it).first.second == i) {
                    std::cout << (*it).second << " ";
                    if (count < dist) ++it;
                    count++;
                }
                else {
                    std::cout << "0 ";
                }
            }
            std::cout << std::endl;
        }
    }
}

/// @brief Performs matrix-vector product.
/// @param m An object of type COOmap representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
std::vector<T> by_vector_dynamic(COOmap<T, S> const& m,
                                 std::vector<T> const& v) {
    std::vector<T> result(m.rows);
    auto mapit = m.matrixptr->begin();

    while (mapit != m.matrixptr->end()) {
        result[mapit->first.first] += mapit->second * v[mapit->first.second];
        ++mapit;
    }

    return std::move(result);
}

}  // namespace algebra

#endif
