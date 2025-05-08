#ifndef COOIMPL_HPP
#define COOIMPL_HPP

#include <cassert>
#include <complex>
#include <fstream>
#include <iostream>
#include <memory>

#include "COO.hpp"
#include "Comparators.hpp"
#include "Concepts.hpp"

using namespace comparators;
namespace algebra {

/// @brief Constructs a COO matrix from index-value pairs.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam B Boolean constant to indicate wether the matrix' size was given
/// as input.
/// @param indexes The container of index pairs.
/// @param values The container of matrix values.
/// @details This constructor sorts the input index-value pairs based on the
/// storage order and ensures that the indices are unique and within bounds. It
/// uses `sort_permutation` and `apply_permutation` to reorder the data.
template <NumericOrComplex T, StorageOrder S>
template <bool B>
COO<T, S>::COO(std::bool_constant<B>, SizetPairContainer auto const& indexes,
               NumericContainer auto const& values) {
#ifdef DEBUG
    size_t index_size =
        static_cast<size_t>(std::distance(indexes.begin(), indexes.end()));
    size_t values_size =
        static_cast<size_t>(std::distance(values.begin(), values.end()));
    assert(index_size == values_size &&
           "Error in COO constructor: sizes don't match.\n");
#endif

    std::vector<T> tempvalues(values.begin(), values.end());
    std::vector<std::pair<size_t, size_t>> tempindexes(indexes.begin(),
                                                       indexes.end());

    comparator = Comparator<S>{};
    std::vector<size_t> p = sort_permutation(tempindexes, comparator);

    apply_permutation(tempindexes, p);
    apply_permutation(tempvalues, p);

    indexptr = std::make_unique<indexlist>();
    valuesptr = std::make_unique<valueslist>();

    size_t max_index_r = 0, max_index_c = 0;
    auto it1 = tempindexes.rbegin();
    auto it2 = tempvalues.rbegin();

    while (it1 != tempindexes.rend() && it2 != tempvalues.rend()) {
#ifdef DEBUG
        if constexpr (B) {
            assert(
                (*it1).first < this->rows && (*it1).second < this->columns &&
                "Error in COO constructor: indexes out of bounds (too big).\n");
        }
        assert((*it1).first >= 0 && (*it1).second >= 0 &&
               "Error in COO constructor: indexes out of bounds (negative).\n");

        if (it1 != tempindexes.rbegin()) {
            auto compare = it1;
            --compare;
            assert((*it1) != (*compare) &&
                   "Error in COO constructor: redefinition of the same element "
                   "(equal indexes).\n");
        }
#endif

        if constexpr (!B) {
            max_index_r = std::max(max_index_r, (*it1).first);
            max_index_c = std::max(max_index_c, (*it1).second);
        }

        indexptr->push_front(*it1);
        valuesptr->push_front(*it2);
        it1++;
        it2++;
    }

    if constexpr (!B) this->resize(max_index_r + 1, max_index_c + 1);
}

/// @brief Constructs a COO matrix from a map of index-value pairs.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam B Boolean constant to indicate wether the matrix' size was given
/// as input.
/// @param m The map of index-value pairs.
/// @details This constructor extracts keys and values from the input map, sorts
/// them based on the storage order, and stores them in the internal data
/// structures.
template <NumericOrComplex T, StorageOrder S>
template <bool B>
COO<T, S>::COO(std::bool_constant<B>,
               std::map<std::pair<size_t, size_t>, T> const& m) {
    std::vector<T> tempvalues;
    std::vector<std::pair<size_t, size_t>> tempindexes;

    tempvalues.reserve(m.size());
    tempindexes.reserve(m.size());

    std::transform(m.begin(), m.end(), std::back_inserter(tempvalues),
                   [](auto const& p) { return p.second; });

    std::transform(m.begin(), m.end(), std::back_inserter(tempindexes),
                   [](auto const& p) { return p.first; });

    comparator = Comparator<S>{};
    std::vector<size_t> p = sort_permutation(tempindexes, comparator);

    apply_permutation(tempindexes, p);
    apply_permutation(tempvalues, p);

    indexptr = std::make_unique<indexlist>();
    valuesptr = std::make_unique<valueslist>();

    size_t max_index_r = 0, max_index_c = 0;
    auto it1 = tempindexes.rbegin();
    auto it2 = tempvalues.rbegin();

    while (it1 != tempindexes.rend() && it2 != tempvalues.rend()) {
#ifdef DEBUG
        if constexpr (B) {
            assert(
                (*it1).first < this->rows && (*it1).second < this->columns &&
                "Error in COO constructor: indexes out of bounds (too big).\n");
        }
        assert((*it1).first >= 0 && (*it1).second >= 0 &&
               "Error in COO constructor: indexes out of bounds (negative).\n");
#endif

        if constexpr (!B) {
            max_index_r = std::max(max_index_r, (*it1).first);
            max_index_c = std::max(max_index_c, (*it1).second);
        }

        indexptr->push_front(*it1);
        valuesptr->push_front(*it2);
        it1++;
        it2++;
    }

    if constexpr (!B) this->resize(max_index_r + 1, max_index_c + 1);
}

/// @brief Constructs a COO matrix by reading data from a file.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param file_name The name of the file to read the matrix data from.
/// @details This constructor reads a `Matrix Market` file, parses its contents,
/// and initializes the matrix. It validates the file format and ensures that
/// the indices are within bounds.
template <NumericOrComplex T, StorageOrder S>
COO<T, S>::COO(std::string& file_name) {
    std::ifstream file(file_name);
#ifdef DEBUG
    assert(file && "Error in COO constructor: cannot open file.");
#endif

    std::string line;
    std::getline(file, line);

#ifdef DEBUG
    assert(line.substr(0, 14) == "%%MatrixMarket" &&
           "Error in COO constructor: invalid Matrix Market file.");
#endif

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] != '%') break;
    }

    size_t nnz;
    std::istringstream dims(line);
    dims >> this->rows >> this->columns >> nnz;

    std::vector<T> tempvalues;
    std::vector<std::pair<size_t, size_t>> tempindexes;

    tempvalues.reserve(nnz);
    tempindexes.reserve(nnz);

    size_t row, col;
    T value;

    for (size_t i = 0; i < nnz && std::getline(file, line); ++i) {
        std::istringstream dims(line);
        dims >> row >> col >> value;
        tempvalues.push_back(value);
        tempindexes.push_back({row - 1, col - 1});
    }

    comparator = Comparator<S>{};
    std::vector<size_t> p = sort_permutation(tempindexes, comparator);

    apply_permutation(tempindexes, p);
    apply_permutation(tempvalues, p);

    indexptr = std::make_unique<indexlist>();
    valuesptr = std::make_unique<valueslist>();

    auto it1 = tempindexes.rbegin();
    auto it2 = tempvalues.rbegin();

    while (it1 != tempindexes.rend() && it2 != tempvalues.rend()) {
#ifdef DEBUG
        assert((*it1).first < this->rows && (*it1).second < this->columns &&
               "Error in COO constructor: indexes out of bounds (too big).\n");
        assert((*it1).first >= 0 && (*it1).second >= 0 &&
               "Error in COO constructor: indexes out of bounds (negative).\n");

        if (it1 != tempindexes.rbegin()) {
            auto compare = it1;
            --compare;
            assert((*it1) != (*compare) &&
                   "Error in COO constructor: redefinition of the same element "
                   "(equal indexes).\n");
        }
#endif

        indexptr->push_front(*it1);
        valuesptr->push_front(*it2);
        it1++;
        it2++;
    }
}

/// @brief Finds the value at the specified position (read-only).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return The value at the specified position.
/// @details This function iterates through the index list to find the matching
/// position. If no match is found, it returns 0.
template <NumericOrComplex T, StorageOrder S>
T COO<T, S>::find_dynamic_const(size_t i, size_t j) const {
    auto it1 = indexptr->begin();
    auto it2 = valuesptr->begin();

    while (it1 != indexptr->end()) {
        if ((*it1).first == i && (*it1).second == j) {
            return *it2;
        }
        else if (!comparator(*it1, std::make_pair(i, j)))
            break;

        it1++;
        it2++;
    }

    return 0;
}

/// @brief Finds the value at the specified position (read-write).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return A reference to the value at the specified position.
/// @details If the position does not exist, this function inserts a new entry
/// into the data structure and returns a reference to it.
template <NumericOrComplex T, StorageOrder S>
T& COO<T, S>::find_dynamic(size_t i, size_t j) {
    auto it1 = indexptr->begin();
    auto next1 = indexptr->before_begin();

    auto it2 = valuesptr->begin();
    auto next2 = valuesptr->before_begin();

    while (it1 != indexptr->end()) {
        if ((*it1).first == i && (*it1).second == j) {
            return *it2;
        }
        else if (!comparator(*it1, std::make_pair(i, j))) {
            indexptr->insert_after(next1, std::make_pair(i, j));
            it2 = valuesptr->insert_after(next2, T{});
            return *it2;
        }

        it1++;
        next1++;
        it2++;
        next2++;
    }

    indexptr->insert_after(next1, std::make_pair(i, j));
    it2 = valuesptr->insert_after(next2, T{});
    return *it2;
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
void COO<T, S>::compress_from_dynamic(
    bool first_it, std::unique_ptr<ivtuple> const& tupleptr) {
    auto static valuesit = (*valuesptr).begin();
    auto static indexit = (*indexptr).begin();

    if (first_it) {
        valuesit = (*valuesptr).begin();
        indexit = (*indexptr).begin();
    }

    (*tupleptr) =
        std::make_tuple((*indexit).first, (*indexit).second, *valuesit);

    ++valuesit;
    ++indexit;
}

/// @brief Takes a triplet and inserts it into the data structure.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param first_it Flags the first iteration.
/// @param tupleptr Pointer to a tuple representing a matrix element.
/// @details This function takes a tuple and inserts its contents into the
/// dynamic data structure.
template <NumericOrComplex T, StorageOrder S>
void COO<T, S>::uncompress_from_triplets(
    bool first_it, std::unique_ptr<ivtuple> const& tupleptr) {
    auto static lastind = indexptr->before_begin();
    auto static lastval = valuesptr->before_begin();

    if (first_it) {
        lastind = indexptr->before_begin();
        lastval = valuesptr->before_begin();
    }

    lastind = indexptr->insert_after(
        lastind,
        std::make_pair(std::get<0>(*tupleptr), std::get<1>(*tupleptr)));
    lastval = valuesptr->insert_after(lastval, std::get<2>(*tupleptr));
}

/// @brief Gets the number of non-zero elements.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @return The number of non-zero elements.
/// @details This function calculates the size of the values list to determine
/// the number of non-zero elements.
template <NumericOrComplex T, StorageOrder S>
size_t COO<T, S>::get_num_elements_dynamic() const {
    return static_cast<size_t>(
        std::distance((*valuesptr).begin(), (*valuesptr).end()));
}

/// @brief Initializes the dynamic storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function allocates memory for the index and values lists.
template <NumericOrComplex T, StorageOrder S>
void COO<T, S>::initialize_dynamic() {
    indexptr = std::make_unique<indexlist>();
    valuesptr = std::make_unique<valueslist>();
}

/// @brief Releases the dynamic storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function deallocates the memory used by the index and values
/// lists.
template <NumericOrComplex T, StorageOrder S>
void COO<T, S>::release_dynamic() {
    indexptr.reset();
    valuesptr.reset();
}

/// @brief Gets the next line in the matrix.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The current row index.
/// @return A pair representing the next line's row and column indexes.
/// @details This function calculates the next line based on the storage order.
template <NumericOrComplex T, StorageOrder S>
std::pair<size_t, size_t> COO<T, S>::next_line(size_t i) const {
    if constexpr (S == rowMajor) {
        return std::make_pair(i + 1, 0);
    }
    else {
        return std::make_pair(0, i + 1);
    }
}

/// @brief Computes the norm of the matrix.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam N The type of norm to compute (Infinity, One, or Frobenius).
/// @return The computed norm value.
/// @details This function calculates the specified norm by iterating through
/// the values list and applying the appropriate formula.
template <NumericOrComplex T, StorageOrder S>
template <NormType N>
double COO<T, S>::norm_dynamic() const {
    if constexpr ((N == Infinity && S == rowMajor) ||
                  (N == One && S == columnMajor)) {
        double res = 0.0;
        double sum = 0.0;
        int i = 0;

        auto it2 = (*valuesptr).begin();
        for (auto it = (*indexptr).begin(); it != (*indexptr).end(); ++it) {
            if (comparator((*it), next_line(i))) {
                sum += std::abs(*it2);
            }
            else {
                res = std::max(sum, res);
                i++;
                sum = std::abs(*it2);
            }
            it2++;
        }
        res = std::max(sum, res);
        return res;
    }

    else if constexpr (N == One && S == rowMajor) {
        std::vector<double> partial_res(this->columns, 0);

        auto it2 = (*valuesptr).begin();
        for (auto it = (*indexptr).begin(); it != (*indexptr).end(); ++it) {
            partial_res[(*it).second] += std::abs(*it2);
            it2++;
        }
        return *std::max_element(partial_res.begin(), partial_res.end());
    }
    else if constexpr (N == Infinity && S == columnMajor) {
        std::vector<double> partial_res(this->rows, 0);

        auto it2 = (*valuesptr).begin();
        for (auto it = (*indexptr).begin(); it != (*indexptr).end(); ++it) {
            partial_res[(*it).first] += std::abs(*it2);
            it2++;
        }
        return *std::max_element(partial_res.begin(), partial_res.end());
    }
    else {
        double sum = 0.0;

        for (auto it = (*valuesptr).begin(); it != (*valuesptr).end(); ++it) {
            sum += std::abs(*it) * std::abs(*it);
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
/// @details This function searches for the specified position in the index list
/// and removes the corresponding entry from both the index and values lists.
template <NumericOrComplex T, StorageOrder S>
bool COO<T, S>::remove_dynamic(size_t i, size_t j) {
    auto prev1 = (*indexptr).before_begin();
    auto prev2 = (*valuesptr).before_begin();

    for (auto it1 = (*indexptr).begin(); it1 != (*indexptr).end(); ++it1) {
        if ((*it1).first == i && (*it1).second == j) {
            (*indexptr).erase_after(prev1);
            (*valuesptr).erase_after(prev2);
            return true;
        }
        ++prev1;
        ++prev2;
    }
    return false;
}

/// @brief Prints the matrix in dynamic format to the standard output.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function iterates through the matrix and prints its elements
/// in a human-readable format. It handles both row-major and column-major
/// storage.
template <NumericOrComplex T, StorageOrder S>
void COO<T, S>::print_dynamic() const {
    auto it = (*indexptr).begin();
    auto it2 = (*valuesptr).begin();
    int count = 1;
    auto dist = std::distance((*indexptr).begin(), (*indexptr).end());

    if constexpr (S == rowMajor) {
        for (size_t i = 0; i < (this)->rows; ++i) {
            for (size_t j = 0; j < (this)->columns; ++j) {
                if ((*it).first == i && (*it).second == j) {
                    std::cout << (*it2) << " ";
                    if (count < dist) {
                        ++it;
                        ++it2;
                        count++;
                    }
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
        for (size_t i = 0; i < (this)->columns; ++i) {
            for (size_t j = 0; j < (this)->rows; ++j) {
                if ((*it).first == j && (*it).second == i) {
                    std::cout << (*it2) << " ";
                    if (count < dist) {
                        ++it;
                        ++it2;
                        count++;
                    }
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
/// @param m An object of type COO representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
std::vector<T> by_vector_dynamic(COO<T, S> const& m, std::vector<T> const& v) {
    std::vector<T> result(m.rows);

    auto indexit = m.indexptr->begin();
    auto valuesit = m.valuesptr->begin();

    while (indexit != m.indexptr->end()) {
        result[indexit->first] += *valuesit * v[indexit->second];
        ++indexit;
        ++valuesit;
    }

    return std::move(result);
}

}  // namespace algebra

#endif
