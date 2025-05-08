#ifndef YALEIMP_HPP
#define YALEIMP_HPP

#include <cassert>
#include <iostream>

#include "Comparators.hpp"
#include "YALE.hpp"

using namespace comparators;
namespace algebra {

/// @brief Constructs a YALE matrix from compressed format data.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam B Boolean constant to indicate wether the matrix' size was given
/// as input.
/// @param out The outer index container.
/// @param in The inner index container.
/// @param val The container of matrix values.
/// @details This constructor initializes the YALE matrix by iterating through
/// the provided containers. It validates the indices and populates the internal
/// data structures.
template <NumericOrComplex T, StorageOrder S>
template <bool B>
YALE<T, S>::YALE(std::bool_constant<B>, SizetContainer auto const& out,
                 SizetContainer auto const& in,
                 NumericContainer auto const& val) {
    size_t inner_index_size =
        static_cast<size_t>(std::distance(in.begin(), in.end()));
    size_t outer_index_size =
        static_cast<size_t>(std::distance(out.begin(), out.end()));

#ifdef DEBUG
    size_t values_size =
        static_cast<size_t>(std::distance(val.begin(), val.end()));
    assert(outer_index_size == values_size &&
           "Error in YALE constructor: sizes don't match.\n");

    if constexpr (B && S == rowMajor) {
        assert(inner_index_size == (this->rows + 1) &&
               "Error in YALE constructor: sizes don't match.\n");
    }
    else if constexpr (B && S == columnMajor) {
        assert(inner_index_size == (this->columns + 1) &&
               "Error in YALE constructor: sizes don't match.\n");
    }
#endif

    comparator = Comparator<S>{};

    innerindex_ptr = std::make_unique<std::vector<size_t>>();
    outerindex_ptr = std::make_unique<std::vector<size_t>>();
    values_ptr = std::make_unique<std::vector<T>>();

    (*innerindex_ptr).reserve(inner_index_size);
    (*outerindex_ptr).reserve(outer_index_size);
    (*values_ptr).reserve(outer_index_size);

    size_t max_size_outer = 0;

    auto val_it = val.begin();
    auto outer_it = out.begin();
    auto inner_it = in.begin();

    while (val_it != val.end()) {
#ifdef DEBUG
        if constexpr (B && S == rowMajor) {
            assert(*outer_it < this->columns &&
                   "Error in YALE constructor: outer index out of bounds (too "
                   "big).\n");
        }
        else if constexpr (B && S == columnMajor) {
            assert(*outer_it < this->rows &&
                   "Error in YALE constructor: outer index out of bounds (too "
                   "big).\n");
        }

        assert(*inner_it <= outer_index_size &&
               "Error in YALE constructor: inner index out of bounds (too "
               "big).\n");
        assert(
            *outer_it >= 0 && *inner_it >= 0 &&
            "Error in YALE constructor: indexes out of bounds (negative).\n");

        if (outer_it != out.begin()) {
            auto prev_outer = std::prev(outer_it);
            auto prev_inner = std::prev(inner_it);

            assert((*outer_it > *prev_outer ||
                    outer_it - out.begin() >= *prev_inner) &&
                   "Error in YALE constructor: redefinition of the same "
                   "element (equal or misordered indexes).\n");
        }
#endif

        if constexpr (!B) max_size_outer = std::max(max_size_outer, *outer_it);

        outerindex_ptr->push_back(*outer_it);
        values_ptr->push_back(*val_it);
        ++val_it;
        ++outer_it;

        while (inner_it != in.end() && outer_it - out.begin() >= *inner_it) {
            innerindex_ptr->push_back(*inner_it);
            ++inner_it;
        }
    }

    if constexpr (!B && S == rowMajor)
        this->resize(inner_index_size - 1, max_size_outer + 1);
    else if constexpr (!B && S == columnMajor)
        this->resize(max_size_outer + 1, inner_index_size - 1);
}

/// @brief Finds the value at the specified position (read-only).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return The value at the specified position.
/// @details This function uses binary search to locate the value in the
/// compressed storage format. If the value is not found, it returns a
/// default-constructed value.
template <NumericOrComplex T, StorageOrder S>
T YALE<T, S>::find_compressed_const(size_t i, size_t j) const {
    auto [in, out] = inner_outer(i, j);
    int cumulated = (*innerindex_ptr)[in];
    int current_line = (*innerindex_ptr)[in + 1] - cumulated;

    auto lower = std::lower_bound(
        (*outerindex_ptr).begin() + cumulated,
        (*outerindex_ptr).begin() + current_line + cumulated - 1, out);

    if (*lower == out) {
        return (*values_ptr)[lower - (*outerindex_ptr).begin()];
    }

    return T{};
}

/// @brief Finds the value at the specified position (read-write).
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return A reference to the value at the specified position.
/// @details If the value does not exist, this function inserts a new entry
/// into the compressed storage format and returns a reference to it.
template <NumericOrComplex T, StorageOrder S>
T& YALE<T, S>::find_compressed(size_t i, size_t j) {
    auto [in, out] = inner_outer(i, j);
    size_t cumulated = (*innerindex_ptr)[in];
    size_t current_line = (*innerindex_ptr)[in + 1] - cumulated;

    auto lower = std::lower_bound(
        (*outerindex_ptr).begin() + cumulated,
        (*outerindex_ptr).begin() + current_line + cumulated - 1, out);

    if (*lower == out) {
        return (*values_ptr)[lower - (*outerindex_ptr).begin()];
    }
    else {
        typename std::vector<T>::iterator ref;
        auto diff = lower - (*outerindex_ptr).begin();

        if (*lower > out) {
            ref = (*values_ptr).insert((*values_ptr).begin() + diff, 0);
            (*outerindex_ptr).insert((*outerindex_ptr).begin() + diff, out);
        }
        else {
            ref = (*values_ptr).insert((*values_ptr).begin() + diff + 1, 0);
            (*outerindex_ptr).insert((*outerindex_ptr).begin() + diff + 1, out);
        }

        for (auto i = 1 + in; i < rows + 1; ++i) {
            (*innerindex_ptr)[i]++;
        }

        return *ref;
    }
}

/// @brief In the process of compressing the matrix sends the next tuple.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param first_it Flags the first iteration.
/// @param tupleptr Pointer to a tuple representing a matrix element.
/// @details This function removes the specified element from the compressed
/// storage format and adjusts the indices accordingly.
template <NumericOrComplex T, StorageOrder S>
void YALE<T, S>::uncompress_from_compressed(
    bool first_it, std::unique_ptr<ivtuple> const& tupleptr) {
    auto static innerit = (*innerindex_ptr).begin();
    auto static outerit = (*outerindex_ptr).begin();
    auto static valuesit = (*values_ptr).begin();

    if (first_it) {
        innerit = (*innerindex_ptr).begin();
        outerit = (*outerindex_ptr).begin();
        valuesit = (*values_ptr).begin();
    }

    auto next = innerit;
    ++next;

    if constexpr (S == rowMajor) {
        (*tupleptr) = std::make_tuple(innerit - (*innerindex_ptr).begin(),
                                      *outerit, *valuesit);
    }
    else {
        (*tupleptr) = std::make_tuple(
            *outerit, innerit - (*innerindex_ptr).begin(), *valuesit);
    }

    ++outerit;
    ++valuesit;

    if (valuesit - (*values_ptr).begin() >= *next) {
        ++innerit;
    }
}

/// @brief Takes a triplet and inserts it into the data structure.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param last_it Flags the last iteration.
/// @param tupleptr Pointer to a tuple representing a matrix element.
/// @details This function processes the triplets and updates the compressed
/// storage format. If it's the last iteration, it finalizes the inner index
/// array.
template <NumericOrComplex T, StorageOrder S>
void YALE<T, S>::compress_from_triplets(
    bool last_it, std::unique_ptr<ivtuple> const& tupleptr) {
    if constexpr (S == rowMajor) {
        (*innerindex_ptr)[std::get<0>(*tupleptr) + 1]++;
        (*outerindex_ptr).push_back(std::get<1>(*tupleptr));
        (*values_ptr).push_back(std::get<2>(*tupleptr));
    }
    else {
        (*innerindex_ptr)[std::get<1>(*tupleptr) + 1]++;
        (*outerindex_ptr).push_back(std::get<0>(*tupleptr));
        (*values_ptr).push_back(std::get<2>(*tupleptr));
    }

    if (last_it) {
        auto prev = (*innerindex_ptr).begin();
        auto next = prev;
        ++next;
        for (; next != (*innerindex_ptr).end(); ++next) {
            *next += *prev;
            ++prev;
        }
    }
}

/// @brief Gets the number of non-zero elements.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @return The number of non-zero elements.
/// @details This function returns the size of the values array, which
/// represents the number of non-zero elements in the matrix.
template <NumericOrComplex T, StorageOrder S>
size_t YALE<T, S>::get_num_elements_compressed() const {
    return (*values_ptr).size();
};

/// @brief Initializes the compressed storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param insize The size of the inner index container.
/// @param outsize The size of the outer index container.
/// @param valsize The size of the values container.
/// @details This function allocates memory for the inner index, outer index,
/// and values containers and resizes them to the specified sizes.
template <NumericOrComplex T, StorageOrder S>
void YALE<T, S>::initialize_compressed(size_t insize, size_t outsize,
                                       size_t valsize) {
    innerindex_ptr = std::make_unique<indexvec>();
    outerindex_ptr = std::make_unique<indexvec>();
    values_ptr = std::make_unique<valuesvec>();

    (*innerindex_ptr).resize(insize);
    (*outerindex_ptr).reserve(outsize);
    (*values_ptr).reserve(valsize);
}

/// @brief Releases the compressed storage format.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function deallocates the memory used by the inner index,
/// outer index, and values containers.
template <NumericOrComplex T, StorageOrder S>
void YALE<T, S>::release_compressed() {
    innerindex_ptr.reset();
    outerindex_ptr.reset();
    values_ptr.reset();
}

/// @brief Computes the norm of the matrix.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @tparam N The type of norm to compute (Infinity, One, or Frobenius).
/// @return The computed norm value.
/// @details This function calculates the specified norm by iterating through
/// the compressed storage format and applying the appropriate formula.
template <NumericOrComplex T, StorageOrder S>
template <NormType N>
double YALE<T, S>::norm_compressed() const {
    if constexpr (N == One && S == rowMajor) {
        std::vector<double> par(this->columns, 0);

        for (size_t i = 0; i < (*values_ptr).size(); ++i) {
            par[(*outerindex_ptr)[i]] += std::abs((*values_ptr)[i]);
        }
        return *std::max_element(par.begin(), par.end());
    }

    else if constexpr (N == Infinity && S == columnMajor) {
        std::vector<double> par(this->rows, 0);

        for (size_t i = 0; i < (*values_ptr).size(); ++i) {
            par[(*outerindex_ptr)[i]] += std::abs((*values_ptr)[i]);
        }
        return *std::max_element(par.begin(), par.end());
    }

    else if constexpr ((N == Infinity && S == rowMajor) ||
                       (N == One && S == columnMajor)) {
        double res = 0.0;
        double sum = 0.0;
        size_t cumulated;
        auto it = (*values_ptr).begin();

        for (size_t k = 1; k < (*innerindex_ptr).size(); ++k) {
            cumulated = (*innerindex_ptr)[k];
            for (auto it = (*values_ptr).begin() + (*innerindex_ptr)[k - 1];
                 it != (*values_ptr).begin() + cumulated; ++it) {
                sum += std::abs(*it);
            }

            res = std::max(res, sum);
            sum = 0.0;
        }
        return res;
    }

    else {
        double sum = 0.0;
        for (auto const& el : *values_ptr) {
            sum += std::abs(el) * std::abs(el);
        }
        return std::sqrt(sum);
    }
}

/// @brief Computes the inner and outer indices for a given position.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return A pair representing the inner and outer indexes.
/// @details This function determines the inner and outer indexes based on
/// the storage order of the matrix.
template <NumericOrComplex T, StorageOrder S>
std::pair<size_t, size_t> YALE<T, S>::inner_outer(size_t i, size_t j) const {
    if constexpr (S == rowMajor) {
        return {i, j};
    }
    else {
        return {j, i};
    }
}

/// @brief Removes the element at the specified position.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @param i The row index.
/// @param j The column index.
/// @return True if the element was removed, false otherwise.
/// @details This function removes the specified element from the compressed
/// storage format and adjusts the indexes accordingly.
template <NumericOrComplex T, StorageOrder S>
bool YALE<T, S>::remove_compressed(size_t i, size_t j) {
    auto [inner_index_to_find, outer_index_to_find] = inner_outer(i, j);

    auto inner_indices = *(this->innerindex_ptr);
    auto outer_indices = *(this->outerindex_ptr);
    auto values = *(this->values_ptr);

    size_t row_start, row_end;

    if constexpr (S == rowMajor) {
        row_start = inner_indices[i];
        row_end = inner_indices[i + 1];
    }
    else {
        row_start = inner_indices[j];
        row_end = inner_indices[j + 1];
    }

    for (size_t k = row_start; k < row_end; ++k) {
        size_t current_outer_index = outer_indices[k];
        if (current_outer_index == outer_index_to_find) {
            values_ptr->erase(values_ptr->begin() + k);
            outerindex_ptr->erase(outerindex_ptr->begin() + k);

            for (size_t l = (S == rowMajor ? i + 1 : j + 1);
                 l < innerindex_ptr->size(); ++l) {
                (*innerindex_ptr)[l]--;
            }
            return true;
        }
    }

    return false;
}

/// @brief Prints the matrix in compressed format to the standard output.
/// @tparam T The type of the matrix elements.
/// @tparam S The storage order (row-major or column-major).
/// @details This function iterates through the compressed storage format and
/// prints values, outer indices, and inner indices.
template <NumericOrComplex T, StorageOrder S>
void YALE<T, S>::print_compressed() const {
    std::cout << "Values: ";
    for (const auto& el : *values_ptr) {
        std::cout << el << " ";
    }
    std::cout << std::endl;
    std::cout << "Outer indexes: ";
    for (auto const& el : *outerindex_ptr) {
        std::cout << el << " ";
    }
    std::cout << std::endl;

    std::cout << "Inner indexes: ";
    for (auto const& el : *innerindex_ptr) {
        std::cout << el << " ";
    }
    std::cout << std::endl;
}

/// @brief Performs matrix-vector product.
/// @param m An object of type YALE representing the matrix, i.e. the lhs.
/// @param v The vector, i.e. the rhs.
/// @tparam T The type of the matrix elements (numeric or complex).
/// @tparam S The storage order (row-major or column-major).
template <NumericOrComplex T, StorageOrder S>
std::vector<T> by_vector_compressed(YALE<T, S> const& m,
                                    std::vector<T> const& v) {
    std::vector<T> result(m.rows, T{});

    const auto inner_indices = *(m.innerindex_ptr);
    const auto outer_indices = *(m.outerindex_ptr);
    const auto values = *(m.values_ptr);
    const size_t num_rows = m.rows;

    for (size_t i = 0; i < num_rows; ++i) {
        size_t start = inner_indices[i];
        size_t end = inner_indices[i + 1];

        for (size_t j = start; j < end; ++j) {
            size_t outer_index = outer_indices[j];
            T value = values[j];

            if constexpr (S == rowMajor) {
                result[i] += value * v[outer_index];
            }
            else {
                result[outer_index] += value * v[i];
            }
        }
    }

    return result;
}

}  // namespace algebra

#endif
