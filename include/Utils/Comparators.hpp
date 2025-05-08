#ifndef COMPARATORS_HPP
#define COMPARATORS_HPP

#include <algorithm>
#include <numeric>
#include <vector>

namespace comparators {

/// @brief Enum representing the type of matrix norm.
enum NormType { One, Infinity, Frobenius };

/// @brief Enum representing the storage order of a matrix.
enum StorageOrder { rowMajor, columnMajor };

/// @brief Comparator for comparing matrix indices based on the storage order.
/// @tparam S The storage order (row-major or column-major).
template <StorageOrder S>
struct Comparator {
    /// @brief Compares two index pairs based on the storage order.
    /// @param a The first index pair.
    /// @param b The second index pair.
    /// @return True if a is `less than` b based on the storage order, false
    /// otherwise.
    bool operator()(auto const& a, auto const& b) const {
        if constexpr (S == rowMajor) {
            if (a.first != b.first) {
                return a.first < b.first;
            }
            else {
                return a.second < b.second;
            }
        }
        else {
            if (a.second != b.second) {
                return a.second < b.second;
            }
            else {
                return a.first < b.first;
            }
        }
    }
};

/// @brief Computes the permutation vector to sort a vector based on a
/// comparator.
/// @tparam T The type of the elements in the vector.
/// @tparam S The storage order (row-major or column-major).
/// @param vec The vector to be sorted.
/// @param comparator The comparator to use for sorting.
/// @return A vector of indices representing the permutation.
template <class T, StorageOrder S>
std::vector<std::size_t> sort_permutation(std::vector<T> const& vec,
                                          Comparator<S> const& comparator) {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);

    std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) {
        return comparator(vec[i], vec[j]);
    });
    return p;
}

/// @brief Applies a permutation to a vector.
/// @tparam T The type of the elements in the vector.
/// @param vec The vector to permute.
/// @param p The permutation vector.
template <typename T>
void apply_permutation(std::vector<T>& vec, std::vector<std::size_t> const& p) {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (done[i]) {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j) {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

}  // namespace comparators
#endif
