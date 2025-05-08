#ifndef DIMENSIONS_HPP
#define DIMENSIONS_HPP

#include <cstddef>

namespace algebra {
/// @brief Represents the dimensions of a matrix.
class Dimensions {
   public:
    /// @brief Gets the number of rows in the matrix.
    /// @return The number of rows.
    size_t get_rows() const;

    /// @brief Gets the number of columns in the matrix.
    /// @return The number of columns.
    size_t get_columns() const;

    /// @brief Resizes the matrix to the specified number of rows and columns.
    /// @param r The new number of rows.
    /// @param c The new number of columns.
    void resize(size_t const& r, size_t const& c);

   protected:
    /// @brief Default constructor for the Dimensions class.
    Dimensions() = default;

    /// @brief Constructs a Dimensions object with the specified rows and
    /// columns.
    /// @param r The number of rows.
    /// @param c The number of columns.
    Dimensions(size_t r, size_t c);

    size_t rows;     ///< The number of rows in the matrix.
    size_t columns;  ///< The number of columns in the matrix.
};

}  // namespace algebra
#endif
