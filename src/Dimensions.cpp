#include "Dimensions.hpp"

namespace algebra {

/// @brief Constructs a Dimensions object with the specified number of rows and
/// columns.
/// @param r The number of rows.
/// @param c The number of columns.
Dimensions::Dimensions(size_t r, size_t c) : rows{r}, columns{c} {};

/// @brief Gets the number of rows in the matrix.
/// @return The number of rows.
size_t Dimensions::get_rows() const { return rows; };

/// @brief Gets the number of columns in the matrix.
/// @return The number of columns.
size_t Dimensions::get_columns() const { return columns; };

/// @brief Resizes the matrix to the specified number of rows and columns.
/// @param r The new number of rows.
/// @param c The new number of columns.
void Dimensions::resize(size_t const& r, size_t const& c) {
    rows = r;
    columns = c;
}

}  // namespace algebra
