#ifndef CONCEPTS_HPP
#define CONCEPTS_HPP

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>
#include <utility>

namespace algebra {

/// @brief Concept to check if a type is convertible to `std::size_t`.
/// @tparam T The type to check.
template <typename T>
concept ConvertibleToSizeT = std::is_convertible_v<T, std::size_t>;

/// @brief Concept to check if a type is numeric or a complex number.
/// @tparam T The type to check.
template <typename T>
concept NumericOrComplex =
    std::is_arithmetic_v<T> ||
    std::is_same_v<T, std::complex<typename T::value_type>>;

/// @brief Concept to check if a type is a pair of `size_t`.
/// @tparam T The type to check.
template <typename T>
concept SizetPair = requires {
    requires std::same_as<
        std::pair<size_t, size_t>,
        std::pair<std::decay_t<decltype(std::declval<T>().first)>,
                  std::decay_t<decltype(std::declval<T>().second)>>>;
};

/// @brief Concept to check if a type is a valid container.
/// @tparam T The type to check.
template <typename T>
concept ValidContainer = requires(T t) {
    typename T::value_type;
    typename T::iterator;
    { t.begin() } -> std::input_iterator;
    { t.end() } -> std::input_iterator;
    requires !std::is_same_v<std::remove_cvref_t<T>,
                             std::string> &&  // Exclude std::string
                 !std::is_same_v<std::remove_cvref_t<T>,
                                 std::wstring> &&  // Exclude std::wstring
                 !std::is_same_v<std::remove_cvref_t<T>,
                                 std::u16string> &&  // Exclude std::u16string
                 !std::is_same_v<std::remove_cvref_t<T>,
                                 std::u32string>;  // Exclude std::u32string
};

/// @brief Concept to check if a type is a container of numeric or complex
/// values.
/// @tparam T The type to check.
template <typename T>
concept NumericContainer =
    ValidContainer<T> && NumericOrComplex<typename T::value_type>;

/// @brief Concept to check if a type is a container of `size_t` values.
/// @tparam T The type to check.
template <typename T>
concept SizetContainer =
    ValidContainer<T> &&
    std::is_same_v<typename T::value_type,
                   std::decay_t<decltype(std::declval<size_t>())>>;

/// @brief Concept to check if a type is a container of `size_t` pairs.
/// @tparam T The type to check.
template <typename T>
concept SizetPairContainer =
    ValidContainer<T> && SizetPair<typename T::value_type>;

}  // namespace algebra
#endif
