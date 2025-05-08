#ifndef TESTS_HPP
#define TESTS_HPP

#include <iostream>

#include "Concepts.hpp"

void run_tests();
void test_concepts();
void test_constructors();
void test_norms();
void test_compress_uncompress();
void test_remove();
void test_matrixmarket_constructors();
void test_matrixvector();
void test_complex();
void test_dotproduct_timing();

template <typename T>
void testSizetPair() {
    if constexpr (algebra::SizetPair<T>) {
        std::cout << __PRETTY_FUNCTION__
                  << " satisfies the \"SizetPair\" concept." << std::endl;
    }
    else {
        std::cout << __PRETTY_FUNCTION__
                  << " does NOT satisfy the \"SizetPair\" concept."
                  << std::endl;
    }
}

template <typename T>
void testNumericOrComplex() {
    if constexpr (algebra::NumericOrComplex<T>) {
        std::cout << __PRETTY_FUNCTION__
                  << " satisfies the \"NumericOrComplex\" concept."
                  << std::endl;
    }
    else {
        std::cout << __PRETTY_FUNCTION__
                  << " does NOT satisfy the \"NumericOrComplex\" concept."
                  << std::endl;
    }
}

template <typename T>
void testSizetPairContainer() {
    if constexpr (algebra::SizetPairContainer<T>) {
        std::cout << __PRETTY_FUNCTION__
                  << " satisfies the \"SizetPairContainer\" concept."
                  << std::endl;
    }
    else {
        std::cout << __PRETTY_FUNCTION__
                  << " does NOT satisfy the \"SizetPairContainer\" "
                     "concept."
                  << std::endl;
    }
}

template <typename T>
void testNumericContainer() {
    if constexpr (algebra::NumericContainer<T>) {
        std::cout << __PRETTY_FUNCTION__
                  << " satisfies the \"NumericContainer\" concept."
                  << std::endl;
    }
    else {
        std::cout << __PRETTY_FUNCTION__
                  << " does NOT satisfy the \"NumericContainer\" "
                     "concept."
                  << std::endl;
    }
}

#endif
