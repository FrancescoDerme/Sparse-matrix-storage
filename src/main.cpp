#include <chrono>
#include <random>

#include "COOImpl.hpp"
#include "COOmap.hpp"
#include "MatrixImpl.hpp"
#include "YALEImpl.hpp"

#ifdef TEST
#include "tests.hpp"
#endif

using namespace algebra;

std::vector<double> fillRandomVector(std::vector<double>& vec, double min,
                                     double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(min, max);

    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] = dist(gen);
    }

    return vec;
}

int main() {
#ifdef TEST
    run_tests();
#else
    std::cout << "TESTING EFFICIENCY OF THE MATRIX-VECTOR PRODUCT" << std::endl;

    std::string s{"matrix.mtx"};
    Matrix<double, YALE, COO, columnMajor> m(s);
    std::vector<double> vec(m.get_columns());
    fillRandomVector(vec, -100.0, 100.0);

    std::cout << std::setprecision(20);

    std::chrono::duration<double, std::milli> duration1;
    for (size_t i = 0; i < 1000; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> res1 = m * vec;
        auto end = std::chrono::high_resolution_clock::now();
        duration1 += end - start;
    }

    std::cout
        << "Matrix-vector product (1000 iteartions) in uncompressed state "
           "took:\t"
        << duration1.count() << " ms\n";

    m.compress();

    std::chrono::duration<double, std::milli> duration2;
    for (size_t i = 0; i < 1000; ++i) {
        auto start2 = std::chrono::high_resolution_clock::now();
        std::vector<double> res2 = m * vec;
        auto end2 = std::chrono::high_resolution_clock::now();
        duration2 += end2 - start2;
    }

    std::cout << "Matrix-vector product (1000 iterations) in compressed state "
                 "took:  \t"
              << duration2.count() << " ms\n";

    std::cout << std::setprecision(6);
    std::cout << std::endl;
#endif
}
