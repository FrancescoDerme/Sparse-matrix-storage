#include "tests.hpp"

#include <COOImpl.hpp>
#include <COOmapImpl.hpp>
#include <MatrixImpl.hpp>
#include <YALEImpl.hpp>
#include <chrono>
#include <complex>
#include <forward_list>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "Matrix.hpp"

void run_tests() {
    test_concepts();
    test_constructors();
    test_norms();
    test_compress_uncompress();
    test_remove();
    test_matrixmarket_constructors();
    test_matrixvector();
    test_complex();
    test_dotproduct_timing();
}

void test_concepts() {
    std::cout << "TESTING CONCEPTS" << std::endl;
    using indexpair = std::pair<size_t, size_t>;

    testSizetPair<indexpair>();
    testSizetPair<std::pair<const size_t, size_t>>();
    testSizetPair<std::pair<const size_t, const size_t>>();
    testSizetPair<std::pair<int, int>>();
    std::cout << std::endl;

    testNumericOrComplex<int>();
    testNumericOrComplex<std::complex<double>>();
    testNumericOrComplex<std::string>();
    std::cout << std::endl;

    testSizetPairContainer<std::vector<indexpair>>();
    testSizetPairContainer<std::forward_list<indexpair>>();
    testSizetPairContainer<std::set<indexpair>>();
    testSizetPairContainer<std::map<indexpair, double>>();
    std::cout << std::endl;

    testNumericContainer<std::vector<double>>();
    std::cout << std::endl;
}

void test_constructors() {
    std::cout << "TESTING CONSTRUCTORS" << std::endl;
    std::cout << "Phase 1: testing how constructors behave if given dimensions "
                 "or not."
              << std::endl;
    using namespace algebra;
    std::vector<std::pair<size_t, size_t>> ind{{0, 0}, {12, 16}};
    std::vector<double> val{1.2, -3.7};

    Matrix<double, YALE, COO, columnMajor> m(UseDynamic{}, ind, val);
    std::cout << "Expected rows: 13,\trows: " << m.get_rows() << std::endl;
    std::cout << "Expected columns: 17,\tcolumns: " << m.get_columns()
              << std::endl;

    Matrix<double, YALE, COO, columnMajor> m1(UseDynamic{}, 20, 20, ind, val);
    std::cout << "Expected rows: 20,\trows: " << m1.get_rows() << std::endl;
    std::cout << "Expected columns: 20,\tcolumns: " << m1.get_columns()
              << std::endl;

    Matrix<double, YALE, COOmap, columnMajor> m2(UseDynamic{}, ind, val);
    std::cout << "Expected rows: 13,\trows: " << m2.get_rows() << std::endl;
    std::cout << "Expected columns: 17,\tcolumns: " << m2.get_columns()
              << std::endl;

    Matrix<double, YALE, COOmap, columnMajor> m3(UseDynamic{}, 20, 20, ind,
                                                 val);
    std::cout << "Expected rows: 20,\trows: " << m3.get_rows() << std::endl;
    std::cout << "Expected columns: 20,\tcolumns: " << m3.get_columns()
              << std::endl;

    std::map<std::pair<size_t, size_t>, int> val1;
    val1[{0, 0}] = 10;
    val1[{13, 14}] = -6;

    Matrix<int, YALE, COO, columnMajor> m4(UseDynamic{}, val1);
    std::cout << "Expected rows: 14,\trows: " << m4.get_rows() << std::endl;
    std::cout << "Expected columns: 15,\tcolumns: " << m4.get_columns()
              << std::endl;

    Matrix<int, YALE, COO, columnMajor> m5(UseDynamic{}, 194, 2077, val1);
    std::cout << "Expected rows: 194,\trows: " << m5.get_rows() << std::endl;
    std::cout << "Expected columns: 2077,\tcolumns: " << m5.get_columns()
              << std::endl;

    Matrix<int, YALE, COOmap, columnMajor> m6(UseDynamic{}, val1);
    std::cout << "Expected rows: 14,\trows: " << m6.get_rows() << std::endl;
    std::cout << "Expected columns: 15,\tcolumns: " << m6.get_columns()
              << std::endl;

    Matrix<int, YALE, COOmap, columnMajor> m7(UseDynamic{}, 15, 15, val1);
    std::cout << "Expected rows: 15,\trows: " << m7.get_rows() << std::endl;
    std::cout << "Expected columns: 15,\tcolumns: " << m7.get_columns()
              << std::endl;

    std::vector<size_t> y1{0, 3, 0, 0, 2, 2, 3};
    std::vector<size_t> y2{0, 2, 3, 5, 7};
    std::vector<int> y3{10, 9, 7, 2, 2, 6, 12};

    Matrix<int, YALE, COO, columnMajor> m8(UseCompressed{}, y1, y2, y3);
    std::cout << "Expected rows: 4,\trows: " << m8.get_rows() << std::endl;
    std::cout << "Expected columns: 4,\tcolumns: " << m8.get_columns()
              << std::endl;

    Matrix<int, YALE, COO, columnMajor> m9(UseCompressed{}, 4, 4, y1, y2, y3);
    std::cout << "Expected rows: 4,\trows: " << m9.get_rows() << std::endl;
    std::cout << "Expected columns: 4,\tcolumns: " << m9.get_columns()
              << std::endl;
    std::cout << std::endl;

    std::cout << "Phase 2: testing edge cases." << std::endl;

    std::vector<size_t> y4{0, 1};
    std::vector<size_t> y5{0, 1, 1, 2};
    std::vector<int> y6{1, 1};

    Matrix<int, YALE, COO, rowMajor> m10(UseCompressed{}, y4, y5, y6);

    std::cout << "Expected values: 1, 1" << std::endl;
    std::cout << "Expected outer indexes: 0, 1" << std::endl;
    std::cout << "Expected inner indexes: 0, 1, 1, 2" << std::endl;
    m10.print();

    std::vector<size_t> y7{0, 1};
    std::vector<size_t> y8{0, 1, 1, 2};
    std::vector<int> y9{1, 1};

    Matrix<int, YALE, COO, columnMajor> m11(UseCompressed{}, y7, y8, y9);

    std::cout << "Expected values: 1, 1" << std::endl;
    std::cout << "Expected outer indexes: 0, 1" << std::endl;
    std::cout << "Expected inner indexes: 0, 1, 1, 2" << std::endl;
    m11.print();

    std::cout << std::endl;
}

void test_norms() {
    std::cout << "TESTING NORMS" << std::endl;
    using namespace algebra;
    double expected = 10.0;
    std::vector<std::pair<size_t, size_t>> ind{
        {0, 0}, {2, 2}, {0, 3}, {1, 1}, {1, 3}};
    std::vector<double> val{8, 4, -2, -3, 4};
    Matrix<double, YALE, COO, rowMajor> m(UseDynamic{}, ind, val);

    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m.norm<Infinity>() << std::endl;

    std::vector<size_t> c1{0, 3, 0, 0, 2, 2, 3};
    std::vector<size_t> r1{0, 2, 3, 5, 7};
    std::vector<int> v1{10, 9, 7, 2, 2, 6, 12};

    expected = 21.0;
    Matrix<int, YALE, COO, columnMajor> m1(UseCompressed{}, c1, r1, v1);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m1.norm<Infinity>() << std::endl;

    expected = 19.0;
    Matrix<int, YALE, COO, columnMajor> m3(UseCompressed{}, c1, r1, v1);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m3.norm<One>() << std::endl;

    expected = 20.445;
    Matrix<int, YALE, COO, columnMajor> m5(UseCompressed{}, c1, r1, v1);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m5.norm<Frobenius>() << std::endl;
    std::cout << std::endl;

    std::vector<size_t> c{0, 1, 1, 3, 2, 3, 4, 5};
    std::vector<size_t> r{0, 2, 4, 7, 8};
    std::vector<int> v{1, 2, 3, 4, 5, 6, 7, 8};

    expected = 18.0;
    Matrix<int, YALE, COO, rowMajor> m6(UseCompressed{}, c, r, v);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m6.norm<Infinity>() << std::endl;

    expected = 10.0;
    Matrix<int, YALE, COO, rowMajor> m7(UseCompressed{}, c, r, v);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m7.norm<One>() << std::endl;

    expected = 14.2829;
    Matrix<int, YALE, COO, rowMajor> m8(UseCompressed{}, c, r, v);
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m8.norm<Frobenius>() << std::endl;
    std::cout << std::endl;
}

void test_compress_uncompress() {
    std::cout << "TESTING COMPRESS / UNCOMPRESS METHODS" << std::endl;
    using namespace algebra;
    std::vector<std::pair<size_t, size_t>> i{{0, 0}, {0, 1}, {1, 0}};
    std::vector<double> v{0, 1, 2};
    Matrix<double, YALE, COOmap, columnMajor> m(UseDynamic{}, i, v);

    std::cout << "Before compressing: " << std::endl;
    m.print();
    m.compress();
    std::cout << "Afer compressing and before uncompressing again: "
              << std::endl;
    m.print();
    m.uncompress();
    std::cout << "Afer uncompressing and before compressing again: "
              << std::endl;
    m.print();
    m.compress();
    std::cout << "Afer compressing and before uncompressing again: "
              << std::endl;
    m.print();
    m.uncompress();
    std::cout << "Afer uncompressing and before compressing again: "
              << std::endl;
    m.print();
    m.compress();
    std::cout << "Afer compressing and before uncompressing again: "
              << std::endl;
    m.print();
    m.uncompress();
    std::cout << "Afer uncompressing the values are: " << std::endl;
    m.print();
    std::cout << std::endl;
}

void test_remove() {
    std::cout << "TESTING REMOVE METHODS" << std::endl;
    using namespace algebra;
    std::vector<std::pair<size_t, size_t>> i{
        {0, 0}, {0, 1}, {1, 0}, {10, 21}, {10, 22}, {10, 23}, {54, 36}};
    std::vector<double> v{3, 1, 2, 173, 174, 175, 994};
    Matrix<double, YALE, COOmap, rowMajor> m(UseDynamic{}, i, v);

    std::cout << "Before removing element (0, 1):" << std::endl;
    m.print();
    std::cout << "After removing element (0, 1):" << std::endl;
    m.remove(0, 1);
    m.print();
    std::cout << "Compressing..." << std::endl;
    m.compress();
    m.remove(1, 0);
    std::cout << "After removing element (1, 0):" << std::endl;
    m.print();
    m.remove(10, 22);
    std::cout << "After removing element (10, 22):" << std::endl;
    m.print();
    std::cout << "After removing element (54, 36):" << std::endl;
    m.remove(54, 36);
    m.print();
    std::cout << std::endl;
}

void test_matrixmarket_constructors() {
    std::cout << "TESTING MATRIX-MARKET-BASED CONSTRUCTOTS" << std::endl;
    using namespace algebra;

    std::string s{"matrix2.txt"};
    Matrix<double, YALE, COO, columnMajor> m(s);
    std::cout << "Expecting the same matrix defined in the file \"matrix.txt\":"
              << std::endl;
    m.print();

    std::cout << std::endl;

    Matrix<double, YALE, COOmap, columnMajor> m1(s);
    std::cout << "Expecting the same matrix defined in the file \"matrix.txt\":"
              << std::endl;
    m1.print();

    std::cout << std::endl;
}

void test_matrixvector() {
    std::cout << "TESTING MATRIX-VECTOR PRODUCT" << std::endl;
    using namespace algebra;
    std::vector<std::pair<size_t, size_t>> i{{0, 0}, {0, 1}, {1, 0}};
    std::vector<double> v{1, 2, 3};
    Matrix<double, YALE, COO, rowMajor> m{UseDynamic{}, i, v};
    Matrix<double, YALE, COOmap, columnMajor> m2{UseDynamic{}, i, v};
    std::vector<double> v2{1, 2};

    std::cout << "Expected result: 5, 3:" << std::endl;
    std::vector<double> v3 = m * v2;
    for (auto const& el : v3) std::cout << el << " ";
    std::cout << std::endl;

    std::cout << "Expected result: 5, 3:" << std::endl;
    std::vector<double> v4 = m2 * v2;
    for (auto const& el : v4) std::cout << el << " ";
    std::cout << std::endl;

    std::cout << "Compressing..." << std::endl;
    m.compress();
    m2.compress();

    std::cout << "Expected result: 5, 3:" << std::endl;
    v3 = m * v2;
    for (auto const& el : v3) std::cout << el << " ";
    std::cout << std::endl;

    std::cout << "Expected result: 5, 3:" << std::endl;
    v4 = m2 * v2;
    for (auto const& el : v4) std::cout << el << " ";
    std::cout << std::endl;

    std::cout << std::endl;
}

void test_complex() {
    std::cout << "TESTING WITH COMPLEX NUMBERS" << std::endl;
    using namespace algebra;
    double expected;

    std::vector<std::pair<size_t, size_t>> i{{0, 0}, {0, 1}, {1, 0}};
    std::vector<std::complex<double>> v{{1, 2}, {2, 3}, {1, 1}};

    Matrix<std::complex<double>, YALE, COO, rowMajor> m{UseDynamic{}, i, v};
    expected = 5.8416;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m.norm<Infinity>() << std::endl;
    expected = 3.6503;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m.norm<One>() << std::endl;
    expected = 4.4721;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m.norm<Frobenius>() << std::endl;

    Matrix<std::complex<double>, YALE, COOmap, rowMajor> m1{UseDynamic{}, i, v};
    expected = 5.8416;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m1.norm<Infinity>() << std::endl;
    expected = 3.6503;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m1.norm<One>() << std::endl;
    expected = 4.4721;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m1.norm<Frobenius>() << std::endl;

    std::vector<size_t> c{0, 1, 0};
    std::vector<size_t> r{0, 2, 3};
    std::vector<std::complex<double>> v1{{1, 2}, {2, 3}, {1, 1}};
    Matrix<std::complex<double>, YALE, COO, rowMajor> m2{UseCompressed{}, c, r,
                                                         v1};
    expected = 5.8416;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m2.norm<Infinity>() << std::endl;
    expected = 3.6503;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m2.norm<One>() << std::endl;
    expected = 4.4721;
    std::cout << "Expected norm: " << expected
              << ",\tcomputed norm: " << m2.norm<Frobenius>() << std::endl;

    std::cout << "Expected result: (-3,13), (0,2),\tcomputed result: ";
    std::vector<std::complex<double>> v2{{1, 1}, {2, 2}};
    std::vector<std::complex<double>> v4;
    v4 = m2 * v2;
    for (auto const& el : v4) std::cout << el << " ";
    std::cout << std::endl;

    std::cout << std::endl;
}

std::vector<double> fillRandomVector1(std::vector<double>& vec, double min,
                                      double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(min, max);

    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] = dist(gen);
    }

    return vec;
}

void test_dotproduct_timing() {
    std::cout << "TESTING EFFICIENCY OF THE MATRIX-VECTOR PRODUCT" << std::endl;
    using namespace algebra;

    std::string s{"matrix.mtx"};
    Matrix<double, YALE, COO, columnMajor> m(s);
    std::vector<double> vec(m.get_columns());
    fillRandomVector1(vec, -100.0, 100.0);

    std::cout << std::setprecision(20);

    std::chrono::duration<double, std::milli> duration1;
    for (size_t i = 0; i < 1000; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> res1 = m * vec;
        auto end = std::chrono::high_resolution_clock::now();
        duration1 += end - start;
    }

    std::cout
        << "Matrix-vector product (1000 ietartions) in uncompressed state "
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
}
