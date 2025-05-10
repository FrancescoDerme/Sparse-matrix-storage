# PACS challenge 2: Sparse Matrix storage

- [Functionality](#functionality)
- [Getting started](#getting-started)
- [Storage methods](#storage-methods)
- [Implementation](#implementation)
  - [Inheritance hierarchy](#inheritance-hierarchy)
  - [Compression and decompression via SRP](#compression-and-decompression-via-srp)
  - [Concepts](#concepts)
- [Formatting](#formatting)

## Functionality
The code implements a templated `Matrix` class capable of storing data in both compressed and uncompressed formats, it allows switching between the two and performing matrix-vector multiplication among many other functionalities. All of it comes with a **very fast implementation** as well as a clean and **modular interface** that easily allows to expand the scope of the project or to insert it in a bigger codebase without breaking a sweat.

- To build the project, simply type _make_ in the repository where you've cloned it. Running the program with _./executable_ will then perform a thousand matrix-vector multiplications between [this matrix](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html) and a randomly generated vector in both compressed and uncompressed states and print the runtime of both.
- To test a broader range of functionalities, compile with _make test_. Running the program will then perform tests on concepts, constructors, norm methods, compress/uncompress methods, remove methods, reading-from-file functionality, matrix-vector multiplications and complex-valued matrices.
- Finally, compiling with _make debug_ will enable many assertions throughout the code that, while disabled by default for efficiency concerns, make it safer to run; indeed if something is not working properly try compiling with this option to see if there's an error in the input or in the sequence of operations or if the code is actually broken.
- _make clean_ and _make doc_ options are available to do what they claim.

> ⚠️ **Warning**: This project must be compiled with the **g++** compiler. On Windows and Linux, there should be no issues since the makefile is configured accordingly. On macOS, the makefile forces the use of g++ only if it was installed via **Homebrew**. If g++ was installed differently, update line 4 of the makefile with the correct path.

## Getting started
To build an object of type `Matrix` you need to decide what kind of values to store (`std::complex` are admissible), which compressed and uncompressed storage formats to use and which storage ordering you prefer. The project comes with `COO` and `COOmap` uncompressed storage formats, `YALE` compressed storage format and `rowMajor` and `columnMajor` orderings.

Once you've elected the ingredients to brew your first matrix it's time to call a constructor and there are plenty, just remember to pass `UseCompressed{}` or `UseDynamic{}` as the first argument to select the state in which the matrix is going to be built. Then you can pass the matrix's dimensions or not, in the second case they will be automatically inferred, and finally you have to pass the actual data. To build a `COO` or `COOmap` matrix you can pass either a `std::map` or a couple of data structures which satisfy the `SizetPairContainer` and `NumericContainer` [concepts](#concepts), to build a `YALE` matrix you have to pass a triplet of data structures of which the first two satisfy the `SizetContainer` concept and the last satisfies the `NumericContainer` concept.

It's really easier than it sounds, here we provide a simple example.

```cpp
std::vector<std::pair<size_t, size_t>> ind{{0, 0}, {12, 16}};
std::vector<double> val{1.2, -3.7};
Matrix<double, YALE, COO, columnMajor> m(UseDynamic{}, ind, val);
m.compress();
m.print();
```

`COO` and `COOmap` matrices also implement a special constructor that takes in the name of a file in *Matrix-market format* to read the matrix from as the only argument.

## Storage methods
`COO` and `COOmap` are the provided uncompressed storage types, `YALE` is the provided compressed one. All of them work with both `rowMajor` and `columnMajor` orderings and new storage methods are quite easy to add if one knows what he's doing. Internally, `COO` uses a couple of `std::forward_list`s, `COOmap` a `std::map` and `YALE` uses three `std::vector`s; such choices were made in careful consideration of the tradeoffs between computational complexity, memory load and programmer time, the latter never having the upper hand.

## Implementation
The code is **doxygen-documented**, the doxygen documentation is the best place to learn more about the inner workings of the code before diving in the source files. What is useful to bring to the reader's attention from the beginning are a couple of details including how the inheritance hierarchy, compression-decompresion mechanism and concepts system work.

Plenty methods are available to interact with your matrices and everything works the same in every format, always. You can remove or add new elements, compress and uncompress your matrix, ask for it's dimensions or the number of non-zero entries and much more, again, look at the docs if you crave for more.

### Inheritance hierarchy
Compressed and dynamic (or uncompressed) matrices are template classes that inherit publicly and **virtually** from the base class `Dimensions` that stores a matrix's number of rows and columns as well as basic methods to access them. Here, the word "virtual" is key because it ensures that when the `Matrix` template class inherits publicly from a compressed and a dynamic matrix it only receives one copy of `Dimensions`, thus eluding the old [diamond problem](https://www.geeksforgeeks.org/diamond-problem-in-cpp/).

### Compression and decompression via SRP
The compression and decompression mechanism follows the **single responsability principle** according to which each entity shall be soely and completely responsible for its data. In particular, each [storage method](#storage-methods) only implements a method to take in an row-column-value triplet and insert it into its data structures and a method to send a triplet to the caller.

Let's make an example: after building a matrix with `Matrix<double, YALE, COO, columnMajor> m(UseDynamic{}, ind, val)` as in the example above, a call to `compress()` will initialize `YALE`'s datastructures, then call N times `COO`'s method `compress_from_dynamic()` and `YALE`'s `compress_from_triplets()` where N is the number of elements currently stored in the matrix. The first gets a triplet from `COO`'s datastructures and the seconds puts the retrived triplets in `YALE`'s datastructures in the right place, finally `COO`'s datastructures are freed. Similarly, `uncompress()` will make use of `YALE`'s `uncompress_from_compressed()` and of `COO`'s `uncompress_from_triplets()`.

Everything is coordinated by the `Matrix` class which is also in charge of initializing and freeing datastructures at the right time. Since the triplets are retrived and inserted in an ordely manner by using `static` iterators, these procedures are practically as efficient as a single ad-hoc method could be, with the difference that the *send triplet* and *receive triplet* methods are **completely agnostic** to which storage method they're interacting with: `COO` doesn't need to know that his triplets are being sent to `YALE` to send them, `YALE` doesn't need to know that the tripets are being sent by `COO` to receive them. If one was to impelment another storage method, the currently-present methods wouldn't break nor would the programmer need to write new methods for old storage types, the new storage method would only need to implement its own *send triplet* and *receive triplet* methods.

### Concepts
Concepts allow to restrict the templated classes and methods to types for which they make sense in a transparent way. The project defines `ConvertibleToSizeT`, `NumericOrComplex`, `SizetPair`, `ValidContainer`, `NumericContainer`, `SizetContainer` and `SizetPairContainer` concepts one on top of the other. The templates are then able to ask specifically for what they need among these options rendering the code legible, robust and efficient.

## Formatting
The code has been evenly formatted to Google's standards with the _clang_format_ formatter with the following settings.

```lua
--style={BasedOnStyle: google, IndentWidth: 4, ColumnLimit : 80, BreakBeforeBraces: Custom, BraceWrapping: {AfterFunction: false, BeforeElse: true}}"
```
