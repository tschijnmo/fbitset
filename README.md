https://circleci.com/gh/tschijnmo/fbitset.svg?style=shield

# fbitset
A header-only C++ bitset library focusing on performance and simplicity

Bitsets can be quite necessary for a lot of combinatorial problems.  However,
the STL `bitset` requires a static size given at the compile-time.  At the same
time, the `dynamic_bitset` from the boost library always put the bits in a
separate `vector`, which could be cache-hostile and introduce significant
overhead.

fbitset aims to be a minimal library attempting to merge the benefits of the
STL `bitset` and boost `dynamic_bitset`.  When the number of bits is less than
a value set at compile time, the bits will be stored right inside the object.
External memory will be allocated only when the number of bits cannot fit.  For
cases where the bits are stored internally, very little overhead is introduced
relative to the STL `bitset`, which is ensured by comparison of the assembly
output from G++ and Clang++.

This library is designed for the ease of development for complex combinatorial
problems.  The interface is not intended to parallel `bitset` or
`dynamic_bitset` exactly.

For best performance and code clarity, latest C++17 language features are
liberally used.

