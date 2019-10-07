[![Build Status](https://circleci.com/gh/tschijnmo/fbitset.svg?style=shield)](https://circleci.com/gh/tschijnmo/fbitset)

# fbitset
A highly optimized and tunable C++ bitset library

Bitsets can be quite necessary for a lot of combinatorial problems, and can
also be used for bitmap-style of indexes for database management systems.
However, the STL `bitset` requires a static size given at the compile-time.  At
the same time, the `dynamic_bitset` from the boost library always put the bits
in a separate `vector`, which could be cache-hostile and introduce significant
overhead for the external memory allocation.

fbitset aims to be a minimalistic library attempting to merge the benefits of
the STL `bitset` and boost `dynamic_bitset`, in the same way as the
short-string optimization commonly found for STL `string`s.  When the number of
bits is less than a value set at compile time, the bits will be stored right
inside the object.  External memory will be allocated only when the number of
bits cannot fit.  For instance, on common 64bit platforms, normally an STL
`vector` would take 3x64bit of space.  So `fbitset` could store the bits
inside, without any external memory allocation, when the size is less than
3x64, and could also automatically switch to external storage with dynamic
allocation when more is needed.

When the bits are stored internally, significant care has been taken to ensure
that little overhead is introduced relative to the STL `bitset`.  By using
latest C++17 compile-time features and template programming, it is ensured that
loops can be correctly unrolled for internal storage, which can be verified by
the assembly output from G++ and Clang++.

Since this library is designed for use cases with shear demand on performance,
in addition to making sure the unrolling of the loops for internal storage,
some very portable GCC compiler intrinsics are also used.  This could utilize
native instructions for some bitwise operations when they are available, like
the `popcnt`, `bsf`, and `bsr` x86-64 instructions.

In addition to being highly optimized, by modern C++ template programming, the
bit set is also very tunable.  Things like the internal integral type used to
store the bits, the number of integers allowed inside, and the container for
external storage, can all be customized with template parameters.  At the same
time, the default setting should be sufficient for all but the most special use
cases.
