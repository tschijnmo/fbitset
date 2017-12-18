/** Tiny library for benchmarking the loop over indices of set bits.
 *
 * In this tiny library, we have a vector of indices, which can be stored
 * either as a vector of indices (as they are here), or into a bit set.  Then
 * we can benchmark the performance of these two ways of storage.  Note that
 * here all the indices are less than 128.
 */

#ifndef FBITSET_TESTLOOP_H
#define FBITSET_TESTLOOP_H

#include <fbitset.hpp>
#include <vector>

/** The dummy function that does nothing.
 */

void sink(fbitset::Size);

/** The indices to be set.
 *
 * All of them should be less than 128.
 */

extern const std::vector<fbitset::Size> idxes;

/** The number of repeat for the benchmark.
 */

constexpr size_t N_REPEATS = 100000000;

#endif
