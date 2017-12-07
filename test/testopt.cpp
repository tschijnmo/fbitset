/** Some small cases for optimization checking.
 *
 * This is a tiny phony library primarily to instantiate some of the templates
 * in fbitset to make sure that the compiler has properly optimized the code.
 */

#include <cstdint>
#include <vector>

#include <fbitset.hpp>

using namespace fbitset;

//
// The primary types interested in.
//

using One64_no_ext = Fbitset<1, uint64_t, No_ext>;

using Two32_no_ext = Fbitset<2, uint32_t, No_ext>;

using Two32_ext = Fbitset<2, uint32_t, std::vector<uint32_t>>;

//
// Simple union of two bitsets.
//
// Union is a binary operation with very good coverage on the core code.  When
// things are properly optimized, the resulted code can be very similar to
// hand-written specialized code most of the times.
//

// For a single 64-bit integer with externals disabled, the resulted assembly
// should be essentially the same as the code for two 64-bit integers.

One64_no_ext union_one64_no_ext(One64_no_ext& o1, One64_no_ext& o2)
{
    return o1 | o2;
}

// The reference code for the union of two 64-bit integers.

uint64_t union_one64_ref(uint64_t& o1, uint64_t& o2) { return o1 | o2; }

Two32_no_ext union_two32_no_ext(Two32_no_ext& o1, Two32_no_ext& o2)
{
    return o1 | o2;
}

Two32_ext union_two32_ext(Two32_ext& o1, Two32_ext& o2) { return o1 | o2; }
