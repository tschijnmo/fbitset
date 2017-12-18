/** Benchmark for looping over indices stored in bit sets.
 */

#include <fbitset.hpp>
#include <vector>

#include <testloop.hpp>

using namespace fbitset;
using Bits = Fbitset<2, size_t, No_ext>;

int main()
{
    for (size_t i = 0; i < N_REPEATS; ++i) {
        Bits bits(128);
        for (auto i : idxes) {
            bits.set(i);
        }

        for (auto it = bits.begin(); it; ++it) {
            sink(*it);
        }
    }

    return 0;
}
