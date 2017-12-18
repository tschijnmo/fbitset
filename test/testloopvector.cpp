/** Benchmark for looping over indices stored in a vector.
 *
 * This benchmark mainly serves as a comparison for looping over indices stored
 * in bit sets.
 */

#include <fbitset.hpp>
#include <vector>

#include <testloop.hpp>

using namespace fbitset;

int main()
{
    for (size_t i = 0; i < N_REPEATS; ++i) {
        // Make a local copy.
        auto vec = idxes;

        for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
            sink(*it);
        }
    }

    return 0;
}
