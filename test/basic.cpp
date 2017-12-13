/** Tests for the basic behaviour of Fbitsets.
 *
 * Due to the structure of (compile-time and run-time) branching inside the
 * implementation, the behaviour should be tested for four different cases,
 *
 * - No external, single limb,
 * - No external, multiple limb,
 * - Allow external, using some limbs,
 * - Allow external, using all limbs,
 * - Allow external, using external.
 *
 */

#include <cstdint>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <catch.hpp>

#include <fbitset.hpp>

using namespace fbitset;

TEST_CASE("Fbitset has basic behaviour")
{
    constexpr Size N_BITS = 64;
    Fbitset<1, uint64_t, No_ext> ne_1l(N_BITS);
    Fbitset<2, uint32_t, No_ext> ne_ml(N_BITS);
    Fbitset<2, uint64_t, std::vector<uint64_t>> ae_sl(N_BITS);
    Fbitset<2, uint32_t, std::vector<uint32_t>> ae_al(N_BITS);
    Fbitset<1, uint32_t, std::vector<uint32_t>> ae_e(N_BITS);

    auto run_on_all = [&](auto act) {
        act(ne_1l);
        act(ne_ml);
        act(ae_sl);
        act(ae_al);
        act(ae_e);
    };

    SECTION("are initialized to all false by default")
    {
        run_on_all([&](const auto& i) {
            CHECK(i.size() >= N_BITS);
            for (Size j = 0; j < N_BITS; ++j) {
                CHECK_FALSE(i[j]);
            }
            CHECK(i.count() == 0);
        });
    }

    SECTION("can be initialized to all true")
    {
        run_on_all([&](const auto& i) {
            using Curr = std::decay_t<decltype(i)>;

            // This could give very good coverage of different branches.
            for (Size n_set : { 0, 16, 32, 48, 64 }) {
                Curr curr(n_set, true);
                CHECK(curr.size() >= n_set);
                for (Size j = 0; j < curr.size(); ++j) {
                    if (j < n_set) {
                        CHECK(curr[j]);
                    } else {
                        CHECK_FALSE(curr[j]);
                    }
                }
                CHECK(curr.count() == n_set);
            }
        });
    }

    SECTION("have correct copy/move constructor/assignment")
    {
        constexpr Size STRIDE = 5;

        auto check_bits = [&](const auto& bits) {
            for (Size i = 0; i < N_BITS; ++i) {
                if (i % STRIDE == 0) {
                    CHECK(bits[i]);
                } else {
                    CHECK_FALSE(bits[i]);
                }
            }
        };

        run_on_all([&](auto& orig) {

            using Bitset = std::remove_reference_t<decltype(orig)>;
            std::hash<Bitset> hash;

            // Setting of bits.
            for (Size i = 0; i < N_BITS; ++i) {
                if (i % STRIDE == 0) {
                    orig.set(i);
                }
            }
            check_bits(orig);

            // Copy initialization.
            Bitset copied(orig);
            check_bits(copied);
            CHECK(copied == orig);
            CHECK(hash(copied) == hash(orig));

            // Move initialization.
            Bitset moved(std::move(copied));
            check_bits(moved);
            CHECK(moved == orig);
            CHECK(hash(moved) == hash(orig));

            // Copy assignment.
            Bitset copy_assigned(N_BITS);
            CHECK(copy_assigned != orig);
            CHECK(hash(copy_assigned) != hash(orig));
            copy_assigned = orig;
            check_bits(copy_assigned);
            CHECK(copy_assigned == orig);
            CHECK(hash(copy_assigned) == hash(orig));

            // Move assignment.
            Bitset move_assigned(N_BITS);
            CHECK(move_assigned != orig);
            CHECK(hash(move_assigned) != hash(orig));
            move_assigned = std::move(copy_assigned);
            check_bits(move_assigned);
            CHECK(move_assigned == orig);
            CHECK(hash(move_assigned) == hash(orig));
        });
    }

    SECTION("given number of lower bits can be set to true and cleared")
    {
        run_on_all([](auto& inp) {
            for (Size i = 0; i <= N_BITS; i += 16) {
                inp.set(63);
                inp.set_all(i);

                for (Size j = 0; j < N_BITS; ++j) {
                    if (j < i || j == 63) {
                        CHECK(inp[j]);
                    } else {
                        CHECK_FALSE(inp[j]);
                    }
                }

                inp.clear();
                for (Size j = 0; j < N_BITS; ++j) {
                    CHECK_FALSE(inp[j]);
                }
            }
        });
    }

    SECTION("bits can be flipped")
    {
        run_on_all([&](auto& inp) {
            constexpr Size S1 = 3;
            constexpr Size S2 = 5;
            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S1 == 0) {
                    inp.set(i);
                }
            }

            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S2 == 0) {
                    inp.flip(i);
                }
            }

            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S1 == 0) {
                    if (i % S2 == 0) {
                        CHECK_FALSE(inp[i]);
                    } else {
                        CHECK(inp[i]);
                    }
                } else if (i % S2 == 0) {
                    CHECK(inp[i]);
                } else {
                    CHECK_FALSE(inp[i]);
                }
            }
        });
    }

    SECTION("has correct union and intersection operation")
    {
        run_on_all([&](auto i) {
            using Bitset = decltype(i);

            constexpr Size S1 = 3;
            constexpr Size S2 = 5;
            Bitset bits1(N_BITS);
            Bitset bits2(N_BITS);

            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S1 == 0) {
                    bits1.set(i);
                }
                if (i % S2 == 0) {
                    bits2.set(i);
                }
            }

            auto union_bits = bits1 | bits2;
            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S1 == 0 || i % S2 == 0) {
                    CHECK(union_bits[i]);
                } else {
                    CHECK_FALSE(union_bits[i]);
                }
            }

            auto inters_bits = bits1 & bits2;
            for (Size i = 0; i < N_BITS; ++i) {
                if (i % S1 == 0 && i % S2 == 0) {
                    CHECK(inters_bits[i]);
                } else {
                    CHECK_FALSE(inters_bits[i]);
                }
            }
        });
    }

    SECTION("can find last set bit correctly")
    {
        run_on_all([&](auto& i) {
            i.set(2);
            i.set(60);

            CHECK(i.count() == 2);

            auto high = i.find_last();
            CHECK(high == 60);
            i.flip(60);

            auto next_high = i.find_last();
            CHECK(next_high == 2);
            i.flip(2);

            auto fin = i.find_last();
            CHECK(fin == -1);
        });
    }
}
