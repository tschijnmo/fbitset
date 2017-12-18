/** The main header file for fbitset.
 *
 * In order to use the library, including this header file is sufficient.
 * Everything resides in the `fbitset` name space.
 */

#ifndef FBITSET_FBITSET_H
#define FBITSET_FBITSET_H

#include <array>
#include <cassert>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

namespace fbitset {

/** The placeholder class to disable the usage of external containers.
 */

struct No_ext {
};

/** The integral type used for bit/limb indices and sizes.
 *
 * Here we use int for performance and compatibility with the `digits` field of
 * numeric_limits.
 */

using Size = int;

namespace internal {

    /** Wrapper over the external container.
     *
     * The given type must be a C++ SequenceContainer.
     */

    template <typename T> struct Ext {
        /** The container.
         */

        T limbs;

        static constexpr bool allow_ext = true;

        /** Default constructor.
         *
         * This constructs the container by default.
         */

        Ext()
            : limbs{}
        {
        }

        /** Test if the container is holding any bits.
         */

        explicit operator bool() const noexcept { return !limbs.empty(); }
    };

    /** Holder for external container when it is disabled.
     *
     * No extra space is taken any more for this case.
     */

    template <> struct Ext<No_ext> {
        static constexpr bool allow_ext = false;
    };

    //
    // Misc bit operations
    //
    // TODO: Make these operations cross-platform by conditional compilation of
    // intrinsic functions and possibly a fallback mode.
    //

    /** Counts the number of leading zeros.
     *
     * The input cannot be zero.
     */

    inline Size clz(unsigned int x) { return __builtin_clz(x); }
    inline Size clz(unsigned long x) { return __builtin_clzl(x); }
    inline Size clz(unsigned long long x) { return __builtin_clzll(x); }

    /** Counts the number of trailing zeros.
     *
     * The input cannot be zero.
     */

    inline Size ctz(unsigned int x) { return __builtin_ctz(x); }
    inline Size ctz(unsigned long x) { return __builtin_ctzl(x); }
    inline Size ctz(unsigned long long x) { return __builtin_ctzll(x); }

    /** Finds the index of the first set bit.
     *
     * The input cannot be zero.
     */

    template <typename T> Size fls(T x)
    {
        assert(x != 0);
        return std::numeric_limits<T>::digits - clz(x) - 1;
    }

    /** Counts the number of set bits.
     */

    inline Size popcount(unsigned int x) { return __builtin_popcount(x); }
    inline Size popcount(unsigned long x) { return __builtin_popcountl(x); }
    inline Size popcount(unsigned long long x)
    {
        return __builtin_popcountll(x);
    }
}

/** A set of bits.
 *
 * The data type aims at a compromise between efficiency and generality towards
 * different problems.  It is basically a mixture of the bitset container in
 * STL and the dynamic_bitset in boost.  A number of limbs can be given at the
 * compile-time.  When the number of bits needed is not greater than what can
 * be hold be the limbs, it will be stored right into the limbs.  Otherwise,
 * memory is going to be allocated on the heap.
 *
 * This class is mostly for the convenience of solving combinatorial problems.
 * So the interface might not be fully compatible with standard bitset or boost
 * dynamic_bitset.  For similar reasons, all binary operations, including
 * equality comparison, can only be performed on bit sets of the same type
 * (compile-time check) and the same size (run-time check).
 *
 * This class also aims to be highly configurable.  Internally, the limbs are
 * stored in little-endian format either inside the object or in external
 * containers.
 *
 * @tparam N The number of limbs allowed inside the object directly (not the
 * number of bits).
 *
 * @tparam L The actual integral type to be used for the limbs.
 *
 * @tparam C The sequence container to be used for the limbs when the limbs
 * inside the object cannot take all the bits.  Or `No_ext` can also be given
 * to completely disable the usage of external container, where any attempt to
 * create bit sets that does not fit will cause assertion error.
 *
 */

template <Size N, typename L = unsigned long long, typename E = std::vector<L>>
class Fbitset {
public:
    // Check the sensibility of the given types.

    static_assert(std::numeric_limits<L>::is_integer
            && !std::numeric_limits<L>::is_signed,
        "Limb needs to be unsigned integral type");

    /** The number of limbs allowed inside.
     */

    static constexpr Size N_LIMBS = N;

    /** The type used for the limbs.
     */

    using Limb = L;

    /** The number of bits hold by a single limb.
     */

    static constexpr Size LIMB_BITS = std::numeric_limits<L>::digits;

    /** The maximum number of bits able to be held in-place.
     */

    static constexpr Size MAX_BITS = N * LIMB_BITS;

    /** The external container type.
     */

    using Ext_limbs = E;

    /** Iterators for iterator over the indices of the set bits.
     *
     * Note that this iterator class does not satisfy the C++ iterator concept.
     * Rather than having a sentinel value to compare for testing the end, here
     * we can explicitly evaluate the truth value of the iterator to see if it
     * still has a values.
     */

    class const_iterator {
    public:
        /** Constructs the iterator for a bit set.
         *
         * This constructor constructs the iterator pointing to the first set
         * bit in the given bit set.  Normally the `begin` method can be used
         * instead.
         */

        const_iterator(const Fbitset& fbitset)
            : fbitset_{ fbitset }
            , curr_limb_{ 0 }
        {
            get_next();
        }

        explicit operator bool() const noexcept
        {
            return curr_limb_ < fbitset_.get_n_limbs();
        }

        Size operator*() const noexcept { return curr_; }

        const_iterator& operator++()
        {
            get_next();
            return *this;
        }

    private:
        void get_next()
        {
            curr_ = exec_limbs(&fbitset_, [this](auto& limbs) -> Size {
                while (*this && limbs[curr_limb_] == 0) {
                    ++curr_limb_;
                }

                if (*this) {
                    auto& limb = limbs[curr_limb_];
                    assert(limb != 0);
                    Size curr_idx = internal::ctz(limb);

                    Limb mask = Limb(1) << curr_idx;
                    assert((limb & mask) != 0);
                    limb ^= mask;

                    return curr_idx + curr_limb_ * LIMB_BITS;
                } else {
                    return -1;
                }
            });
        }

        /** A copy of the original bit set.
         *
         * The set bits in this copy are going to be gradually toppled.
         */

        Fbitset fbitset_;

        /** The current bit index.
         */

        Size curr_;

        /** The current limb index.
         */

        Size curr_limb_;
    };

    /** Initializes to an all-false bit set of a given size.
     *
     * @param size The number of bits that need to be held.
     */

    Fbitset(Size size, bool set_true = false)
    {
        if (size <= MAX_BITS) {
            for (Size i = 0; i < N; ++i) {
                limbs_[i] = 0;
            }
        } else {
            if constexpr (allow_ext) {
                Size n_limbs = (size + LIMB_BITS - 1) / LIMB_BITS;
                ext_.limbs.assign(n_limbs, 0);
            } else {
                assert(0);
            }
        }

        if (set_true) {
            set_all(size);
        }
    }

    // All Special functions are from the default implementation.
    //
    // The default gang-of-five should work just fine.

    //
    // General operations.
    //

    /** Gets the maximum number of bits that can be hold.
     *
     * Note that this size is not necessarily the size of bits given to the
     * constructor, but rather the number of bits that *can be* hold by the bit
     * set.
     */

    Size size() const noexcept { return get_n_limbs() * LIMB_BITS; }

    /** Makes equality comparison.
     */

    friend bool operator==(const Fbitset& o1, const Fbitset& o2)
    {
        return exec_limbs(&o1, &o2,
            [](const auto& i, const auto& j) -> bool { return i == j; });
    }

    friend bool operator!=(const Fbitset& o1, const Fbitset& o2)
    {
        return !(o1 == o2);
    }

    /** Computes the hash of a subset.
     *
     * Note that we assume that only bit sets of the same size are going to be
     * compared.  As a result, the size is *not* put into consideration for
     * this hash function.
     */

    size_t hash() const
    {
        return exec_limbs(this, [this](const auto& limbs) -> size_t {
            return hash(limbs.cbegin(), get_n_limbs());
        });
    }

    /** Sets a given bit.
     */

    void set(Size idx)
    {
        get_limb(idx) |= get_mask(idx);
        return;
    }

    /** Sets all given number of lower bits.
     *
     * The lower `num` bits will be all toppled to true.
     */

    void set_all(Size num)
    {
        assert(num <= size());

        exec_limbs(this, [num](auto& limbs) {
            Size idx = 0;
            auto remain = num;
            while (remain >= LIMB_BITS) {
                limbs[idx] = std::numeric_limits<Limb>::max();
                ++idx;
                remain -= LIMB_BITS;
            }
            if (remain > 0) {
                limbs[idx] |= (static_cast<Limb>(1) << remain) - 1;
            }
        });
    }

    /** Clears all bits inside.
     */

    void clear() noexcept
    {
        exec_limbs(this, [this](auto& limbs) {
            for (Size i = 0; i < get_n_limbs(); ++i) {
                limbs[i] = 0;
            }
        });
    }

    /** Test a given bit.
     *
     * Note that different from STL bitset, here the result is a pr-value of
     * boolean type, rather than a proxy for the actual bit.
     */

    bool operator[](Size idx) const
    {
        return (get_limb(idx) & get_mask(idx)) != 0;
    }

    /** Flips a given bit.
     */

    void flip(Size idx)
    {
        get_limb(idx) ^= get_mask(idx);
        return;
    }

    //
    // Arithmetic operations
    //

    /** Computes the bitwise or (union).
     */

    Fbitset& operator|=(const Fbitset& other)
    {
        zip_limbs(other, [](Limb& i, Limb j) { i |= j; });
        return *this;
    }

    Fbitset operator|(const Fbitset& other) const
    {
        Fbitset res(*this);
        res |= other;
        return res;
    }

    /** Computes the bitwise and (intersection).
     */

    Fbitset& operator&=(const Fbitset& other)
    {
        zip_limbs(other, [](Limb& i, Limb j) { i &= j; });
        return *this;
    }

    Fbitset operator&(const Fbitset& other) const
    {
        Fbitset res(*this);
        res &= other;
        return res;
    }

    //
    // Misc bit operations.
    //

    /** Finds the index of the last (highest) set bit.
     */

    Size find_last() const
    {
        for (Size i = get_n_limbs(); i != 0; --i) {
            const auto& limb = get_limb_lidx(i - 1);
            if (limb != 0) {
                Size idx = internal::fls(limb);
                return idx + LIMB_BITS * (i - 1);
            }
        }
        return -1;
    }

    /** Counts the number of all set bits inside the set.
     */

    Size count() const noexcept
    {
        return exec_limbs(this, [this](const auto& limbs) -> Size {
            Size res = 0;
            for (Size i = 0; i < get_n_limbs(); ++i) {
                res += internal::popcount(limbs[i]);
            }
            return res;
        });
    }

private:
    //
    // Internal core functions.
    //
    // Generic utilities.
    //

    /** Gets the limb index for a given bit index.
     */

    static Size get_lidx(Size idx)
    {
        if constexpr (!allow_ext && N_LIMBS == 1) {
            assert(idx < LIMB_BITS);
            return 0;
        } else {
            return idx / LIMB_BITS;
        }
    }

    /** Get the number of limbs.
     *
     * Note that for bits fit inside the object, the result will always be all
     * the limbs hold by the object.
     */

    Size get_n_limbs() const
    {
        if constexpr (!allow_ext) {
            return N_LIMBS;
        } else {
            if (ext_) {
                return ext_.limbs.size();
            } else {
                return N_LIMBS;
            }
        }
    }

    /** Gets the limb mask for a given bit index.
     */

    static Limb get_mask(Size idx)
    {
        Limb res = 1;
        res <<= idx % LIMB_BITS;
        return res;
    }

    /** Gets the limb for a particular bit index.
     *
     * This method aims at as little overhead over the raw use of the
     * underlying integral type as possible.
     */

    const Limb& get_limb(Size idx) const
    {
        assert(idx < size());
        return get_limb_lidx(get_lidx(idx));
    }

    Limb& get_limb(Size idx)
    {
        assert(idx < size());
        return get_limb_lidx(get_lidx(idx));
    }

    /** Gets the limb at a particular limb index.
     */

    const Limb& get_limb_lidx(Size lidx) const
    {
        assert(lidx < get_n_limbs());

        return exec_limbs(this,
            [lidx](auto& limbs) -> decltype(auto) { return limbs[lidx]; });
    }

    Limb& get_limb_lidx(Size lidx)
    {
        return const_cast<Limb&>(
            static_cast<const Fbitset*>(this)->get_limb_lidx(lidx));
    }

    //
    // Core visiting functions.
    //
    // These functions aims to abstract away from the actual data layout of the
    // bit sets.  Usually a generic callable needs to be given, which can treat
    // both `Limbs` arguments or `Ext_limbs` arguments.  For looping over
    // limbs, normally we loop an index from zero up to `get_n_limbs()`.  It
    // has been shown that both g++ and clang++ optimize this loop better than
    // the begin/end iterator pair paradigm.
    //

    /** Takes an unary action on the limbs.
     */

    template <typename B, typename U>
    static decltype(auto) exec_limbs(B* o, U act)
    {
        static_assert(std::is_same_v<std::decay_t<B>, Fbitset>);

        if constexpr (!allow_ext) {
            return act(o->limbs_);
        } else {
            if (o->ext_) {
                return act(o->ext_.limbs);
            } else {
                return act(o->limbs_);
            }
        }
    }

    /** Takes a binary operation on the limbs.
     *
     * The two bit sets need to be of the same type and have the same size.
     */

    template <typename B, typename C, typename U>
    static decltype(auto) exec_limbs(B* o1, C* o2, U act)
    {
        static_assert(std::is_same_v<std::decay_t<B>, Fbitset>);
        static_assert(std::is_same_v<std::decay_t<C>, Fbitset>);
        assert(o1->size() == o2->size());

        if constexpr (!allow_ext) {
            return act(o1->limbs_, o2->limbs_);
        } else {
            if (o1->ext_) {
                assert(o2->ext_);
                return act(o1->ext_.limbs, o2->ext_.limbs);
            } else {
                assert(!o2->ext_);
                return act(o1->limbs_, o2->limbs_);
            }
        }
    }

    /** Apply the given callable to matching pairs of limbs.
     */

    template <typename T> void zip_limbs(const Fbitset& o, T act)
    {
        exec_limbs(
            this, &o, [&act, this](auto& self, const auto& other) -> void {
                for (Size i = 0; i < get_n_limbs(); ++i) {
                    act(self[i], other[i]);
                }
            });

        return;
    }

    //
    // Hash related.
    //

    /** Combines the given hash values.
     *
     * The algorithm is adapted from the boost hash library.
     */

    static inline void combine_hash(size_t& seed, size_t value)
    {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    template <typename It> static size_t hash(It first, Size n)
    {
        size_t curr = 0;
        std::hash<Limb> hasher{};
        for (Size i = 0; i < n; ++i) {
            combine_hash(curr, hasher(*first));
            ++first;
        }
        return curr;
    }

    //
    // Internal data fields.
    //

    /** The container for in-place limbs.
     */

    using Limbs = std::array<Limb, N>;

    /** The actual limbs stored in-place.
     */

    Limbs limbs_;

    // The fall-back external container.

    using Ext = internal::Ext<E>;

    static constexpr bool allow_ext = Ext::allow_ext;

    Ext ext_{};
};
}

// Inject the hash function into the std namespace.

namespace std {

/** Computes the hash of a subset.
 */

template <fbitset::Size N, typename L, typename E>
struct hash<fbitset::Fbitset<N, L, E>> {
    size_t operator()(const fbitset::Fbitset<N, L, E>& o) const
    {
        return o.hash();
    }
};
}

#endif
