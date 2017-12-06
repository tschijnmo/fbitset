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

        T cont;

        static constexpr bool allow_ext = true;

        /** Default constructor.
         *
         * This constructs the container by default.
         */

        Ext()
            : cont{}
        {
        }

        /** Test if the container is holding any bits.
         */

        explicit operator bool() const { return cont.size() > 0; }
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

    /** Finds the index of the first set bit.
     *
     * The input cannot be zero.
     */

    template <typename T> Size fls(T x)
    {
        assert(x != 0);
        return std::numeric_limits<T>::digits - clz(x) - 1;
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

    using Ext_container = E;

    /** Initializes to an all-false bit set of a given size.
     */

    Fbitset(Size size)
        : size_(size)
    {
        if (size <= MAX_BITS) {
            for (Size i = 0; i < N; ++i) {
                limbs_[i] = 0;
            }
        } else {
            if constexpr (allow_ext) {
                Size n_limbs = (size + LIMB_BITS - 1) / LIMB_BITS;
                ext_.cont.assign(n_limbs, 0);
            } else {
                assert(0);
            }
        }
    }

    // All Special functions are from the default implementation.
    //
    // The default gang-of-five should work just fine.

    //
    // General operations.
    //

    /** Gets the size of the bit set.
     */

    Size size() const { return size_; }

    /** Makes equality comparison.
     */

    friend bool operator==(const Fbitset& o1, const Fbitset& o2)
    {
        assert(o1.size() == o2.size());

        if constexpr (!allow_ext) {
            return o1.limbs_ == o2.limbs_;
        } else {
            if (o1.ext_) {
                return o1.ext_.cont == o2.ext_.cont;
            } else {
                return o1.limbs_ == o2.limbs_;
            }
        }
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
        if constexpr (!allow_ext) {
            return hash(limbs_.cbegin(), N);
        } else {
            if (ext_) {
                return hash(ext_.cont.cbegin(), ext_.cont.size());
            } else {
                return hash(limbs_.cbegin(), N);
            }
        }
    }

    /** Sets a given bit.
     */

    void set(Size idx)
    {
        get_limb(idx) |= get_mask(idx);
        return;
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

private:
    //
    // Internal core functions.
    //
    // Generic utilities.
    //

    /** Gets the limb index for a given bit index.
     */

    static Size get_lidx(Size idx) { return idx / LIMB_BITS; }

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
                return ext_.cont.size();
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
        assert(idx < size_);
        return get_limb_lidx(get_lidx(idx));
    }

    Limb& get_limb(Size idx)
    {
        assert(idx < size_);
        return get_limb_lidx(get_lidx(idx));
    }

    /** Gets the limb at a particular limb index.
     */

    const Limb& get_limb_lidx(Size lidx) const
    {
        assert(lidx < get_n_limbs());

        if constexpr (!allow_ext) {
            if constexpr (N_LIMBS == 1) {
                assert(lidx == 0);
                return limbs_[0];
            } else {
                return limbs_[lidx];
            }
        } else {
            if (ext_) {
                return ext_.cont[lidx];
            } else {
                return limbs_[lidx];
            }
        }
    }

    Limb& get_limb_lidx(Size lidx)
    {
        return const_cast<Limb&>(
            static_cast<const Fbitset*>(this)->get_limb_lidx(lidx));
    }

    /** Apply the given callable to matching pairs of limbs.
     */

    template <typename T> void zip_limbs(const Fbitset& other, T act)
    {
        if constexpr (!allow_ext) {
            for (Size i = 0; i < N_LIMBS; ++i) {
                act(limbs_[i], other.limbs_[i]);
            }
        } else {
            if (ext_) {
                auto i = ext_.cont.begin();
                auto j = other.ext_.cont.cbegin();
                while (i != ext_.cont.end()) {
                    act(*i, *j);
                    ++i;
                    ++j;
                }
            } else {
                for (Size i = 0; i < N_LIMBS; ++i) {
                    act(limbs_[i], other.limbs_[i]);
                }
            }
        }
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

    /** The number of bits.
     */

    Size size_;

    /** The actual limbs stored in-place.
     */

    std::array<Limb, N> limbs_;

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
