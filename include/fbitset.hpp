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
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace fbitset {

/** The placeholder class to disable the usage of external containers.
 */
struct No_ext {
};

/** The integral type used for bit/limb indices and sizes.
 *
 * Here we just use the native size_t.
 */
using Size = size_t;

namespace internal {
    /** Utility for checking if the container is `No_ext`.
     */
    template <typename T> constexpr bool is_no_ext = false;
    template <> constexpr bool is_no_ext<No_ext> = true;

    //
    // Misc bit operations
    //
    // TODO: Make these operations cross-platform by conditional compilation of
    // intrinsic functions and possibly a fallback mode.
    //
    // TODO: Investigate a SIMD-based solution to these problems for
    // SSE/AVX/Neon instructions.
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

/** The core base class.
 *
 * This base class handles the basic data storage and resource management.
 * Derived class can transparently access the limbs via the access methods.
 */
template <Size N, typename L, typename E> class Fbitset_base {
public:
    /** The type used for the limbs.
     */
    using Limb = L;

    /** The external container type.
     */
    using Ext_limbs = E;

    /** The number of limbs allowed inside.
     */
    static constexpr Size N_LIMBS = N ? N : sizeof(Ext_limbs) / sizeof(Limb);

    /** The number of bits hold by a single limb.
     */
    static constexpr Size LIMB_BITS = std::numeric_limits<L>::digits;

    /** The maximum number of bits able to be held in-place.
     */
    static constexpr Size MAX_BITS = N_LIMBS * LIMB_BITS;

    /** The internal container type.
     */
    using Int_limbs = std::array<Limb, N_LIMBS>;

    // Check the sensibility of the given types.
    static_assert(std::numeric_limits<L>::is_integer
            && !std::numeric_limits<L>::is_signed,
        "Limb needs to be unsigned integral type");

    // static_assert(sizeof(Int_limbs) >= sizeof(Ext_limbs),
    //     "Too small number of limbs allowed inside, wasting space");

    /** If external storage is not allowed at all.
     *
     * This value will be used with `if constexpr` liberally for performance
     * boosting when we are sure that the data are inside.  Then all branching
     * about the decision if the bits is inside/outside can be skipped.
     */
    static constexpr bool NO_EXT = internal::is_no_ext<Ext_limbs>;

    /** The basic constructor.
     */
    Fbitset_base(Size size)
        : size_(size)
    {
        if constexpr (NO_EXT) {
            assert(size <= MAX_BITS);
            limbs_.int_.fill(0);
        } else {
            if (is_inplace()) {
                limbs_.int_.fill(0);
            } else {
                new (&limbs_.ext_) Ext_limbs(n_limbs(), 0);
            }
        }
    }

    /** The number of stored bits.
     */
    inline Size size() const { return size_; }

    /** Are the bits stored inplace in the object?
     */
    inline bool is_inplace() const
    {
        if constexpr (NO_EXT) {
            return true;
        } else {
            return size_ <= MAX_BITS;
        }
    }

    /** The number of limbs needed to store the bits.
     */
    inline Size n_limbs() const { return (size() + LIMB_BITS - 1) / LIMB_BITS; }

    /** Default constructor.
     */
    Fbitset_base()
        : Fbitset_base(0)
    {
    }

    /** Copy constructor.
     */
    Fbitset_base(const Fbitset_base& o)
        : size_(o.size())
    {
        init(o);
    }

    /** Move constructor.
     */
    Fbitset_base(Fbitset_base&& o)
        : size_(o.size())
    {
        init(std::move(o));
    }

    /** Destructor.
     */
    ~Fbitset_base()
    {
        if constexpr (NO_EXT)
            return;
        if (!is_inplace()) {
            limbs_.ext_.~E();
        }
    }

    /** Copy assignment.
     */
    Fbitset_base& operator=(const Fbitset_base& o)
    {
        if (this == &o)
            return *this;

        this->~Fbitset_base();
        init(o);
        return *this;
    }

    /** Move assignment.
     */
    Fbitset_base& operator=(Fbitset_base&& o)
    {
        if (this == &o)
            return *this;

        this->~Fbitset_base();
        init(std::move(o));
        return *this;
    }

protected:
    //
    // Core visiting functions.
    //
    // These functions aims to abstract away the actual data layout of the bit
    // sets, so that subclass can access them transparently.  Note that they
    // are generally written as template functions to work automatically with
    // both const and non-const this, while preserving the const-correctness.
    // Generally pointers are used, instead of references, since the argument
    // is most likely `this`.
    //

    /** Gets a pointer to the limbs.
     *
     * This might be the most intuitive way to access the limbs, but note that
     * it is almost surely not the most efficient way to loop limbs over,
     * specially when the limbs are stored in-place.  See functions like
     * `exec_limbs`, which can make it much easier for the compilers to unroll
     * the loop for in-place storage.  This function is more for random access
     * patterns.
     */
    template <typename B> static auto limbs(B* o)
    {
        static_assert(std::is_base_of_v<Fbitset_base, std::decay_t<B>>);
        auto& first_int = o->limbs_.int_[0];
        if constexpr (NO_EXT) {
            return &first_int;
        } else {
            return &(o->is_inplace() ? first_int : o->limbs_.ext_[0]);
        }
    }

    /** Takes an unary action on the limbs.
     *
     * @param act a callable which can be called for both `Int_limbs` arguments
     * or `Ext_limbs` arguments.  For looping over limbs, range-based loop or
     * using the `size()` of the argument is recommended.  `n_limbs()` may
     * necessitate additional computing of the number of limbs and branching.
     * Since generally the in-place limbs can be stored in a single cache line,
     * when there are trailing zero limbs, the additional processing of the
     * suffix generally incurs little overhead.  For instance, to compute the
     * or of two limbs stored in place, with external storage disabled, g++
     * 9.2.0 gives
     *
     * ```asm
     * 	movl	(%rsi), %ecx
     * 	movdqu	8(%rsi), %xmm0
     * 	movq	%rdi, %rax
     * 	movl	%ecx, (%rdi)
     * 	movq	8(%rsi), %rdi
     * 	movups	%xmm0, 8(%rax)
     * 	orq	8(%rdx), %rdi
     * 	movq	16(%rsi), %rcx
     * 	movq	%rdi, 8(%rax)
     * 	orq	16(%rdx), %rcx
     * 	movq	%rcx, 16(%rax)
     * 	ret
     * ```
     *
     * compared with the same operation with `n_limbs()`,
     *
     * ```asm
     * movl	(%rsi), %ecx
     * movdqu	8(%rsi), %xmm0
     * movq	%rdi, %rax
     * movl	%ecx, (%rdi)
     * movups	%xmm0, 8(%rdi)
     * testl	%ecx, %ecx
     * jle	.L18
     * addl	$63, %ecx
     * movq	8(%rdx), %rsi
     * orq	%rsi, 8(%rdi)
     * sarl	$6, %ecx
     * cmpl	$1, %ecx
     * je	.L18
     * movq	16(%rdx), %rsi
     * orq	%rsi, 16(%rdi)
     * cmpl	$2, %ecx
     * je	.L18
     * movq	24(%rdx), %rsi
     * orq	%rsi, 24(%rdi)
     * cmpl	$3, %ecx
     * je	.L18
     * movl	$3, %esi
     * .L20:
     * movq	8(%rdx,%rsi,8), %rdi
     * orq	%rdi, 8(%rax,%rsi,8)
     * addq	$1, %rsi
     * cmpl	%esi, %ecx
     * jg	.L20
     * .L18:
     * ret
     * ```
     *
     * Clearly, for example, when the second limb is zero because of a size
     * less than 65, the additional cost to skip the processing of the second
     * limb dynamically may not be worthy of the effort.
     */
    template <typename B, typename U>
    static decltype(auto) exec_limbs(B* o, U act)
    {
        static_assert(std::is_base_of_v<Fbitset_base, std::decay_t<B>>);

        if constexpr (NO_EXT) {
            return act(o->limbs_.int_);
        } else {
            if (o->is_inplace()) {
                return act(o->limbs_.int_);
            } else {
                return act(o->limbs_.ext_);
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
        static_assert(std::is_base_of_v<Fbitset_base, std::decay_t<B>>);
        static_assert(std::is_base_of_v<Fbitset_base, std::decay_t<C>>);
        assert(o1->size() == o2->size());

        if constexpr (NO_EXT) {
            return act(o1->limbs_.int_, o2->limbs_.int_);
        } else {
            if (o1->is_inplace()) {
                assert(o2->is_inplace());
                return act(o1->limbs_.int_, o2->limbs_.int_);
            } else {
                assert(!o2->is_inplace());
                return act(o1->limbs_.ext_, o2->limbs_.ext_);
            }
        }
    }

    /** Apply the given callable to matching pairs of limbs.
     */
    template <typename T> void zip_limbs(const Fbitset_base& o, T act)
    {
        exec_limbs(
            this, &o, [&act, this](auto& self, const auto& other) -> void {
                for (Size i = 0; i < self.size(); ++i) {
                    act(self[i], other[i]);
                }
            });

        return;
    }

private:
    /** The number of bits inside.
     */
    Size size_;

    /** The core container.
     *
     * The constructor and destructor does nothing, the resource management is
     * handled at Fbitset_base.
     */
    union Limbs_ {
        Int_limbs int_;
        Ext_limbs ext_;

        Limbs_() {}

        ~Limbs_() {}
    } limbs_;

    /** Initialize the current bit set from another.
     *
     * The T type must be the same Fbitset_base type.  Copy will be performed
     * for in-place storage.  For external storage, the content will be copied
     * or moved according to the value category of `o`.
     */
    template <typename T> void init(T&& o)
    {
        static_assert(std::is_same_v<Fbitset_base, std::decay_t<T>>,
            "The same Fbitset type is expected");

        if (o.is_inplace()) {
            limbs_.int_ = o.limbs_.int_;
        } else {
            if constexpr (std::is_lvalue_reference_v<T>) {
                // o is an L-value.
                new (&limbs_.ext_) E(o.limbs_.ext_);
            } else {
                new (&limbs_.ext_) E(std::move(o.limbs_.ext_));
            }
        }
    }
};

/** A highly optimized bit set.
 *
 * This data type aims at a compromise between efficiency and generality
 * towards different problems.  It is basically a mixture of the bitset
 * container in STL and the dynamic_bitset in boost.  A size can be given at
 * the compile-time.  When the number of bits needed is no greater than it,
 * it will be stored right into the container without external memory
 * allocation. Otherwise, bits will be stored in a given container type,
 * which allocates external memory.  This is essentially short-string
 * optimization applied to bit sets.
 *
 * In addition to the short-string optimization, native machine instructions
 * for some bit operations, like the x86-64 `lzcnt`, will be utilized for deep
 * optimization.
 *
 * This class is mostly for the convenience of solving combinatorial
 * problems. So the interface is not fully compatible with standard bitset
 * or boost dynamic_bitset.  For similar reasons, all binary operations,
 * including equality comparison, can only be performed on bit sets of the
 * same type (compile-time check) and the same size (run-time check).
 *
 * This class also aims to be highly configurable.  Internally, the bits are
 * organized into limbs, which are stored in a little-endian format either
 * inside the object or in external containers.  The details of data storage
 * and resource management are all defined in the class `Fbitset_base`.
 *
 * @tparam N The number of limbs allowed inside the object directly (not the
 * number of bits).  The total size of the limbs should not be less than the
 * total size of the external container for maximum space efficiency.  By
 * the default value of 0, limbs will be stored inside if the place for the
 * external storage can accommodate them, which is common for short-string
 * optimization.
 *
 * @tparam L The actual integral type to be used for the limbs.
 *
 * @tparam C The sequence container to be used for the limbs when the limbs
 * inside the object cannot take all the bits.  Or `No_ext` can also be
 * given to completely disable the usage of external container, where any
 * attempt to create bit sets that does not fit will cause assertion error.
 *
 */
template <Size N = 0, typename L = unsigned long long,
    typename E = std::vector<L>>
class Fbitset : public Fbitset_base<N, L, E> {
public:
    /** Convenient name for the base.
     */
    using Base = Fbitset_base<N, L, E>;

    // Forward some static things from the base class.
    using Limb = typename Base::Limb;
    static constexpr auto N_LIMBS = Base::N_LIMBS;
    static constexpr auto LIMB_BITS = Base::LIMB_BITS;
    static constexpr auto MAX_BITS = Base::MAX_BITS;
    static constexpr auto NO_EXT = Base::NO_EXT;

    /** Iterators for iterator over the indices of the set bits.
     *
     * Note that this iterator class does not satisfy the C++ iterator concept.
     * Rather than having a sentinel value to compare for testing the end, here
     * we can explicitly evaluate the truth value of the iterator to see if it
     * still has a values, similar to the Java-style `hasNext()`.
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
            : curr_{ Base::limbs(&fbitset) }
            , last_{ curr_ + fbitset.n_limbs() }
        {
            if (*this) {
                limb_ = *curr_;
            }
            get_next();
        }

        /** Is the iterator still pointing to a valid element.
         */
        explicit operator bool() const noexcept { return curr_ < last_; }

        /** Gets the index of the current index of the set bit.
         *
         * PR-value of the Size type will be directly returned.
         */
        Size operator*() const noexcept
        {
            assert(*this && limb_ != 0);
            Size curr_idx = internal::ctz(limb_);
            return curr_idx + base_;
        }

        /** Advances the iterator.
         */
        const_iterator& operator++() noexcept
        {
            assert(*this && limb_ != 0);
            Size curr_idx = internal::ctz(limb_);
            Limb mask = Limb(1) << curr_idx;
            assert((limb_ & mask) != 0);
            limb_ ^= mask;
            get_next();
            return *this;
        }

    private:
        /** Move the current pointer to the next non-zero limb.
         *
         * The copy of the current limb and base will also be updated to keep
         * the invariant.
         */
        void get_next() noexcept
        {
            while (curr_ < last_ && limb_ == 0) {
                ++curr_;
                limb_ = *curr_;
                base_ += LIMB_BITS;
            }
        }

        /** A pointer to the current limb to loop over.
         */
        const Limb* curr_;

        /** The sentinel for looping the limbs over.
         */
        const Limb* last_;

        /** The base to which limb-internal indices will be added to.
         *
         * This will be incremented when we loop over to the next limb.
         */
        Size base_{ 0 };

        /** A copy of the current limb.
         *
         * It always starts with being a copy of the limb pointed to by
         * `curr_`.  Bits in this copy will gradually be toppled to zero during
         * the loop.
         */
        Limb limb_;
    };

    /** Initializes to an all-false bit set of a given size.
     *
     * @param size The number of bits that need to be held.
     */
    Fbitset(Size size, bool set_true = false)
        : Base(size)
    {
        if (set_true) {
            set_all(size);
        }
    }

    /** Constructs the bit set with the bits at the given indices set.
     *
     * The indices should be given as a pair of iterators giving the indices.
     */
    template <typename It>
    Fbitset(Size size, It first_idx, It last_idx)
        : Fbitset(size)
    {
        for (; first_idx != last_idx; ++first_idx) {
            set(*first_idx);
        }
    }

    /** Constructs a bit set with the set-bit indices given in a list.
     *
     * This constructor can be useful for creating bit sets as a simple
     * literal.
     */
    Fbitset(Size size, std::initializer_list<Size> idxes)
        : Fbitset(size, idxes.begin(), idxes.end())
    {
    }

    // All Special functions are from the base class.  The default gang-of-zero
    // should work just fine.

    //
    // General operations.
    //

    /** Makes equality comparison.
     *
     * Only available for bit sets of the *same* type.  Note that two bit sets
     * are equal only if they were constructed to be of the same size, or they
     * shall be considered unequal even if the same set of bits are set.
     */
    friend bool operator==(const Fbitset& o1, const Fbitset& o2) noexcept
    {
        return o1.size() == o2.size()
            && Base::exec_limbs(&o1, &o2,
                [](const auto& i, const auto& j) -> bool { return i == j; });
    }

    friend bool operator!=(const Fbitset& o1, const Fbitset& o2) noexcept
    {
        return !(o1 == o2);
    }

    /** Computes the hash of a subset.
     *
     * To be consistent with the equality comparison, here, two bit sets will
     * hash differently if they are constructed to have different sizes.
     */
    size_t hash() const noexcept
    {
        return Base::exec_limbs(this, [this](const auto& limbs) -> size_t {
            size_t curr = this->size();
            std::hash<Limb> hasher{};
            for (auto i : limbs) {
                // The algorithm is adapted from the boost hash library.
                curr ^= hasher(i) + 0x9e3779b9 + (curr << 6) + (curr >> 2);
            }
            return curr;
        });
    }

    /** Sets a given bit.
     */
    void set(Size idx) noexcept
    {
        get_limb(idx) |= get_mask(idx);
        return;
    }

    /** Sets all given number of lower bits.
     *
     * The lower `num` bits will be all toppled to true.
     */
    void set_all(Size num) noexcept
    {
        assert(num <= this->size());

        Base::exec_limbs(this, [num](auto& limbs) {
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
        Base::exec_limbs(this, [this](auto& limbs) {
            for (auto& i : limbs) {
                i = 0;
            }
        });
    }

    /** Test a given bit.
     *
     * Note that different from STL bitset, here the result is a pr-value of
     * boolean type, rather than a proxy for the actual bit.
     */
    bool operator[](Size idx) const noexcept
    {
        return (get_limb(idx) & get_mask(idx)) != 0;
    }

    /** Flips a given bit.
     */
    void flip(Size idx) noexcept
    {
        get_limb(idx) ^= get_mask(idx);
        return;
    }

    //
    // Arithmetic operations
    //

    /** Computes the bitwise or (union).
     */
    Fbitset& operator|=(const Fbitset& other) noexcept
    {
        Base::zip_limbs(other, [](Limb& i, Limb j) { i |= j; });
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
    Fbitset& operator&=(const Fbitset& other) noexcept
    {
        Base::zip_limbs(other, [](Limb& i, Limb j) { i &= j; });
        return *this;
    }

    Fbitset operator&(const Fbitset& other) const
    {
        Fbitset res(*this);
        res &= other;
        return res;
    }

    /** Computes the bitwise xor (disjunctive union).
     */
    Fbitset& operator^=(const Fbitset& other) noexcept
    {
        Base::zip_limbs(other, [](Limb& i, Limb j) { i ^= j; });
        return *this;
    }

    Fbitset operator^(const Fbitset& other) const
    {
        Fbitset res(*this);
        res ^= other;
        return res;
    }

    //
    // Misc bit operations.
    //

    /** Finds the index of the last (highest) set bit.
     */
    Size find_last() const noexcept
    {
        // For this, we start with n_limbs to skip some zeroes.
        for (Size i = this->n_limbs(); i != 0; --i) {
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
        return Base::exec_limbs(this, [this](const auto& limbs) -> Size {
            Size res = 0;
            for (auto& i : limbs) {
                res += internal::popcount(i);
            }
            return res;
        });
    }

    //
    // Bits as container.
    //

    /** Gets the iterator for the indices of the set bits.
     *
     * The resulted iterator is of type `Fbitset::const_iterator` class. This
     * makes the bit set similar to an iterable container for the indices of
     * the set bits.
     */
    const_iterator begin() const { return const_iterator(*this); }

private:
    //
    // Internal core functions.
    //
    // Generic utilities.
    //

    /** Gets the limb index for a given bit index.
     */
    static Size get_lidx(Size idx) noexcept
    {
        if constexpr (NO_EXT && N_LIMBS == 1) {
            return 0;
        } else {
            return idx / LIMB_BITS;
        }
    }

    /** Gets the limb mask for a given bit index.
     */
    static Limb get_mask(Size idx) noexcept
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
    const Limb& get_limb(Size idx) const noexcept
    {
        assert(idx < this->size());
        return get_limb_lidx(get_lidx(idx));
    }

    Limb& get_limb(Size idx) noexcept
    {
        assert(idx < this->size());
        return get_limb_lidx(get_lidx(idx));
    }

    /** Gets the limb at a particular limb index.
     */
    const Limb& get_limb_lidx(Size lidx) const noexcept
    {
        assert(lidx < this->n_limbs());
        return Base::limbs(this)[lidx];
    }

    Limb& get_limb_lidx(Size lidx) noexcept
    {
        assert(lidx < this->n_limbs());
        return Base::limbs(this)[lidx];
    }
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
