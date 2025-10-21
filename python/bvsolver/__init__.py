from __future__ import annotations
from collections.abc import Generator
import random
from typing import Generic, TypeVar, overload
import bvsolver._lib

# Solve equations on bit vectors,
# e.g. 4 bits, a to d, and we know that:
# a ^ 1 == 0
# a ^ b ^ 1 == 0
# a ^ c == 0
# a ^ d ^ == 0
# We can solve: [a, b, c, d] == [1, 0, 1, 0]
# The equations can only have ^, no other operations
# To represent values computed from the unknown bits,
# Use integer to represent each bit (including the constant 1):
# 1 -> 1, a -> 2, b -> 4, c -> 8, d -> 16
# then, a ^ b -> 2 ^ 4
# So each intermediate bit is mapped to an integer
# Each bit vector corresponds to a integer,
# where each bit correponds to a bit of solver


class BitVector:
    # from LSB to MSB
    _bits: list[int]

    def __init__(self, bits: list[int]) -> None:
        self._bits = bits.copy()

    def __repr__(self) -> str:
        bits = ", ".join([str(bit) for bit in self._bits])
        return f"<BitVector [{bits}]>"

    def __xor__(self, other: BitVector | int) -> BitVector:
        res: list[int] = []
        if isinstance(other, BitVector):
            new_len = max(len(self._bits), len(other._bits))
            res = [0] * new_len
            for i, b in enumerate(self._bits):
                res[i] ^= b
            for i, b in enumerate(other._bits):
                res[i] ^= b
        else:
            assert other.bit_length() <= len(self._bits)
            for i, b in enumerate(self._bits):
                res.append(b ^ (1 if other & (1 << i) != 0 else 0))
        return BitVector(res)

    def __and__(self, other: int) -> BitVector:
        res: list[int] = []
        for i in range(min(len(self._bits), other.bit_length())):
            res.append(self._bits[i] if other & (1 << i) != 0 else 0)
        return BitVector(res)

    def __mul__(self, other: int) -> BitVector:
        # only one bit is supported
        assert len(self._bits) == 1
        res: list[int] = []
        for i in range(other.bit_length()):
            res.append(self._bits[0] if other & (1 << i) != 0 else 0)
        return BitVector(res)

    def __rshift__(self, other: int) -> BitVector:
        assert 0 <= other <= len(self._bits)
        return BitVector(self._bits[other:])

    def __lshift__(self, other: int) -> BitVector:
        return BitVector([0] * other + self._bits)

    def get(self, solution: int) -> int:
        res = 0
        for i, b in enumerate(self._bits):
            # odd number of 1s
            if (solution & b).bit_count() & 1:
                res += 1 << i
        return res


class Solver:
    _bitvectors: list[BitVector]
    _size: int

    def __init__(self, sizes: list[int]) -> None:
        """Create a solver to solve some BitVector, the width of them are passed as argument"""

        # assign integer to each bit of each bitvector
        i = 1  # reserved for constant 1
        self._bitvectors = []
        for size in sizes:
            self._bitvectors.append(BitVector([1 << (i + j) for j in range(size)]))
            i += size
        self._size = i - 1

    def __repr__(self) -> str:
        return f"<Solver bitvectors={self._bitvectors}, size={self._size}>"

    def bitvectors(self) -> list[BitVector]:
        return self._bitvectors

    def solve(self, zeros: list[BitVector], use_python=False) -> Generator[int]:
        # flatten
        flatten: list[int] = []
        for zero in zeros:
            flatten += zero._bits

        if use_python:
            # gauss elimination
            done = 0
            # for each bit, record its row with the bit set to 1
            mapping = [0] * self._size
            for i in range(len(flatten)):
                e = flatten[i]
                if e == 1:
                    # no solution
                    return None
                # elminate with existing rows in reduced row echelon form
                added = False
                e_no_low = e ^ (e & 1)
                while e_no_low != 0:
                    lowbit = e_no_low & (-e_no_low)
                    # minus 2: 0b10 -> mapping index 0
                    bit = lowbit.bit_length() - 2
                    if mapping[bit] != 0:
                        e ^= mapping[bit]
                        e_no_low = e ^ (e & 1)
                    else:
                        mapping[bit] = e
                        done += 1
                        added = True
                        break
                if not added and e != 0:
                    # eliminate to non-zero, no solutions
                    return None

            # find free variables
            free = []
            for i in range(self._size):
                if mapping[i] == 0:
                    free.append(i)

            # convert to row reduced echelon form and compute solution if all free variables are zero
            base_res = 1  # constant term
            for i in range(self._size):
                for j in range(i + 1, self._size):
                    if mapping[i] & (1 << (j + 1)) != 0:
                        # eliminate
                        mapping[i] ^= mapping[j]
            for i in range(self._size):
                # after elimination, compute solution if free variables are zero
                base_res |= (mapping[i] & 1) << (i + 1)

            # now all free variables are found
            assert done + len(free) == self._size

            # compute contribution of each free variable to solution
            contribution = []
            for bit in free:
                c = 0
                for i in range(self._size):
                    if (mapping[i] & (1 << (bit + 1))) != 0 or i == bit:
                        c |= 1 << (i + 1)
                contribution.append(c)
        else:
            base_res, contribution = bvsolver._lib.solve(self._size, flatten)
            if base_res == -1:
                # no solution
                return None

        # we got 2 ** len(contribution) solutions
        # recover all solutions
        num_sols = 1 << len(contribution)
        sols = 0
        while sols < num_sols:
            # compute solution
            # assign free variables
            res = base_res
            for i in range(len(contribution)):
                if sols & (1 << i) != 0:
                    res ^= contribution[i]

            yield res
            sols += 1
        return None


@overload
def xor_reduce(bv: BitVector) -> BitVector: ...
@overload
def xor_reduce(bv: int) -> int: ...
def xor_reduce(bv):
    if isinstance(bv, BitVector):
        res = 0
        for bit in bv._bits:
            res ^= bit
        return BitVector([res])
    else:
        return bv.bit_count() & 1


T = TypeVar("T", BitVector, int)


class FibonacciLFSR(Generic[T]):
    _width: int
    _poly: int
    _state: T

    def __init__(self, width: int, poly: int, state: T) -> None:
        self._width = width
        self._poly = poly
        self._state = state

    def gen(self) -> T:
        res = self._state & 1
        # shift
        self._state = (self._state >> 1) ^ (
            xor_reduce(self._state & self._poly) << (self._width - 1)
        )
        return res


class MersenneTwister(Generic[T]):
    _index: int
    _state: list[T]
    _w: int
    _n: int
    _m: int
    _r: int
    _a: int
    _u: int
    _d: int
    _s: int
    _b: int
    _t: int
    _c: int
    _l: int
    _f: int

    def __init__(
        self,
        state: list[T],
        index: int,
        w: int,
        n: int,
        m: int,
        r: int,
        a: int,
        u: int,
        d: int,
        s: int,
        b: int,
        t: int,
        c: int,
        l: int,
        f: int,
    ) -> None:
        assert len(state) == n
        self._state = state.copy()
        self._index = index
        self._w = w
        self._n = n
        self._m = m
        self._r = r
        self._a = a
        self._u = u
        self._d = d
        self._s = s
        self._b = b
        self._t = t
        self._c = c
        self._l = l
        self._f = f

    def gen(self) -> T:
        # see genrand_uint32 from _randommodule.c or operator()() from random.tcc
        if self._index >= self._n:
            mask = (1 << self._w) - 1
            upper_mask = (mask << self._r) & mask
            lower_mask = mask ^ upper_mask
            for kk in range(self._n):
                y = (self._state[kk] & upper_mask) ^ (
                    self._state[(kk + 1) % self._n] & lower_mask
                )
                # mag01[y & 0x1U]: if y & 0x1U == 0, then 0; else MATRIX_A
                mag01 = (y & 1) * self._a
                self._state[kk] = (
                    self._state[(kk + self._m) % self._n] ^ (y >> 1) ^ mag01
                )
            self._index = 0

        y = self._state[self._index]
        self._index += 1
        y = y ^ (y >> self._u) & self._d
        y = y ^ (y << self._s) & self._b
        y = y ^ (y << self._t) & self._c
        y = y ^ (y >> self._l)
        return y

    def seed(self, seed: int) -> None:
        # see seed(result_type __sd) in random.tcc, init_genrand in _randomodule.c
        mask = (1 << self._w) - 1
        self._state[0] = seed & mask
        for i in range(1, self._n):
            # mt[mti] =
            #   (1812433253U * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
            self._state[i] = (
                self._f * (self._state[i - 1] ^ (self._state[i - 1] >> (self._w - 2)))
                + i
            ) & mask
        self._index = self._n


class MT19937(MersenneTwister[T]):
    def __init__(self, state: list[T], index: int = 624) -> None:
        MersenneTwister.__init__(
            self,
            state,
            index,
            32,
            624,
            397,
            31,
            0x9908B0DF,
            11,
            0xFFFFFFFF,
            7,
            0x9D2C5680,
            15,
            0xEFC60000,
            18,
            1812433253,
        )


class MT19937_64(MersenneTwister[T]):
    def __init__(self, state: list[T], index: int = 312) -> None:
        MersenneTwister.__init__(
            self,
            state,
            index,
            64,
            312,
            156,
            31,
            0xB5026F5AA96619E9,
            29,
            0x5555555555555555,
            17,
            0x71D67FFFEDA60000,
            37,
            0xFFF7EEE000000000,
            43,
            6364136223846793005,
        )


class CPythonRandom(MT19937[T]):
    def __init__(self, state: list[T], index: int = 624) -> None:
        MT19937.__init__(
            self,
            state,
            index,
        )

    def genrand_uint32(self) -> T:
        # genrand_uint32 from _randommodule.c
        return self.gen()

    def getrandbits(self, bits) -> T:
        # _random_Random_getrandbits_impl in _randommodule.c
        if isinstance(self._state[0], BitVector):
            res = BitVector([0])
        else:
            res = 0

        if bits == 0:
            return res
        elif bits <= 32:
            return self.genrand_uint32() >> (32 - bits)

        words = (bits - 1) // 32 + 1
        for i in range(words):
            r = self.genrand_uint32()
            if bits < 32:
                r = r >> (32 - bits)
            res = res ^ (r << (32 * i))
            bits -= 32
        return res

    def to_cpython_random(self) -> random.Random:
        res = random.Random(0)
        res.setstate((3, (*self._state, self._index), None))
        return res

    @classmethod
    def from_cpython_random(cls, rng: random.Random) -> CPythonRandom:
        state = rng.getstate()
        res = CPythonRandom(list(state[1][:-1]), state[1][-1])
        return res

    @classmethod
    def recover_seed(cls, rng: random.Random) -> int:
        state = rng.getstate()
        N = 624
        assert state[1][-1] == N

        # learned from https://github.com/PKU-GeekGame/geekgame-3rd/blob/279868d814fd8371075e06cf705730c54a3e2a13/official_writeup/prob08-cookie/README.md
        # and https://github.com/jailctf/challenges-2024/blob/master/stupid-crypto-chall/solve/solve.py
        # see init_by_array in _randommodule.c
        mt = list(state[1][:-1])
        assert mt[0] == 0x80000000
        # assume key_length == N

        # compute init_genrand(self, 19650218U)
        ma = [19650218]
        for i in range(1, N):
            # mt[mti] =
            #   (1812433253U * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
            ma.append((1812433253 * (ma[i - 1] ^ (ma[i - 1] >> 30)) + i) & 0xFFFFFFFF)

        # for the second loop, reverse the process
        i = 2
        mt[0] = mt[N - 1]
        for k in range(1, N):
            # if (i>=N) { mt[0] = mt[N-1]; i=1; }
            if i == 1:
                i = N
            # i++;
            # mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941U))
            #      - (uint32_t)i; /* non linear */
            i -= 1
            mt[i] += i
            mt[i] ^= (mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941
            mt[i] &= 0xFFFFFFFF

        # for the first loop, reverse the process
        init_key = []
        i = 2
        j = N
        mt[0] = mt[N - 1]
        for k in range(N):
            # if (i>=N) { mt[0] = mt[N-1]; i=1; }
            if i == 1:
                i = N
            # i++; j++;
            j -= 1
            i -= 1
            # mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525U))
            #      + init_key[j] + (uint32_t)j; /* non linear */
            m = mt[i] - j
            # what is mt[i] before this assignment?
            t = mt[i] if k == 0 else ma[i]
            m -= t ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525)
            if i == 1:
                # handle mt[0] = mt[N-1];
                # recover mt[0] to its initial value
                mt[0] = ma[0]
            # got key
            init_key.append(m & 0xFFFFFFFF)

        # see random_seed in _randommodule.c
        # 32 bit chunks from the right
        seed = 0
        for key in init_key:
            seed = (seed << 32) | key

        return seed
