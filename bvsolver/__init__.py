from __future__ import annotations
from collections.abc import Generator
import functools
import random
from typing import Generic, TypeVar, overload

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
        assert other.bit_length() <= len(self._bits)
        for i in range(other.bit_length()):
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

    def solve(self, zeros: list[BitVector]) -> Generator[int]:
        # flatten
        flatten: list[int] = []
        for zero in zeros:
            flatten += zero._bits

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
        base_res = 1 # constant term
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

        # we got 2 ** len(free) solutions
        # recover all solutions
        num_sols = 1 << len(free)
        sols = 0
        while sols < num_sols:
            # compute solution
            # assign free variables
            res = base_res
            for i in range(len(free)):
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


# some constants from _randommodule.c
N = 624
M = 397
MATRIX_A = 0x9908B0DF
UPPER_MASK = 0x80000000
LOWER_MASK = 0x7FFFFFFF


class CPythonRandom(Generic[T]):
    _index: int
    _state: list[T]  # 624 32-bit numbers

    def __init__(self, state: list[T], index: int) -> None:
        assert len(state) == N
        self._state = state.copy()
        self._index = index

    def genrand_uint32(self) -> T:
        # see genrand_uint32 from _randommodule.c
        if self._index >= N:
            for kk in range(N):
                y = (self._state[kk] & UPPER_MASK) ^ (
                    self._state[(kk + 1) % N] & LOWER_MASK
                )
                # mag01[y & 0x1U]: if y & 0x1U == 0, then 0; else MATRIX_A
                mag01 = (y & 1) * MATRIX_A
                self._state[kk] = self._state[(kk + M) % N] ^ (y >> 1) ^ mag01
            self._index = 0

        y = self._state[self._index]
        self._index += 1
        y = y ^ (y >> 11)
        y = y ^ (y << 7) & 0x9D2C5680
        y = y ^ (y << 15) & 0xEFC60000
        y = y ^ (y >> 18)
        return y

    def getrandbits(self, bits) -> T:
        # _random_Random_getrandbits_impl in _randommodule.c
        if bits <= 32:
            return self.genrand_uint32() >> (32 - bits)
        else:
            assert False

    def to_cpython_random(self) -> random.Random:
        res = random.Random(0)
        res.setstate((3, (*self._state, self._index), None))
        return res

    @classmethod
    def from_cpython_random(cls, rng: random.Random) -> CPythonRandom:
        state = rng.getstate()
        res = CPythonRandom(list(state[1][:-1]), state[1][-1])
        return res
