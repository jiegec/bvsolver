from __future__ import annotations
from collections.abc import Generator
import functools
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
        equations: list[int] = []
        for zero in zeros:
            equations += zero._bits

        # gauss elimination
        free = []  # free variables
        done = 0
        for i in range(1, self._size + 1):
            # find one entry with the bit set
            for j in range(done, len(equations)):
                if equations[j] & (1 << i) != 0:
                    # found
                    # eliminate others
                    for k in range(0, len(equations)):
                        if j != k and equations[k] & (1 << i) != 0:
                            equations[k] ^= equations[j]

                    # swap equations[j] and equations[done]
                    equations[j], equations[done] = equations[done], equations[j]
                    done += 1

                    break
        # now in reduced row echelon form
        # find more free variables
        # some free variables do not appear in any equation
        appearing = 0
        for e in equations:
            if e == 1:
                # no solution
                return None
            # ignore the constant term
            if e & 1 == 1:
                e ^= 1
            appearing |= e
            count = e.bit_count()
            if count >= 2:
                # correlated free variables
                # add all except the last one
                e &= e - 1
                while e != 0:
                    free.append(e & -e)
                    e &= e - 1

        if appearing.bit_count() != self._size:
            # count not appearing
            for i in range(1, self._size + 1):
                if appearing & (1 << i) == 0:
                    free.append(1 << i)

        # now all free variables are found
        print(equations, done, free)
        assert done + len(free) == self._size

        # we got 2 ** len(free) solutions
        # recover all solutions
        num_sols = 1 << len(free)
        sols = 0
        while sols < num_sols:
            # compute solution
            # assign free variables
            res = 0
            for i in range(len(free)):
                if sols & (1 << i) != 0:
                    res |= free[i]

            # solve other variables given the known free variables
            # recover solution from equations
            for e in equations:
                # find low bit in e
                e_drop_1 = e ^ 1 if e & 1 == 1 else e
                lowbit = e_drop_1 & (-e_drop_1)
                # e.g. e=0b101, and res=0b000, then lowbit=0b100 must be set
                # e.g. e=0b11101, and res=0b01000, then lowbit=0b100 must not be set
                if (e_drop_1 & res).bit_count() & 1 != e & 1:
                    # the low bit must be set to make equation equal
                    res |= lowbit
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

    def __init__(self, state: list[T]) -> None:
        assert len(state) == N
        self._state = state
        self._index = 0  # the initial index does not matter

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
        y = y ^ (y << 18)
        return y

    def getrandbits(self, bits) -> T:
        # TODO
        assert bits == 32
        return self.genrand_uint32()
