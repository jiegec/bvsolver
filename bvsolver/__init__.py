from __future__ import annotations
from collections.abc import Generator
import functools

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
    _bits: list[int]

    def __init__(self, bits: list[int]) -> None:
        self._bits = bits.copy()

    def __repr__(self) -> str:
        bits = ", ".join([str(bit) for bit in self._bits])
        return f"<BitVector [{bits}]>"

    def __xor__(self, other: BitVector | int) -> BitVector:
        res: list[int] = []
        if isinstance(other, BitVector):
            assert len(self._bits) == len(other._bits)
            for l, r in zip(self._bits, other._bits):
                res.append(l ^ r)
        else:
            assert other.bit_length() <= len(self._bits)
            for i, b in enumerate(self._bits):
                res.append(b ^ (1 if other & (1 << i) != 0 else 0))
        return BitVector(res)

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
        missing = []
        for i in range(1, self._size + 1):
            # find one entry with the bit set
            found = False
            for j in range(len(equations)):
                if equations[j] & (1 << i) != 0:
                    # found
                    # eliminate others
                    for k in range(j + 1, len(equations)):
                        if equations[k] & (1 << i) != 0:
                            equations[k] ^= equations[j]

                    found = True
                    break
            if not found:
                missing.append(i)

        # recover solution from equations
        res = 0
        for e in equations:
            count = e.bit_count()
            if count == 2:
                # the bit is one, one of them is the bit 0
                assert e & 1 == 1
                res |= e ^ 1
            elif count == 1:
                if e == 1:
                    # no solution
                    return None
                # otherwise, the bit is zero, do nothing
                pass
            else:
                assert False

        # recover all solutions
        num_sols = 1 << len(missing)
        sols = 0
        base_res = res
        while sols < num_sols:
            yield res

            # next solution
            sols += 1
            res = base_res
            for i in range(len(missing)):
                if sols & (1 << i) != 0:
                    res |= 1 << missing[i]
        return None
