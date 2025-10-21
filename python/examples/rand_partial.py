import bvsolver
import random

# assume we got some random numbers from the remote
# 2025 qiangwangbei quals crypto ezran
known = []
for i in range(3200):
    r1 = random.getrandbits(8)
    r2 = random.getrandbits(16)
    x = (pow(r1, 2 * i, 257) & 0xFF) ^ r2
    known.append(x)
rest = [random.getrandbits(32) for i in range(1000)]

# create equations
solver = bvsolver.Solver([32] * 624)
state = solver.bitvectors()

rng = bvsolver.CPythonRandom(
    state,
    # newly created Random instance has its index equal to 624
    624,
)
# newly created Random instance has its first element equal to 0x80000000
# see init_by_array in _randommodule.c
zeros = [state[0] ^ 0x80000000]
for value in known:
    rng.getrandbits(8)
    zeros.append((rng.getrandbits(16) >> 8) ^ (value >> 8))

# find solutions
for sol in solver.solve(zeros):
    print("Found solution")
    state_recover = [s.get(sol) for s in state]

    # recreate random number generator
    rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
    for i in range(len(known)):
        r1 = rng_recover.getrandbits(8)
        r2 = rng_recover.getrandbits(16)
        x = (pow(r1, 2 * i, 257) & 0xFF) ^ r2
        assert x == known[i]

    for i in range(len(rest)):
        assert rng_recover.getrandbits(32) == rest[i]
