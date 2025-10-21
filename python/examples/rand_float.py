import bvsolver
import random

# assume we got some random numbers from the remote
# https://github.com/qxxxb/ctf/tree/master/2021/zh3r0_ctf/real_mersenne
all_nums = [random.random() for i in range(10000)]
known = all_nums[:623]

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
    # see _random_Random_random_impl in _randommodule.c
    s = int(value * 9007199254740992.0)
    a = s // 67108864
    b = s % 67108864
    zeros.append((rng.getrandbits(32) >> 5) ^ a)
    zeros.append((rng.getrandbits(32) >> 6) ^ b)
# find solutions
for sol in solver.solve(zeros):
    print("Found solution")
    state_recover = [s.get(sol) for s in state]

    # recreate random number generator
    rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
    for i in range(len(all_nums)):
        assert rng_recover.random() == all_nums[i]
